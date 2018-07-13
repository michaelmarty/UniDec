import numpy as np
import scipy.optimize as opt

from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


def make_mass_list(massdat, arrayin, psfun, startarray, *args):
    """
    For parameters in arrayin, make mass peaks to be fit to massdat.
    Arrayin must be P x 3 array of (mid, fwhm, area).
    :param massdat: Data to be fit
    :param arrayin: Parameters defining peaks
    :param psfun: Peak shape function integer code
    :param startarray: Starting array (not currently used but can be used to restrict some guesses)
    :param args: Extra arguments (such as nonorm to prevent normalization)
    :return:
    """
    length = len(massdat)
    output = np.zeros(length)
    num = len(arrayin.flatten()) / 3
    array = np.reshape(arrayin, (int(num), 3))

    for i in range(0, len(array)):
        output += ud.make_peak_shape(massdat[:, 0], psfun, array[i, 1], array[i, 0], norm_area=True) * array[i, 2]
    if np.amax(output) != 0 and "nonorm" in args:
        output = output / np.amax(output) * np.amax(massdat[:, 1])
        # if "smallguess" in args:
        #    if np.any(array[:,1])>3:
        #        output=output*0
    if np.any(array[:, 1:] < 0):
        output *= 0.1
    # if (np.any(array[:,1]>startarray[:,1]*100000) or np.any(np.ravel(array)<0)):
    #    output=0*output
    return output


def error_function(array, massdat, psfun, startarray, *args):
    """
    Error function for least_squares_minimize
    :param array: Array of test parameters
    :param massdat: Mass data to be fit
    :param psfun: Peak shape function integer code
    :param startarray: Starting array (not currently used but can be used to restrict some guesses)
    :param args: Extra arguments for make_mass_list
    :return: Squared errors array
    """
    error = (massdat[:, 1] - make_mass_list(massdat, array, psfun, startarray, *args)) ** 2
    return error


def least_squares_minimize(massdat, array, psfun, *args):
    """
    Perform least squares minimization of peaks defined in array to massdat.
    :param massdat: Data to fit
    :param array: Array of parameters for defining peaks that will be fit
    :param psfun: Peak shape function integer code
    :param args: Extra arguments for make_mass_list
    :return: Best fit of values in array
    """
    fit = opt.leastsq(error_function, np.ravel(array), args=(massdat, psfun, array, args))[0]
    num = len(array.flatten()) / 3
    fit = np.reshape(fit, (int(num), 3))
    return fit


class MassFitter:
    """
    Class for fitting zero-charge mass spectrum to overlapping peaks.
    """

    def __init__(self, massdat, guessarray, psfun, *args):
        """
        Initialize values and set initial guesses.

        Guess array can be fed in several things. If there are five columns, it expects a P x 5 grid.

        [mass, mass fwhm, charge (not used), charge std (not used), intensity]

        If there are only two columns, it expects a P x 2 grid

        [mass, intensity]

        :param massdat: Mass data to be fit
        :param guessarray: Array of initial peak parameters
        :param psfun: Peak shape function integer code
        :param args: Arguments passed to fits and guesses.
        :return: fitdat
        """
        self.massdat = massdat
        self.finarray = guessarray
        self.psfun = psfun
        self.fit = None
        self.fitdat = None
        try:
            self.initguess = np.array([[self.finarray[i, 0], self.finarray[i, 1] * 1,
                                        self.finarray[i, 4] * self.finarray[i, 1] * 1 * np.sqrt(2 * np.pi)] for i in
                                       range(0, len(self.finarray))])
        except IndexError:
            if "smallguess" in args:
                self.initguess = np.array(
                    [[self.finarray[i, 0], 0.5, self.finarray[i, 1]] for i in range(0, len(self.finarray))])
            elif "microguess" in args:
                self.initguess = np.array(
                    [[self.finarray[i, 0], 0.1, self.finarray[i, 1]] for i in range(0, len(self.finarray))])
                self.initguess[len(self.initguess)-1,1]=10
            else:
                self.initguess = np.array(
                    [[self.finarray[i, 0], 500., self.finarray[i, 1]] for i in range(0, len(self.finarray))])
        print("Inital Guess: ", self.initguess)

    def perform_fit(self, *args):
        """
        Run least squares fitting
        :param args: Arguments passed to fitting
        :return: fitdat, fit (fit to data and fit parameters in a P x 3 array of (mid, fwhm, area))
        """
        self.fit = least_squares_minimize(self.massdat, self.initguess, self.psfun, *args)
        self.fitdat = make_mass_list(self.massdat, self.fit, self.psfun, self.finarray, *args)
        if "sort" in args:
            self.fit = self.fit[self.fit[:, 0].argsort()]
        return self.fitdat, self.fit

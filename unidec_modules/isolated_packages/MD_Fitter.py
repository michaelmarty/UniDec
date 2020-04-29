import numpy as np
import os
import unidec_modules.unidectools as ud
import scipy.optimize as opt


def make_mass_list(massdat, arrayin, psfun, posarray, nparam=3, maxshift=None, *args):
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
    num = len(arrayin.flatten()) / nparam
    array = np.reshape(arrayin, (int(num), nparam))

    for i in np.arange(0, len(array)):
        output += ud.make_peak_shape(massdat[:, 0], psfun, array[i, 1], array[i, 0], norm_area=True) * array[i, 2]
    if np.amax(output) != 0:
        output = output / np.amax(output) * np.amax(massdat[:, 1])
    if np.any(array[:, 1:] < 0):
        output *= 0.1
    if maxshift is not None:
        # if np.abs(shift)>maxshift:
        #    output *=0.1
        if np.any(np.abs(array[:, 0] - posarray) > maxshift):
            output *= 0.1
    return output


def error_function(array, massdat, psfun, posarray, nparam, maxshift, *args):
    error = (massdat[:, 1] - make_mass_list(massdat, array, psfun, posarray, nparam, maxshift, *args)) ** 2
    return error

def error_function2(array, massdat, psfun, posarray, nparam, maxshift, *args):
    error = (massdat[:, 1] - make_mass_list(massdat, array, psfun, posarray, nparam, maxshift, *args)) ** 2
    error = np.sum(error)#/np.sum(massdat[:,1])
    #num = len(array.flatten()) / nparam
    #array = np.reshape(array, (int(num), nparam))
    #error2 = (array[:,0]-posarray)**2
    #error2 = np.sum(error2)/np.sum(posarray)
    return error


def least_squares_minimize(massdat, array, psfun, posarray, nparam=3, maxshift=None, *args):
    """
    Perform least squares minimization of peaks defined in array to massdat.
    :param massdat: Data to fit
    :param array: Array of parameters for defining peaks that will be fit
    :param psfun: Peak shape function integer code
    :param args: Extra arguments for make_mass_list
    :return: Best fit of values in array
    """

    #fit = opt.leastsq(error_function, np.ravel(array), args=(massdat, psfun, posarray, nparam, maxshift, args))[0]
    fit = opt.minimize(error_function2, np.ravel(array), args=(massdat, psfun, posarray, nparam, maxshift, args), method="trust-constr").x


    num = len(array.flatten()) / nparam
    fit = np.reshape(fit, (int(num), nparam))
    return fit


def MD_Fitter(data, mds=None, maxshift=0.1, widthguess=0.05, shiftguess=None, plot=False):
    data[:, 1] = data[:, 1] - np.amin(data[:, 1])
    data[:, 1] = data[:, 1] / np.amax(data[:, 1])
    peaks = ud.peakdetect(data, threshold=0.1)
    #print(peaks, mds)
    mds = np.array(mds)

    if shiftguess is None:
        a, b = np.meshgrid(peaks[:, 0], mds)
        minarg = np.argmin(np.abs(a - b))
        minarg = np.unravel_index(minarg, a.shape)
        shiftguess = a - b
        shiftguess = shiftguess[minarg]
        if np.abs(shiftguess) > maxshift:
            print("Shift Guess too High:", shiftguess, maxshift)
            shiftguess = 0
        else:
            print("Guess for shift:", shiftguess)

    widths = [0.75, 1, 1.25]
    #widths=[1]
    error = 100000000000
    finalfit = None
    finalfitdat = None
    #shifts = [0, shiftguess, -1*shiftguess]
    shifts=[shiftguess]
    for w in widths:
        for s in shifts:
            widths = np.ones(len(mds)) * w * widthguess

            heights = []
            for m in mds:
                height = ud.data_extract(data, m + s, 0)
                heights.append(height)
            heights = np.array(heights)

            try:
                guessarray = np.transpose([mds + s, widths, heights]).ravel()
            except Exception as e:
                print(e)
                print(mds, s, widths, heights)
                guessarray = np.transpose([mds, widths, heights]).ravel()

            fit = least_squares_minimize(data, guessarray, 0, mds, nparam=3, maxshift=maxshift)
            fitdat = make_mass_list(data, fit, 0, mds, nparam=3)

            newerror = np.sum((data[:, 1] - fitdat) ** 2)
            if newerror < error:
                error = newerror
                finalfit = fit
                finalfitdat = fitdat
    fit = finalfit
    fitdat = finalfitdat
    print(fit, error)

    if plot:
        import matplotlib.pyplot as plt
        plt.subplot(121)
        plt.plot(data[:, 0], data[:, 1], color="k")
        plt.plot(data[:, 0], fitdat, linestyle="--", color="r")
        for i, m in enumerate(mds):
            fdat = make_mass_list(data, fit[i], 0, mds, nparam=3)
            plt.plot(data[:,0], fdat*fit[i,2], label=str(i))
        plt.legend()
        plt.subplot(122)
        plt.bar(np.arange(len(mds)), fit[:, -1])
        plt.show()

    return fit, fitdat


if __name__ == "__main__":
    dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\Peptide_Paper2\\LactoferricinB\\Mass Defect Extracts"
    file = "20191216_DMPG_LactoferricinBwZero_3.0_1D_Mass_Defects.txt"
    file = "20191217_DMPG_LactoferricinB_3.0_1D_Mass_Defects.txt"
    file = "20191219_DMPG_LactoferricinB_9.0_1D_Mass_Defects.txt"
    file = "20191219_DMPG_LactoferricinB_12.0_1D_Mass_Defects.txt"
    path = os.path.join(dir, file)

    data = np.loadtxt(path)

    centermode = 1
    refmass = 667
    basemass = 44088
    omass = 3123.9
    range = [0, 12]

    indexes = np.arange(range[0], range[1])
    masses = basemass + omass * indexes
    kmasses = masses / float(refmass)
    if centermode == 1:
        nominalkmass = np.floor(kmasses)
    else:
        nominalkmass = np.round(kmasses)
    mds = kmasses - nominalkmass

    MD_Fitter(data, mds, widthguess=0.05, maxshift=0.05, shiftguess=-0.02, plot=True)

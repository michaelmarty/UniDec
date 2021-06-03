import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
from scipy import special
import scipy.optimize as opt


def poisson(x, mu, A):
    return stats.poisson.pmf(x, mu) * A


def poisson_fit(xvals, yvals):
    xvals = np.array(xvals)
    Aguess = np.amax(yvals)
    mguess = xvals[np.argmax(yvals)]
    guess = [mguess, Aguess]
    fits = curve_fit(poisson, xvals, yvals, p0=guess, maxfev=1000000)[0]
    fitdat = poisson(xvals, fits[0], fits[1])
    return fits, fitdat


def binomial(x, n, p, A):
    return stats.binom.pmf(x, n, p) * A


def binomial_fit(xvals, yvals):
    xvals = np.array(xvals)
    Aguess = np.amax(yvals)
    nguess = np.amax(xvals)
    pguess = float(xvals[np.argmax(yvals)]) / float(nguess)
    guess = [nguess, pguess, Aguess]
    print(guess)
    fits = curve_fit(binomial, xvals, yvals, p0=guess, maxfev=1000000)[0]
    fitdat = binomial(xvals, fits[0], fits[1], fits[2])
    return fits, fitdat


def ndis_std(x, mid, sig, a=1, norm_area=False):
    """
    Normal Gaussian function normalized to the max of 1.
    :param x: x values
    :param mid: Mean of Gaussian
    :param sig: Standard Deviation
    :param a: Maximum amplitude (default is 1)
    :param norm_area: Boolean, Whether to normalize so that the area under the distribution is 1.
    :return: Gaussian distribution at x values
    """
    if norm_area:
        a *= 1 / (sig * np.sqrt(2 * np.pi))
    return a * np.exp(-(x - mid) * (x - mid) / (2.0 * sig * sig))


def ndis(x, mid, fwhm, **kwargs):
    """
    Gaussian function normalized to a max of 1.

    Note: x and mid are interchangable. At least one should be a single float. The other may be an array.
    :param x: x values
    :param mid: Mean
    :param fwhm: Full Width at Half Max (2.35482 * standard deviation)
    :param kwargs: Allows norm_area flag to be passed
    :return: Gaussian distribution at x values
    """
    sig = fwhm / 2.35482
    return ndis_std(x, mid, sig, **kwargs)


def ndis_fit(x, s, m, a, b):
    """
    Function for fitting normal distribution to peak.
    Adds a background to ndis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return ndis(x, m, s, a=a, norm_area=True) + b


def ldis(x, mid, fwhm, a=1, norm_area=False):
    """
    Lorentzian function normalized to a max of 1.
    Note: x and mid are interchangable. At least one should be a single float. The other may be an array.
    :param x: x values
    :param mid: Mean
    :param fwhm: Full Width at Half Max
    :param a: Amplitude (default is 1)
    :param norm_area: Boolean, Whether to normalize so that the area under the distribution is 1.
    :return: Lorentzian distribution at x values
    """
    if norm_area:
        a *= ((1 / np.pi) * (fwhm / 2.))
    else:
        a *= ((fwhm / 2.0) * (fwhm / 2.0))
    return a / ((x - mid) * (x - mid) + (fwhm / 2.0) * (fwhm / 2.0))


def ldis_fit(x, s, m, a, b):
    """
    Function for fitting Lorentzian distribution to peak.
    Adds a background to ldis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return ldis(x, m, s, a=a, norm_area=True) + b


def splitdis(x, mid, fwhm, a=1, norm_area=False):
    """
    Split Gaussain/Lorentzian function normalized to a max of 1.

    Gaussian < mid. Lorentzian > mid.

    :param mid: Mid point (point of peak intensity)
    :param x: x value or values
    :param fwhm: Full Width at Half Max
    :return: Split Gaussian/Lorentzian distribution at x value
    """
    sig2 = fwhm / (2 * np.sqrt(2 * np.log(2)))
    if norm_area:
        a1 = a * ((1 / np.pi) / (fwhm / 2.)) / 0.83723895067
        a2 = a * 2. / (fwhm * np.pi) / 0.83723895067
    else:
        a1 = a
        a2 = a
    try:
        if mid < x:
            return ldis(x, mid, fwhm, a=a1)
        else:
            return ndis_std(x, mid, sig2, a=a2)
    except ValueError:
        output = np.zeros_like(x)
        output[x > mid] = ldis(x[x > mid], mid, fwhm, a=a1)
        output[x <= mid] = ndis_std(x[x <= mid], mid, sig2, a=a2)
        return output


def splitdis_fit(x, s, m, a, b):
    """
    Function for fitting Split G/L distribution to peak.
    Adds a background to splitdis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return splitdis(x, m, s, a=a, norm_area=True) + b


def voigt(x, mu=0, sigma=1, gamma=1, amp=1, background=0):
    """\
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))
    """
    if sigma == 0:
        return ldis(x, mu, gamma * 2., amp) + background
    elif gamma == 0:
        return ndis_std(x, mu, sigma, amp) + background
    else:
        z = (x - mu + 1j * gamma) / (sigma * np.sqrt(2))
        w = special.wofz(z)
        v = w.real / (sigma * np.sqrt(2 * np.pi))
        v *= (amp / np.amax(v))
        v += background
    return v


def exp_decay(x, gamma=1, a=1, background=0):
    return np.exp(-gamma * x) * a + background


def exp_fit(xvals, yvals, gguess=None, aguess=None, bguess=None):
    """
    Simple exponential fitting function.
    """
    xvals = np.array(xvals)
    if gguess is None:
        gguess = 2. / np.mean(xvals)
    if aguess is None:
        aguess = np.amax(yvals)
    if bguess is None:
        guess = [gguess, aguess]
    else:
        guess = [gguess, aguess, bguess]
    fits = curve_fit(exp_decay, xvals, yvals, p0=guess, maxfev=1000000)[0]
    if bguess is None:
        fitdat = exp_decay(xvals, gamma=fits[0], a=fits[1])  # , background=fits[2])
    else:
        fitdat = exp_decay(xvals, gamma=fits[0], a=fits[1], background=fits[2])
    return fits, fitdat


def lin_fit(xvals, yvals):
    xvals = np.array(xvals)
    outputs = stats.linregress(xvals, yvals)
    slope, intercept, rvalue, pvalue, slope_std_error = outputs
    fitdat = xvals * slope + intercept
    return np.array([slope, intercept]), fitdat


def logistic(x, mid, slope, aspan, ymin):
    """Logistic function
    """
    y = (aspan / (1 + np.exp(-slope * (x - mid)))) + ymin
    return y


def sig_fit(xvals, yvals):
    xvals = np.array(xvals)
    yguess = np.amin(yvals)
    Aguess = np.amax(yvals) - np.amin(yvals)
    mguess = np.mean(xvals)
    if Aguess != 0:
        sguess = (yvals[len(yvals) - 1] - yvals[0]) / (xvals[len(xvals) - 1] - xvals[0]) / Aguess
    else:
        Aguess = 1
        sguess = 1
    guess = [mguess, sguess, Aguess, yguess]
    fits = curve_fit(logistic, xvals, yvals, p0=guess, maxfev=10000000)[0]
    fitdat = logistic(xvals, fits[0], fits[1], fits[2], fits[3])
    return fits, fitdat


def gaussfit(xvals, yvals, mguess=None, sguess=0.1, aguess=None, cleanup=True):
    """
    Simple gaussian fitting function.
    :param xvals: X values
    :param yvals: Y values
    :param mguess: Guess for midpoint
    :param sguess: Guess for standard deviation
    :param aguess: Guess for amplitude
    :param cleanup: Boolean Flag, if True, will clean up data to remove background and normalize
    :return: Gaussian fit parameters [mid, sig, amp]
    """
    if cleanup:
        yvals -= np.amin(yvals)
        yvals /= np.amax(yvals)
    if mguess is None:
        mguess = xvals[np.argmax(yvals)]
    if aguess is None:
        aguess = np.amax(yvals)
    guess = [mguess, sguess, aguess]
    fits = curve_fit(ndis_std, xvals, yvals, p0=guess, maxfev=1000000)[0]
    print(fits)
    return fits


def psfit(x, s, m, a=1, b=0, psfun=0):
    """
    Make peak shape from fit
    :param x: x values
    :param s: fwhm
    :param m: max position
    :param a: amplitude
    :param b: background
    :param psfun: peak shape function integer code
    :return: peak shape fit data
    """
    if psfun == 0:
        return ndis_fit(x, s, m, a, b)
    elif psfun == 1:
        return ldis_fit(x, s, m, a, b)
    elif psfun == 2:
        return splitdis_fit(x, s, m, a, b)


def voigt_fit(xvals, yvals, mguess=0, sguess=0.1, gguess=0, aguess=0, bguess=0):
    """

    """
    guess = [mguess, sguess, gguess, aguess, bguess]
    popt, pcov = curve_fit(voigt, xvals, yvals, p0=guess, maxfev=1000000)
    fitdat = voigt(xvals, popt[0], popt[1], popt[2], popt[3], popt[4])
    return popt, np.sqrt(np.diag(pcov)), fitdat


def fit_peak(xvals, yvals, psfun, midguess, fwhmguess, aguess, bguess):
    """
    Fit peak from xvals and yvals data to defined peak shape function.
    :param xvals: x values of data
    :param yvals: y values of data
    :param psfun: peak shape function integer code
    :param midguess: midpoint guess
    :param fwhmguess: fwhm guess
    :param aguess: amplitude guess
    :param bguess: background guess
    :return: popt, perr, fitdat (optimized parameters [fwhm, mid, a, b], std error of parameters, fit to data)
    """
    guess = [fwhmguess, midguess, aguess, bguess]
    if psfun == 0:
        popt, pcov = curve_fit(ndis_fit, xvals, yvals, p0=guess)
    elif psfun == 1:
        popt, pcov = curve_fit(ldis_fit, xvals, yvals, p0=guess)
    elif psfun == 2:
        popt, pcov = curve_fit(splitdis_fit, xvals, yvals, p0=guess)
    else:
        popt = guess
        pcov = np.ones((len(guess), len(guess)))
        print("Failed")

    fitdat = psfit(xvals, popt[0], popt[1], popt[2], popt[3], psfun)
    try:
        out2 = np.sqrt(np.diag(pcov))
    except:
        out2 = 0
    return popt, out2, fitdat


def weighted_std_2(values, weights):
    """
    Calculate weighted standard deviation.
    :param values: Values
    :param weights: Weights
    :return: Weighted standard deviation.
    """
    average = np.average(values, weights=weights)
    variance = np.average((np.array(values) - average) ** 2, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)


def isolated_peak_fit(xvals, yvals, psfun, **kwargs):
    """
    Fit an isolated peak to the peak shape model.
    :param xvals: x values of data
    :param yvals: y values of data
    :param psfun: peak shape function integer code
    :param kwargs: keywords (unused)
    :return: fit_array, fit_data (4 x 2 array of (fit, error) for [fwhm,mid,amp,background],fit to data)
    """
    midguess = xvals[np.argmax(yvals)]
    bguess = np.amin(yvals)
    sigguess = weighted_std_2(xvals, yvals - bguess) * 1
    # Two rounds to guess at area
    if psfun < 3:
        testdat = psfit(xvals, sigguess, midguess, 1, bguess, psfun)
        aguess = np.amax(yvals) / np.amax(testdat)
        testdat = psfit(xvals, sigguess, midguess, aguess, bguess, psfun)
        aguess = aguess * np.amax(yvals) / np.amax(testdat)
    else:
        testdat = psfit(xvals, sigguess, midguess, 1, bguess, 0)
        aguess = np.amax(yvals) / np.amax(testdat)
        testdat = psfit(xvals, sigguess, midguess, aguess, bguess, 0)
        aguess = aguess * np.amax(yvals) / np.amax(testdat)
    # Fit it
    if psfun < 3:
        fit, err, fitdat = fit_peak(xvals, yvals, psfun, midguess, sigguess, aguess, bguess)
    else:
        fit, err, fitdat = voigt_fit(xvals, yvals, midguess, sigguess, 0, aguess, bguess)
    return np.transpose([fit, err]), fitdat


def poly_fit(datatop, degree=1, weights=None):
    x = datatop[:, 0]
    y = datatop[:, 1]
    fit = np.polyfit(x, y, degree, w=weights)
    p = np.poly1d(fit)
    fitdat = p(x)
    return fit, fitdat


def poisson_mono_di(datatop):
    bool1 = datatop[:, 0] % 2 == 0
    odds = datatop[np.logical_not(bool1)]
    fit, fitdata = poisson_fit(odds[:, 0], odds[:, 1])
    fitdat = poisson(datatop[:, 0], fit[0], fit[1])
    evens = datatop[bool1]
    evens[:, 1] = evens[:, 1] - fitdat[bool1]
    fit2, fitdata2 = poisson_fit(evens[:, 0], evens[:, 1])
    fitdatatot = fitdat + poisson(datatop[:, 0], fit2[0], fit2[1]) * bool1
    print(fit[1], fit2[1])
    return fit, fit2, fitdatatot


def multipoisson(array, datatop, oarray, background=False):
    l = len(oarray)
    fitdat = np.zeros_like(datatop[:, 1])
    integrals = np.empty_like(oarray).astype(np.float)
    integrals2 = np.empty_like(oarray).astype(np.float)
    for i in range(0, l):
        oligomer = oarray[i]
        mu = array[i * 2]
        A = array[i * 2 + 1]
        bool1 = datatop[:, 0] % oligomer == 0
        fdat = poisson(datatop[:, 0], mu, A) * bool1
        fitdat += fdat
        integrals[i] = np.sum(fdat)
        integrals2[i] = np.sum(fdat * datatop[:, 0])
    if background:
        fitdat += array[-1]
    return fitdat, integrals, integrals2


def mpinit(datatop, oarray, background=False):
    l = len(oarray)
    array = []
    for i in range(0, l):
        oligomer = oarray[i]
        bool1 = datatop[:, 0] % oligomer == 0
        data = datatop[bool1]
        if len(data) > 3:
            fit, fitdat = poisson_fit(data[:, 0], data[:, 1])
            array.append(fit[0])
            array.append(fit[1])
        else:
            array.append(np.amin(data[:, 0]))
            array.append(np.amax(data[:, 1]))
    if background:
        array.append(0)
    return array


def mperror(array, datatop, oarray, background):
    bool1 = array < -1
    array[bool1] = 0
    error = (datatop[:, 1] - multipoisson(array, datatop, oarray, background)[0]) ** 2.
    return error


def complex_poisson(datatop, oarray=[1, 2], background=False):
    array = mpinit(datatop, oarray, background)
    print(array)
    fit = opt.leastsq(mperror, array, args=(datatop, oarray, background))[0]
    fitdat, integrals, integrals2 = multipoisson(fit, datatop, oarray, background)
    return fit, fitdat, integrals, integrals2


if __name__ == "__main__":
    path = "C:\\UniDecPastedSpectra\PastedSpectrum_2017_Dec_11_09_02_49_unidecfiles\Extract_total_2D_Extract.txt"
    path = "C:\\UniDecPastedSpectra\PastedSpectrum_2017_Dec_11_11_30_45_unidecfiles\Extract_total_2D_Extract.txt"
    data = np.loadtxt(path)[1:18]
    data[:, 1] -= np.amin(data[:, 1])

    import matplotlib.pyplot as plt

    oarray = [1, 2, 3, 4, 5, 6, 8]
    fits, fitdat, integrals, integrals2 = complex_poisson(data, oarray=oarray, background=False)
    print(fits)
    print(integrals, integrals2)
    amps = integrals2  # fits[1::2]

    plt.bar(oarray, amps / np.amax(amps))
    plt.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]))
    plt.plot(data[:, 0], fitdat / np.amax(data[:, 1]))
    plt.show()

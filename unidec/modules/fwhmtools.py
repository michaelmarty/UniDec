import numpy as np
from numba import njit
import unidec.tools as ud
from copy import deepcopy

@njit(fastmath=True)
def calc_sigma(data):
    """
    https://arxiv.org/pdf/1907.07241
    """
    s = 0
    x1, y1 = data[0]
    for d in data[1:]:
        x2, y2 = d
        s += (y1 + y2) * (x2 - x1) / 2.0
        x1, y1 = x2, y2

    # diffs = np.diff(data[:,0])
    # ys = data[:-1,1]
    # ys2 = data[1:,1]
    # ys = (ys + ys2) / 2.0  # Average the y values
    # sum = np.sum(ys * diffs)
    return s/(np.sqrt(2 * np.pi) * np.amax(data[:, 1]))


@njit(fastmath=True)
def calc_ab(data, c):
    b1 = data[:,1] > np.amax(data[:,1]) * 0.25  # Remove points with very low intensity
    data = data[b1]  # Remove zero intensity points
    ys = data[:, 1]**2
    m11 = np.sum(ys)
    m12 = np.sum(data[:,0] * ys)
    m22 = np.sum(data[:,0]**2 * ys)

    M = np.array([[m11, m12], [m12, m22]])

    logdat = np.log(data[:, 1])

    v1 = np.sum(ys * logdat) - c * m22
    v2 = np.sum(data[:,0] * ys * logdat) - c * np.sum(data[:,0] **3 * ys)

    v = np.array([v1, v2])

    ab = np.linalg.solve(M, v)
    return ab

@njit(fastmath=True)
def calc_gauss(data):
    """
    Calculate the Gaussian fit parameters for the data.
    https://arxiv.org/pdf/1907.07241
    :param data: Data array with m/z and intensity
    :return: Gaussian fit parameters (amplitude, mean, sigma)
    """
    # mean = np.sum(data[:, 0] * data[:, 1]) / np.sum(data[:, 1])
    sigma = calc_sigma(data)
    # amplitude = np.amax(data[:, 1])
    c = -1 / (2 * sigma ** 2)
    a, b = calc_ab(data, c)
    mean = -b / (2 * c)
    amplitude = np.exp(a - (b ** 2 / (4*c)))
    return mean, sigma, amplitude

@njit(fastmath=True)
def check_tol(x, y, ppmtol):
    if ppmtol is None:
        return True
    diff = np.abs(x - y)
    diffcut = x * ppmtol * 1e-6
    if diff < diffcut:
        return True
    else:
        return False

@njit(fastmath=True)
def linear_interpolation(d1, d2, y):
    """
    Perform linear interpolation between two points.
    :param d1: First data point (x1, y1)
    :param d2: Second data point (x2, y2)
    :param y: y-coordinate where to interpolate
    :return: interpolated x value at y
    """
    x1, y1 = d1
    x2, y2 = d2
    if y1 == y2:
        return (x1 + x2) / 2
    # Calculate the slope (m) and intercept (b) of the line
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    # Calculate the x value at the given y
    x = (y - b) / m
    return x



@njit(fastmath=True)
def find_fwhm(data, index, wfactor=4, maxppmtol=1000, maxasymmetry=4):
    """
    Find the full width at half maximum (FWHM) of a peak in the data.
    :param data: Data array with m/z and intensity
    :param index: Index of the peak in the data
    :param wfactor: Width factor to adjust the amount of indexes included
    :param maxppmtol: Maximum ppm tolerance for peak matching
    :param maxasymmetry: Maximum asymmetry allowed for the peak
    :return: FWHM value
    """
    mz = data[index, 0]
    intensity = data[index, 1]

    fwhmdefault = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

    if intensity == 0:
        return fwhmdefault

    half_max = intensity / 2.0

    # Find left side
    left_index = index
    while left_index > 0 and data[left_index, 1] >= half_max and check_tol(data[left_index, 0], mz, maxppmtol):
        left_index -= 1

    # Find right side
    right_index = index
    while right_index < len(data) - 1 and data[right_index, 1] >= half_max and check_tol(data[right_index, 0], mz, maxppmtol):
        right_index += 1

    leftmz = linear_interpolation(data[left_index], data[left_index + 1], half_max)
    rightmz = linear_interpolation(data[right_index - 1], data[right_index], half_max)
    if rightmz < mz:
        # If the right mz is less than the peak mz, we cannot calculate FWHM properly
        return fwhmdefault
    if leftmz > mz:
        # If the left mz is greater than the peak mz, we cannot calculate FWHM properly
        return fwhmdefault

    # print("Left MZ:", leftmz, "Right MZ:", rightmz, "Peak MZ:", mz)
    leftwidth = mz - leftmz
    rightwidth = rightmz - mz
    fwhm = leftwidth + rightwidth

    centroid = np.sum(data[left_index:right_index + 1, 0] * data[left_index:right_index + 1, 1]) / np.sum(data[left_index:right_index + 1, 1])

    indexwidth = right_index - left_index

    # left_index = ud.nearest(data[:, 0], mz-wfactor*leftwidth)
    # right_index = ud.nearest(data[:, 0], mz+wfactor*rightwidth)
    left_index = index - int(wfactor) * (index-left_index)
    right_index = index + int(wfactor) * (right_index-index)

    if leftwidth == 0 or rightwidth == 0:
        # If either width is zero, we cannot calculate FWHM properly
        return fwhmdefault

    if left_index < 0:
        left_index = 0
    if right_index >= len(data):
        right_index = len(data) - 1

    # If leftwidth is far apart from rightwidth, we can assume the peak is not symmetric
    if maxasymmetry is not None:
        if leftwidth / rightwidth > maxasymmetry or rightwidth / leftwidth > maxasymmetry:
            return fwhmdefault
    sigma = calc_sigma(data[left_index:right_index + 1])

    return [fwhm, leftwidth, rightwidth, left_index, right_index, sigma, centroid, intensity, indexwidth]

@njit(fastmath=True)
def fast_fwhm(data, peaks, sort=False, wfactor=4, maxppmtol=1000):
    """
    Fast FWHM calculation for a set of peaks.
    :param data: Data array with m/z and intensity
    :param peaks: List of peak m/z values
    :return: FWHM values for each peak
    """
    if sort:
        data = data[np.argsort(data[:, 0])]
        peaks = peaks[np.argsort(peaks[:, 0])]
    # Otherwise assume sorted

    # matchedindexes = np.array([ud.nearest(data[:, 0], i) for i in peaks[:, 0]])
    matchedindexes = match_indexes(data[:,0], peaks[:,0])
    fwhmdefault = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    fwhmlist = []
    for i, p in enumerate(peaks):
        if p[1] == 0:
            fwhmlist.append(fwhmdefault)
            continue

        index = matchedindexes[i]
        fwhm = find_fwhm(data, index, wfactor=wfactor, maxppmtol=maxppmtol)

        fwhmlist.append(fwhm)

    return np.array(fwhmlist)


@njit(fastmath=True)
def match_indexes(array, targetarray):
    """
    Assuming that both arrays are sorted, find the indexes of the elements in targetarray that match the elements in array.
    """
    indexes = []
    index = 0
    for t in targetarray:
        while index < len(array)-1 and array[index] < t:
            index += 1

        i1 = array[index]
        i2 = array[index - 1] if index > 0 else 0

        d1 = np.abs(i1 - t)
        d2 = np.abs(i2 - t)
        if d1 < d2:
            indexes.append(index)
        else:
            indexes.append(index - 1)

    return np.array(indexes)

@njit(fastmath=True)
def remove_peaks_from_fwhmlist(data, fwhmlist):
    """
    Remove peaks from the data based on FWHM values.
    :param data: Data array with m/z and intensity
    :param fwhmlist: List of FWHM values for each peak
    :return: Data with peaks removed
    """
    b1 = np.ones(len(data), dtype=np.bool_)
    for i, fwhm in enumerate(fwhmlist):
        if fwhm[0] == -1:
            continue
        s = int(fwhm[3])
        e = int(fwhm[4])
        b1[s:e] = 0
    # print("Removing Peaks, N Peaks:", np.sum(~b1), "Total Data Points:", len(data))
    return data[b1]

@njit(fastmath=True)
def remove_peaks(data, isodist, wfactor=6, maxppmtol=1000):
    print("Data Length:", len(data), "Isodist Length:", len(isodist))
    fwhms = fast_fwhm(data, isodist, wfactor=wfactor, maxppmtol=maxppmtol)
    print("FWHM List:", len(fwhms))
    # print(fwhms)
    leftoverdata = remove_peaks_from_fwhmlist(data, fwhms)
    # print("Leftover Data Length:", len(leftoverdata))
    return leftoverdata

@njit(fastmath=True)
def ndis_std(x: np.ndarray, mid: float, sig: float, a: float = 1.0) -> np.array:
    """
    Normal Gaussian function normalized to the max of 1.
    :param x: x values
    :param mid: Mean of Gaussian
    :param sig: Standard Deviation
    :param a: Maximum amplitude (default is 1)
    :return: Gaussian distribution at x values
    """
    # x = np.array(x).astype(float)
    return a * np.exp(-(x - mid) * (x - mid) / (2.0 * sig * sig))

# @njit(fastmath=True)
def create_peaks(data, isodist=None, fwhmlist=None):
    if fwhmlist is None:
        fwhmlist = fast_fwhm(data, isodist)
    if fwhmlist is None and isodist is None:
        raise ValueError("Either isodist or fwhmlist must be provided.")

    xvals = data[:,0]
    yvals = np.zeros_like(xvals)
    for i, fwhm in enumerate(fwhmlist):
        if fwhm[0] <= 0:
            continue
        # Generate the peak using the isodist and FWHM
        # peakmz = isodist[i, 0]
        # peakint = isodist[i, 1]

        # Pull the sigma and start/end indexes from the FWHM list
        sigma = fwhm[5]
        mean = fwhm[6]
        amp = fwhm[7]
        startindex = int(fwhm[3])
        endindex = int(fwhm[4])

        xvalues = xvals[startindex:endindex + 1]
        yvalues = ndis_std(xvalues, float(mean), float(sigma), float(amp))

        yvals[startindex:endindex + 1] += yvalues

    return np.column_stack((xvals, yvals))

def create_simdata(data, fwhmlist):
    xvals = data[:,0]
    simdata = []
    for i, fwhm in enumerate(fwhmlist):
        yvals = np.zeros_like(xvals)
        if fwhm[0] <= 0:
            continue
        # Pull the sigma and start/end indexes from the FWHM list
        sigma = fwhm[5]
        mean = fwhm[6]
        amp = fwhm[7]
        startindex = int(fwhm[3])
        endindex = int(fwhm[4])

        xvalues = xvals[startindex:endindex + 1]
        yvalues = ndis_std(xvalues, float(mean), float(sigma), float(amp))

        yvals[startindex:endindex + 1] += yvalues
        simdata.append(yvals)

    return simdata

def intensity_decon(data, dlist, n=10):
    # dmax = np.amax(data[:, 1])
    simsum = np.sum(dlist, axis=0)
    b1 = data[:, 1] < 0
    data[b1, 1] = 0
    for i in range(n):
        ratios = ud.safedivide(data[:,1], simsum)
        for j in range(len(dlist)):
            dsum = np.sum(dlist[j])
            if dsum <= 0:
                continue
            avgratio = np.sum(ratios * dlist[j]) / dsum
            dlist[j] = dlist[j] * avgratio
        simsum = np.sum(dlist, axis=0)

    return simsum

def interpolate_minima(data, nslices=10, weights=None):
    ldata = len(data)
    slicepoints = np.linspace(0, ldata, nslices+1)
    minpoints = []
    for i in range(nslices):
        s = slicepoints[i]
        e = slicepoints[i+1]
        slicedata = data[int(s):int(e)]
        if weights is None:
            mincent = np.argmin(slicedata[:,1])
            minval = slicedata[mincent]
        else:
            minyval = np.average(slicedata[:,1], weights=weights[int(s):int(e)])
            minxval = np.average(slicedata[:,0], weights=weights[int(s):int(e)])
            minval = np.array([minxval, minyval])
        minpoints.append(minval)

    minpoints = np.array(minpoints)
    # Now interpolate between minima
    newy = np.interp(data[:,0], minpoints[:,0], minpoints[:,1])
    return newy

def interpolate_fwhms(data, fwhmlist):
    x = fwhmlist[:,6]
    y = fwhmlist[:,0]
    newy = np.interp(data[:,0], x, y)
    return newy

def variable_pw(data, window=5):
    # Find FWHM at highest intensity peak
    peakindex = np.argmax(data[:, 1])
    peakfwhm = find_fwhm(data, peakindex, maxppmtol=None)
    indexwidthintensity = peakfwhm[4] - peakfwhm[3]
    print("Index Width at Intensity Peak:", indexwidthintensity)

    # Use these max values for the background subtraction
    data = deepcopy(data)
    data = ud.datacompsub(data, indexwidthintensity)
    # Remove data below 0
    data[:, 1][data[:, 1] < 0] = 0

    peaks = ud.peakdetect(data, window=window)
    fwhms = fast_fwhm(data, peaks, maxppmtol=None)

    # Remove peaks with invalid FWHM
    fwhmsvalid = fwhms[:,0] > 0
    peaks = peaks[fwhmsvalid]
    fwhms = fwhms[fwhmsvalid]

    # Remove peaks with very small FWHM double the min line
    minyline = interpolate_minima(np.transpose([fwhms[:,6], fwhms[:,7]]), nslices=10, weights=None)
    b1 = fwhms[:,7] > 1.25*minyline
    peaks = peaks[b1]
    fwhms = fwhms[b1]

    # # Deconvolve peaks
    # peakdata = create_simdata(data, fwhmlist=fwhms)
    # idecon = intensity_decon(data, peakdata, n=10)

    flist = interpolate_fwhms(data, fwhms)

    return peaks, fwhms, data, flist


if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt
    from unidec.IsoDec.plots import cplot
    from unidec.engine import UniDec
    from unidec.modules.unidecwrapper import run_unidec_core
    import time

    os.chdir("C:\\Data\\UniDecTest\\")

    data = np.loadtxt("ADH.txt")
    # data = np.loadtxt("BSA.txt")
    data[:,1] = data[:,1] - np.amin(data[:,1])
    data[:,1] = data[:,1] / np.amax(data[:,1])

    peaks, fwhms, procdata, flist = variable_pw(data)

    # Max fwhm
    print("Max FWHM:", np.amax(flist))

    print(len(peaks), "Peaks Found")
    print(len(fwhms), "FWHMs Calculated")

    eng = UniDec()
    eng.pass_data_in(procdata)
    # flist = flist * 0 + 10
    eng.data.fwhmlist = flist
    eng.config.psig = 1
    eng.config.beta = 0
    eng.config.mzsig = 16
    eng.config.psfun = 2
    eng.config.minratio = 0


    eng.config.variablepw = 0
    ts1 = time.perf_counter()
    run_unidec_core(eng)

    te1 = time.perf_counter()
    m1 = deepcopy(eng.data.massdat)

    eng.config.variablepw = 1
    # eng.config.minratio = 0.1
    ts2 = time.perf_counter()
    run_unidec_core(eng)
    print("Time1:", te1 - ts1)
    print("Time2:", time.perf_counter()-ts2)
    m2 = deepcopy(eng.data.massdat)

    plt.subplot(121)
    plt.plot(procdata[:, 0], procdata[:, 1], 'k')
    plt.subplot(122)
    plt.plot(m1[:,0], m1[:,1], 'k')
    plt.plot(m2[:,0], m2[:,1], 'g')
    plt.show()

    exit()

    plt.subplot(121)
    # plt.plot(data[:, 0], data[:, 1])
    plt.plot(procdata[:, 0], procdata[:, 1], 'b')
    # plt.plot(procdata[:, 0], idecon, 'r')
    # Plot the max FWHM as vertical lines from 0 to 1
    # argmax = np.argmax(fwhms[:, 0])
    # maxfwhm = fwhms[argmax, 0]
    # peakmz = peaks[argmax, 0]
    # plt.vlines(peakmz - maxfwhm / 2, 0, 1, colors='r', linestyles='dashed')
    # plt.vlines(peakmz + maxfwhm / 2, 0, 1, colors='r', linestyles='dashed')
    # plt.vlines(peakmz, 0, 1, colors='g', linestyles='dashed')
    #
    # # Plot the max sigma as vertical lines from 0 to 1
    # argsigma = np.argmax(fwhms[:, 5])
    # maxsigma = fwhms[argsigma, 5]
    # peakmz2 = peaks[argsigma, 0]
    # plt.vlines(peakmz2 - maxsigma * 2.35482 / 2, 0, 1, colors='b', linestyles='dashed')
    # plt.vlines(peakmz2 + maxsigma * 2.35482 / 2, 0, 1, colors='b', linestyles='dashed')
    # plt.vlines(peakmz2, 0, 1, colors='g', linestyles='dashed')
    #
    #
    cplot(peaks, color='g')
    plt.subplot(122)
    # plt.plot(fwhms[:, 6], fwhms[:, 5], 'r.')
    plt.plot(fwhms[:, 6], fwhms[:,4]-fwhms[:,3], 'b.')

    newy = interpolate_minima(np.transpose([fwhms[:,6], fwhms[:,4]-fwhms[:,3]]), nslices=10, weights=None)
    plt.plot(fwhms[:,6], newy, 'g-')

    plt.plot(data[:,0], flist)

    plt.show()




    # os.chdir("C:\\Data\\Yuri\\")
    #
    # data = np.loadtxt("processed.txt")
    #
    # isodist = np.loadtxt("isodist.txt")
    # print(len(data), len(isodist))
    # st = time.perf_counter()
    # leftoverdata = remove_peaks(data, isodist)
    # gdat = create_peaks(data, isodist)
    #
    # print("Removed Peaks, Time:", time.perf_counter() - st)
    # plt.plot(data[:, 0], data[:, 1], label='Processed Data')
    # plt.plot(leftoverdata[:, 0], leftoverdata[:, 1], label='Leftover Data', linestyle='--')
    # plt.plot(gdat[:, 0], gdat[:, 1], label='Generated Peaks', linestyle=':')
    # plt.xlabel('m/z')
    # plt.ylabel('Intensity')
    # plt.show()
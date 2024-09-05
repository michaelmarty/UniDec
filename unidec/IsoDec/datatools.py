import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.tools as ud

import time
import matchms
from copy import deepcopy
from numba import njit
import numba as nb

# from bisect import bisect_left

'''
@njit(fastmath=True)
def bisect_left(a, x):
    """Similar to bisect.bisect_left(), from the built-in library."""
    M = len(a)
    for i in range(M):
        if a[i] >= x:
            return i
    return M'''


@njit(fastmath=True)
def bisect_left(a, x):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(i, x) will
    insert just before the leftmost x already there.
    """

    lo = 0
    hi = len(a)

    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid

    return lo


@njit(fastmath=True)
def fastnearest(array, target):
    """
    In a sorted array, quickly find the position of the element closest to the target.
    :param array: Array
    :param target: Value
    :return: np.argmin(np.abs(array - target))
    """
    i = int(bisect_left(array, target))
    if i <= 0:
        return 0
    elif i >= len(array) - 1:
        return len(array) - 1
    if np.abs(array[i] - target) > np.abs(array[i + 1] - target):
        i += 1
    elif np.abs(array[i] - target) > np.abs(array[i - 1] - target):
        i -= 1
    return int(i)

@njit(fastmath=True)
def calculate_cosinesimilarity(dist1, dist2):
    ab = 0
    a2 = 0
    b2 = 0
    for i in range(len(dist1)):
        ab += dist1[i] * dist2[i]
        a2 += dist1[i] ** 2
        b2 += dist2[i] ** 2
    return ab / (a2 ** 0.5 * b2 ** 0.5)

@njit(fastmath=True)
def calculate_cosinesimilarity2(isodist, isomatches, centroids, centmatches):
    ab = 0
    a2 = 0
    b2 = 0

    match_index = 0

    for i in range(len(isodist)):
        if i in isomatches:
            ab += isodist[i][1] * centroids[centmatches[match_index]][1]
            a2 += isodist[i][1] ** 2
            b2 += centroids[centmatches[match_index]][1] ** 2
            match_index += 1
        else:
            ab += isodist[i][1] * 0
            a2 += isodist[i][1] ** 2
            b2 += 0
    if ab == 0 or a2 == 0 or b2 == 0:
        return 0
    return ab / (a2 ** 0.5 * b2 ** 0.5)

@njit(fastmath=True)
def fastwithinppmtol(array, target, ppmtol):
    result = []
    nearest_idx = fastnearest(array, target)

    if ud.within_ppm(array[nearest_idx], target, ppmtol):
        result.append(nearest_idx)
    else:
        return result

    adding_upper = True
    adding_lower = True

    current_upper = nearest_idx + 1
    current_lower = nearest_idx - 1

    while adding_upper or adding_lower:
        if adding_upper:
            if current_upper >= len(array):
                adding_upper = False
            elif ud.within_ppm(array[current_upper], target, ppmtol):
                result.append(current_upper)
                current_upper += 1
            else:
                adding_upper = False

        if adding_lower:
            if current_lower < 0:
                adding_lower = False
            elif ud.within_ppm(array[current_lower], target, ppmtol):
                result.append(current_lower)
                current_lower -= 1
            else:
                adding_lower = False

    return result

@njit(fastmath=True)
def fastwithin_abstol(array,target, tol):
    result = []
    nearest_idx = fastnearest(array, target)

    if np.abs(array[nearest_idx] - target) <= tol:
        result.append(nearest_idx)
    else:
        return result

    adding_upper = True
    adding_lower = True

    current_upper = nearest_idx + 1
    current_lower = nearest_idx - 1

    while adding_upper or adding_lower:
        if adding_upper:
            if current_upper >= len(array):
                adding_upper = False
            elif np.abs(array[current_upper] - target) <= tol:
                result.append(current_upper)
                current_upper += 1
            else:
                adding_upper = False

        if adding_lower:
            if current_lower < 0:
                adding_lower = False
            elif np.abs(array[current_lower] - target) <= tol:
                result.append(current_lower)
                current_lower -= 1
            else:
                adding_lower = False

    return result

@njit(fastmath=True)
def fastpeakdetect(data, window:int=10, threshold=0.0, ppm=None, norm=True):
    """
    Simple peak detection algorithm.

    Detects a peak if a given data point is a local maximum within plus or minus config.peakwindow.
    Peaks must also be above a threshold of config.peakthresh * max_data_intensity.

    The mass and intensity of peaks meeting these criteria are output as a P x 2 array.

    :param data: Mass data array (N x 2) (mass intensity)
    :param window: Tolerance window of the x values
    :param threshold: Threshold of the y values
    :param ppm: Tolerance window in ppm
    :param norm: Whether to normalize the data before peak detection
    :return: Array of peaks positions and intensities (P x 2) (mass intensity)
    """
    peaks = []
    length = len(data)
    shape = np.shape(data)
    if length == 0 or shape[1] != 2:
        return np.array(peaks)

    if norm:
        maxval = np.amax(data[:, 1])
    else:
        maxval = 1
    for i in range(0, length):
        if data[i, 1] > maxval * threshold:
            if ppm is not None:
                ptmass = data[i, 0]
                newwin = ppm * 1e-6 * ptmass
                start = fastnearest(data[:, 0], ptmass - newwin)
                end = fastnearest(data[:, 0], ptmass + newwin)
            else:
                start = i - window
                end = i + window

                start = int(start)
                end = int(end) + 1

                if start < 0:
                    start = 0
                if end > length:
                    end = length

            testmax = np.amax(data[start:end, 1])
            if data[i, 1] == testmax and np.all(data[i, 1] != data[start:i, 1]):
                peaks.append([data[i, 0], data[i, 1]])

    return np.array(peaks)


@njit(fastmath=True)
def fastcalc_FWHM(peak, data):
    index = fastnearest(data[:, 0], peak)
    int = data[index, 1]
    leftwidth = 0
    rightwidth = 0
    counter = 1
    leftfound = False
    rightfound = False
    while rightfound is False or leftfound is False:
        if leftfound is False:
            if index - counter < 0:
                leftfound = True
            elif data[index - counter, 1] <= int / 2.:
                leftfound = True
                leftwidth += 1
            else:
                leftwidth += 1
        if rightfound is False:
            if index + counter >= len(data):
                rightfound = True
            elif data[index + counter, 1] <= int / 2.:
                rightfound = True
                rightwidth += 1
            else:
                rightwidth += 1
        counter += 1

    indexstart = index - leftwidth
    indexend = index + rightwidth
    if indexstart < 0:
        indexstart = 0
    if indexend >= len(data):
        indexend = len(data) - 1

    FWHM = data[indexend, 0] - data[indexstart, 0]
    return FWHM, [data[indexstart, 0], data[indexend, 0]]


def get_noise(data, n=20):
    """
    Get the noise level of the data.
    :param data: 2D numpy array of data
    :return: float, noise level
    """
    ydat = data[:, 1]
    # take n random data points from the data
    noise = np.concatenate([ydat[:n], ydat[-n:]])
    # return the standard deviation of the noise
    return np.std(noise)


def remove_noise_cdata(data, localmin=100, factor=1.5, mode="median"):
    """
    Remove noise from the data.
    :param data: 2D numpy array of data
    :param localmin: int, number of data points local width to take for min calcs
    :return: data with noise removed
    """
    if len(data) < localmin:
        return data

    ydat = data[:, 1]

    # find local minima
    if mode == "median":
        localmins = [np.median(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    elif mode == "min":
        localmins = [np.amin(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    elif mode == "mean":
        localmins = [np.mean(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    else:
        localmins = [np.median(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    # Extend the last bit at the end value to be the same length as data
    localmins = np.concatenate([localmins, np.full(len(data) - len(localmins), localmins[-1])])
    # smooth
    noiselevel = np.convolve(localmins, np.ones(localmin*1) / localmin*1, mode="same")
    # Subtract the noise level from the data
    newdata = data[:, 1] - noiselevel * factor

    data = data[newdata > 0]
    # return the standard deviation of the noise
    return data


def remove_noise_cdata2(data, min_snr=2):
    ydata = data[:, 1]
    noise = np.median(ydata)
    #plt.hist(ydata, bins=1000)
    #plt.show()
    max_signal = data[np.argmax(ydata), 1]
    print("Max SNR:", max_signal / noise)
    min_signal = noise * min_snr
    data = data[ydata > min_signal]
    return data


def get_top_peak_mz(data):
    """
    Get the m/z value of the top peak in the data.
    :param data: 2D numpy array of data
    :return: float, m/z value of the top peak
    """
    # get the index of the maximum value in the data
    maxindex = np.argmax(data[:, 1])
    # return the m/z value of the peak
    return data[maxindex, 0]


@njit(fastmath=True)
def get_fwhm_peak(data, peakmz):
    """
    Get the full width half max of the peak.
    :param data: 2D numpy array of data
    :param peakmz: float, m/z value of the peak
    :return: float, full width half max of the peak
    """
    fwhm, interval = fastcalc_FWHM(peakmz, data)
    return fwhm


@njit(fastmath=True)
def get_centroid(data, peakmz, fwhm=1):
    """
    Get the centroid of the peak.
    :param data: 2D numpy array of data
    :param peakmz: float, m/z value of the peak
    :return: float, centroid of the peak
    """
    b1 = data[:, 0] > peakmz - fwhm
    b2 = data[:, 0] < peakmz + fwhm
    b = b1 & b2
    if np.sum(b) == 0:
        return peakmz

    d = data[b]
    return np.sum(d[:, 0] * d[:, 1]) / np.sum(d[:, 1])


@njit(fastmath=True)
def get_all_centroids(data, window=5, threshold=0.0001):
    if len(data) < 3:
        return np.empty((0, 2))

    fwhm = get_fwhm_peak(data, data[np.argmax(data[:, 1]), 0])
    peaks = fastpeakdetect(data, window=window, threshold=threshold)

    p0 = [get_centroid(data, p[0], fwhm) for p in peaks]
    outpeaks = np.empty((len(p0), 2))
    outpeaks[:, 0] = p0
    outpeaks[:, 1] = peaks[:, 1]

    return outpeaks


# @njit(fastmath=True)
def get_centroids(data, peakmz, mzwindow=None):
    if mzwindow is None:
        mzwindow = [-1.5, 3.5]
    chopped = data[(data[:, 0] > peakmz + mzwindow[0]) & (data[:, 0] < peakmz + mzwindow[1])]

    if len(chopped) < 3:
        return np.array([]), np.array([])

    fwhm = get_fwhm_peak(chopped, peakmz)

    peaks = fastpeakdetect(chopped, window=5, threshold=0.01)
    # print(peaks)

    for p in peaks:
        d = chopped[(chopped[:, 0] > p[0] - fwhm) & (chopped[:, 0] < p[0] + fwhm)]
        if len(d) == 0:
            continue
        p[0] = np.sum(d[:, 0] * d[:, 1]) / np.sum(d[:, 1])
    # print(peaks)
    return peaks, chopped


def calc_match_simp(centroids, isodist):
    spectrum1 = matchms.Spectrum(mz=centroids[:, 0], intensities=centroids[:, 1], metadata={"precursor_mz": 1})
    spectrum2 = matchms.Spectrum(mz=isodist[:, 0], intensities=isodist[:, 1], metadata={"precursor_mz": 1})
    cosine_greedy = matchms.similarity.CosineGreedy(tolerance=0.01)
    score = cosine_greedy.pair(spectrum1, spectrum2)
    print(score)
    return score


def isotope_finder(data, mzwindow=1.5):
    nl = get_noise(data)
    # Chop data below noise level
    peaks = fastpeakdetect(data, window=50, threshold=nl / np.amax(data[:, 1]))
    # sort peaks
    peaks = np.array(sorted(peaks, key=lambda x: x[1], reverse=True))
    return peaks


def simp_charge(centroids, silent=False):
    """
    Simple charge prediction based on the spacing between peaks. Picks the largest peak and looks at the spacing
    between the two nearest peaks. Takes 1/avg spacing as the charge.
    :param centroids: Centroid data, 2D numpy array as [m/z, intensity]
    :param silent: Whether to print the charge state, default False
    :return: Charge state as int
    """
    diffs = np.diff(centroids[:, 0])
    maxindex = np.argmax(centroids[:, 1])
    start = maxindex - 1
    end = maxindex + 2
    if start < 0:
        start = 0
    if end > len(diffs):
        end = len(diffs)
    centraldiffs = diffs[start:end]
    if len(centraldiffs) == 0:
        return 0
    avg = np.mean(centraldiffs)
    if avg == 0:
        return 0
    charge = 1 / avg
    charge = round(charge)
    if not silent:
        print("Simple Prediction:", charge)
    return charge

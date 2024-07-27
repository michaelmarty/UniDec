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


@njit(fastmath=True)
def bisect_left(a, x):
    """Similar to bisect.bisect_left(), from the built-in library."""
    M = a.size
    for i in range(M):
        if a[i] >= x:
            return i
    return M


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
    if np.abs(array[i] - target) > np.abs(array[i - 1] - target):
        i -= 1
    return int(i)


@njit(fastmath=True)
def fastpeakdetect(data, config=None, window=10, threshold=0.0, ppm=None, norm=True):
    """
    Simple peak detection algorithm.

    Detects a peak if a given data point is a local maximum within plus or minus config.peakwindow.
    Peaks must also be above a threshold of config.peakthresh * max_data_intensity.

    The mass and intensity of peaks meeting these criteria are output as a P x 2 array.

    :param data: Mass data array (N x 2) (mass intensity)
    :param config: UniDecConfig object
    :param window: Tolerance window of the x values
    :param threshold: Threshold of the y values
    :param ppm: Tolerance window in ppm
    :param norm: Whether to normalize the data before peak detection
    :return: Array of peaks positions and intensities (P x 2) (mass intensity)
    """
    if config is not None:
        window = config.peakwindow / config.massbins
        threshold = config.peakthresh
        norm = config.normthresh

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
    charge = 1 / avg
    charge = round(charge)
    if not silent:
        print("Simple Prediction:", charge)
    return charge

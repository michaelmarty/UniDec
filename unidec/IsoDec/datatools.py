import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.tools as ud
from unidec.modules.isotopetools import *
import time
import matchms
from copy import deepcopy


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


def get_fwhm_peak(data, peakmz):
    """
    Get the full width half max of the peak.
    :param data: 2D numpy array of data
    :param peakmz: float, m/z value of the peak
    :return: float, full width half max of the peak
    """
    fwhm, interval = ud.calc_FWHM(peakmz, data)
    return fwhm


def get_centroids(data, peakmz, mzwindow=1.5):
    if type(mzwindow) is float:
        mzwindow = [-mzwindow, mzwindow]
    chopped = data[(data[:, 0] > peakmz + mzwindow[0]) & (data[:, 0] < peakmz + mzwindow[1])]

    fwhm = get_fwhm_peak(chopped, peakmz)

    peaks = ud.peakdetect(chopped, window=3, threshold=0.01)
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
    peaks = ud.peakdetect(data, window=50, threshold=nl / np.amax(data[:, 1]))
    # sort peaks
    peaks = np.array(sorted(peaks, key=lambda x: x[1], reverse=True))
    return peaks


def create_isodist(peakmz, charge, data):
    mass = peakmz * float(charge)
    isodist = calc_averagine_isotope_dist(mass, charge=charge, crop=True)
    isodist[:, 1] *= np.amax(data[:, 1]) / np.amax(isodist[:, 1])
    # shift isodist so that maxes are aligned with data
    isodist[:, 0] = isodist[:, 0] + peakmz - isodist[np.argmax(isodist[:, 1]), 0]
    return isodist


def simp_charge(centroids, silent=False):
    diffs = np.diff(centroids[:, 0])
    maxindex = np.argmax(centroids[:, 1])
    start = maxindex - 1
    end = maxindex + 2
    if start < 0:
        start = 0
    if end > len(diffs):
        end = len(diffs)
    centraldiffs = diffs[start:end]
    avg = np.mean(centraldiffs)
    charge = 1 / avg
    charge = round(charge)
    if not silent:
        print("Simple Prediction:", charge)
    return charge

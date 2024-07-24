import numpy as np
import time
import matchms
from copy import deepcopy
from unidec.IsoDec.datatools import simp_charge, create_isodist
import matplotlib.pyplot as plt
from numba import njit

mass_diff_c = 1.0033


def cplot(centroids, color='r', factor=1, base=0, mask=None, mfactor=-1, mcolor="g", z=0, zcolor="b", zfactor=1):
    """
    Simple script to plot centroids
    :param centroids: Centroid array with m/z in first column and intensity in second
    :param color: Color
    :param factor: Mutiplicative factor for intensity. -1 will set below the axis
    :param base: Base of the lines. Default is 0. Can be adjusted to shift up or down.
    :return: None
    """
    plt.hlines(0, np.amin(centroids[:, 0]), np.amax(centroids[:, 0]), color="k")

    if mask is not None:
        if len(centroids) > len(mask):
            mask = np.append(mask, np.zeros(len(centroids) - len(mask)))
        else:
            mask = mask[:len(centroids)]
        for c in centroids[mask.astype(bool)]:
            plt.vlines(c[0], base, base + mfactor * c[1], color=mcolor)

    if z != 0:
        isodist = create_isodist(centroids[np.argmax(centroids[:, 1]), 0], z, centroids)
        for c in isodist:
            plt.vlines(c[0], base, base + zfactor * c[1], color=zcolor, linewidth=3)

    for c in centroids:
        plt.vlines(c[0], base, base + factor * c[1], color=color)

@njit(fastmath=True)
def optimize_shift(centroids, isodist, z, tol=0.01, maxshift=3):
    shiftrange = np.arange(-maxshift, maxshift + 1)
    matchedindexes, isomatches = match_peaks(centroids, isodist)
    #mc = centroids[matchedindexes]
    #mi = isodist[isomatches]

    mc = [centroids[i, 1] for i in matchedindexes]
    mi = [isodist[i, 1] for i in isomatches]

    mc = np.array(mc)
    mi = np.array(mi)

    #overlaps = []
    bestshift = -1
    sum = -1
    for i, shift in enumerate(shiftrange):
        overlap = mc * np.roll(mi, shift)
        s = np.sum(overlap)
        #overlaps.append(s)

        if s > sum:
            sum = s
            bestshift = shift
    #bestshift = shiftrange[np.argmax(overlaps)]

    isodist[:, 0] = isodist[:, 0] + bestshift * mass_diff_c / z
    return isodist

@njit(fastmath=True)
def match_peaks(centroids, isodist, tol=0.01):
    matchingpeaks = matchms.similarity.spectrum_similarity_functions.find_matches(centroids[:, 0],
                                                                                  isodist[:, 0], tol)
    matchedindexes = [match[0] for match in matchingpeaks]
    isomatches = [match[1] for match in matchingpeaks]
    return matchedindexes, isomatches


def match_charge(centroids, peakmz, charge_range=[1, 50], tol=0.01, threshold=0.85, nextbest=0.8, baseline=0.1,
                 silent=True):
    startime = time.perf_counter()
    topcent = deepcopy(centroids)

    b1 = centroids[:, 1] > baseline * np.amax(centroids[:, 1])
    centroids = centroids[b1]
    if len(centroids) < 3:
        return 0, [], [], []

    target = matchms.Spectrum(mz=centroids[:, 0], intensities=centroids[:, 1], metadata={"precursor_mz": peakmz})
    cosine_greedy = matchms.similarity.CosineGreedy(tolerance=tol)

    ztab = np.arange(charge_range[0], charge_range[1])
    charge_guess = simp_charge(centroids, silent=True)
    # sort so that things close to charge_guess come first in ztab
    ztab = np.array(sorted(ztab, key=lambda x: np.abs(x - charge_guess)))

    for i, z in enumerate(ztab):
        isodist = create_isodist(peakmz, z, centroids)

        matchindexes, isomatches = match_peaks(topcent, isodist, tol)
        if len(matchindexes) < 3:
            continue

        spectrum2 = matchms.Spectrum(mz=isodist[:, 0], intensities=isodist[:, 1],
                                     metadata={"charge": z, "precursor_mz": peakmz})

        score = cosine_greedy.pair(target, spectrum2)['score']

        if score > threshold:

            if not silent:
                print("Matched:", z, score)
                endtime = time.perf_counter()
                print("Time:", endtime - startime)
            return z, isodist, matchindexes, isomatches
    if not silent:
        print("No Match")
        endtime = time.perf_counter()
        print("Time:", endtime - startime)
    return 0, [], [], []

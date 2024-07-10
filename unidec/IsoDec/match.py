import numpy as np
import time
import matchms
from copy import deepcopy
from IsoDec.datatools import simp_charge, create_isodist
import matplotlib.pyplot as plt

mass_diff_c = 1.0033


def cplot(centroids, color='r', factor=1):
    for c in centroids:
        plt.vlines(c[0], 0, factor * c[1], color=color)


def optimize_shift(centroids, isodist, z, tol=0.01, maxshift=3):
    shiftrange = np.arange(-maxshift, maxshift + 1)
    matchedindexes, isomatches = match_peaks(centroids, isodist)
    mc = centroids[matchedindexes]
    mi = isodist[isomatches]

    overlaps = []
    for i, shift in enumerate(shiftrange):
        overlap = mc[:, 1] * np.roll(mi[:, 1], shift)
        overlaps.append(np.sum(overlap))
    bestshift = shiftrange[np.argmax(overlaps)]

    isodist[:, 0] = isodist[:, 0] + bestshift * mass_diff_c/z
    return isodist


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
        spectrum2 = matchms.Spectrum(mz=isodist[:, 0], intensities=isodist[:, 1],
                                     metadata={"charge": z, "precursor_mz": peakmz})

        score = cosine_greedy.pair(target, spectrum2)['score']

        if score > threshold:
            matchindexes, isomatches = match_peaks(topcent, isodist, tol)

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

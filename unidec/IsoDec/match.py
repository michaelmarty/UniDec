import numpy as np
import time
import matchms
from copy import deepcopy
from unidec.IsoDec.datatools import simp_charge
import matplotlib.pyplot as plt
from numba import njit
from typing import List, Tuple
import matplotlib as mpl
from unidec.modules.isotopetools import *
import pickle as pkl


# mass_diff_c = 1.0033
class MatchedCollection:
    """
    Class for collecting matched peaks
    """

    def __init__(self):
        self.peaks = []
        self.colormap = mpl.colormaps.get_cmap("tab10")

    def add_peak(self, peak):
        """
        Add peak to collection
        :param peak: Add peak to collection
        :return:
        """
        newcolor = self.colormap(len(self.peaks) % 10)
        peak.color = newcolor
        self.peaks.append(peak)

    def save_pks(self, filename="peaks.pkl"):
        with open(filename, "wb") as f:
            pkl.dump(self.peaks, f)
            print(f"Saved {len(self.peaks)} peaks to {filename}")

    def load_pks(self, filename="peaks.pkl"):
        with open(filename, "rb") as f:
            self.peaks = pkl.load(f)
            print(f"Loaded {len(self.peaks)} peaks from {filename}")


class MatchedPeak:
    """
    Matched peak object for collecting data on peaks with matched distributions
    """

    def __init__(self, centroids, isodist, z, mz, matchedindexes=None, isomatches=None):
        self.mz = mz
        self.z = z
        self.centroids = centroids
        self.isodist = isodist
        self.matchedcentroids = None
        self.matchedisodist = None
        self.matchedindexes = None
        self.isomatches = None
        self.mask = None
        self.color = "g"
        self.scan = -1
        self.massdist = None
        self.monoiso = -1
        if matchedindexes is not None:
            self.matchedindexes = matchedindexes
            self.matchedcentroids = centroids[matchedindexes]
            self.mask = np.zeros(len(centroids))
            self.mask[matchedindexes] = 1
        if isomatches is not None:
            self.isomatches = isomatches
            self.matchedisodist = isodist[isomatches]


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
def create_isodist(peakmz, charge, data, adductmass=1.007276467):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    mass = (peakmz - adductmass) * charge
    isodist = fast_calc_averagine_isotope_dist(mass, charge=charge)
    isodist[:, 1] *= np.amax(data[:, 1])
    # shift isodist so that maxes are aligned with data
    mzshift = peakmz - isodist[np.argmax(isodist[:, 1]), 0]
    isodist[:, 0] = isodist[:, 0] + mzshift
    return isodist


@njit(fastmath=True)
def create_isodist_full(peakmz, charge, data, adductmass=1.007276467):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    mass = (peakmz - adductmass) * charge
    isodist, massdist = fast_calc_averagine_isotope_dist_dualoutput(mass, charge=charge)
    isodist[:, 1] *= np.amax(data[:, 1])
    massdist[:, 1] *= np.amax(data[:, 1])
    # shift isodist so that maxes are aligned with data
    mzshift = peakmz - isodist[np.argmax(isodist[:, 1]), 0]
    isodist[:, 0] = isodist[:, 0] + mzshift

    massshift = mzshift * charge
    monoiso = mass + massshift
    massdist[:, 0] = massdist[:, 0] + massshift
    return isodist, massdist, monoiso


@njit(fastmath=True)
def optimize_shift(centroids, z, tol=0.01, maxshift=2):
    # Limit max shifts if necessary
    if z < 3:
        maxshift = 1
    elif z < 10:
        maxshift = 2
    else:
        maxshift = maxshift

    peakmz = centroids[np.argmax(centroids[:, 1]), 0]
    isodist, massdist, monoiso = create_isodist_full(peakmz, z, centroids)

    matchedindexes, isomatches = match_peaks(centroids, isodist)
    # mc = centroids[matchedindexes]
    # mi = isodist[isomatches]

    mc = [centroids[i, 1] for i in matchedindexes]
    mi = [isodist[i, 1] for i in isomatches]

    mc = np.array(mc)
    mi = np.array(mi)

    if maxshift > len(mc):
        maxshift = len(mc) - 1

    shiftrange = np.arange(-maxshift, maxshift + 1)

    # overlaps = []
    bestshift = -1
    sum = -1
    for i, shift in enumerate(shiftrange):
        overlap = mc * np.roll(mi, shift)
        s = np.sum(overlap)
        # overlaps.append(s)
        if s > sum:
            sum = s
            bestshift = shift
    # bestshift = shiftrange[np.argmax(overlaps)]

    # Correct masses based on shift
    shiftmass = bestshift * mass_diff_c
    monoiso = monoiso + shiftmass
    massdist[:, 0] = massdist[:, 0] + shiftmass
    # Correct m/z based on shift
    shiftmz = shiftmass / z
    isodist[:, 0] = isodist[:, 0] + shiftmz
    peakmz = peakmz + shiftmz
    # Match it again
    matchedindexes, isomatches = match_peaks(centroids, isodist, tol=tol)

    return isodist, matchedindexes, isomatches, peakmz, monoiso, massdist


@njit(fastmath=True)
def find_matches(spec1_mz: np.ndarray, spec2_mz: np.ndarray,
                 tolerance: float, shift: float = 0) -> Tuple[List[int], List[int]]:
    """Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from
    low to high m/z).

    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_mz:
        Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    tolerance
        Peaks will be considered a match when <= tolerance appart.
    shift
        Shift peaks of second spectra by shift. The default is 0.

    Returns
    -------
    matches
        List containing entries of type (idx1, idx2).

    """
    lowest_idx = 0
    m1 = []
    m2 = []
    for peak1_idx in range(spec1_mz.shape[0]):
        mz = spec1_mz[peak1_idx]
        low_bound = mz - tolerance
        high_bound = mz + tolerance
        for peak2_idx in range(lowest_idx, spec2_mz.shape[0]):
            mz2 = spec2_mz[peak2_idx] + shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx + 1
            else:
                m1.append(peak1_idx)
                m2.append(peak2_idx)
    return m1, m2


@njit(fastmath=True)
def match_peaks(centroids: np.array, isodist: np.array, tol: float = 0.01) -> Tuple[List[int], List[int]]:
    """
    matchingpeaks = matchms.similarity.spectrum_similarity_functions.find_matches(centroids[:, 0],
                                                                                  isodist[:, 0], tol)
    matchedindexes = [match[0] for match in matchingpeaks]
    isomatches = [match[1] for match in matchingpeaks]"""
    matchedindexes, isomatches = find_matches(centroids[:, 0], isodist[:, 0], tol)
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


if __name__ == "__main__":
    isodist = [[860.51031724, 66796.046875],
               [861.51031724, 31534.41372969],
               [862.51031724, 8754.51400417],
               [863.51031724, 1790.77682158]]

    centroids = [[859.29762565, 1947.690552],
                 [859.89066066, 3896.894775],
                 [860.03604576, 2534.360352],
                 [860.18160474, 3057.42041],
                 [860.3255286, 2427.568115],
                 [860.51031724, 66796.046875],
                 [861.10237412, 1194.873047],
                 [861.31358998, 950.150452],
                 [861.51228606, 30179.460938],
                 [861.85051775, 1136.555054],
                 [862.04103128, 2694.974609],
                 [862.0847034, 1198.278931],
                 [862.51676745, 9102.647461],
                 [862.74771343, 2822.634277],
                 [863.03157199, 2646.25415],
                 [863.43871169, 1349.041504],
                 [863.75527503, 5128.32373],
                 [863.97688287, 1786.120972]]

    optimize_shift(np.array(centroids), np.array(isodist), 1)

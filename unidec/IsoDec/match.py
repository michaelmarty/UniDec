import numpy as np
import time
import matchms
from copy import deepcopy
from unidec.IsoDec.datatools import simp_charge, fastnearest, fastwithinppmtol, calculate_cosinesimilarity
import matplotlib.pyplot as plt
from numba import njit
from typing import List, Tuple
import numba as nb
import matplotlib as mpl
from unidec.modules.isotopetools import *
import pickle as pkl
import unidec.tools as ud
import unidec.IsoDec.msalign_export as msalign
import math
import pandas as pd


# mass_diff_c = 1.0033
class MatchedCollection:
    """
    Class for collecting matched peaks
    """

    def __init__(self):
        self.peaks = []
        self.masses = []
        self.monoisos = np.array([])
        self.colormap = mpl.colormaps.get_cmap("tab10")

    def __str__(self):
        outstring = ""
        for p in self.peaks:
            outstring += str(p)
        return outstring

    def __len__(self):
        return len(self.peaks)

    def __getitem__(self, item):
        return self.peaks[item]

    def add_peak(self, peak):
        """
        Add peak to collection
        :param peak: Add peak to collection
        :return:
        """
        newcolor = self.colormap(len(self.peaks) % 10)
        peak.color = newcolor
        self.peaks.append(peak)

    def add_peaks(self, peaks):
        self.peaks = peaks
        return self

    def save_pks(self, filename="peaks.pkl"):
        with open(filename, "wb") as f:
            pkl.dump(self.peaks, f)
            print(f"Saved {len(self.peaks)} peaks to {filename}")

    def load_pks(self, filename="peaks.pkl"):
        with open(filename, "rb") as f:
            self.peaks = pkl.load(f)
            print(f"Loaded {len(self.peaks)} peaks from {filename}")
        return self

    def add_pk_to_masses(self, pk, ppmtol):
        """
        Checks if an existing mass matches to this peak, if so adds it to that mass, otherwise creates a new mass
        The list of masses is constantly kept in order of monoisotopic mass.
        """
        if len(self.masses) == 0:
            self.monoisos = np.append(self.monoisos, pk.monoiso)
            self.masses.append(MatchedMass(pk))


        else:
            idx = fastnearest(self.monoisos, pk.monoiso)
            nearest_mass = self.monoisos[idx]
            if pk.scan - self.masses[idx].scans[len(self.masses[idx].scans) - 1] <= 100 and ud.within_ppm(
                    self.masses[idx].monoiso, pk.monoiso, ppmtol):
                self.masses[idx].scans = np.append(self.masses[idx].scans, pk.scan)
                if pk.matchedintensity is not None:
                    if pk.matchedintensity > self.masses[idx].maxintensity:
                        self.masses[idx].maxintensity = pk.matchedintensity
                        self.masses[idx].maxscan = pk.scan
                        self.masses[idx].maxrt = pk.rt

                if not np.isin(self.masses[idx].zs, pk.z).any():
                    self.masses[idx].zs = np.append(self.masses[idx].zs, pk.z)
                    self.masses[idx].mzs = np.append(self.masses[idx].mzs, pk.mz)

            else:
                if pk.monoiso > nearest_mass:
                    idx += 1
                if idx == len(self.monoisos):
                    self.monoisos = np.append(self.monoisos, pk.monoiso)
                    self.masses.append(MatchedMass(pk))
                else:
                    self.monoisos = np.insert(self.monoisos, idx, pk.monoiso)
                    self.masses.insert(idx, MatchedMass(pk))

    def export_prosightlite(self, filename="prosight.txt"):
        with open(filename, "w") as f:
            for p in self.masses:
                f.write(str(p.monoiso) + "\n")

    def merge_missed_monoisotopics(self, ppm_tolerance=20, max_mm=1):
        mass_diff_c = 1.0033
        to_remove = []
        for i1 in range(len(self.peaks)):
            for i2 in range(len(self.peaks)):
                if i1 == i2 or self.peaks[i1].z != self.peaks[i2].z:
                    continue
                if abs(self.peaks[i1].monoiso - self.peaks[i2].monoiso) < max_mm * 1.0033 * 1.1:
                    mm_count = int(round(self.peaks[i1].monoiso - self.peaks[i2].monoiso))
                    if self.peaks[i1].matchedintensity > self.peaks[i2].matchedintensity:
                        adj_mass = self.peaks[i2].monoiso + mm_count * mass_diff_c
                        if ud.within_ppm(self.peaks[i1].monoiso, adj_mass, ppm_tolerance):
                            to_remove.append(i2)
                    else:
                        adj_mass = self.peaks[i1].monoiso - mm_count * mass_diff_c
                        if ud.within_ppm(self.peaks[i2].monoiso, adj_mass, ppm_tolerance):
                            to_remove.append(i1)
        self.peaks = [self.peaks[i] for i in range(len(self.peaks)) if i not in to_remove]
        return

    def copy_to_string(self):
        outstring = ""
        for p in self.peaks:
            outstring += str(p.monoiso) + "\n"
        return outstring

    def export_msalign(self, reader, filename="export.msalign", act_type="HCD", max_precursors=None):
        print("Exporting to", filename, "with activation type", act_type, "N:", len(self.peaks))
        ms1_scan_dict = msalign.sort_by_scan_order(self, 1)
        ms2_scan_dict = msalign.sort_by_scan_order(self, 2)

        msalign.write_ms1_msalign(ms1_scan_dict, filename)
        msalign.write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, filename, max_precursors=max_precursors,
                                  act_type=act_type)

    def to_df(self):
        """
        Convert a MatchedCollection object to a pandas dataframe
        :param pks: MatchedCollection object
        :return: Pandas dataframe
        """
        data = []
        for p in self.peaks:
            d = {"Charge": p.z, "Most Abundant m/z": p.mz, "Monoisotopic Mass": p.monoiso, "Scan": p.scan,
                 "Most Abundant Mass": p.peakmass, "Abundance": p.matchedintensity}
            data.append(d)
        df = pd.DataFrame(data)
        return df

    def export_tsv(self, filename="export.tsv"):
        df = self.to_df()
        df.to_csv(filename, sep="\t", index=False)


class MatchedMass:
    """
    Matched mass object for collecting data on MatchedPeaks with matched masses.
    """

    def __init__(self, pk):
        self.monoiso = pk.monoiso
        self.scans = np.array([pk.scan])
        self.maxintensity = pk.matchedintensity
        self.maxscan = pk.scan
        self.maxrt = pk.rt
        self.mzs = np.array([pk.mz])
        self.zs = np.array([pk.z])


'''@nb.experimental.jitclass([("centroids", nb.types.Array(nb.float64, 2, "C")),
                           ("isodist", nb.types.Array(nb.float64, 2, "C")),
                           ("matchedcentroids", nb.types.Array(nb.float64, 2, "C")),
                           ("matchedisodist", nb.types.Array(nb.float64, 2, "C")),
                           ("matchedindexes", nb.types.Array(nb.int64, 1, "C")),
                           ("isomatches",  nb.types.Array(nb.int64, 1, "C")),
                           ("massdist", nb.types.Array(nb.float64, 2, "C")),
                           ])'''


class MatchedPeak:
    """
    Matched peak object for collecting data on peaks with matched distributions
    """
    mz: float
    z: int
    centroids: nb.optional(np.ndarray)
    isodist: nb.optional(np.ndarray)
    matchedintensity: nb.optional(float)
    matchedcentroids: nb.optional(np.ndarray)
    matchedisodist: nb.optional(np.ndarray)
    matchedindexes: nb.optional(np.ndarray)
    isomatches: nb.optional(np.ndarray)
    color: str
    scan: int
    rt: float
    ms_order: int
    massdist: nb.optional(np.ndarray)
    monoiso: float
    peakmass: float
    avgmass: float
    startindex: int
    endindex: int
    matchedion: str

    def __init__(self, z, mz, centroids=None, isodist=None, matchedindexes=None, isomatches=None):
        self.mz = mz
        self.z = z
        self.centroids = centroids
        self.isodist = isodist
        self.matchedintensity = None
        self.matchedcentroids = None
        self.matchedisodist = None
        self.matchedindexes = None
        self.isomatches = None
        self.color = "g"
        self.scan = -1
        self.rt = -1
        self.ms_order = -1
        self.massdist = None
        self.monoiso = -1
        self.peakmass = -1
        self.avgmass = -1
        self.startindex = -1
        self.endindex = -1
        self.matchedion = ""
        if matchedindexes is not None:
            self.matchedindexes = matchedindexes
            self.matchedcentroids = centroids[np.array(matchedindexes)]
            yvals = np.array([c[1] for c in centroids[np.array(matchedindexes)]])
            self.matchedintensity = np.sum(yvals)
        if isomatches is not None:
            self.isomatches = isomatches
            self.matchedisodist = isodist[np.array(isomatches)]

    def __str__(self):
        return f"MatchedPeak: mz={self.mz}, z={self.z}, monoiso={self.monoiso}\n"


def df_to_matchedcollection(df, monoiso="Monoisotopic Mass", peakmz="Most Abundant m/z", peakmass="Most Abundant Mass",
                            scan="Scan", z="Charge", intensity="Abundance", ion="Ion"):
    """
    Convert a pandas dataframe to a MatchedCollection object
    :param df: Pandas dataframe with columns mz, intensity, z, scan, rt, monoiso, peakmass, avgmass
    :return: MatchedCollection object
    """
    mc = MatchedCollection()
    for i, row in df.iterrows():
        pk = MatchedPeak(row[z], row[peakmz])
        pk.scan = row[scan]
        pk.monoiso = row[monoiso]
        pk.peakmass = row[peakmass]
        pk.matchedion = row[ion]
        if pk.matchedion == math.nan:
            pk.matchedion = None
        isodist = create_isodist2(pk.monoiso, pk.z, row[intensity])
        pk.isodist = isodist
        mc.add_peak(pk)
    print("Loaded", len(mc.peaks), "peaks from dataframe")
    return mc


def peak_mz_z_df_to_matchedcollection(df, data=None):
    """
    Convert a pandas dataframe of peak mzs (col1) and zs (col2) to a MatchedCollection object
    :param df: Pandas dataframe with columns mz, z
    :param data: Optional data to match to as 2D numpy array [m/z, intensity]
    :return: MatchedCollection object
    """
    mc = MatchedCollection()
    for i, row in df.iterrows():
        pk = MatchedPeak(z=int(row['z']), mz=row['peakmz'])
        pk.monoiso = (pk.mz - 1.007276467) * pk.z
        mc.add_peak(pk)
        isodist = fast_calc_averagine_isotope_dist(pk.monoiso, pk.z)

        if data is not None:
            # Find nearest peak in data
            mz = isodist[np.argmax(isodist[:, 1]), 0]
            startindex = fastnearest(data[:, 0], mz)
            # endindex = fastnearest(data[:, 0], mz + window[1])

            dmax = data[startindex, 1]


        else:
            dmax = 10000

        for i in range(len(isodist)):
            isodist[i, 1] *= dmax
        pk.mz = isodist[np.argmax(isodist[:, 1]), 0]
        pk.isodist = isodist
    print("Loaded", len(mc.peaks), "peaks from dataframe")
    return mc


def compare_matchedcollections(coll1, coll2, ppmtol=50, objecttocompare="monoiso"):
    masses1 = []
    masses2 = []
    if objecttocompare == "monoiso":
        coll1.peaks = sorted(coll1.peaks, key=lambda x: x.monoiso)
        coll2.peaks = sorted(coll2.peaks, key=lambda x: x.monoiso)
        masses1 = np.array([p.monoiso for p in coll1.peaks])
        masses2 = np.array([p.monoiso for p in coll2.peaks])
    elif objecttocompare == "peakmass":
        coll1.peaks = sorted(coll1.peaks, key=lambda x: x.peakmass)
        coll2.peaks = sorted(coll2.peaks, key=lambda x: x.peakmass)
        masses2 = np.array([p.peakmass for p in coll2.peaks])
        masses1 = np.array([p.peakmass for p in coll1.peaks])
    elif objecttocompare == "peakmz":
        coll1.peaks = sorted(coll1.peaks, key=lambda x: x.mz)
        coll2.peaks = sorted(coll2.peaks, key=lambda x: x.mz)
        masses2 = np.array([p.mz for p in coll2.peaks])
        masses1 = np.array([p.mz for p in coll1.peaks])

    unique1 = []
    shared = []
    coll2_matched_indices = []

    for i in range(len(coll1.peaks)):
        foundmatch = False

        within_tol = fastwithinppmtol(masses2, masses1[i], ppmtol)
        if len(within_tol) > 0:
            for j in within_tol:
                if coll1.peaks[i].z == coll2.peaks[j].z:
                    coll2_matched_indices.append(j)
                    foundmatch = True

        if foundmatch:
            shared.append(coll1.peaks[i])
        else:
            unique1.append(coll1.peaks[i])

    unique2 = [coll2.peaks[i] for i in range(len(coll2.peaks)) if i not in coll2_matched_indices]

    shared = MatchedCollection().add_peaks(shared)
    unique1 = MatchedCollection().add_peaks(unique1)
    unique2 = MatchedCollection().add_peaks(unique2)

    return shared, unique1, unique2


def get_unique_matchedions(coll1, coll2):
    shared = []
    unique1 = []
    unique2 = []

    coll2_matchedindices = []
    for p1 in range(len(coll1.peaks)):
        if not pd.isnull(coll1.peaks[p1].matchedion):
            matched = False
            for p2 in range(len(coll2.peaks)):
                if (not pd.isnull(coll2.peaks[p2].matchedion) and
                        coll1.peaks[p1].matchedion == coll2.peaks[p2].matchedion and
                        coll1.peaks[p1].z == coll2.peaks[p2].z):
                    coll2_matchedindices.append(p2)
                    matched = True

            if matched:
                shared.append(coll1.peaks[p1])
            else:
                unique1.append(coll1.peaks[p1])

    for p2 in range(len(coll2.peaks)):
        if p2 not in coll2_matchedindices and not pd.isnull(coll2.peaks[p2].matchedion):
            unique2.append(coll2.peaks[p2])

    return shared, unique1, unique2


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
def create_isodist2(monoiso, charge, maxval, adductmass=1.007276467):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    isodist = fast_calc_averagine_isotope_dist(monoiso, charge=charge)
    isodist[:, 1] *= maxval
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


#@njit(fastmath=True)
def optimize_shift2(config, centroids: np.ndarray, z, peakmz, tol=0.01, maxshift=2, gamma=1):
    peaks = []

    min_score_diff = 0.05

    # Limit max shifts if necessary
    if z < 3:
        maxshift = 1
    elif z < 6:
        maxshift = 2
    else:
        maxshift = maxshift

    # peakmz = centroids[np.argmax(centroids[:, 1]), 0]
    isodist, massdist, monoiso = create_isodist_full(peakmz, z, centroids)

    # matchedindexes, isomatches = match_peaks(centroids, isodist)

    cent_intensities = find_matched_intensities(centroids[:, 0], centroids[:, 1], isodist[:, 0], maxshift, tolerance=5,
                                                z=z)

    shiftrange = np.arange(-maxshift, maxshift + 1)
    # shifts_scores = [[0, 0] for i in range(len(shiftrange))]
    shifts_scores = np.zeros((len(shiftrange), 2))
    # overlaps = []
    bestshift = -1
    sum = -1
    meanratio = 1
    for i, shift in enumerate(shiftrange):
        s = calculate_cosinesimilarity(cent_intensities ** gamma, isodist[:, 1] ** gamma, shift, maxshift)
        shifts_scores[i, 0] = shift
        shifts_scores[i, 1] = s
        #print(z, peakmz, shift, s)

        if s > sum:
            sum = s
            bestshift = shift

    max_score = max([x[1] for x in shifts_scores])
    if max_score < config.css_thresh:
        return peaks

    accepted_shifts = [x for x in shifts_scores if x[1] > max_score - min_score_diff]

    for shift in accepted_shifts:
        shiftmass = shift[0] * mass_diff_c
        monoiso_new = monoiso + shiftmass
        massdist_new = massdist.copy()
        massdist_new[:, 0] = massdist[:, 0] + shiftmass
        # Correct m/z based on shift
        shiftmz = shiftmass / z
        isodist_new = isodist.copy()
        isodist_new[:, 0] = isodist[:, 0] + shiftmz
        peakmz_new = peakmz + shiftmz
        # Match it again
        matchedindexes, isomatches = match_peaks(centroids, isodist_new, tol=tol)

        matchediso = isodist[np.array(isomatches).astype(int)]

        areaper = np.sum(matchediso[:, 1]) / np.sum(isodist[:, 1])

        minpeaks = config.minpeaks
        if z == 1:
            if len(matchedindexes) == 2:
                if isomatches[0] == 0 and isomatches[1] == 1:
                    int1 = centroids[matchedindexes[0], 1]
                    int2 = centroids[matchedindexes[1], 1]
                    if int1 == 0:
                        ratio = 0
                    else:
                        ratio = int2 / int1
                    if config.plusoneintwindow[0] < ratio < config.plusoneintwindow[1]:
                        config.minpeaks = 2
                        areaper = 1

        passing = False
        if areaper >= config.minmatchper or shift[1] >= config.css_thresh:
            passing = True

        if len(matchedindexes) >= config.minpeaks and passing:
            m = MatchedPeak(int(z), float(peakmz_new), centroids, isodist_new, matchedindexes, isomatches)
            m.scan = config.activescan
            m.rt = config.activescanrt
            m.ms_order = config.activescanorder
            m.monoiso = monoiso_new
            m.massdist = massdist_new
            m.avgmass = np.average(massdist_new[:, 0], weights=massdist_new[:, 1])
            m.peakmass = np.sum(massdist_new[:, 0] * massdist_new[:, 1]) / np.sum(massdist_new[:, 1])
            peaks.append(m)

    return peaks


# @njit(fastmath=True)
def _shift(centroids, z, peakmz, tol=0.01, maxshift=2, gamma=0.5):
    # Limit max shifts if necessary
    if z < 3:
        maxshift = 1
    elif z < 6:
        maxshift = 2
    else:
        maxshift = maxshift

    # peakmz = centroids[np.argmax(centroids[:, 1]), 0]
    isodist, massdist, monoiso = create_isodist_full(peakmz, z, centroids)

    # matchedindexes, isomatches = match_peaks(centroids, isodist)

    cent_intensities = find_matched_intensities(centroids[:, 0], centroids[:, 1], isodist[:, 0], maxshift, tolerance=5,
                                                z=z)

    # mc = centroids[matchedindexes]
    # mi = isodist[isomatches]

    # mc = [centroids[i, 1] for i in matchedindexes]
    # mi = [isodist[i, 1] for i in isomatches]

    # mc = np.array(mc)
    # mi = np.array(mi)

    # Renormalize so that the sums are equal
    # s1 = np.sum(mc)
    # s2 = np.sum(mi)
    # if s2 == 0:
    #   return isodist, matchedindexes, isomatches, 0, 0, massdist
    # mi *= s1 / s2
    # isodist[:, 1] *= s1 / s2

    # if maxshift > len(mc):
    #   maxshift = len(mc) - 1

    shiftrange = np.arange(-maxshift, maxshift + 1)

    # overlaps = []
    bestshift = -1
    sum = -1
    meanratio = 1
    for i, shift in enumerate(shiftrange):
        # TODO: This is actually unsafe, if fast. If there are gaps, it will roll over them
        # However, I'm not sure how to do it better without sacrificing a decent amount of speed
        # Gonna fix it in the c code though
        s = calculate_cosinesimilarity(cent_intensities ** gamma, isodist[:, 1] ** gamma, shift, maxshift)
        # rmsd = calculate_RMSD(cent_intensities**gamma, isodist[:, 1]**gamma, shift, maxshift)
        # s = 1 / rmsd
        # print(shift, s)
        # overlap = mc ** gamma * roll ** gamma
        # s = np.sum(overlap)
        # s = calculate_cosinesimilarity(mc*gamma, roll*gamma)

        # calculate rmsd between roll and mc
        # s = 1 / np.sqrt(np.sum((mc - roll) ** 2))

        # overlaps.append(s)
        if s > sum:
            sum = s
            bestshift = shift
            # meanratio = np.average(mc / roll, weights=mc)
    # bestshift = shiftrange[np.argmax(overlaps)]
    # isodist[:, 1] *= meanratio

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
    # print(bestshift, sum, isodist[0])
    return isodist, matchedindexes, isomatches, peakmz, monoiso, massdist


@njit(fastmath=True)
def calculate_cosinesimilarity(cent_intensities, iso_intensities, shift: int, max_shift: int):
    ab = 0
    a2 = 0
    b2 = 0

    for i in range(len(iso_intensities)):
        a_val = cent_intensities[i + max_shift + shift]
        b_val = iso_intensities[i]
        ab += a_val * b_val
        a2 += a_val ** 2
        b2 += b_val ** 2
    if ab == 0 or a2 == 0 or b2 == 0:
        return 0
    return ab / (math.sqrt(a2) * math.sqrt(b2))


def calculate_RMSD(cent_intensities, iso_intensities, shift, max_shift):
    sse = 0
    for i in range(len(iso_intensities)):
        sse += (cent_intensities[i + max_shift + shift] - iso_intensities[i]) ** 2
    return math.sqrt(sse / len(iso_intensities))


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
    ppm_tol = 5
    lowest_idx = 0
    m1 = []
    m2 = []
    for peak1_idx in range(spec1_mz.shape[0]):
        mz = spec1_mz[peak1_idx]
        low_bound = mz - (mz * ppm_tol * 1e-6)
        high_bound = mz + (mz * ppm_tol * 1e-6)
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
def find_matched_intensities(spec1_mz: np.ndarray, spec1_intensity: np.ndarray, spec2_mz: np.ndarray,
                             max_shift: int, tolerance: float, z: int) -> List[float]:
    """
    Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from low to high m/z).

    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_mz:
        Theoretical isotope peak m/z values as numpy array. Peak mz values must be ordered.
    tolerance
        Peaks will be considered a match when within [tolerance] ppm of the theoretical value.
    shift
        Shift peaks of second spectra by shift. The default is 0.

    Returns
    -------
    matches
        List containing entries of type (centroid intensity)

    """
    mass_diff_c = 1.0033

    query_mzs = [float(0) for i in range(len(spec2_mz) + 2 * max_shift)]

    cent_intensities = np.zeros(len(query_mzs))

    mono_idx = max_shift

    for i in range(max_shift + 1):
        if i == 0:
            continue
        else:
            query_mzs[max_shift + len(spec2_mz) - 1 + i] = spec2_mz[len(spec2_mz) - 1] + (i * mass_diff_c) / z
            query_mzs[mono_idx - i] = spec2_mz[0] - (i * mass_diff_c) / z

    for i in range(len(spec2_mz)):
        query_mzs[i + max_shift] = spec2_mz[i]

    for i in range(len(query_mzs)):
        mz = query_mzs[i]
        low_bound = mz - (mz * tolerance * 1e-6)
        high_bound = mz + (mz * tolerance * 1e-6)
        for j in range(len(spec1_mz)):
            mz2 = spec1_mz[j]
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                continue
            else:
                cent_intensities[i] += spec1_intensity[j]

    return cent_intensities


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

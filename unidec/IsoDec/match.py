
import numpy as np

from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol_withnearest, fastwithin_abstol
from typing import List, Tuple
from numba import njit
import matplotlib as mpl
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist, mass_diff_c, fast_calc_averagine_isotope_dist_dualoutput
import pickle as pkl
import unidec.tools as ud
import unidec.IsoDec.msalign_export as msalign
import math
import pandas as pd


class IsoDecConfig: # (udstruct.UniDecConfig):
    def __init__(self):
        """
        Configuration class for IsoDec Engine. Holds the key parameters for the data processing and deconvolution.
        """
        super().__init__()
        self.filepath = ""
        self.batch_size = 32
        self.test_batch_size = 2048
        self.current_KD_round = 0
        self.activescan = -1
        self.activescanrt = -1
        self.activescanorder = -1
        self.meanpeakspacing_thresh = 0.01  ##This needs to be optimized.
        self.isolation_width = 3.0
        self.adjust_css = True

        self.adductmass = 1.007276467
        self.mass_diff_c = 1.0033

        self.verbose = False

        self.peakwindow = 80

        self.phaseres = 8

        self.matchtol = 5  # In ppm
        self.minpeaks = 3
        self.peakthresh = 0.0001
        self.css_thresh = 0.7
        self.maxshift = 3  # This will get overwritten for smaller z, where it's dangerous to have more than 1 or 2
        self.mzwindow = [-1.05, 2.05]
        self.plusoneintwindow = [0.1, 0.6]
        self.knockdown_rounds = 5
        self.min_score_diff = 0.1
        self.minareacovered = 0.20
        self.minusoneaszero = True
        self.isotopethreshold = 0.01
        self.datathreshold = 0.05
        self.zscore_threshold = 0.95

    def set_scan_info(self, s, reader=None):
        """
        Sets the active scan info
        :param s: The current scan
        :param reader: The reader object
        :return: None
        """
        self.activescan = s
        if reader is not None:
            self.activescanrt = reader.get_scan_time(s)
            self.activescanorder = reader.get_ms_order(s)
            polarity = reader.polarity
            if polarity == "Negative":
                self.adductmass = -1.007276467
            else:
                self.adductmass = 1.007276467



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

    def get_z_dist(self):
        return np.array([p.z for p in self.peaks])

    def get_z(self, item):
        return self.peaks[item].z

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

    def group_peaks(self, peaks):
        #Groups peaks with identical intensities
        grouped_peaks = []
        grouped_peak_intensities = np.array([])
        for p in peaks:
            if len(grouped_peaks) == 0:
                grouped_peaks.append(p)
                grouped_peak_intensities = np.append(grouped_peak_intensities, p.mz)
                continue
            else:
                nearest_idx = fastnearest(grouped_peak_intensities, p.matchedintensity)
                intensity_diff = np.abs(p.matchedintensity - grouped_peak_intensities[nearest_idx])
                if intensity_diff < 0.001 and p.z == grouped_peaks[nearest_idx].z:
                    grouped_peaks[nearest_idx].monoisos.append(p.monoiso)
                else:
                    #Insert into the grouped peaks/grouped_mzs at the proper position
                    if p.matchedintensity > grouped_peak_intensities[nearest_idx]:
                        nearest_idx += 1
                    grouped_peaks.insert(nearest_idx, p)
                    grouped_peak_intensities = np.insert(grouped_peak_intensities, nearest_idx, p.matchedintensity)
        return grouped_peaks


    def add_pk_to_masses(self, pk, ppmtol):
        """
        Checks if an existing mass matches to this peak, if so adds it to that mass, otherwise creates a new mass
        The list of masses is constantly kept in order of monoisotopic mass.
        """
        if len(self.masses) == 0:
            self.monoisos = np.append(self.monoisos, pk.monoiso)
            self.masses.append(MatchedMass(pk))

        else:
            #need to update this code to allow for overlapped monoisos.
            idx = fastnearest(self.monoisos, pk.monoiso)
            nearest_mass = self.monoisos[idx]
            if pk.scan - self.masses[idx].scans[len(self.masses[idx].scans) - 1] <= 100 and ud.within_ppm(
                    self.masses[idx].monoiso, pk.monoiso, ppmtol):
                self.masses[idx].scans = np.append(self.masses[idx].scans, pk.scan)
                self.masses[idx].totalintensity += pk.matchedintensity

                if pk.matchedintensity is not None:
                    if pk.matchedintensity > self.masses[idx].maxintensity:
                        self.masses[idx].maxintensity = pk.matchedintensity
                        self.masses[idx].maxscan = pk.scan
                        self.masses[idx].maxrt = pk.rt

                if not np.isin(self.masses[idx].zs, pk.z).any():
                    self.masses[idx].zs = np.append(self.masses[idx].zs, pk.z)
                    self.masses[idx].mzs = np.append(self.masses[idx].mzs, pk.mz)
                    self.masses[idx].mzints = np.append(self.masses[idx].mzints, pk.peakint)
                    self.masses[idx].isodists = np.vstack([self.masses[idx].isodists, pk.isodist])

            else:
                if pk.monoiso > nearest_mass:
                    idx += 1
                if idx == len(self.monoisos):
                    self.monoisos = np.append(self.monoisos, pk.monoiso)
                    self.masses.append(MatchedMass(pk))
                else:
                    self.monoisos = np.insert(self.monoisos, idx, pk.monoiso)
                    self.masses.insert(idx, MatchedMass(pk))

    def add_pk_to_masses2(self, pk, ppmtol, rt_tol=2, maxshift=3):
        """
        Checks if an existing mass matches to this peak, if so adds it to that mass, otherwise creates a new mass
        The list of masses is constantly kept in order of monoisotopic mass.
        """
        if len(self.masses) == 0:
            self.monoisos = np.append(self.monoisos, pk.monoiso)
            self.masses.append(MatchedMass(pk))

        else:
            #Get the index of the nearest mass
            idx, indices = fastwithin_abstol_withnearest(self.monoisos, pk.monoiso, maxshift*1.25)
            nearest_mass = self.monoisos[idx]
            matched_indices = []
            if len(indices) > 0:
                for i in indices:
                    #Filter out the indices of peaks that are too far away in RT
                    if np.abs(pk.rt - self.masses[i].maxrt) > rt_tol:
                        continue
                    for m in self.masses[i].monoisos:
                        if ud.within_ppm(m, pk.monoiso, ppmtol):
                            matched_indices.append(i)
                            break
                        if i in matched_indices:
                            break

                if len(matched_indices) == 0:
                    if pk.monoiso > nearest_mass:
                        idx += 1
                    if idx == len(self.monoisos):
                        self.monoisos = np.append(self.monoisos, pk.monoiso)
                        self.masses.append(MatchedMass(pk))
                    else:
                        self.monoisos = np.insert(self.monoisos, idx, pk.monoiso)
                        self.masses.insert(idx, MatchedMass(pk))

                elif len(matched_indices) == 1:
                    #Add the peak to the matched mass
                    self.masses[matched_indices[0]].scans = np.append(self.masses[matched_indices[0]].scans, pk.scan)
                    self.masses[idx].totalintensity += pk.matchedintensity
                    if pk.matchedintensity is not None:
                        if pk.matchedintensity > self.masses[matched_indices[0]].maxintensity:
                            self.masses[matched_indices[0]].maxintensity = pk.matchedintensity
                            self.masses[matched_indices[0]].maxscan = pk.scan
                            self.masses[matched_indices[0]].maxrt = pk.rt
                    if pk.z not in self.masses[matched_indices[0]].zs:
                        self.masses[matched_indices[0]].zs = np.append(self.masses[matched_indices[0]].zs, pk.z)
                        self.masses[matched_indices[0]].mzs = np.append(self.masses[matched_indices[0]].mzs, pk.mz)
                else:
                    #The new peak matches to multiple, select the one with the closest RT
                    rt_diffs = np.abs([pk.rt - self.masses[i].maxrt for i in matched_indices])
                    closest = np.argmin(rt_diffs)

                    self.masses[matched_indices[closest]].scans = np.append(self.masses[matched_indices[closest]].scans, pk.scan)
                    self.masses[idx].totalintensity += pk.matchedintensity
                    if pk.matchedintensity is not None:
                        if pk.matchedintensity > self.masses[matched_indices[closest]].maxintensity:
                            self.masses[matched_indices[closest]].maxintensity = pk.matchedintensity
                            self.masses[matched_indices[closest]].maxscan = pk.scan
                            self.masses[matched_indices[closest]].maxrt = pk.rt
                    if pk.z not in self.masses[matched_indices[closest]].zs:
                        self.masses[matched_indices[closest]].zs = np.append(self.masses[matched_indices[closest]].zs, pk.z)
                        self.masses[matched_indices[closest]].mzs = np.append(self.masses[matched_indices[closest]].mzs, pk.mz)


            else:
                #Find where to insert a new mass
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

        # TODO: Fix masses as well with this
        masses_to_remove = []
        for i1 in range(len(self.masses)):
            if i1 in masses_to_remove:
                continue
            for i2 in range(len(self.masses)):
                i2_matched = False
                if i2 in masses_to_remove or i1 == i2:
                    continue
                for mono1 in self.masses[i1].monoisos:
                    for mono2 in self.masses[i2].monoisos:
                        if ud.within_ppm(mono1, mono2, ppm_tolerance):
                            #Remove the mass at i2 and add unmatched monoisos to the mass at i1
                            for mono in self.masses[i2].monoisos:
                                if any(ud.within_ppm(x, mono, ppm_tolerance) for x in self.masses[i1].monoisos):
                                    continue
                                else:
                                    self.masses[i1].monoisos = np.append(self.masses[i1].monoisos, mono)
                            #Add charge states from i2 that don't exist in i1
                            for z in self.masses[i2].zs:
                                if z not in self.masses[i1].zs:
                                    self.masses[i1].zs = np.append(self.masses[i1].zs, z)
                            masses_to_remove.append(i2)
                            i2_matched = True
                            break
                        if i2_matched:
                            break
                    if i2_matched:
                        break
                if i2_matched:
                    break
        self.masses = [self.masses[i] for i in range(len(self.masses)) if i not in masses_to_remove]
        return

    def filter(self, minval, maxval, type="mz"):
        if type == "mz":
            self.peaks = [p for p in self.peaks if minval < p.mz < maxval]
        elif type == "monoiso":
            self.peaks = [p for p in self.peaks if minval < p.monoiso < maxval]
        else:
            print("Type not recognized")
        return self

    def copy_to_string(self):
        outstring = ""
        for p in self.peaks:
            for m in p.monoisos:
                outstring += str(m) + "\n"
        return outstring

    def export_msalign(self, config, reader, filename="export.msalign", act_type="HCD", max_precursors=None):
        print("Exporting to", filename, "with activation type", act_type, "N:", len(self.peaks))
        ms1_scan_dict = msalign.sort_by_scan_order(self, 1)
        ms1_features = 0
        for k, v in ms1_scan_dict.items():
            ms1_features += len(v)
        print("MS1 Features:", ms1_features)
        ms2_scan_dict = msalign.sort_by_scan_order(self, 2)
        ms2_features = 0
        for k, v in ms2_scan_dict.items():
            ms2_features += len(v)
        print("MS2 Features:", ms2_features)

        msalign.write_ms1_msalign(ms1_scan_dict, ms2_scan_dict, filename, config)
        msalign.write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, filename, config, max_precursors=max_precursors, act_type=act_type)

    def to_df(self):
        """
        Convert a MatchedCollection object to a pandas dataframe
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
        print("Exporting to", filename)
        df = self.to_df()
        df.to_csv(filename, sep="\t", index=False)

    def to_mass_spectrum(self, binsize=0.1):
        if len(self.peaks) == 0:
            print("Zero Sized Peak Array")
            return np.array([])
        minmass = np.min([p.monoiso for p in self.peaks]) - 10 * binsize
        maxmass = np.max([p.monoiso for p in self.peaks]) + 100 * binsize

        # Round minmass and maxmass to nearest binsize
        minmass = np.floor(minmass / binsize) * binsize
        maxmass = np.ceil(maxmass / binsize) * binsize

        masses = np.arange(minmass, maxmass, binsize)

        intensities = np.zeros(len(masses))
        for p in self.peaks:
            idx = round((p.monoiso - minmass) / binsize)
            intensities[idx] += p.matchedintensity

        return np.transpose([masses, intensities])


class MatchedMass:
    """
    Matched mass object for collecting data on MatchedPeaks with matched masses.
    """

    def __init__(self, pk):
        self.monoiso = pk.monoiso
        self.monoisos = pk.monoisos
        self.scans = np.array([pk.scan])
        self.maxintensity = pk.matchedintensity
        self.maxscan = pk.scan
        self.maxrt = pk.rt
        self.mzs = np.array([pk.mz])
        self.zs = np.array([pk.z])
        self.totalintensity = pk.matchedintensity
        self.mzints = np.array([pk.peakint])
        self.isodists = pk.isodist



class MatchedPeak:
    """
    Matched peak object for collecting data on peaks with matched distributions
    """
    def __init__(self, z, mz, centroids=None, isodist=None, matchedindexes=None, isomatches=None):
        self.mz = mz
        self.z = z
        self.centroids = centroids
        self.isodist = isodist
        self.matchedintensity = 0
        self.peakint = 0
        self.matchedcentroids = None
        self.matchedisodist = None
        self.matchedindexes = None
        self.isomatches = None
        self.color = "g"
        self.scan = -1
        self.rt = -1
        self.ms_order = -1
        self.massdist = None
        self.monoisos = []
        self.peakmass = -1
        self.avgmass = -1
        self.startindex = -1
        self.endindex = -1
        self.monoiso = -1
        self.bestshift = -1
        self.acceptedshifts = None
        self.matchedions = []
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


def read_msalign_to_matchedcollection(file, data=None, mz_type="monoiso"):
    """
    Read an msalign file to a MatchedCollection object
    :param file: Path to msalign file
    :return: MatchedCollection object
    """

    mc = MatchedCollection()

    current_scan = 1
    current_rt = 0
    scan_peaks = []
    with open(file, "r") as f:
        for line in f:
            # print(line)
            if line[0] == "#":
                continue
            elif "BEGIN IONS" in line:
                continue
            elif "END IONS" in line:
                grouped_peaks = mc.group_peaks(scan_peaks)
                mc.peaks.extend(grouped_peaks)
                scan_peaks = []
                continue
            elif "SCANS" in line:
                current_scan = int(line.split("=")[1])
                continue
            elif "RETENTION_TIME" in line:
                current_rt = float(line.split("=")[1])/60
                continue
            elif "=" not in line and len(line) > 1:
                split = line.split("\t")
                mz = (float(split[0]) + (1.007276467 * int(split[2]))) / int(split[2])
                pk = MatchedPeak(int(split[2]), mz)
                pk.monoisos = [float(split[0])]
                pk.scan = current_scan
                pk.monoiso = float(split[0])
                pk.matchedintensity = float(split[1])
                pk.rt = current_rt
                scan_peaks.append(pk)

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
                if mz_type == "monoiso":
                    pk.mz = isodist[np.argmax(isodist[:, 1]), 0]
                pk.isodist = isodist
    print("Loaded", len(mc.peaks), "peaks from msalign file")
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
        pk.monoisos = [pk.monoiso]
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

def compare_matchedcollections(coll1, coll2, ppmtol=50, objecttocompare="monoisos", maxshift=3, ignorescan=False):
    unique1 = []
    shared = []
    unique2 = []

    diffs = []

    if not ignorescan:
        scans = np.unique([p.scan for p in coll1.peaks])
        scans2 = np.unique([p.scan for p in coll2.peaks])
        scans = np.union1d(scans, scans2)
    else:
        scans = [0]
    for s in scans:
        if not ignorescan:
            # Pull out peaks for each individual scan.
            coll1_sub = MatchedCollection().add_peaks([p for p in coll1.peaks if p.scan == s])
            coll2_sub = MatchedCollection().add_peaks([p for p in coll2.peaks if p.scan == s])
        else:
            coll1_sub = coll1
            coll2_sub = coll2

        # Sort the peak lists and the masses by monoisotopic mass
        coll1_sub.peaks = sorted(coll1_sub.peaks, key=lambda x: x.monoiso)
        coll2_sub.peaks = sorted(coll2_sub.peaks, key=lambda x: x.monoiso)
        masses1 = np.array([p.monoiso for p in coll1_sub.peaks])
        masses2 = np.array([p.monoiso for p in coll2_sub.peaks])
        # print(masses1, masses2)
        coll2_matched_indices = []

        for i in range(len(coll1_sub.peaks)):
            foundmatch = False
            within_tol = fastwithin_abstol(masses2, masses1[i], int(math.ceil(maxshift + 0.5)))
            if len(within_tol) > 0:
                for j in within_tol:
                    if coll1_sub.peaks[i].z == coll2_sub.peaks[j].z:
                        # now check if a monoiso is within the ppmtol
                        for m1 in coll1_sub.peaks[i].monoisos:
                            for m2 in coll2_sub.peaks[j].monoisos:
                                if ud.within_ppm(m1, m2, ppmtol):
                                    coll2_matched_indices.append(j)
                                    foundmatch = True
                                    diffs.append((m1 - m2) / m1 * 1e6)

                                else:
                                    if maxshift > 0:
                                        diff = np.abs(m1 - m2)
                                        # print(m1, m2, diff)
                                        if diff < maxshift * 1.1:
                                            mm_count = int(round(m1 - m2))
                                            newm2 = m2 + mm_count * 1.0033
                                            if ud.within_ppm(m1, newm2, ppmtol):
                                                coll2_matched_indices.append(j)
                                                foundmatch = True
                                                diffs.append((m1 - newm2) / m1 * 1e6)

            if foundmatch:
                shared.append(coll1_sub.peaks[i])
            else:
                unique1.append(coll1_sub.peaks[i])

        coll2_unique_matchedinds = list(set(coll2_matched_indices))
        for i in range(len(coll2_sub.peaks)):
            if i not in coll2_matched_indices:
                unique2.append(coll2_sub.peaks[i])

    shared = MatchedCollection().add_peaks(shared)
    unique1 = MatchedCollection().add_peaks(unique1)
    unique2 = MatchedCollection().add_peaks(unique2)

    print("Avg PPM Diff:", np.mean(diffs))

    return shared, unique1, unique2

def compare_matchedmasses(coll1, coll2, ppmtol=20, maxshift=3, rt_tol = 2, f_shared_zs=0.5):
    #make a list of the masses within the list of matched masses contained by each collection
    monoisos1 = np.array([m.monoiso for m in coll1.masses])
    monoisos2 = np.array([m.monoiso for m in coll2.masses])

    sum_intensity1 = np.sum(np.array([m.totalintensity for m in coll1.masses]))
    sum_intensity2 = np.sum(np.array([m.totalintensity for m in coll2.masses]))

    #now order the monoisos and apply that ordering to the masses
    sorted_indices1 = np.argsort(monoisos1)
    sorted_indices2 = np.argsort(monoisos2)
    monoisos1 = monoisos1[sorted_indices1]
    monoisos2 = monoisos2[sorted_indices2]
    coll1.masses = [coll1.masses[i] for i in sorted_indices1]
    coll2.masses = [coll2.masses[i] for i in sorted_indices2]

    matched_coll1 = []
    matched_coll2 = []

    unique_coll1 = []
    unique_coll2 = []

    #Initially we'll store the indices of the matches/unique masses
    for i in range(len(monoisos1)):
        foundmatch = False
        within_tol = fastwithin_abstol(monoisos2, monoisos1[i], int(math.ceil(maxshift + 0.5)))
        if len(within_tol) > 0:
            #Order the within_tol by descending matched intensity
            within_tol = sorted(within_tol, key=lambda x: coll2.masses[x].totalintensity, reverse=True)
            for j in within_tol:
                #ensure that the retention times are close enough
                if abs(coll1.masses[i].maxrt - coll2.masses[j].maxrt) < rt_tol:
                    #ensure that the monoisos match within the ppmtol allowing the maxshift
                    possible_monoisos = [monoisos1[i] + x * 1.0033 for x in range(-maxshift, maxshift + 1)]
                    if any([ud.within_ppm(m, monoisos2[j], ppmtol) for m in possible_monoisos]):
                        #ensure that the z values overlap enough
                        z_overlap = np.intersect1d(coll1.masses[i].zs, coll2.masses[j].zs)
                        min_frac = len(z_overlap) / max(len(coll1.masses[i].zs), len(coll2.masses[j].zs))
                        if min_frac >= f_shared_zs:
                            matched_coll1.append(i)
                            foundmatch = True
                            break
        if not foundmatch:
            unique_coll1.append(i)

    for i in range(len(monoisos2)):
        if i not in matched_coll2:
            foundmatch = False
            within_tol = fastwithin_abstol(monoisos1, monoisos2[i], int(math.ceil(maxshift + 0.5)))
            if len(within_tol) > 0:
                for j in within_tol:
                    #ensure that the retention times are close enough
                    if abs(coll1.masses[j].maxrt - coll2.masses[i].maxrt) < rt_tol:
                        #ensure that the monoisos match within the ppmtol allowing the maxshift
                        possible_monoisos = [monoisos2[i] + x * 1.0033 for x in range(-maxshift, maxshift + 1)]
                        if any([ud.within_ppm(m, monoisos1[j], ppmtol) for m in possible_monoisos]):
                            #ensure that the z values overlap enough
                            z_overlap = np.intersect1d(coll1.masses[j].zs, coll2.masses[i].zs)
                            min_frac = len(z_overlap) / max(len(coll1.masses[j].zs), len(coll2.masses[i].zs))
                            if min_frac >= f_shared_zs:
                                matched_coll2.append(i)
                                foundmatch = True
                                break
            if not foundmatch:
                unique_coll2.append(i)


    #Now we'll convert the indices to the actual masses
    matched_coll1 = [coll1.masses[i] for i in matched_coll1]
    matched_coll2 = [coll2.masses[i] for i in matched_coll2]
    unique_coll1 = [coll1.masses[i] for i in unique_coll1]
    unique_coll2 = [coll2.masses[i] for i in unique_coll2]

    matched_intensity1 = np.sum(np.array([m.totalintensity for m in matched_coll1]))
    matched_intensity2 = np.sum(np.array([m.totalintensity for m in matched_coll2]))



    percent_matched1 = (len(matched_coll1) / len(coll1.masses)) * 100
    percent_matched2 = (len(matched_coll2) / len(coll2.masses)) * 100
    percent_intensity1 = (matched_intensity1 / sum_intensity1) * 100
    percent_intensity2 = (matched_intensity2 / sum_intensity2) * 100
    print("% Masses Matched Collection 1:", (len(matched_coll1) / len(coll1.masses)) * 100)
    print("% Masses Matched Collection 2:", (len(matched_coll2) / len(coll2.masses)) * 100)
    print("% Matched Intensity Collection 1:", (matched_intensity1 / sum_intensity1) * 100)
    print("% Matched Intensity Collection 2:", (matched_intensity2 / sum_intensity2) * 100)

    return percent_matched1, percent_matched2, percent_intensity1, percent_intensity2, len(matched_coll1), len(matched_coll2), len(monoisos1), len(monoisos2)


def compare_annotated(l1, l2, ppmtol, maxshift):
    unique_annotated = []
    unique_experimental = []
    shared = []

    annotated_mzs = []
    annotated_charges = []
    annotated_indices = []
    for i in range(len(l1.peaks)):
        annotated_mzs.append(l1[i].mz)
        annotated_charges.append(l1[i].z)
        annotated_indices.append(i)
        currentshift = 1
        while currentshift <= maxshift:
            negative_mz = l1[i].mz - (1 / l1[i].z) * currentshift
            annotated_mzs.append(negative_mz)
            annotated_charges.append(l1[i].z)
            annotated_indices.append(i)
            positive_mz = l1[i].mz + (1 / l1[i].z) * currentshift
            annotated_mzs.append(positive_mz)
            annotated_charges.append(l1[i].z)
            annotated_indices.append(i)
            currentshift += 1

    #Sort annotated mzs and apply that sort to the other lists
    annotated_mzs, annotated_charges, annotated_indices = zip(*sorted(zip(annotated_mzs, annotated_charges, annotated_indices)))

    matched_annotated_indices = []
    for i in range(len(l2.peaks)):
        foundmatch = False
        curr_mz = l2[i].mz
        curr_tol = (curr_mz / 1e6) * ppmtol
        matched_indices = fastwithin_abstol(np.array(annotated_mzs), curr_mz, curr_tol)
        for index in matched_indices:
            matched_charge = annotated_charges[index]
            if matched_charge == l2[i].z:
                shared.append(l2[i])
                matched_annotated_indices.append(annotated_indices[index])
                foundmatch = True
                break
        if not foundmatch:
            unique_experimental.append(l2[i])

    unique_annotated = MatchedCollection().add_peaks([l1[i] for i in range(len(l1.peaks)) if i not in matched_annotated_indices])
    shared = MatchedCollection().add_peaks(shared)
    unique_experimental = MatchedCollection().add_peaks(unique_experimental)

    return shared, unique_annotated, unique_experimental

def is_close(mz1, mz2, tolerance):
    return abs(mz1 - mz2) <= tolerance


def compare_matched_ions(coll1, coll2, other_alg=None):
    shared1 = []
    unique1 = []

    shared2 = []
    unique2 = []

    for p1 in coll1.peaks:
        matched = False
        for p2 in coll2.peaks:
            if set(p1.matchedions) & set(p2.matchedions):
                if p1 not in shared1:
                    shared1.append(p1)
                matched = True
        if not matched and len(p1.matchedions) > 0:
            unique1.append(p1)

    for p2 in coll2.peaks:
        matched = False
        for p1 in coll1.peaks:
            if set(p1.matchedions) & set(p2.matchedions):
                if p2 not in shared2:
                    shared2.append(p2)
                matched = True
        if not matched and len(p2.matchedions) > 0:
            unique2.append(p2)

    other_name = "Other"
    if other_alg is not None:
        other_name = other_alg

    print("% Shared IsoDec:", len(shared1) / (len(shared1) + len(unique1)) * 100)
    print("% Shared", other_name, ":", len(shared2) / (len(shared2) + len(unique2)) * 100)
    print("% Unique IsoDec:", len(unique1) / (len(shared1) + len(unique1)) * 100)
    print("% Unique", other_name, ":", len(unique2) / (len(shared2) + len(unique2)) * 100)

    print("Number shared IsoDec:", len(shared1))
    print("Number shared", other_name, ":", len(shared2))
    print("Number unique IsoDec:", len(unique1))
    print("Number unique", other_name, ":", len(unique2))

    #Create matched collections of shared and unique peaks
    shared1 = MatchedCollection().add_peaks(shared1)
    unique1 = MatchedCollection().add_peaks(unique1)
    shared2 = MatchedCollection().add_peaks(shared2)
    unique2 = MatchedCollection().add_peaks(unique2)


    return shared1, unique1, shared2, unique2

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


#@njit(fastmath=True)
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
    charge = float(charge)
    isodist = fast_calc_averagine_isotope_dist(monoiso, charge=charge)
    isodist[:, 1] *= maxval
    return isodist


#@njit(fastmath=True)
def create_isodist_full(peakmz, charge, data, adductmass=1.007276467, isotopethresh: float = 0.01):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    mass = (peakmz - adductmass) * charge
    isodist, massdist = fast_calc_averagine_isotope_dist_dualoutput(mass, charge=charge, adductmass=adductmass,
                                                                    isotopethresh=isotopethresh)
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
def get_accepted_shifts(cent_intensities, isodist, maxshift, min_score_diff, css_thresh, minusoneaszero=True):
    shiftrange = np.arange(-maxshift, maxshift + 1)
    # shifts_scores = [[0, 0] for i in range(len(shiftrange))]
    shifts_scores = np.zeros((len(shiftrange), 2))
    # overlaps = []
    bestshift = -1
    sum = -1
    meanratio = 1
    for i, shift in enumerate(shiftrange):
        s = calculate_cosinesimilarity(cent_intensities, isodist[:, 1], shift, maxshift,
                                       minusoneareaszero=minusoneaszero)
        shifts_scores[i, 0] = shift
        shifts_scores[i, 1] = s
        if s > sum:
            sum = s
            bestshift = shift
    max_score = np.amax(shifts_scores[:, 1])
    if max_score < css_thresh:
        return None

    accepted_shifts = [x for x in shifts_scores if x[1] > max_score - min_score_diff]

    return accepted_shifts


@njit(fastmath=True)
def make_shifted_peak(shift: int, shiftscore: float, monoiso: float, massdist: np.ndarray, isodist: np.ndarray,
                      peakmz: float, z: int,
                      centroids: np.ndarray, matchtol: float,
                      minpeaks: int, p1low: float,
                      p1high: float, css_thresh: float, minareacovered: float, verbose=True):
    b1 = isodist[:, 1] > 0
    shiftmass = float(shift) * mass_diff_c
    monoiso_new = monoiso + shiftmass
    massdist_new = massdist[b1]
    massdist_new[:, 0] = massdist_new[:, 0] + shiftmass
    # Correct m/z based on shift
    shiftmz = shiftmass / float(z)
    isodist_new = isodist[b1]
    isodist_new[:, 0] = isodist_new[:, 0] + shiftmz
    peakmz_new = peakmz

    # Match it again
    matchedindexes, isomatches = find_matches(centroids, isodist_new, matchtol)
    matchediso = isodist_new[np.array(isomatches)]
    # if verbose:
    #    print("Find Matches:", matchedindexes, isomatches, matchtol, z)

    if z == 1:
        if len(matchedindexes) == 2:
            if isomatches[0] == 0 and isomatches[1] == 1:
                int1 = centroids[matchedindexes[0], 1]
                int2 = centroids[matchedindexes[1], 1]
                if int1 == 0:
                    ratio = 0
                else:
                    ratio = int2 / int1
                if p1low < ratio < p1high:
                    minpeaks = 2

    areacovered = np.sum(matchediso[:, 1]) / np.sum(centroids[:, 1])
    # print("Area Covered:", areacovered)
    matchedcentroids = centroids[np.array(matchedindexes)]
    # Find the top three most intense peaks in matched iso
    topthreeiso = np.sort(matchedcentroids[:, 1])[::-1][:minpeaks]
    # Find the top three most intense peaks in centroids
    topthreecent = np.sort(centroids[:, 1])[::-1][:minpeaks]

    # Check if the top three peaks in the matched iso are the same as the top three peaks in the centroids
    if np.array_equal(topthreeiso, topthreecent):
        topthree = True
    else:
        topthree = False

    if verbose:
        print("Matched Peaks:", len(matchedindexes), "Shift Score:", shiftscore, "Area Covered:",
              areacovered, "Top Three:", topthree)
    if len(matchedindexes) >= minpeaks and (
            shiftscore >= css_thresh) and (areacovered > minareacovered or topthree):
        return peakmz_new, isodist_new, matchedindexes, isomatches, monoiso_new, massdist_new
    else:
        if verbose:
            print("Failed Peak:", len(matchedindexes), shiftscore, areacovered, topthree)
            # Determine which of the tests it failed
            if len(matchedindexes) < minpeaks:
                print("Failed Min Peaks")
            if not shiftscore >= css_thresh:
                print("Failed CSS")
            if not (areacovered > minareacovered or topthree):
                if areacovered < minareacovered:
                    print("Failed Min Area Covered")
                if not topthree:
                    print("Failed Top Three")
        return None, None, None, None, None, None


# @njit(fastmath=True)
def optimize_shift2(config, centroids: np.ndarray, z, peakmz):
    peaks = []

    # Limit max shifts if necessary
    if z < 3:
        maxshift = 1
    elif z < 6:
        maxshift = 2
    else:
        maxshift = config.maxshift

    isodist, massdist, monoiso = create_isodist_full(peakmz, z, centroids, adductmass=config.adductmass,
                                                     isotopethresh=config.isotopethreshold)

    cent_intensities = find_matched_intensities(centroids[:, 0], centroids[:, 1], isodist[:, 0], maxshift,
                                                tolerance=config.matchtol, z=z, peakmz=peakmz)

    norm_factor = max(cent_intensities) / max(isodist[:, 1])
    isodist[:, 1] *= norm_factor

    accepted_shifts = get_accepted_shifts(cent_intensities, isodist, maxshift, config.min_score_diff,
                                          config.css_thresh, config.minusoneaszero)

    if config.verbose:
        print("Accepted Shifts:", accepted_shifts)
    if accepted_shifts is None:
        return None

    # The plan is to report a set of possible monoisotopics, but only include the centroids, isodist, massdist, and matches for the best one.
    m = MatchedPeak(int(z), float(peakmz), centroids, isodist, None, None)

    bestshift = accepted_shifts[np.argmax([x[1] for x in accepted_shifts])]
    m.bestshift = bestshift
    m.acceptedshifts = np.array(accepted_shifts)

    peakmz_new, isodist_new, matchedindexes, isomatches, monoiso_new, massdist_new = make_shifted_peak(
        bestshift[0], bestshift[1], monoiso, massdist, isodist, peakmz, z, centroids, config.matchtol, config.minpeaks,
        config.plusoneintwindow[0], config.plusoneintwindow[1], config.css_thresh,
        config.minareacovered,
        config.verbose)


    if peakmz_new is None:
        return None

    m.monoiso = monoiso_new
    m.scan = config.activescan
    m.rt = config.activescanrt
    m.ms_order = config.activescanorder
    m.isodist = isodist_new
    m.matchedindexes = matchedindexes
    m.matchedintensity = np.sum(isodist_new[:, 1])
    m.isomatches = isomatches
    m.massdist = massdist_new
    m.monoisos = [monoiso + accepted_shifts[i][0] * mass_diff_c for i in range(len(accepted_shifts))]
    m.avgmass = np.average(massdist_new[:, 0], weights=massdist_new[:, 1])
    m.peakmass = np.sum(massdist_new[:, 0] * massdist_new[:, 1]) / np.sum(massdist_new[:, 1])
    if config.verbose:
        print("Accepted Peak:", m)
    return m


# We already have this function in datatools
@njit(fastmath=True)
def calculate_cosinesimilarity(cent_intensities, iso_intensities, shift: int, max_shift: int,
                               minusoneareaszero: bool = True):
    ab = 0
    a2 = 0
    b2 = 0

    if minusoneareaszero:
        a_val = cent_intensities[max_shift + shift - 1]
        b_val = 0
        ab += a_val * b_val
        a2 += a_val ** 2
        b2 += b_val ** 2

    for i in range(len(iso_intensities)):
        a_val = cent_intensities[i + max_shift + shift]
        b_val = iso_intensities[i]
        ab += a_val * b_val
        a2 += a_val ** 2
        b2 += b_val ** 2
    if ab == 0 or a2 == 0 or b2 == 0:
        return 0
    return ab / (math.sqrt(a2) * math.sqrt(b2))


@njit(fastmath=True)
def find_matches(spec1: np.ndarray, spec2: np.ndarray,
                 tolerance: float) -> Tuple[List[int], List[int]]:
    """Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from
    low to high m/z).

    Parameters
    ----------
    spec1:
        Spectrum peak m/z and int values as numpy array. Peak mz values must be ordered.
    spec2:
        Isotope distribution peak m/z and int values as numpy array. Peak mz values must be ordered.
    tolerance
        Peaks will be considered a match when <= tolerance appart in ppm.
    Returns
    -------
    matches
        List containing entries of type (idx1, idx2).

    """
    lowest_idx = 0
    m1 = []
    m2 = []
    diff = spec1[0, 0] * tolerance * 1e-6

    for iso_idx in range(spec2.shape[0]):
        mz = spec2[iso_idx, 0]
        low_bound = mz - diff
        high_bound = mz + diff

        topint = 0
        topindex = -1

        for peak_idx in range(lowest_idx, spec1.shape[0]):
            mz2 = spec1[peak_idx, 0]
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak_idx + 1
            else:
                newint = spec1[peak_idx, 1]
                if newint > topint:
                    topint = newint
                    topindex = peak_idx

        if topindex != -1:
            m1.append(topindex)
            m2.append(iso_idx)

    return m1, m2


@njit(fastmath=True)
def find_matched_intensities(spec1_mz: np.ndarray, spec1_intensity: np.ndarray, spec2_mz: np.ndarray,
                             max_shift: int, tolerance: float, z: int, peakmz: float) -> List[float]:
    """
    Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from low to high m/z).

    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_mz:
        Theoretical isotope peak m/z values as numpy array. Peak mz values must be ordered.
    tolerance:
        Peaks will be considered a match when within [tolerance] ppm of the theoretical value.

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

    diff = (peakmz * tolerance * 1e-6)
    for i in range(len(query_mzs)):
        mz = query_mzs[i]
        low_bound = mz - diff
        high_bound = mz + diff
        for j in range(len(spec1_mz)):
            mz2 = spec1_mz[j]
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                continue
            else:
                value = spec1_intensity[j]
                if value > cent_intensities[i]:
                    cent_intensities[i] = spec1_intensity[j]

    return cent_intensities


def remove_noise_peaks(pks, noiselevel):
    """
    Remove peaks below a certain noise level
    :param pks: List of MatchedPeaks
    :param noiselevel: Noise level to remove
    :return: List of MatchedPeaks
    """
    newcollection = MatchedCollection()
    filteredpeaks = [p for p in pks if p.matchedintensity > noiselevel]
    newcollection.add_peaks(filteredpeaks)
    return newcollection


# Will generate list of matchedpeak objects from text file
# Just specify the desired field with charge
# Text file should be in the format: field(monoiso or mz) " " charge
def read_manual_annotations(path=None, delimiter=' ', data=None):
    if path is None:
        path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt"
    mc = MatchedCollection()
    peaks = []
    with open(path, "r") as file:
        for line in file:
            if line.strip():
                row = line.strip().split(delimiter, 1)
                peak = round(float(row[0]), 2)
                currCharge = row[1].strip()
                z = MatchedPeak(mz=peak, z=int(currCharge))
                peakmass = (peak - 1.007276467) * int(currCharge)
                monoiso = get_estimated_monoiso(peakmass)
                isodist = fast_calc_averagine_isotope_dist(monoiso, currCharge)
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

                z.isodist = isodist

                mc.add_peak(z)
                peaks.append(z)
    return mc


def get_estimated_monoiso(peakmass):
    """
    Estimates the monoisotopic mass from the peak mass.
    This is a lazy approximation, but is good enough for most purposes.
    --JGP
    Args:
        peakmass: Most abundant isotopologue mass
    Returns:
        estimated monoisotopic mass
    """
    most_intense_iso = (int)(0.0006*peakmass + 0.4074)
    return peakmass - (most_intense_iso * 1.0033)


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

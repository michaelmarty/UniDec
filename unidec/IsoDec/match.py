import numpy as np
from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol_withnearest
from typing import List, Tuple
from numba import njit
import matplotlib as mpl
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist, mass_diff_c, \
    fast_calc_averagine_isotope_dist_dualoutput
import pickle as pkl
import unidec.tools as ud
import unidec.IsoDec.msalign_export as msalign
import math
from copy import deepcopy
import pandas as pd
from unidec.modules.fwhmtools import fast_fwhm, ndis_std, intensity_decon
from unidec.modules.unidecstructure import IsoDecConfig


class MatchedCollection:
    """
    Class for collecting matched peaks
    """

    def __init__(self):
        self.peaks = []
        self.masses = []
        self.monoisos = np.array([])
        self.precursors = np.array([])
        self.colormap = mpl.colormaps.get_cmap("tab10")
        self.simdata = []
        self.sumsim = []

    def __str__(self):
        outstring = ""
        for p in self.peaks:
            outstring += str(p) + "\n"
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
        # Groups peaks with identical intensities
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
                    # Insert into the grouped peaks/grouped_mzs at the proper position
                    if p.matchedintensity > grouped_peak_intensities[nearest_idx]:
                        nearest_idx += 1
                    grouped_peaks.insert(nearest_idx, p)
                    grouped_peak_intensities = np.insert(grouped_peak_intensities, nearest_idx, p.matchedintensity)
        return grouped_peaks

    def add_pk_to_masses(self, pk, config=None, rt_tol=None, scan_tol=100):
        """
        Checks if an existing mass matches to this peak, if so adds it to that mass, otherwise creates a new mass
        The list of masses is constantly kept in order of monoisotopic mass.
        """
        if config is None:
            config = IsoDecConfig()

        if len(self.masses) == 0:
            self.monoisos = np.append(self.monoisos, pk.monoiso)
            self.masses.append(MatchedMass(pk))
            return

        # Get the index of the nearest mass
        idx, indices = fastwithin_abstol_withnearest(self.monoisos, pk.monoiso, config.maxshift * 1.25)
        nearest_mass = self.monoisos[idx]
        matched_indices = []

        # Loop through the indices and check for a match
        for i in indices:
            massobj = self.masses[i]
            # Check massobj for match
            check = massobj.check_if_match(pk, config, scan_tol=scan_tol, rt_tol=rt_tol)
            if check:
                matched_indices.append(i)
            pass

        # If one match is found, use it
        if len(matched_indices) == 1:
            matchedindex = matched_indices[0]
            self.masses[matchedindex].merge_in_pk(pk)
        elif len(matched_indices) > 1:
            # The new peak matches to multiple, select the one with the closest RT
            rt_diffs = np.abs([pk.rt - self.masses[i].apexrt for i in matched_indices])
            closest = np.argmin(rt_diffs)
            matchedindex = matched_indices[closest]
            self.masses[matchedindex].merge_in_pk(pk)
        else:
            # If no matches are found, create a new mass
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

    def merge_missed_monoisotopics(self, ppm_tolerance=20, max_mm=1, mass_diff_c=1.0033):
        to_remove = []
        for i1 in range(len(self.peaks)):
            for i2 in range(len(self.peaks)):
                if i1 == i2 or self.peaks[i1].z != self.peaks[i2].z:
                    continue
                if abs(self.peaks[i1].monoiso - self.peaks[i2].monoiso) < max_mm * mass_diff_c * 1.1:
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
                            # Remove the mass at i2 and add unmatched monoisos to the mass at i1
                            for mono in self.masses[i2].monoisos:
                                if any(ud.within_ppm(x, mono, ppm_tolerance) for x in self.masses[i1].monoisos):
                                    continue
                                else:
                                    self.masses[i1].monoisos = np.append(self.masses[i1].monoisos, mono)
                            # Add charge states from i2 that don't exist in i1
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
        ms2_scan_dict = msalign.sort_by_scan_order(self, 2)
        ms2_features = 0
        for k, v in ms2_scan_dict.items():
            ms2_features += len(v)
        msalign.write_ms1_msalign(ms1_scan_dict, ms2_scan_dict, filename, config)
        msalign.write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, filename, config,
                                  max_precursors=max_precursors, act_type=act_type,
                                  report_multiple_monoisos=config.report_multiple_monoisos,
                                  write_noprec_scans=config.write_scans_without_precs)

    def to_df(self, avg=False, report_multiple_monoisos=True):
        """
        Convert a MatchedCollection object to a pandas dataframe
        :return: Pandas dataframe
        """
        col = "Average Mass" if avg else "Monoisotopic Mass"

        data = []
        for p in self.peaks:
            try:
                avg_mono_delta = p.avgmass - p.monoiso
            except:
                avg_mono_delta = 0
            if report_multiple_monoisos:
                for mono in p.monoisos:
                    d = {"Charge": p.z, "Most Abundant m/z": p.mz, col: mono + avg_mono_delta if avg else mono,
                         "Scan": p.scan,
                         "Most Abundant Mass": p.peakmass, "Abundance": p.matchedintensity}
                    data.append(d)
            else:
                d = {"Charge": p.z, "Most Abundant m/z": p.mz, col: p.avgmass if avg else p.monoiso, "Scan": p.scan,
                     "Most Abundant Mass": p.peakmass, "Abundance": p.matchedintensity}
                data.append(d)
        df = pd.DataFrame(data)
        return df

    def export_tsv(self, filename="export.tsv", avg=False,
                   report_multiple_monoisos=True):
        print("Exporting to", filename)
        df = self.to_df(avg, report_multiple_monoisos=report_multiple_monoisos)
        df.to_csv(filename, sep="\t", index=False)

    def to_mass_spectrum(self, binsize=0.1):
        if len(self.masses) == 0:
            print("Zero Sized Peak Array")
            return np.array([])
        minmass = np.min([p.monoiso for p in self.masses]) - 10 * binsize
        maxmass = np.max([p.monoiso for p in self.masses]) + 100 * binsize

        # Round minmass and maxmass to nearest binsize
        minmass = np.floor(minmass / binsize) * binsize
        maxmass = np.ceil(maxmass / binsize) * binsize

        masses = np.arange(minmass, maxmass, binsize)

        intensities = np.zeros(len(masses))
        for p in self.masses:
            idx = round((p.monoiso - minmass) / binsize)
            intensities[idx] += p.apexintensity

        return np.transpose([masses, intensities])

    def to_merged_isodists(self):
        isodistall = []
        for p in self.peaks:
            if p.isodist is not None:
                for d in p.isodist:
                    isodistall.append(d)
        # to numpy array
        isodistall = np.array(isodistall)
        if len(isodistall) == 0:
            print("Zero Sized Isodist Array")
            return np.array([])
        # Sort by mass
        isodistall = isodistall[np.argsort(isodistall[:, 0])]
        return isodistall

    def fill_fwhm(self, profiledata, sort=False, wfactor=4, maxppmtol=1000):
        for p in self.peaks:
            p.set_fwhms(profiledata, sort=sort, wfactor=wfactor, maxppmtol=maxppmtol)

    def gen_peaks(self, profiledata):
        newpeaks = []
        for p in self.peaks:
            newpeak = p.gen_peak(profiledata)
            if newpeak is not None:
                newpeaks.append(newpeak)
            if np.isnan(newpeak).any():
                print("NaN in simdata")
                print(p)
        self.simdata = np.array(newpeaks)
        self.sumsim = self.sum_peaks_data(profiledata)
        huskleft = self.remove_peaks_from_data(profiledata)
        return newpeaks, self.sumsim, huskleft

    def sum_peaks_data(self, data):
        firstsum = np.sum(self.simdata, axis=0)
        sum = intensity_decon(data, self.simdata, n=10)
        error1 = np.sum((data[:, 1] - firstsum) ** 2) / np.sum(data[:, 1] ** 2)
        error2 = np.sum((data[:, 1] - sum) ** 2) / np.sum(data[:, 1] ** 2)
        print("Sum of Squares Error: Simple Sum:", error1, "Deconvoluted Sum:", error2)
        return sum

    def remove_peaks_from_data(self, data, percentcutoff=0.25):
        huskleft = data[:, 1] - self.sumsim
        huskleft[huskleft < 0] = 0

        if percentcutoff > 0:
            pointerror = ud.safedivide(huskleft, data[:, 1])
            huskleft[pointerror < percentcutoff] = 0

        return huskleft


class MatchedMass:
    """
    Matched mass object for collecting data on MatchedPeaks with matched masses.
    """

    def __init__(self, pk):
        self.monoiso = pk.monoiso
        self.monoisos = pk.monoisos
        self.scans = np.array([pk.scan])
        self.apexintensity = pk.peakint
        self.maxscan = pk.scan
        self.maxrt = pk.rt
        self.minscan = pk.scan
        self.minrt = pk.rt
        self.apexscan = pk.scan
        self.apexrt = pk.rt
        self.mzs = np.array([pk.mz])
        self.zs = np.array([pk.z])
        self.totalintensity = pk.matchedintensity
        self.mzints = np.array([pk.peakint])
        self.scan_intensities = {pk.scan: pk.matchedintensity}
        self.isodists = pk.isodist
        self.totalpeaks = 1
        self.avgmass = pk.avgmass
        self.decon_centroids = pk.calc_mass_dists()
        self.massdist = pk.massdist

        # List of MatchedPeak objects
        self.clusters = [pk]

    def __str__(self):
        return f"Mass: {self.monoiso:.4f}, Charge States: {self.zs}, Total Intensity: {self.totalintensity}, Total Peaks: {self.totalpeaks}"

    def merge_in_pk(self, newpk):
        # Add to list of peaks
        self.clusters.append(newpk)

        # Update Apex intensity and scans
        if newpk.peakint is not None:
            if newpk.peakint > self.apexintensity:
                self.apexintensity = newpk.peakint
                self.apexscan = newpk.scan
                self.apexrt = newpk.rt

        # Update min and max scans if necessary
        if newpk.scan > self.maxscan:
            self.maxscan = newpk.scan
            self.maxrt = newpk.rt
        if newpk.scan < self.minscan:
            self.minscan = newpk.scan
            self.minrt = newpk.rt

        # Update mz, z, mzints, isodists if new z state
        if not np.isin(self.zs, newpk.z).any():
            self.zs = np.append(self.zs, newpk.z)
            self.mzs = np.append(self.mzs, newpk.mz)
            self.mzints = np.append(self.mzints, newpk.peakint)
            self.isodists = np.vstack([self.isodists, newpk.isodist])

        # Update self.decon_centroids
        newdist = newpk.calc_mass_dists()
        if newdist is not None:
            if self.decon_centroids is not None:
                self.decon_centroids = np.vstack([self.decon_centroids, newdist])
            else:
                self.decon_centroids = newdist
        self.decon_centroids = merge_decon_centroids(self.decon_centroids, ppm_tol=newpk.config.matchtol)

        # Update monoisotopic masses if necessary
        for m in newpk.monoisos:
            closeindex = fastnearest(self.monoisos, m)
            # Check if within ppm tolerance of the closest monoisotopic mass
            if ud.within_ppm(m, self.monoisos[closeindex], newpk.config.matchtol):
                intensity = newpk.matchedintensity if newpk.matchedintensity is not None else 1
                # Average the monoisotopic masses weighted by intensity of the cluster
                self.monoisos[closeindex] = self.monoisos[closeindex] * self.totalintensity + m * intensity / (
                        self.totalintensity + intensity)
            else:
                self.monoisos = np.append(self.monoisos, m)

        # Update monoisotopic mass
        if ud.within_ppm(self.monoiso, newpk.monoiso, newpk.config.matchtol):
            # If monoisotopic mass is exactly the same, average weighted by intensity
            intensity = newpk.matchedintensity if newpk.matchedintensity is not None else 1
            self.monoiso = (self.monoiso * self.totalintensity + newpk.monoiso * intensity) / (
                    self.totalintensity + intensity)
        else:
            # If monoisotopics masses don't match, compare each with the deconvolved centroids and keep the one with the best CSS
            # There's a faster way to do this where you don't match each, but I'm not up for it today
            # print("Warning: Merging peaks with different monoisotopic masses:", self.monoiso, newpk.monoiso)
            css_new = calc_css_from_data(self.decon_centroids, newpk.massdist, ppm_tol=newpk.config.matchtol)
            css_old = calc_css_from_data(self.decon_centroids, self.massdist, ppm_tol=newpk.config.matchtol)
            if css_new > css_old:
                self.monoiso = newpk.monoiso
                self.massdist = newpk.massdist
            # else:
            #     print("Warning: Merging peaks with different monoisotopic masses:", self.monoiso, newpk.monoiso,)
            # Do nothing if peak is not better

        # Fit Massdist to decon_centroids, this is somewhat crude, but I think it's probably good enough
        self.massdist = merge_massdist(np.array(self.massdist), self.decon_centroids, newpk.config.matchtol)

        # Add to total intensity
        if newpk.scan not in self.scans:
            self.scans = np.append(self.scans, newpk.scan)
        if newpk.matchedintensity is not None:
            # Add the intensity to the scan intensity dictionary
            if self.scan_intensities.get(newpk.scan) is not None:
                self.scan_intensities[newpk.scan] += newpk.matchedintensity
            else:
                self.scan_intensities[newpk.scan] = newpk.matchedintensity
            self.totalintensity += newpk.matchedintensity
        else:
            self.totalintensity += 1
        self.totalpeaks += 1

    def check_if_match(self, pk, config, scan_tol=100, rt_tol=None):
        nearest_mass = self.monoiso
        # Close check to see if the masses are anywhere near each other
        close_check = abs(pk.monoiso - nearest_mass) <= config.maxshift * config.mass_diff_c * 1.1
        if not close_check:
            return False

        # Check if the scan is close enough
        scan_check = abs(pk.scan - self.scans[len(self.scans) - 1]) <= scan_tol
        # If RT tolerance is set, check if within RT tolerance
        if rt_tol is not None and self.maxrt > 0 and self.minrt > 0 and pk.rt > 0:
            rt_check = abs(pk.rt - self.maxrt) <= rt_tol or abs(pk.rt - self.minrt) <= rt_tol
            scan_check = scan_check and rt_check
        if not scan_check:
            return False

        # Check if within tolerance with monoisotopic shifts
        ppm_check = within_ppm_plus_mm(nearest_mass, pk.monoiso, config.matchtol, config.maxshift, config.mass_diff_c)
        if not ppm_check:
            return False

        # Check if the isodist matches on each
        isodist_check = isodist_match(self.decon_centroids, pk.massdist, css_threshold=config.css_thresh,
                                      ppm_tol=config.matchtol)
        if not isodist_check:
            # print("Isodist check failed for peak", pk.monoiso, "and mass", self.monoiso)
            return False

        return True


class MatchedPeak:
    """
    Matched peak object for collecting data on peaks with matched distributions
    """

    def __init__(self, z, mz, avgmass=None, centroids=None, isodist=None, matchedindexes=None, isomatches=None,
                 config=None):
        self.mz = mz
        self.z = z
        self.centroids = centroids
        self.decon_centroids = None

        self.isodist = isodist
        self.massdist = None

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

        self.monoiso = -1
        self.monoisos = []

        self.peakmass = -1
        self.startindex = -1
        self.endindex = -1

        self.bestshift = -1
        self.avgmass = avgmass
        self.acceptedshifts = None
        self.matchedions = []

        self.fwhms = []
        self.avgfwhm = 0
        self.fwhmratio = 0
        self.avgsigma = 0

        if matchedindexes is not None:
            self.matchedindexes = matchedindexes
            self.matchedcentroids = centroids[np.array(matchedindexes)]
            yvals = np.array([c[1] for c in centroids[np.array(matchedindexes)]])
            self.matchedintensity = np.sum(yvals)
        if isomatches is not None:
            self.isomatches = isomatches
            self.matchedisodist = isodist[np.array(isomatches)]

        # Keep a copy of config here solely for calc_mass_dists
        self.config = config
        if self.config is None:
            self.config = IsoDecConfig()

    def calc_mass_dists(self):
        if self.centroids is not None:
            self.decon_centroids = np.column_stack(
                (self.centroids[:, 0] * self.z - self.config.adductmass * self.z, self.centroids[:, 1]))
            return self.decon_centroids
        else:
            print("No centroids")
            return None

    def strip_from_data(self, data):
        """
        Strips the matched peak from the data
        :param data: The data to strip from
        :return: None
        """
        minmz = np.amin(self.isodist[:, 0]) - 0.5
        maxmz = np.amax(self.isodist[:, 0]) + 0.5
        data = ud.dataremove(data, minmz, maxmz)
        return data

    def set_fwhms(self, profiledata, sort=False, wfactor=4, maxppmtol=1000):
        if ud.isempty(self.isodist):
            self.fwhms = np.array([])
        else:
            self.fwhms = fast_fwhm(profiledata, self.isodist, sort=sort, wfactor=wfactor, maxppmtol=maxppmtol)

        goodb = self.fwhms[:, 0] > 0
        goodfwhms = self.fwhms[goodb]
        if len(goodfwhms) > 0:
            self.avgfwhm = np.mean(goodfwhms[:, 0])
            self.fwhmratio = np.mean(goodfwhms[:, 1] / goodfwhms[:, 2])
            self.avgsigma = np.mean(goodfwhms[:, 5])

    def gen_peak(self, data):
        if self.isodist is None or len(self.isodist) == 0:
            print("No isodist for peak", self.mz, self.z)
            return None
        if len(self.fwhms) == 0:
            print("No FWHM for peak", self.mz, self.z)
            return None
        if self.avgfwhm == 0 or np.isnan(self.avgfwhm):
            print("No valid FWHM for peak", self.mz, self.z)
            return None

        outdata = np.zeros_like(data)
        for i, fwhm in enumerate(self.fwhms):
            if fwhm[0] <= 0:
                continue
            # Generate the peak using the isodist and FWHM
            peakmz = self.isodist[i, 0]
            peakint = self.isodist[i, 1]
            # fwhm_value = fwhm[0]
            sigma = fwhm[5]
            startindex = int(fwhm[3])
            endindex = int(fwhm[4])
            centroid = fwhm[6]

            if sigma <= 0:
                sigma = self.avgsigma

            xvalues = data[startindex:endindex + 1, 0]
            yvalues = ndis_std(xvalues, centroid, sigma, peakint)

            outdata[startindex:endindex + 1, 1] += yvalues

        return outdata[:, 1]

    def __str__(self):
        return f"MatchedPeak: mz={self.mz}, z={self.z}, monoiso={self.monoiso}"

@njit(fastmath=True)
def within_ppm_plus_mm(mass1, mass2, ppm_tol=20, max_mm=1, mass_diff_c=1.0033):
    for mm in range(-max_mm, max_mm + 1):
        adjusted_mass2 = mass2 + mm * mass_diff_c
        if ud.within_ppm(mass1, adjusted_mass2, ppm_tol):
            return True
    return False

@njit(fastmath=True)
def merge_massdist(massdist: np.ndarray, decon_centroids: np.ndarray, matchtol: float):
    max_intensity = 0.0
    diff_at_max = 0.0
    normfactor = 1.0
    for d in massdist:
        if d[1] > max_intensity:
            # Check if it matches anything in self.decon_centroids
            idx = fastnearest(decon_centroids[:, 0], d[0])
            if ud.within_ppm(decon_centroids[idx, 0], d[0], matchtol):
                max_intensity = d[1]
                diff_at_max = abs(decon_centroids[idx, 0] - d[0])
                normfactor = decon_centroids[idx, 1] / d[1] if d[1] != 0 else 1
    if diff_at_max < matchtol * massdist[0,0] / 1e6:
        massdist[:, 0] += diff_at_max
    if max_intensity != 0:
        massdist[:, 1] = massdist[:, 1] * normfactor

    return massdist


@njit(fastmath=True)
def merge_decon_centroids(centroids, ppm_tol=10):
    if centroids is None or len(centroids) == 0:
        return None

    centroids = centroids[np.argsort(centroids[:, 0])]
    merged = np.zeros((len(centroids), 2))
    m_idx = 0

    current_mass = centroids[0, 0]
    current_intensity = centroids[0, 1]
    weighted_sum = current_mass * current_intensity
    total_intensity = current_intensity

    for i in range(1, len(centroids)):
        next_mass = centroids[i, 0]
        next_intensity = centroids[i, 1]
        ppm = abs(current_mass - next_mass) / current_mass * 1e6
        if ppm <= ppm_tol:
            weighted_sum += next_mass * next_intensity
            total_intensity += next_intensity
            current_mass = weighted_sum / total_intensity
            current_intensity = total_intensity
        else:
            merged[m_idx, 0] = current_mass
            merged[m_idx, 1] = current_intensity
            m_idx += 1
            current_mass = next_mass
            current_intensity = next_intensity
            weighted_sum = current_mass * current_intensity
            total_intensity = current_intensity

    merged[m_idx, 0] = current_mass
    merged[m_idx, 1] = current_intensity
    return merged[:m_idx+1]

@njit(fastmath=True)
def isodist_match(isodists1, isodist2, css_threshold=0.7, ppm_tol=20):
    if isodists1 is None or isodist2 is None:
        return False
    if len(isodists1) == 0 or len(isodist2) == 0:
        return False

    # Calc cosine similarity
    css = calc_css_from_data(isodists1, isodist2, ppm_tol=ppm_tol)
    # print("css=", css)
    if css < css_threshold:
        return False
    return True


def remove_multiple_monoisos(mc):
    # Find the objects with the same everything except for monoiso
    mcnew = MatchedCollection()
    duplicate_indices = []
    for i, pk in enumerate(mc.peaks):
        pklist = [pk.scan, pk.peakmass, pk.matchedion, pk.matchedintensity, pk.z]
        if i in duplicate_indices:
            continue  # Skip if this peak has already been identified as a duplicate
        for j, pk2 in enumerate(mc.peaks[i + 1:]):
            pklist2 = [pk2.scan, pk2.peakmass, pk2.matchedion, pk2.matchedintensity, pk2.z]
            # check each list to see if they are the same
            if pklist == pklist2 and pk.monoiso != pk2.monoiso:
                # Add the monoiso to the list of monoisos
                if pk2.monoiso not in pk.monoisos:
                    pk.monoisos = np.append(pk.monoisos, pk2.monoiso)
                # print("Found multiple monoisos for the same peak:", pk.monoisos, duplicate_indices[-3:])
                duplicate_indices.append(j + i + 1)
            else:
                # Add the peak to the new collection
                mcnew.add_peak(pk)
                break

    return mcnew


# @njit(fastmath=True)
def create_isodist(peakmz, charge, data, adductmass=1.007276467):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    if charge == 0:
        mass = peakmz
    else:
        mass = (peakmz - adductmass) * charge
    isodist = fast_calc_averagine_isotope_dist(mass, charge=charge, adductmass=adductmass)
    isodist[:, 1] *= np.amax(data[:, 1])
    # shift isodist so that maxes are aligned with data
    mzshift = peakmz - isodist[np.argmax(isodist[:, 1]), 0]
    isodist[:, 0] = isodist[:, 0] + mzshift
    return isodist


# @njit(fastmath=True)
def create_isodist_full(peakmz, charge, data, adductmass=1.007276467, isotopethresh: float = 0.01):
    """
    Create an isotopic distribution based on the peak m/z and charge state.
    :param peakmz: Peak m/z value as float
    :param charge: Charge state
    :param data: Data to match to as 2D numpy array [m/z, intensity]
    :return: Isotopic distribution as 2D numpy array [m/z, intensity]
    """
    charge = float(charge)
    if charge == 0:
        mass = peakmz
    else:
        mass = (peakmz - adductmass) * charge

    isodist, massdist = fast_calc_averagine_isotope_dist_dualoutput(mass, charge=charge, adductmass=adductmass,
                                                                    isotopethresh=isotopethresh)
    isodist[:, 1] *= np.amax(data[:, 1])
    massdist[:, 1] *= np.amax(data[:, 1])
    # shift isodist so that maxes are aligned with data
    mzshift = peakmz - isodist[np.argmax(isodist[:, 1]), 0]
    isodist[:, 0] = isodist[:, 0] + mzshift

    if charge == 0:
        massshift = mzshift
    else:
        massshift = mzshift * charge
    monoiso = mass + massshift
    massdist[:, 0] = massdist[:, 0] + massshift
    return isodist, massdist, monoiso


# @njit(fastmath=True)
def get_accepted_shifts(centroids, z, peakmz, config):
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

    shiftrange = np.arange(-maxshift, maxshift + 1)
    # shifts_scores = [[0, 0] for i in range(len(shiftrange))]
    shifts_scores = np.zeros((len(shiftrange), 2))
    # overlaps = []
    bestshift = -1
    sum = -1
    meanratio = 1
    for i, shift in enumerate(shiftrange):
        s = calculate_cosinesimilarity(cent_intensities, isodist[:, 1], shift, maxshift,
                                       minusoneaszero=config.minusoneaszero)
        shifts_scores[i, 0] = shift
        shifts_scores[i, 1] = s
        if s > sum:
            sum = s
            bestshift = shift
    max_score = np.amax(shifts_scores[:, 1])
    if max_score < config.css_thresh:
        return None, None, None, None

    accepted_shifts = [x for x in shifts_scores if x[1] > max_score - config.min_score_diff]

    return accepted_shifts, monoiso, massdist, isodist


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
    accepted_shifts, monoiso, massdist, isodist = get_accepted_shifts(centroids, z, peakmz, config)

    if config.verbose:
        print("Accepted Shifts:", accepted_shifts)
    if accepted_shifts is None:
        return None

    # The plan is to report a set of possible monoisotopics, but only include the centroids, isodist, massdist, and matches for the best one.
    m = MatchedPeak(int(z), float(peakmz), centroids=centroids, isodist=isodist, config=config)

    bestshift = accepted_shifts[np.argmax([x[1] for x in accepted_shifts])]
    m.bestshift = bestshift
    m.acceptedshifts = np.array(accepted_shifts)

    peakmz_new, isodist_new, matchedindexes, isomatches, monoiso_new, massdist_new = make_shifted_peak(
        bestshift[0], bestshift[1], monoiso, massdist, isodist, peakmz, z, centroids, config.matchtol, config.minpeaks,
        config.plusoneintwindowlb, config.plusoneintwindowub, config.css_thresh,
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
                               minusoneaszero: bool = True):
    ab = 0
    a2 = 0
    b2 = 0

    if minusoneaszero:
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
    if z == 0:
        z = 1

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
    most_intense_iso = (int)(0.0006 * peakmass + 0.4074)
    return peakmass - (most_intense_iso * 1.0033)

@njit(fastmath=True)
def calc_css_from_data(centroids, isodist, ppm_tol=5):
    """
    Check if the data matches the isotopic distribution.
    :param centroids: Centroid data as numpy array [m/z, intensity]
    :param isodist: Isotopic distribution as numpy array [m/z, intensity]
    :param ccsthresh: Cosine similarity threshold for a match
    :return: True if the data matches the isotopic distribution, False otherwise
    """
    isodist = isodist.copy()

    cent_intensities = find_matched_intensities(centroids[:, 0], centroids[:, 1], isodist[:, 0], 0,
                                                tolerance=ppm_tol, z=1, peakmz=isodist[0, 0])

    norm_factor = max(cent_intensities) / max(isodist[:, 1])
    isodist[:, 1] *= norm_factor

    css = calculate_cosinesimilarity(cent_intensities, isodist[:, 1], 0, 0)
    return css


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

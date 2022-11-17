import numpy as np
import unidec_modules.unidectools as ud
import molmass
from copy import deepcopy
import pandas as pd
import os
import time
import matplotlib.pyplot as plt
import unidec_modules.peakstructure as ps

pd.set_option('mode.chained_assignment', None)


def matcher(array, target, tolerance):
    index = ud.nearest(array, target)
    match = array[index]
    if np.abs(match - target) < tolerance * 1e-6 * target:
        return index
    else:
        return -1


class LipiDecRunner:
    def __init__(self, datapath, libpath, dir, datarange=None):
        self.datapath = datapath
        self.libpath = libpath
        self.directory = dir
        os.chdir(self.directory)
        print("Data File:", self.datapath)
        print("Library File:", self.libpath)

        self.eng = None
        self.results = None

        libext = os.path.splitext(libpath)[1]
        if libext == ".xlsx":
            libdf = pd.read_excel(libpath)
        elif libext == ".npz":
            npz = np.load(libpath, allow_pickle=True)
            libdf = pd.DataFrame.from_dict({item: npz[item] for item in npz.files}, orient='columns')
            try:
                libdf = libdf.drop(columns=["site", "fa-chain"])
            except:
                pass
        self.libdf = libdf

        fastfile = self.datapath[:-4] + ".npz"
        self.polarity = ud.get_polarity(self.datapath)
        if not os.path.isfile(fastfile):
            data = ud.load_mz_file(self.datapath)
            np.savez(fastfile, data=data)
        else:
            data = np.load(fastfile, allow_pickle=True)["data"]
        self.data = data

        if datarange is not None:
            self.data = ud.datachop(self.data, datarange[0], datarange[1])
            # data = ud.datachop(data, 742, 772)
            # data = ud.datachop(data, 870, 876)

        self.eng = LipiDec(self.data, self.libdf)

    def run(self, ignore_list=[], outdir=None):
        st = time.perf_counter()
        self.eng.cleanup_db(mode=self.polarity, ignore_list=ignore_list)
        self.eng.find_matches()
        self.eng.run_decon()
        self.results = self.eng.find_alternates()
        if outdir is None:
            self.results.to_excel(self.datapath[:-4] + "_results.xlsx")
        else:
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            self.results.to_excel(os.path.join(outdir, self.datapath[:-4] + "_results.xlsx"))
        # results = calc_fromula_from_mass(peaks[:,0], tolerance=5)
        print("Done:", time.perf_counter() - st)

    def make_plot(self, show=True):
        plt.figure()
        plt.plot(self.data[:, 0], self.data[:, 1])
        self.eng.plot_hits(plt.gca())
        plt.hlines(self.eng.noiselevel, np.amin(self.data[:, 0]), np.amax(self.data[:, 0]))
        plt.hlines(self.eng.noiselevel * 5, np.amin(self.data[:, 0]), np.amax(self.data[:, 0]))
        plt.xlim(np.amin(self.data[:, 0]), np.amax(self.data[:, 0]))
        plt.ylim(self.eng.noiselevel, np.amax(self.data[:, 1]))
        plt.gca().set_yscale('symlog')
        if show:
            plt.show()


class LipiDec:
    def __init__(self, data, topdf):
        self.topdf = topdf
        self.df = self.topdf
        self.hdf = None
        self.peaks = None
        self.data = data
        self.tolerance = 3
        self.peaktol = self.tolerance
        # self.sort_df()
        self.minisomatches = 2
        self.intensity_error_cutoff = 0.2
        self.coal_tol_mult = 5
        self.coal_rat_cutoff = 10
        self.db_round_decimal = 3


    def sort_df(self):
        self.df.sort_values("Mass")

    def cleanup_db(self, mode="positive", ignore_list=[], drop_duplicates=True):
        if mode is None:
            mode = "positive"
        if mode.lower() == "positive":
            print("Positive Mode")
            b1 = self.topdf["charge"].to_numpy() > 0
            self.topdf = self.topdf.loc[b1]
        else:
            print("Negative Mode")
            b1 = self.topdf["charge"].to_numpy() < 0
            self.topdf = self.topdf.loc[b1]

        for ignore in ignore_list:
            self.topdf = self.topdf.drop(self.topdf[self.topdf["subclassKey"] == ignore].index)

        self.topdf = self.topdf.reset_index()

        if drop_duplicates:
            self.unique_masses, self.unique_indexes, self.unique_inverse = np.unique(
                np.round(self.topdf["Mz"], decimals=self.db_round_decimal), return_index=True, return_inverse=True)
            print("Cutting Duplicate Masses:", len(self.topdf), len(self.unique_masses))
            self.df = self.topdf.loc[self.unique_indexes]
            self.df = self.df.reset_index()
        else:
            self.df = self.topdf

    def find_matches(self, tolerance=None):
        if tolerance is not None:
            self.tolerance = tolerance
        self.find_initial_matches()
        self.correct_peaks()
        self.find_isotope_matches(minmatches=self.minisomatches)

    def get_peaks(self):
        self.noiselevel = ud.noise_level2(self.data[self.data[:, 1] > 0], percent=0.50)
        self.peaks = ud.peakdetect(self.data, threshold=self.noiselevel * 5, ppm=self.peaktol, norm=False)
        print("Peaks:", len(self.peaks))
        # Change peak apex to centroids
        self.pks = ps.Peaks()
        self.pks.add_peaks(self.peaks)
        ud.peaks_error_FWHM(self.pks, self.data)
        #self.peaks[:, 0] = self.pks.centroids
        self.pks.integrate(self.data)

    def find_initial_matches(self):
        self.get_peaks()
        dbmasses = self.df["Mz"].to_numpy()
        matches = np.array([matcher(self.peaks[:, 0], mass, self.tolerance) for mass in dbmasses])
        b1 = matches > -1
        self.hdf = deepcopy(self.df.loc[b1])
        self.hdf["Match"] = matches[b1]
        self.hdf["Height"] = self.peaks[matches[b1], 1]
        self.hdf["Errors"] = self.peaks[matches[b1], 0] - self.hdf["Mz"]

    def correct_peaks(self):
        median_error = np.median(self.hdf["Errors"])
        self.peaks[:, 0] -= median_error
        self.hdf["Errors"] = self.peaks[self.hdf["Match"], 0] - self.hdf["Mz"]
        median_error2 = np.median(self.hdf["Errors"])
        print("Corrected Peaks:", median_error, median_error2)

    def detect_coalescence(self, isomasses, isoints, isomatches, row, index=2):
        newmatch = matcher(self.peaks[:, 0], isomasses[index], self.tolerance * self.coal_tol_mult)
        if newmatch > -1:
            predicted_int = row["Height"] * isoints[index]
            measured_int = self.peaks[newmatch, 1]
            ratio = measured_int / predicted_int
            # print(isomasses[index], self.peaks[newmatch], predicted_int, ratio)
            if ratio > self.coal_rat_cutoff:
                isomatches[index] = newmatch
                isomasses[index] = self.peaks[newmatch, 0]
        return isomatches, isomasses

    def find_isotope_matches(self, minmatches=2):
        isomassarray = []
        isointarray = []
        isomatcharray = []
        isohitarray = []
        for i, row in self.hdf.iterrows():
            formula = row["Formula"]
            adductmass = row["AdductMass"]
            charge = row["charge"]
            formula = formula
            f = molmass.Formula(formula)
            isotopes = np.array(f.spectrum(min_intensity=0.1).dataframe())
            isomasses = (isotopes[:, 0] + adductmass) / np.abs(charge)
            isoints = isotopes[:, 2] / 100.
            isomatches = np.array([matcher(self.peaks[:, 0], m, self.tolerance) for m in isomasses])
            hits = np.sum(isomatches > -1)

            if hits >= minmatches:
                for j in range(2, len(isomatches)):
                    if isomatches[j] == -1:
                        isomatches, isomasses = self.detect_coalescence(isomasses, isoints, isomatches, row, index=j)

                b1 = isomatches > -1
                isomatcharray.append(isomatches[b1])
                isomassarray.append(isomasses[b1])
                isointarray.append(isoints[b1])
                isohitarray.append(True)
            else:
                isohitarray.append(False)

        self.hdf = self.hdf[isohitarray]
        self.hdf["Isomasses"] = isomassarray
        self.hdf["Isoints"] = isointarray
        self.hdf["Isomatch"] = isomatcharray

    def run_decon(self):
        for i in range(0, 10):
            spectrum = self.make_spectrum()
            ratio = ud.safedivide(self.peaks[:, 1], spectrum)
            self.apply_ratio(ratio)

        self.find_good_species()
        self.find_missing_peaks()
        self.find_alternates()

    def make_spectrum(self):
        spectrum = np.zeros_like(self.peaks[:, 1])
        for i, row in self.hdf.iterrows():
            isomatches = np.array(row["Isomatch"]).astype(int)
            isoints = row["Isoints"] * row["Height"]
            spectrum[isomatches] = isoints + spectrum[isomatches]
        error = np.log(np.sum((self.peaks[:, 1] - spectrum) ** 2))
        print("Error", error)
        self.diffs = np.abs(self.peaks[:, 1] - spectrum) / np.array(self.peaks[:, 1])
        self.spectrum = spectrum
        return spectrum

    def apply_ratio(self, ratio):
        heights = []
        for i, row in self.hdf.iterrows():
            isomatches = row["Isomatch"]
            isoints = row["Isoints"] * row["Height"]
            avgrat = np.average(np.sqrt(ratio[isomatches]), weights=isoints)
            # avgrat = np.average(ratio[isomatches], weights=isoints)
            heights.append(row["Height"] * avgrat)
        self.hdf["Height"] = heights

    def find_missing_peaks(self):
        spectrum = self.make_spectrum()
        missing_peaks = []
        matched_peaks = []
        good_peaks = []
        bad_peaks = []
        for i, p in enumerate(self.peaks):
            if spectrum[i] == 0:
                missing_peaks.append(p)
            else:
                matched_peaks.append(p)
            if self.diffs[i] > self.intensity_error_cutoff:
                bad_peaks.append(p)
            else:
                good_peaks.append(p)

        print("Percentage of Missing Peaks:", 100 * len(missing_peaks) / (len(self.peaks)))
        print("Percentage of Good Peaks:", 100 * len(good_peaks) / (len(self.peaks)))
        return missing_peaks, matched_peaks

    def find_good_species(self):
        self.make_spectrum()
        errors = []
        good_things = []
        bad_things = []
        for i, row in self.hdf.iterrows():
            isomatches = row["Isomatch"][:2]
            isoints = row["Isoints"][:2]
            avgerr = np.average(self.diffs[isomatches], weights=isoints ** 2)
            errors.append(avgerr)
            if avgerr < self.intensity_error_cutoff:
                good_things.append(row)
            else:
                bad_things.append(row)
        self.hdf["IntErrors"] = errors
        print("Percentage of good matches: ", 100 * len(good_things) / len(self.hdf))
        self.hdf = self.hdf[self.hdf["IntErrors"] < self.intensity_error_cutoff]
        # print(self.hdf[["SumComp", "subclassKey", "Errors" ]])

    def plot_hits(self, ax):
        spectrum = self.make_spectrum()

        for i, p in enumerate(self.peaks):
            if spectrum[i] == 0:
                ax.plot(p[0], p[1], marker="o", linestyle=" ", color="r")
            elif self.diffs[i] > self.intensity_error_cutoff:
                ax.plot(p[0], p[1], marker="o", linestyle=" ", color="y")
            else:
                # ax.plot(p[0], spectrum[i], marker="o", linestyle=" ", color="y")
                pass

        for i, row in self.hdf.iterrows():
            name = row["subclassKey"] + " " + row["SumComp"]
            isomasses = row["Isomasses"]
            isoints = row["Isoints"] * row["Height"]

            ax.plot(isomasses, isoints, marker="o", color="g", linestyle="--")
            ax.text(isomasses[0], isoints[0] * 1.05, str(name))

    def find_alternates(self):
        mzvals = np.round(self.topdf["Mz"], decimals=self.db_round_decimal)
        self.resultsdf = pd.DataFrame()
        for i, row in self.hdf.iterrows():
            mz = np.round(row["Mz"], decimals=self.db_round_decimal)
            b1 = mzvals == mz
            row["PeakNumber"] = i
            newdf = self.topdf[b1]
            for column in row.keys():
                if column not in newdf.keys():
                    try:
                        newdf[column] = row[column]
                    except:
                        newdf[column] = str(row[column])

            self.resultsdf = pd.concat([self.resultsdf, newdf])
        return self.resultsdf.sort_values("Height", ascending=False)

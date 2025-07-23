import os
import time
import numpy as np
from pathlib import Path
import sys

path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

from unidec.IsoDec.datatools import get_all_centroids, check_spacings, remove_noise_cdata, datacompsub
from unidec.IsoDec.match import MatchedCollection
from unidec.modules.unidecstructure import IsoDecConfig

from copy import deepcopy
from unidec.IsoDec.c_interface import IsoDecWrapper, example
from unidec.IsoDec.plots import plot_pks
from unidec.IsoDec.altdecon import thrash_predict
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import unidec.modules.fwhmtools as fwhmtools

class IsoDecRuntime:
    """
    Main class for IsoDec Engine
    """

    def __init__(self, phaseres=8, verbose=False):
        """
        Initialize the IsoDec Engine
        :param phaseres: Bit depth of the phase encoding. 8 is default.
        """
        self.config = IsoDecConfig()
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.test_centroids = []
        self.training_centroids = []
        self.config.verbose = verbose
        self.version = "1.0.0"
        self.pks = MatchedCollection()
        self.wrapper = IsoDecWrapper()
        self.config.phaseres = phaseres
        self.maxz = 50
        self.reader = None
        self.predmode = 0
        self.showavg = False
        self.analyte_type = None

    def phase_predictor(self, centroids):
        """
        Predict the charge of a peak
        :param centroids: Set of centroid data for a peak with m/z in first column and intensity in second
        :return: Charge state, integer
        """
        z = IsoDecWrapper().predict_charge(centroids)
        return [int(z), 0]

    def thrash_predictor(self, centroids):
        return [thrash_predict(centroids), 0]

    def batch_process_spectrum(self, data, window=5, type=None, threshold=0.0001, centroided=False, refresh=False):
        """
        Process a spectrum and identify the peaks. It first identifies peak cluster, then predicts the charge,
        then checks the peaks. If all is good, it adds them to the MatchedCollection as a MatchedPeak object.

        :param data: Spectrum data, m/z in first column, intensity in second
        :param window: Window for peak selection
        :param threshold: Threshold for peak selection
        :param centroided: Whether the data is already centroided. If not, it will centroid it.
        :return: MatchedCollection of peaks
        """

        if window is None:
            window = self.config.peakwindow
        if threshold is None:
            threshold = self.config.peakthresh

        if centroided:
            centroids = deepcopy(data)
        else:
            newdat = deepcopy(data)
            if self.config.background_subtraction > 0:
                newdat = datacompsub(newdat, self.config.background_subtraction)
            centroids = get_all_centroids(newdat, window=window, threshold=threshold * 0.1)

        med_spacing = check_spacings(centroids)
        if med_spacing <= self.config.meanpeakspacing_thresh:
            if self.config.verbose:
                print("Median Spacing:", med_spacing, "Removing noise.")
            centroids = remove_noise_cdata(centroids, 100, factor=1.5, mode="median")

        if refresh:
            self.pks = MatchedCollection()

        self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.config, type)

        return self.pks

    def pks_to_mass(self, binsize=0.1):
        """
        Convert the MatchedCollection to mass
        :return: None
        """
        return self.pks.to_mass_spectrum(binsize)

    def process_file(self, file, scans=None, check_centroided=True, assume_centroided=False):
        starttime = time.perf_counter()
        self.config.filepath = file
        # Get importer and check it
        reader = ImporterFactory.create_importer(file)
        self.reader = reader
        try:
            print("File:", file, "N Scans:", np.amax(reader.scans))
        except Exception as e:
            print("Could not open:", file)
            print(e)
            return []

        ext = file.split(".")[-1]
        t2 = time.perf_counter()
        if (ext == "raw" or ext == "RAW") and not os.path.isdir(file):
            isThermo = True
        else:
            isThermo = False
        # Loop over all scans
        for s in reader.scans:
            if scans is not None:
                if s not in scans:
                    continue

            # Open the scan and get the spectrum
            try:
                if isThermo:
                    spectrum = reader.grab_centroid_data(s)
                    reader.centroided = True
                else:
                    # mzml and mzxml will auto detect if it is centroided from pymzml
                    spectrum = reader.get_single_scan(s)
                    if check_centroided and not assume_centroided:
                        reader.check_centroided()
                    elif assume_centroided:
                        reader.centroided = True

            except Exception as e:
                print("Error Reading Scan", s, e)
                continue
            # If the spectrum is too short, skip it
            if len(spectrum) < 3:
                continue
            self.config.set_scan_info(s, reader)
            # print("Scan:", s, "Length:", len(spectrum), "Centroided:", reader.centroided)
            self.batch_process_spectrum(spectrum, centroided=reader.centroided)

            if s % 10 == 0:
                print("Scan:", s, "Length:", len(spectrum), "Avg. Time per scan:", (time.perf_counter() - t2) / 10.)
                t2 = time.perf_counter()

        print("Time:", time.perf_counter() - starttime)
        print("N Peaks:", len(self.pks.peaks))

        # self.pks.save_pks()
        return reader

    def export_peaks(self, type="prosightlite", filename=None, reader=None, act_type="HCD", max_precursors=1):
        print("Filename:", filename)
        if filename is None:
            filename = "peaks"

        if reader is None:
            reader = self.reader

        if type == "prosightlite":
            self.pks.export_prosightlite(filename)
        elif type == "msalign":
            if self.showavg:
                print("MSAlign not supported for Avg, defaulting to Monoisotopic")
            self.pks.export_msalign(self.config, reader, filename, act_type=act_type, max_precursors=max_precursors)
        elif type == "pkl":
            self.pks.save_pks()
        elif type == "tsv":
            self.pks.export_tsv(filename, self.showavg, self.config.report_multiple_monoisos)
        else:
            raise ValueError("Unknown Export Type", type)

    def remove_assigned_peaks(self, data):
        """
        Return a processed spectrum where the assigned peaks have been removed.
        """
        st = time.perf_counter()
        self.batch_process_spectrum(data)
        print("IsoDec time:", time.perf_counter() - st)
        isodists = self.pks.to_merged_isodists()
        if len(isodists) == 0:
            print("No peaks found, returning original data.")
            return data, isodists
        data = fwhmtools.remove_peaks(data, isodists, wfactor=6, maxppmtol=100)

        # for p in self.pks.peaks:
        #     # print("Removing:", p)
        #     data = p.strip_from_data(data)
        # data = ud.remove_middle_zeros(data)
        print("Removed Assigned Peaks, N Peaks:", len(self.pks.peaks), "Time:", time.perf_counter() - st)
        return data, isodists

    def batch_centroid(self, files, window=5, threshold=0.00001):
        """
        Batch centroiding of files.
        :param files: List of files to centroid
        :return: List of centroided data
        """
        outfiles = []
        for file in files:
            reader = ImporterFactory.create_importer(file)
            for s in reader.scans:
                spectrum = reader.get_single_scan(s)
                if len(spectrum) < 3:
                    continue
                centroids = get_all_centroids(spectrum, window=window, threshold=threshold)
                basepath = os.path.splitext(file)[0]
                outfile = f"{basepath}_{s}_centroided.txt"
                print("Centroiding file:", file, "Scan:", s, "N Centroids:", len(centroids))
                np.savetxt(outfile, centroids)
                outfiles.append(outfile)
        return outfiles




if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecRuntime()

    # files = ["C:\Data\Yuri\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3.csv",
    #          "C:\Data\Yuri\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms.csv",
    #          "C:\Data\Yuri\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms_NT3p5.csv"]
    # eng.batch_centroid(files)

    eng.analyte_type = "None"
    import platform
    if platform.system() == "Windows":
        file = "C:\\Python\\UniDecDev\\unidec\\bin\\TestSpectra\\test_2.txt"
    elif platform.system() == "Linux":
        file = "/mnt/c/Python/UniDecDev/unidec/bin/TestSpectra/test_2.txt"


    eng.process_file(file, assume_centroided=True)

    exit()

    c = example
    import matplotlib.pyplot as plt
    pks = eng.batch_process_spectrum(c, centroided=True)
    print("Time:", time.perf_counter() - starttime)
    exit()
    plot_pks(pks)

    plt.show()

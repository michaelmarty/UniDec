import os
import time
import numpy as np
from pathlib import Path
import sys



path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

from unidec.IsoDec.datatools import get_all_centroids, check_spacings, remove_noise_cdata
from unidec.IsoDec.match import IsoDecConfig, MatchedCollection

from copy import deepcopy
from unidec.IsoDec.c_interface import IsoDecWrapper, example
from unidec.IsoDec.plots import plot_pks
from unidec.IsoDec.altdecon import thrash_predict
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

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
        self.selection_type = None

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
            centroids = deepcopy(get_all_centroids(data, window=window, threshold=threshold * 0.1))

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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    starttime = time.perf_counter()
    eng = IsoDecRuntime(phaseres=8)
    eng.selection_type = "None"



    file ="C:\\Python\\UniDec3\\unidec\\bin\\TestSpectra\\test_2.txt"
    eng.process_file(file)

    exit()

    c = example
    pks = eng.batch_process_spectrum(c, centroided=True)
    print("Time:", time.perf_counter() - starttime)
    exit()
    plot_pks(pks)

    plt.show()

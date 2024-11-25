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

    def batch_process_spectrum(self, data, window=None, threshold=None, centroided=False, refresh=False):
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

        # TODO: Need a way to test for whether data is centroided already
        if centroided:
            centroids = deepcopy(data)
        else:
            centroids = deepcopy(get_all_centroids(data, window=5, threshold=threshold * 0.1))

        med_spacing = check_spacings(centroids)
        if med_spacing <= self.config.meanpeakspacing_thresh:
            if self.config.verbose:
                print("Median Spacing:", med_spacing, "Removing noise.")
            centroids = remove_noise_cdata(centroids, 100, factor=1.5, mode="median")

        if refresh:
            self.pks = MatchedCollection()

        self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.config)

        return self.pks

    def pks_to_mass(self, binsize=0.1):
        """
        Convert the MatchedCollection to mass
        :return: None
        """
        return self.pks.to_mass_spectrum(binsize)

    def process_file(self, file, scans=None):
        starttime = time.perf_counter()
        self.config.filepath = file
        # Get importer and check it
        reader = ImporterFactory.create_importer(file)
        self.reader = reader
        ext = os.path.splitext(file)[1]
        try:
            print("File:", file, "N Scans:", np.amax(reader.scans))
        except Exception as e:
            print("Could not open:", file)
            return []

        if "centroid" in file:
            centroided = True
            print("Assuming Centroided Data")
        else:
            centroided = False

        t2 = time.perf_counter()
        # Loop over all scans
        for s in reader.scans:
            if scans is not None:
                if s not in scans:
                    continue

            # Open the scan and get the spectrum
            try:
                if ext == ".raw":
                    spectrum = reader.grab_centroid_data(s)
                    centroided = True
                else:
                    spectrum = reader.grab_scan_data(s)
            except Exception as e:
                print("Error Reading Scan", s, e)
                continue
            # If the spectrum is too short, skip it
            if len(spectrum) < 3:
                continue

            self.config.set_scan_info(s, reader)
            # b1 = spectrum[:,1] > 0
            # spectrum = spectrum[b1]
            self.batch_process_spectrum(spectrum, centroided=centroided)

            if s % 10 == 0:
                print("Scan:", s, "Length:", len(spectrum), "Avg. Time per scan:", (time.perf_counter() - t2) / 10.)
                t2 = time.perf_counter()

        print("Time:", time.perf_counter() - starttime)
        print("N Peaks:", len(self.pks.peaks))

        # self.pks.save_pks()
        return reader

    def export_peaks(self, type="prosightlite", filename=None, reader=None, act_type="HCD", max_precursors=None):
        print("Filename:", filename)
        if filename is None:
            filename = "peaks"

        if reader is None:
            reader = self.reader

        if type == "prosightlite":
            self.pks.export_prosightlite(filename)
        elif type == "msalign":
            self.pks.export_msalign(self.config, reader, filename, act_type=act_type, max_precursors=max_precursors)
        elif type == "pkl":
            self.pks.save_pks()
        elif type == "tsv":
            self.pks.export_tsv(filename)
        else:
            raise ValueError("Unknown Export Type", type)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    starttime = time.perf_counter()
    eng = IsoDecRuntime(phaseres=8)

    c = example
    pks = eng.batch_process_spectrum(c, centroided=True)
    print("Time:", time.perf_counter() - starttime)
    exit()
    plot_pks(pks)

    plt.show()

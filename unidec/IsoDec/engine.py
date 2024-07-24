import time

import numpy as np
import torch
from torch.utils.data import DataLoader
from unidec.IsoDec.models import IsoDecClassifier, IsoDecSegmenter, example, IsoDecMixedModel, PhaseModel
from unidec.IsoDec.datatools import create_isodist, fastpeakdetect, get_centroids, get_all_centroids
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import data_dirs, encode_noise, encode_mask, charge_phase_calculator
import os
import unidec.tools as ud
import pickle as pkl
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    mpl.use("WxAgg")
except:
    pass

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


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
        if matchedindexes is not None:
            self.matchedindexes = matchedindexes
            self.matchedcentroids = centroids[matchedindexes]
            self.mask = np.zeros(len(centroids))
            self.mask[matchedindexes] = 1
        if isomatches is not None:
            self.isomatches = isomatches
            self.matchedisodist = isodist[isomatches]


class IsoDecDataset(torch.utils.data.Dataset):
    """
    Dataset class for IsoDec
    """

    def __init__(self, data, mtype=0):
        self.data = data
        self.mtype = mtype

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        out = self.data[idx]
        if self.mtype != 3 and self.mtype != "phase":
            sample = {"emat": out[0], "mask": out[1], "z": out[2]}
        else:
            sample = {"emat": out[0], "z": out[2]}
        return sample


class IsoDecEngine:
    """
    Main class for IsoDec Engine
    """

    def __init__(self, mtype=0):
        self.batch_size = 32
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.pks = MatchedCollection()
        self.classifier = IsoDecClassifier()
        self.segmenter = IsoDecSegmenter()
        self.mixedmodel = IsoDecMixedModel()
        self.phasemodel = PhaseModel()
        self.activemodel = None
        self.mtype = mtype
        self.activescan = -1

    def set_model(self, mtype=0):
        if mtype == 0 or mtype == "classifier":
            self.activemodel = self.classifier
        elif mtype == 1 or mtype == "segmenter":
            self.activemodel = self.segmenter
        elif mtype == 2 or mtype == "mixed":
            self.activemodel = self.mixedmodel
        elif mtype == 3 or mtype == "phase":
            self.activemodel = self.phasemodel
        else:
            raise ValueError("Model type not recognized. Use 'classifier', 'segmenter', or 'mixed'.")
        self.mtype = mtype

    def create_training_dataloader(self, training_path, test_path=None, batch_size=None, noise_percent=0.1, mask=False,
                                   batchsize=None):
        if batchsize is not None:
            self.batch_size = batchsize

        if ".pth" not in training_path:
            training_path = "training_data_" + training_path + ".pth"

        if test_path is None:
            test_path = training_path.replace("training", "test")
        elif ".pth" not in test_path:
            test_path = "test_data_" + test_path + ".pth"

        if batch_size is not None:
            self.batch_size = batch_size

        self.training_data = torch.load(training_path)
        self.test_data = torch.load(test_path)

        if noise_percent > 0:
            nnoise = int(len(self.training_data) * noise_percent)
            print(f"Adding {nnoise} noise samples to training data")
            for i in range(nnoise):
                d = self.training_data[i]
                noise = encode_noise(d, mtype=self.mtype)
                self.training_data.append(noise)

            nnoisetest = int(len(self.test_data) * noise_percent)
            print(f"Adding {nnoisetest} noise samples to test data")
            for i in range(nnoisetest):
                d = self.test_data[i]
                noise = encode_noise(d, mtype=self.mtype)
                self.test_data.append(noise)

        self.training_data = IsoDecDataset(self.training_data, self.mtype)
        self.test_data = IsoDecDataset(self.test_data, self.mtype)

        if mask:
            print("Applying Masks")
            tmod = []
            for s in self.training_data:
                s = encode_mask(s)
                tmod.append(s)
            self.training_data = IsoDecDataset(tmod)

            tmod = []
            for s in self.test_data:
                s = encode_mask(s)
                tmod.append(s)
            self.test_data = IsoDecDataset(tmod)
            print("Masks Done")

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.batch_size, shuffle=True)

        for s in self.train_dataloader:
            X = s["emat"]
            # y = s["mask"]
            z = s["z"]
            print(f"Shape of X [N, C, H, W]: {X.shape}")
            # print(f"Shape of y: {y.shape} {y.dtype}")
            print(f"Shape of z: {z.shape} {z.dtype}")
            break

        for s in self.train_dataloader:
            X = s["emat"]
            # y = s["mask"]
            z = s["z"]
            print(f"Shape of X [N, C, H, W]: {X.shape}")
            # print(f"Shape of y: {y.shape} {y.dtype}")
            print(f"Shape of z: {z.shape} {z.dtype}")
            break

    def train_model(self, epochs=30, save=True, mtype=None):
        if mtype is None:
            mtype = self.activemodel
        self.set_model(mtype)

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.activemodel.train_model(self.train_dataloader)
            self.activemodel.evaluate_model(self.test_dataloader)
        print("Done!")
        if save:
            self.activemodel.save_model()

    def save_bad_data(self, filename="bad_data.pkl", maxbad=50):
        if self.activemodel is None:
            self.set_model(mtype=self.mtype)
        if self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        bad_data = self.activemodel.evaluate_model(self.test_dataloader, maxbad=maxbad)
        with open(filename, "wb") as f:
            pkl.dump(bad_data, f)
            print(f"Saved {len(bad_data)} bad data points to {filename}")

    def single_charge_prediction(self, centroids):
        """
        Predict the charge state of a single set of centroids. No segmentation or fancy business. Just prediction.
        :param centroids: Centroids as 2D numpy array [m/z, intensity], must be at least 3 peaks
        :return: Charge state as int
        """
        if len(centroids) < 3:
            return 0
        z, mask = self.classifier.predict(centroids)
        return int(z)

    def two_stage_charge_prediction(self, centroids):
        """
        Predict the charge state of set of centroids by first segmenting them and then predicting the charge.
        :param centroids: Centroids as 2D numpy array [m/z, intensity], must be at least 3 peaks
        :return: Charge state as int, mask as 1D numpy array of index values from centroids that are peaks
        """
        if len(centroids) < 3:
            return 0, []
        _, maskedindexes = self.segmenter.predict(centroids)

        mcent = centroids[maskedindexes]
        z = self.single_charge_prediction(mcent)
        return int(z), maskedindexes

    def complex_charge_prediction(self, centroids):
        if len(centroids) < 3:
            return 0, []
        z, mask = self.mixedmodel.predict(centroids)
        return int(z), mask

    def phase_calculator(self, centroids):
        z, mask = charge_phase_calculator(centroids)
        return int(z), mask

    def phase_predictor(self, centroids):
        if len(centroids) < 3:
            return 0
        z = self.phasemodel.predict(centroids)
        return int(z)

    def two_stage_phase_prediction(self, centroids):
        """
        Predict the charge state of set of centroids by first segmenting them and then predicting the charge.
        :param centroids: Centroids as 2D numpy array [m/z, intensity], must be at least 3 peaks
        :return: Charge state as int, mask as 1D numpy array of index values from centroids that are peaks
        """
        if len(centroids) < 3:
            return 0, []
        _, maskedindexes = self.segmenter.predict(centroids)

        mcent = centroids[maskedindexes]
        z = self.phase_predictor(mcent)
        return int(z), maskedindexes

    def get_matches(self, centroids, pks=None):
        #starttime = time.perf_counter()

        z = self.phase_predictor(centroids)

        # z2, mask2 = self.complex_charge_prediction(mcent)

        if z == 0:
            return 0, [], []

        peakmz = centroids[np.argmax(centroids[:, 1]), 0]
        isodist = create_isodist(peakmz, z, centroids)
        isodist = optimize_shift(centroids, isodist, z)

        matchedindexes, isomatches = match_peaks(centroids, isodist)

        isoper = len(isomatches) / len(isodist)

        #print("Charge Predictions:", z, "Time:", time.perf_counter() - starttime)

        if len(matchedindexes) > 2 and isoper > 0.5:
            m = MatchedPeak(centroids, isodist, z, peakmz, matchedindexes, isomatches)
            m.scan = self.activescan
            if pks is not None:
                pks.add_peak(m)
            else:
                self.pks.add_peak(m)
            return z, matchedindexes, isomatches
        else:
            #print("No Match", z, len(matchedindexes), isoper)
            return 0, [], []

    def fancy_prediction(self, centroids, pks=None):
        if pks is None:
            pks = MatchedCollection()
        z, matchedindexes, isomatches = self.get_matches(centroids, pks)
        if z == 0:
            return pks

        remaining_centroids = np.delete(centroids, matchedindexes, axis=0)
        if remaining_centroids.shape[0] == 0:
            return pks
        z2, matchedindexes2, isomatches2 = self.get_matches(remaining_centroids, pks)

        return pks

    def process_spectrum(self, data, window=50, threshold=0.0001, mzwindow=[-1.5, 3.5]):
        starttime = time.perf_counter()
        peaks = fastpeakdetect(data, window=window, threshold=threshold)
        # sort peaks
        peaks = np.array(sorted(peaks, key=lambda x: x[1], reverse=True))
        # TODO: Need a way to test for whether data is centroided already
        centroids = get_all_centroids(data, window=5, threshold=threshold * 0.1)
        for i, p in enumerate(peaks):
            peakmz = p[0]

            # Find all centroids in the neighborhood of the peak
            b1 = centroids[:, 0] > peakmz + mzwindow[0]
            b2 = centroids[:, 0] < peakmz + mzwindow[1]
            b = b1 & b2
            if len(centroids[b]) < 3:
                continue

            # Get the centroids around the peak
            z, matchedindexes, isomatches = self.get_matches(centroids[b], self.pks)
            if z == 0:
                continue
            # correct matchedindexes based on b
            startindex = np.where(b)[0][0]
            # print(i, peakmz, startindex, matchedindexes, len(centroids))
            matchedindexes = matchedindexes + startindex

            centroids = np.delete(centroids, matchedindexes, axis=0)
            if centroids.shape[0] == 0:
                break
        print("Time:", time.perf_counter() - starttime)
        return self.pks

    def process_file(self, file):
        starttime = time.perf_counter()
        # Get importer and check it
        reader = ud.get_importer(file)
        try:
            print("File:", file, "N Scans:", np.amax(reader.scans))
        except Exception as e:
            print("Could not open:", file)
            return []

        # Loop over all scans
        for s in reader.scans:
            # Open the scan and get the spectrum
            try:
                spectrum = reader.grab_scan_data(s)
            except:
                print("Error Reading Scan", s)
                continue
            # If the spectrum is too short, skip it
            if len(spectrum) < 3:
                continue
            print("Scan:", s, "Length:", len(spectrum))
            self.activescan = s
            self.process_spectrum(spectrum)
        print("Time:", time.perf_counter() - starttime)
        print("N Peaks:", len(self.pks.peaks))

        self.pks.save_pks()
        return reader

    def plot_pks(self, data, scan=-1, mzwindow=[-1.5, 3.5], show=False):
        plt.plot(data[:, 0], data[:, 1])
        for p in self.pks.peaks:
            if scan == -1 or p.scan==scan:
                color = p.color
                isodist = p.isodist
                cplot(isodist, color=color, factor=-1)
                centroids = p.centroids
                peakmz = p.mz
                cplot(centroids)
        if show:
            plt.show()


def get_charge_nn(centroids):
    """
    Simple function to launch the engine, predict the charge, and return the charge
    :param centroids: Centroid data, 2D numpy array as [m/z, intensity]
    :return: Charge state as int
    """
    engine = IsoDecEngine()
    z = engine.single_charge_prediction(centroids)
    return z


if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecEngine(mtype=3)
    topdirectory = "C:\\Data\\IsoNN\\training"
    if True:
        for i, d in enumerate(data_dirs[4:]):
            directory = os.path.join(topdirectory, d)
            print(i, d, directory)
            os.chdir(directory)

            eng.create_training_dataloader("phase81", mask=False, noise_percent=0.2, batchsize=32)
            eng.train_model(epochs=10, mtype="phase")
            eng.save_bad_data()
        print("Time:", time.perf_counter() - starttime)
        exit()
    c = example
    pks = eng.fancy_prediction(c)
    cplot(c)
    for p in pks.peaks:
        cplot(p.centroids, mask=p.mask, z=p.z, mcolor=p.color, zcolor=p.color)

    # z, p = eng.classifier.predict(c)
    # print(p)
    # cplot(c, mask=p.flatten(), z=z)
    plt.show()

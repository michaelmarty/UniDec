import numpy as np
import torch
from torch.utils.data import DataLoader
from IsoDec.models import IsoDecClassifier, IsoDecSegmenter, example
from IsoDec.datatools import create_isodist
from IsoDec.match import *
import os
import pickle as pkl
import matplotlib.pyplot as plt

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


class MatchedCollection:
    def __init__(self):
        self.peaks = []

    def add_peak(self, peak):
        self.peaks.append(peak)


class MatchedPeak:
    def __init__(self, centroids, isodist, z, mz, matchedindexes=None, isomatches=None):
        self.mz = mz
        self.z = z
        self.centroids = centroids
        self.isodist = isodist
        self.matchedcentroids = None
        self.matchedisodist = None
        self.matchedindexes = None
        self.isomatches = None
        if matchedindexes is not None:
            self.matchedindexes = matchedindexes
            self.matchedcentroids = centroids[matchedindexes]
        if isomatches is not None:
            self.isomatches = isomatches
            self.matchedisodist = isodist[isomatches]


class IsoDecEngine:
    def __init__(self, mtype=0):
        self.mtype = mtype
        self.isomod = self.set_model(self.mtype)
        self.batch_size = 128
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.pks = MatchedCollection()

    def set_model(self, mtype=0):
        self.mtype = mtype
        if mtype == 0:
            self.isomod = IsoDecClassifier()
        elif mtype == 1:
            self.isomod = IsoDecSegmenter()
        else:
            raise ValueError("Model type not recognized")
        return self.isomod

    def create_training_dataloader(self, training_path, test_path, batch_size=None):
        if batch_size is not None:
            self.batch_size = batch_size

        self.training_data = torch.load(training_path)
        self.test_data = torch.load(test_path)

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.batch_size, shuffle=True)

        for X, y in self.train_dataloader:
            print(f"Shape of X [N, C, H, W]: {X.shape}")
            print(f"Shape of y: {y.shape} {y.dtype}")
            break

        for X, y in self.test_dataloader:
            print(f"Shape of X [N, C, H, W]: {X.shape}")
            print(f"Shape of y: {y.shape} {y.dtype}")
            break

    def train_model(self, epochs=30, save=True):
        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.isomod.train_model(self.train_dataloader)
            self.isomod.evaluate_model(self.test_dataloader)
        print("Done!")
        if save:
            self.isomod.save_model()

    def single_charge_prediction(self, centroids):
        z = self.isomod.predict(centroids)
        return z

    def get_matches(self, centroids, pks=None):
        z = self.single_charge_prediction(centroids)
        peakmz = centroids[np.argmax(centroids[:, 1]), 0]
        print(peakmz)
        #isodist = create_isodist(peakmz, z, centroids)
        #isodist = optimize_shift(centroids, isodist, z)

        matchedindexes, isomatches = match_peaks(centroids, isodist)

        isoper = len(isomatches) / len(isodist)
        print(isoper)

        if len(matchedindexes) > 2 and isoper > 0.6:
            m = MatchedPeak(centroids, isodist, z, peakmz, matchedindexes, isomatches)
            if pks is not None:
                pks.add_peak(m)
            else:
                self.pks.add_peak(m)
            return z, matchedindexes, isomatches
        else:
            return 0, [], []

    def fancy_prediction(self, centroids):
        pks = MatchedCollection()
        z, matchedindexes, isomatches = self.get_matches(centroids, pks)
        if z == 0:
            print("No Match")
            return pks

        remaining_centroids = np.delete(centroids, matchedindexes, axis=0)
        if remaining_centroids.shape[0] == 0:
            return pks
        z2, matchedindexes2, isomatches2 = self.get_matches(remaining_centroids, pks)
        if z2 == 0:
            print("No Match")

        return pks


def get_charge_nn(centroids):
    engine = IsoDecEngine()
    z = engine.single_charge_prediction(centroids)
    return z


if __name__ == "__main__":
    os.chdir("C:\\Data\\IsoNN\\multi")
    eng = IsoDecEngine(1)
    eng.create_training_dataloader("training_data_large.pth", "test_data_large.pth")
    eng.train_model(epochs=5)

    c = example
    p = eng.isomod.predict(c)
    indexes = eng.isomod.indexes
    cplot(c)
    colors = ['g', 'b', 'c', 'm', 'y', 'k', 'w']
    for j, vec in enumerate(p):
        v = vec.astype(bool)
        i = indexes[v]
        b1 = np.zeros(len(c))
        b1[i.astype(int)] = 1
        b1 = b1.astype(bool)
        if len(b1) > len(c):
            b1 = b1[:len(c)]

        d = c[b1]
        cplot(d, color=colors[j], factor=-1)
    plt.show()

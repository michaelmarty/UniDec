import time
import warnings
import numpy as np
from itertools import chain
import torch
from torch.utils.data import DataLoader
from unidec.IsoDec.models import example, PhaseModel
from unidec.IsoDec.datatools import fastpeakdetect, get_all_centroids, fastnearest
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import data_dirs, encode_noise, encode_phase_all, small_data_dirs, \
    encode_double
import os
import math
import unidec.tools as ud
import pickle as pkl
import matplotlib.pyplot as plt
from unidec.IsoDec.c_interface import IsoDecWrapper
from unidec.IsoDec.plots import *


class IsoDecDataset(torch.utils.data.Dataset):
    """
    Dataset class for IsoDec
    """

    def __init__(self, emat, z):
        self.emat = [torch.as_tensor(e, dtype=torch.float32) for e in emat]
        self.z = torch.as_tensor(z, dtype=torch.long)

    def __len__(self):
        return len(self.z)

    def __getitem__(self, idx):
        return [self.emat[idx], self.z[idx]]


class IsoDecEngine:
    """
    Main class for IsoDec Engine
    """

    def __init__(self, phaseres=8):
        self.batch_size = 32
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.pks = MatchedCollection()
        self.phaseres = phaseres
        self.phasemodel = PhaseModel()
        if self.phaseres != 8:
            self.phasemodel.dims = [50, self.phaseres]
            self.phasemodel.modelid = 0
        self.activescan = -1
        self.activescanrt = -1
        self.activescanorder = -1
        self.matchtol = 0.01
        self.minpeaks = 3
        self.minmatchper = 0.67
        self.maxshift = 3  # This will get overwritten for smaller z, where it's dangerous to have more than 1 or 2
        self.mzwindow = [-1.5, 2.5]
        self.plusoneintwindow = [0.1, 0.5]
        self.knockdown_rounds = 3
        self.peakwindow = 10
        self.test_centroids = []
        self.training_centroids = []
        self.test_batch_size = 2048

        self.use_wrapper = False

        if self.use_wrapper:
            self.wrapper = IsoDecWrapper()
        else:
            self.wrapper = None

    def add_noise(self, noise_percent):
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        nnoise = int(ltraining * noise_percent)
        print(f"Adding {nnoise} noise samples to training data")
        emats = []
        zs = []
        centroids = []
        for i in range(nnoise):
            index = np.random.randint(0, ltraining - 1)
            centroid = self.training_data[2][index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.phaseres)
            emats.append(emat)
            zs.append(0)
            centroids.append(centroid)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)
        self.training_data[2] = self.training_data[2] + centroids

        nnoisetest = int(ltest * noise_percent)
        print(f"Adding {nnoisetest} noise samples to test data")
        emats = []
        zs = []
        centroids = []
        for i in range(nnoisetest):
            index = np.random.randint(0, ltest - 1)
            centroid = self.test_data[2][index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.phaseres)
            emats.append(emat)
            zs.append(0)
            centroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + centroids

    def add_doubles(self, double_percent):
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        ndouble = int(ltraining * double_percent)
        print(f"Adding {ndouble} double samples in training data")
        emats = []
        zs = []
        centroids = []
        for i in range(ndouble):
            index = np.random.randint(0, ltraining - 1)
            centroid1 = self.training_data[2][index]
            index2 = np.random.randint(0, ltraining - 1)
            centroid2 = self.training_data[2][index2]
            z = self.training_data[1][index]
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.phaseres)
            emats.append(emat)
            zs.append(z)
            centroids.append(centroid)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)
        self.training_data[2] = self.training_data[2] + centroids

        ndoubletest = int(ltest * double_percent)
        print(f"Adding {ndoubletest} double samples in test data")
        emats = []
        zs = []
        centroids = []
        for i in range(ndoubletest):
            index = np.random.randint(0, ltest - 1)
            index2 = np.random.randint(0, ltest - 1)
            centroid1 = self.test_data[2][index]
            centroid2 = self.test_data[2][index2]
            z = self.test_data[1][index]
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.phaseres)
            emats.append(emat)
            zs.append(z)
            centroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + centroids

    def load_training_data(self, training_path, test_path=None, noise_percent=0.1, double_percent=0.1):
        ext = ".npz"
        if ext not in training_path:
            training_path = "training_data_" + training_path + ext

        if test_path is None:
            test_path = training_path.replace("training", "test")
        elif ext not in test_path:
            test_path = "test_data_" + test_path + ext

        td = np.load(training_path, allow_pickle=True)
        self.training_data = [td["emat"], td["z"], list(td["centroids"])]
        td = np.load(test_path, allow_pickle=True)
        self.test_data = [td["emat"], td["z"], list(td["centroids"])]

        if double_percent > 0:
            self.add_doubles(double_percent)

        if noise_percent > 0:
            self.add_noise(noise_percent)

        print("Loaded:", len(self.training_data[0]), "Training Samples")

    def create_training_dataloader(self, training_path, test_path=None, noise_percent=0.1, batchsize=None,
                                   double_percent=0.1):
        if batchsize is not None:
            self.batch_size = batchsize

        self.load_training_data(training_path, test_path=test_path, noise_percent=noise_percent,
                                double_percent=double_percent)

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.test_batch_size, shuffle=False,
                                          pin_memory=False)



    def create_merged_dataloader(self, dirs, training_path, noise_percent=0.1, batchsize=None, double_percent=0.1):
        if batchsize is not None:
            self.batch_size = batchsize

        training_data = []
        test_data = []

        for d in dirs:
            os.chdir(d)
            print(d)
            self.load_training_data(training_path, noise_percent=noise_percent, double_percent=double_percent)
            training_data.append(self.training_data)
            test_data.append(self.test_data)

        self.training_data = [np.concatenate([t[0] for t in training_data], axis=0),
                              np.concatenate([t[1] for t in training_data], axis=0)]
        self.test_data = [np.concatenate([t[0] for t in test_data], axis=0),
                          np.concatenate([t[1] for t in test_data], axis=0)]

        self.training_centroids = list(chain(*[t[2] for t in training_data]))
        self.test_centroids = list(chain(*[t[2] for t in test_data]))

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.test_batch_size, shuffle=False,
                                          pin_memory=True)

        print(f"Training Data Length: {len(self.training_data)}")
        print(f"Test Data Length: {len(self.test_data)}")

        #plot_zdist(self)

    def train_model(self, epochs=30, save=True, lossfn="crossentropy"):
        starttime = time.perf_counter()

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")
        self.phasemodel.get_class_weights(self.train_dataloader)
        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.phasemodel.train_model(self.train_dataloader, lossfn=lossfn)
            self.phasemodel.evaluate_model(self.test_dataloader)

        if save:
            self.phasemodel.save_model()
        print("Done! Time:", time.perf_counter() - starttime)

    def save_bad_data(self, filename="bad_data.pkl", maxbad=50):
        if self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        bad_data = self.phasemodel.evaluate_model(self.test_dataloader, savebad=True)
        output = []
        for b in bad_data:
            i1 = b[0]
            i2 = b[1]
            index = i1 * self.test_batch_size + i2
            if index > len(self.test_centroids):
                print("Index Error:", index, len(self.test_centroids), i1, i2, self.test_batch_size)
                continue
            centroid = self.test_centroids[index]
            output.append([centroid, b[2], b[3]])

        outindexes = np.random.randint(0, len(output), maxbad)
        output = [output[o] for o in outindexes]

        with open(filename, "wb") as f:
            pkl.dump(output, f)
            print(f"Saved {len(output)} bad data points to {filename}")

    def phase_predictor(self, centroids):
        if self.use_wrapper:
            z = IsoDecWrapper().predict_charge(centroids)
        else:
            z = self.phasemodel.predict(centroids)
        return int(z)

    def save_peak(self, centroids, z, peakmz, pks=None):
        isodist, matchedindexes, isomatches, peakmz, monoiso, massdist = optimize_shift(centroids, z, peakmz,
                                                                                        tol=self.matchtol,
                                                                                        maxshift=self.maxshift)

        # isoper = len(isomatches) / len(isodist)
        areaper = np.sum(isodist[isomatches, 1]) / np.sum(isodist[:, 1])
        # print("Charge Predictions:", z, "Time:", time.perf_counter() - starttime)
        # print(peakmz, areaper, len(matchedindexes), z)

        # Elaborate set of conditions to determine if the peak is a 1+ peak and if it's matching 2 peaks ok
        # If so, let it pass with just two matches.
        minpeaks = self.minpeaks
        if z == 1:
            if len(matchedindexes) == 2:
                if isomatches[0] == 0 and isomatches[1] == 1:
                    int1 = centroids[matchedindexes[0], 1]
                    int2 = centroids[matchedindexes[1], 1]
                    ratio = int2 / int1
                    if self.plusoneintwindow[0] < ratio < self.plusoneintwindow[1]:
                        minpeaks = 2
                        areaper = 1

        if len(matchedindexes) >= minpeaks and areaper >= self.minmatchper:
            m = MatchedPeak(z, peakmz, centroids, isodist, matchedindexes, isomatches)
            m.scan = self.activescan
            m.rt = self.activescanrt
            m.ms_order = self.activescanorder
            m.monoiso = monoiso
            m.massdist = massdist
            m.avgmass = np.sum(massdist[:, 0] * massdist[:, 1]) / np.sum(massdist[:, 1])
            m.peakmass = massdist[np.argmax(massdist[:, 1]), 0]
            if pks is not None:
                pks.add_peak(m)
                pks.add_pk_to_masses(m, 10)
            else:
                self.pks.add_peak(m)
                self.pks.add_pk_to_masses(m, 10)
            return matchedindexes
        else:
            return []

    def get_matches(self, centroids, z, peakmz, pks=None):
        # starttime = time.perf_counter()
        if len(centroids) < self.minpeaks:
            return []
        if z == 0 or z > 50:
            return []
        matchedindexes = self.save_peak(centroids, z, peakmz, pks=pks)
        return matchedindexes

    def batch_process_spectrum(self, data, window=None, threshold=0.001, centroided=False):
        starttime = time.perf_counter()
        if window is None:
            window = self.peakwindow

        # TODO: Need a way to test for whether data is centroided already
        if centroided:
            centroids = data
        else:
            centroids = get_all_centroids(data, window=5, threshold=threshold * 0.1)

        if self.use_wrapper:
            self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.activescan)
        else:
            kwindow = window
            threshold = threshold
            for i in range(self.knockdown_rounds):

                if i > 1:
                    kwindow = 2
                peaks = fastpeakdetect(centroids, window=kwindow, threshold=threshold)
                # print("Knockdown:", i, "Peaks:", len(peaks))
                if len(peaks) == 0:
                    break
                emats, peaks, centlist, indexes = encode_phase_all(centroids, peaks, lowmz=self.mzwindow[0],
                                                                   highmz=self.mzwindow[1], phaseres=self.phaseres)
                emats = [torch.as_tensor(e, dtype=torch.float32) for e in emats]
                # emats = torch.as_tensor(emats, dtype=torch.float32).to(self.phasemodel.device)
                data_loader = DataLoader(emats, batch_size=1024, shuffle=False, pin_memory=True)
                preds = self.phasemodel.batch_predict(data_loader)
                knockdown = []
                ngood = 0
                # print(peaks, len(peaks))
                for j, p in enumerate(peaks):
                    z = preds[j]
                    kindex = fastnearest(centroids[:, 0], p[0])
                    if kindex in knockdown:
                        continue
                    if z == 0:
                        knockdown.append(kindex)
                        continue
                    # Get the centroids around the peak
                    matchedindexes = self.get_matches(centlist[j], z, p[0], pks=self.pks)

                    if len(matchedindexes) > 0:
                        ngood += 1
                        # Find matches
                        indval = indexes[j]
                        matchindvals = indval[matchedindexes]
                        # Knock them down
                        knockdown.extend(matchindvals)
                # print("NGood:", ngood)
                if len(knockdown) == 0:
                    break
                knockdown = np.array(knockdown)
                centroids = np.delete(centroids, knockdown, axis=0)

                if len(centroids) < 3:
                    break

            # print("Time:", time.perf_counter() - starttime)
        return self.pks

    def process_file(self, file, scans=None):
        starttime = time.perf_counter()
        # Get importer and check it
        reader = ud.get_importer(file)
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

            self.activescan = s
            self.activescanrt = reader.get_scan_time(s)
            self.activescanorder = reader.get_ms_order(s)
            # b1 = spectrum[:,1] > 0
            # spectrum = spectrum[b1]
            self.batch_process_spectrum(spectrum, centroided=centroided)

            if s % 10 == 0:
                print("Scan:", s, "Length:", len(spectrum), "Avg. Time per scan:", (time.perf_counter() - t2) / 10.)
                t2 = time.perf_counter()

        print("Time:", time.perf_counter() - starttime)
        print("N Peaks:", len(self.pks.peaks))

        self.pks.save_pks()
        return reader



class IsoDecConfig:
    def __init__(self, filepath):
        self.filepath = filepath
        self.batch_size = 32
        self.matchtol = 0.01
        self.minpeaks = 3
        self.minmatchper = 0.67
        self.maxshift = 3  # This will get overwritten for smaller z, where it's dangerous to have more than 1 or 2
        self.mzwindow = [-1.5, 5.5]
        self.plusoneintwindow = [0.1, 0.5]
        self.knockdown_rounds = 3

if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecEngine(phaseres=4)
    topdirectory = "C:\\Data\\IsoNN\\training"

    dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]
    eng.create_merged_dataloader(dirs, "phase4", noise_percent=0.0, batchsize=32, double_percent=0.0)
    #eng.train_model(epochs=3)
    eng.train_model(epochs=3, lossfn="focal")

    # eng.create_merged_dataloader([os.path.join(topdirectory, small_data_dirs[2])], "phase82", noise_percent=0.2,
    #                             batchsize=32, double_percent=0.2)
    # eng.save_bad_data()
    exit()
    c = example
    pks = eng.batch_process_spectrum(c, centroided=True)
    cplot(c)
    for p in pks.peaks:
        print("Charge:", p.z)
        cplot(p.centroids, z=p.z, mcolor=p.color, zcolor=p.color)

    # z, p = eng.classifier.predict(c)
    # print(p)
    # cplot(c, mask=p.flatten(), z=z)
    plt.show()

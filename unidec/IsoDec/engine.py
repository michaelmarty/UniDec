import time
import warnings
import numpy as np
import torch
from torch.utils.data import DataLoader
from unidec.IsoDec.models import example, PhaseModel
from unidec.IsoDec.datatools import fastpeakdetect, get_all_centroids
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import data_dirs, encode_noise, charge_phase_calculator, encode_phase_all, small_data_dirs, \
    encode_double
import os
import unidec.tools as ud
import pickle as pkl
import matplotlib.pyplot as plt

try:
    mpl.use("WxAgg")
except:
    pass


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

    def __init__(self, mtype=3):
        self.batch_size = 32
        self.training_data = None
        self.test_data = None
        self.train_dataloader = None
        self.test_dataloader = None
        self.pks = MatchedCollection()
        # self.classifier = IsoDecClassifier()
        # self.segmenter = IsoDecSegmenter()
        # self.mixedmodel = IsoDecMixedModel()
        self.phasemodel = PhaseModel()
        self.activemodel = self.phasemodel
        self.mtype = mtype
        self.activescan = -1
        self.matchtol = 0.01
        self.minpeaks = 3
        self.minmatchper = 0.5
        self.maxshift = 2  # This will get overwritten for smaller z, where it's dangerous to have more than 1 or 2
        self.mzwindow = [-1.5, 5.5]
        self.knockdown_rounds = 2
        self.test_centroids = []
        self.test_batch_size = 2048

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

    def add_noise(self, noise_percent):
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        nnoise = int(ltraining * noise_percent)
        print(f"Adding {nnoise} noise samples to training data")
        emats = []
        zs = []
        for i in range(nnoise):
            index = i % ltest
            centroid = self.test_centroids[index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]))
            emats.append(emat)
            zs.append(z)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)

        nnoisetest = int(ltest * noise_percent)
        print(f"Adding {nnoisetest} noise samples to test data")
        emats = []
        zs = []
        for i in range(nnoisetest):
            index = i % ltest
            centroid = self.test_centroids[index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]))
            emats.append(emat)
            zs.append(z)
            self.test_centroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)

    def add_doubles(self, double_percent):
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        ndouble = int(ltraining * double_percent)
        print(f"Adding {ndouble} double samples in training data")
        emats = []
        zs = []
        for i in range(ndouble):
            index = i % ltest
            centroid = self.test_centroids[index]
            index2 = np.random.randint(0, ltest - 1)
            centroid2 = self.test_centroids[index2]
            z = self.training_data[1][index]
            emat, centroid = encode_double(centroid, centroid2)
            emats.append(emat)
            zs.append(z)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)

        ndoubletest = int(ltest * double_percent)
        print(f"Adding {ndoubletest} double samples in test data")
        emats = []
        zs = []
        for i in range(ndoubletest):
            index = i % ltest
            centroid = self.test_centroids[index]
            index2 = np.random.randint(0, ltest - 1)
            centroid2 = self.test_centroids[index2]
            z = self.test_data[1][index]
            emat, centroid = encode_double(centroid, centroid2)
            emats.append(emat)
            zs.append(z)
            self.test_centroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)

    def load_training_data(self, training_path, test_path=None, noise_percent=0.1, double_percent=0.1):
        ext = ".npz"
        if ext not in training_path:
            training_path = "training_data_" + training_path + ext

        if test_path is None:
            test_path = training_path.replace("training", "test")
        elif ext not in test_path:
            test_path = "test_data_" + test_path + ext

        td = np.load(training_path)
        self.training_data = [td["emat"], td["z"]]
        td = np.load(test_path, allow_pickle=True)
        self.test_data = [td["emat"], td["z"]]
        test_centroids = list(td["centroids"])
        self.test_centroids.extend(test_centroids)

        if double_percent > 0:
            self.add_doubles(double_percent)

        if noise_percent > 0:
            self.add_noise(noise_percent)

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])
        print("Loaded:", len(self.training_data), "Training Samples")

    def create_training_dataloader(self, training_path, test_path=None, noise_percent=0.1, batchsize=None,
                                   double_percent=0.1):
        if batchsize is not None:
            self.batch_size = batchsize

        self.load_training_data(training_path, test_path=test_path, noise_percent=noise_percent,
                                double_percent=double_percent)

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.test_batch_size, shuffle=False,
                                          pin_memory=True)

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

        self.training_data = torch.utils.data.ConcatDataset(training_data)
        self.test_data = torch.utils.data.ConcatDataset(test_data)

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.test_batch_size, shuffle=False,
                                          pin_memory=True)

        print(f"Training Data Length: {len(self.training_data)}")
        print(f"Test Data Length: {len(self.test_data)}")

    def train_model(self, epochs=30, save=True, mtype=None):
        starttime = time.perf_counter()
        if mtype is not None:
            self.set_model(mtype)

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.activemodel.train_model(self.train_dataloader)
            self.activemodel.evaluate_model(self.test_dataloader)

        if save:
            self.activemodel.save_model()
        print("Done! Time:", time.perf_counter() - starttime)

    def save_bad_data(self, filename="bad_data.pkl", maxbad=50):
        if self.activemodel is None:
            self.set_model(mtype=self.mtype)
        if self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        bad_data = self.activemodel.evaluate_model(self.test_dataloader, savebad=True)
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

    def phase_calculator(self, centroids):
        z, mask = charge_phase_calculator(centroids)
        return int(z), mask

    def phase_predictor(self, centroids):
        z = self.phasemodel.predict(centroids)
        return int(z)

    def save_peak(self, centroids, z, pks=None):

        isodist, matchedindexes, isomatches, peakmz, monoiso, massdist = optimize_shift(centroids, z, tol=self.matchtol,
                                                                                        maxshift=self.maxshift)

        isoper = len(isomatches) / len(isodist)
        # print("Charge Predictions:", z, "Time:", time.perf_counter() - starttime)

        if len(matchedindexes) >= self.minpeaks and isoper > self.minmatchper:
            m = MatchedPeak(centroids, isodist, z, peakmz, matchedindexes, isomatches)
            m.scan = self.activescan
            m.monoiso = monoiso
            m.massdist = massdist
            if pks is not None:
                pks.add_peak(m)
            else:
                self.pks.add_peak(m)
            return matchedindexes
        else:
            # print("No Match", z, len(matchedindexes), isoper)
            return []

    def get_matches(self, centroids, pks=None, z=None):
        # starttime = time.perf_counter()
        if len(centroids) < self.minpeaks:
            return []
        if z is None:
            z = self.phase_predictor(centroids)
        if z == 0 or z > 50:
            return []
        matchedindexes = self.save_peak(centroids, z, pks=pks)
        return matchedindexes

    def batch_process_spectrum(self, data, window=5, threshold=0.0001, centroided=False):
        starttime = time.perf_counter()

        # TODO: Need a way to test for whether data is centroided already
        if centroided:
            centroids = data
        else:
            centroids = get_all_centroids(data, window=5, threshold=threshold * 0.1)

        for i in range(self.knockdown_rounds):
            peaks = fastpeakdetect(centroids, window=window, threshold=threshold)

            emats, peaks, centlist, indexes = encode_phase_all(centroids, peaks, lowmz=self.mzwindow[0],
                                                               highmz=self.mzwindow[1])
            emats = [torch.as_tensor(e, dtype=torch.float32) for e in emats]
            # emats = torch.as_tensor(emats, dtype=torch.float32).to(self.phasemodel.device)
            data_loader = DataLoader(emats, batch_size=1024, shuffle=False, pin_memory=True)
            preds = self.phasemodel.batch_predict(data_loader)

            knockdown = []
            for i, p in enumerate(peaks):
                z = preds[i]
                if z == 0:
                    continue
                # Get the centroids around the peak
                matchedindexes = self.get_matches(centlist[i], self.pks, z=z)
                # Find matches
                indval = indexes[i]
                matchindvals = indval[matchedindexes]
                # Knock them down
                knockdown.extend(matchindvals)
            if len(knockdown) == 0:
                break
            knockdown = np.array(knockdown)
            centroids = np.delete(centroids, knockdown, axis=0)

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
            self.batch_process_spectrum(spectrum, centroided=centroided)

            if s % 10 == 0:
                print("Scan:", s, "Length:", len(spectrum), "Avg. Time per scan:", (time.perf_counter() - t2) / 10.)
                t2 = time.perf_counter()

        print("Time:", time.perf_counter() - starttime)
        print("N Peaks:", len(self.pks.peaks))

        self.pks.save_pks()
        return reader

    def plot_pks(self, data, scan=-1, show=False):
        plt.subplot(121)
        plt.plot(data[:, 0], data[:, 1])

        for p in self.pks.peaks:
            if scan == -1 or p.scan == scan:
                color = p.color
                isodist = p.isodist
                plt.subplot(121)
                cplot(isodist, color=color, factor=-1)
                centroids = p.centroids
                peakmz = p.mz
                cplot(centroids)
                plt.subplot(122)
                massdist = p.massdist
                cplot(massdist, color=color)
        if show:
            plt.show()


if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecEngine()
    topdirectory = "C:\\Data\\IsoNN\\training"

    dirs = [os.path.join(topdirectory, d) for d in data_dirs]
    eng.create_merged_dataloader(dirs, "phase82", noise_percent=0.2, batchsize=32, double_percent=0.2)
    eng.train_model(epochs=30)
    eng.save_bad_data()

    exit()
    c = example
    pks = eng.batch_process_spectrum(c, centroided=True)
    cplot(c)
    for p in pks.peaks:
        cplot(p.centroids, mask=p.mask, z=p.z, mcolor=p.color, zcolor=p.color)

    # z, p = eng.classifier.predict(c)
    # print(p)
    # cplot(c, mask=p.flatten(), z=z)
    plt.show()

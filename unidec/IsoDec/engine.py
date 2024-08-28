import time
import warnings
import numpy as np
from itertools import chain
import torch
from torch.utils.data import DataLoader
from unidec.IsoDec.models import example, PhaseModel
from unidec.IsoDec.datatools import fastpeakdetect, get_all_centroids, fastnearest, calculate_cosinesimilarity2
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import data_dirs, encode_noise, encode_phase_all, small_data_dirs, \
    encode_double
import os
import math
from copy import deepcopy
import unidec.tools as ud
import pickle as pkl
import matplotlib.pyplot as plt
from unidec.IsoDec.c_interface import IsoDecWrapper
from unidec.IsoDec.plots import *
import platform
import numba as nb


class IsoDecDataset(torch.utils.data.Dataset):
    """
    Dataset class for IsoDec
    """

    def __init__(self, emat, z):
        """
        Initialize the dataset
        :param emat: List of encoded matrices
        :param z: List of charge assignments
        """
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

        self.pks = MatchedCollection()

        self.phaseres = phaseres
        self.phasemodel = PhaseModel()
        self.maxz = 50
        self.phasemodel.dims = [self.maxz, self.phaseres]
        if self.phaseres == 8:
            self.phasemodel.modelid = 1
        elif self.phaseres == 4:
            self.phasemodel.modelid = 0
        else:
            self.phasemodel.modelid = 2

        self.use_wrapper = True
        if platform.system() == "Linux":
            self.use_wrapper = False

        if self.use_wrapper:
            self.wrapper = IsoDecWrapper()
        else:
            self.wrapper = None

        self.reader = None

    def add_noise(self, noise_percent):
        """
        Add noise to the training and test data
        :param noise_percent: Percent of total data to add as noise
        :return: None
        """
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        nnoise = int(ltraining * noise_percent)
        print(f"Adding {nnoise} noise samples to training data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(nnoise):
            index = np.random.randint(0, ltraining - 1)
            centroid = self.training_data[2][index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.phaseres)
            emats.append(emat)
            zs.append(0)
            tempcentroids.append(centroid)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)
        self.training_data[2] = self.training_data[2] + tempcentroids

        nnoisetest = int(ltest * noise_percent)
        print(f"Adding {nnoisetest} noise samples to test data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(nnoisetest):
            index = np.random.randint(0, ltest - 1)
            centroid = self.test_data[2][index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.phaseres)
            emats.append(emat)
            zs.append(0)
            tempcentroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + tempcentroids

    def add_doubles(self, double_percent):
        """
        Add double peaks to the training and test data
        :param double_percent: Percent of total data to add as double peaks
        :return: None
        """
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        ndouble = int(ltraining * double_percent)
        print(f"Adding {ndouble} double samples in training data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(ndouble):
            index = np.random.randint(0, ltraining - 1)
            centroid1 = self.training_data[2][index]
            index2 = np.random.randint(0, ltraining - 1)
            centroid2 = self.training_data[2][index2]
            z = self.training_data[1][index]
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)
        self.training_data[2] = self.training_data[2] + tempcentroids

        ndoubletest = int(ltest * double_percent)
        print(f"Adding {ndoubletest} double samples in test data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(ndoubletest):
            index = np.random.randint(0, ltest - 1)
            index2 = np.random.randint(0, ltest - 1)
            centroid1 = self.test_data[2][index]
            centroid2 = self.test_data[2][index2]
            z = self.test_data[1][index]
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + tempcentroids

    def load_training_data(self, training_path, test_path=None, noise_percent=0.1, double_percent=0.1):
        """
        Load training data from a file
        :param training_path: Path to the training data file or name of the file tag
        :param test_path: Optional path to the test data file or name of the file tag.
                            If not, will default to same as training_path
        :param noise_percent: The percent of noise to add to the training and test data
        :param double_percent: The percent of double peaks to add to the training and test data
        :return: None
        """
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
        """
        Create the training and test dataloaders from a single file path
        :param training_path: Path to the training data file or name of the file tag
        :param test_path: Optional path to the test data file or name of the file tag. Default is same as training
        :param noise_percent: Percent of noise to add to the training and test data
        :param batchsize: Batch size for training
        :param double_percent: Percent of double peaks to add to the training and test data
        :return:
        """
        if batchsize is not None:
            self.config.batch_size = batchsize

        self.load_training_data(training_path, test_path=test_path, noise_percent=noise_percent,
                                double_percent=double_percent)

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.config.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.config.test_batch_size, shuffle=False,
                                          pin_memory=False)

    def create_merged_dataloader(self, dirs, training_path, noise_percent=0.1, batchsize=None, double_percent=0.1):
        """
        Create a merged dataloader from multiple directories. Looks for common file names and merges them together
        :param dirs: Directories to look in
        :param training_path: File name or tag, fed to load_training_data
        :param noise_percent: Percent of noise to add to the training and test data
        :param batchsize: Batch size for training
        :param double_percent: Percent of double peaks to add to the training and test data
        :return:
        """
        if batchsize is not None:
            self.config.batch_size = batchsize

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

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.config.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.config.test_batch_size, shuffle=False,
                                          pin_memory=True)

        print(f"Training Data Length: {len(self.training_data)}")
        print(f"Test Data Length: {len(self.test_data)}")

        # plot_zdist(self)

    def train_model(self, epochs=30, save=True, lossfn="crossentropy", forcenew=False):
        """
        Train the model
        :param epochs: Number of epochs
        :param save: Whether to save it. Default is True
        :param lossfn: Loss function, default is crossentropy. Options are crossentropy, weightedcrossentropy, focal
        :return: None
        """
        starttime = time.perf_counter()

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")
        self.phasemodel.get_class_weights(self.train_dataloader)
        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.phasemodel.train_model(self.train_dataloader, lossfn=lossfn, forcenew=forcenew)
            self.phasemodel.evaluate_model(self.test_dataloader)

        if save:
            self.phasemodel.save_model()
        print("Done! Time:", time.perf_counter() - starttime)

    def save_bad_data(self, filename="bad_data.pkl", maxbad=50):
        """
        Save bad data to a file. Evaluates the model, collects bad data, and saves it
        :param filename: Filename to save too. Default is bad_data.pkl
        :param maxbad: How many to save, default is 50
        :return: None
        """
        if self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")

        bad_data = self.phasemodel.evaluate_model(self.test_dataloader, savebad=True)
        output = []
        for b in bad_data:
            i1 = b[0]
            i2 = b[1]
            index = i1 * self.config.test_batch_size + i2
            if index > len(self.test_centroids):
                print("Index Error:", index, len(self.test_centroids), i1, i2, self.config.test_batch_size)
                continue
            centroid = self.test_centroids[index]
            output.append([centroid, b[2], b[3]])

        outindexes = np.random.randint(0, len(output), maxbad)
        output = [output[o] for o in outindexes]

        with open(filename, "wb") as f:
            pkl.dump(output, f)
            print(f"Saved {len(output)} bad data points to {filename}")

    def phase_predictor(self, centroids):
        """
        Predict the charge of a peak
        :param centroids: Set of centroid data for a peak with m/z in first column and intensity in second
        :return: Charge state, integer
        """
        if self.use_wrapper:
            z = IsoDecWrapper().predict_charge(centroids)
        else:
            z = self.phasemodel.predict(centroids)
        return int(z)

    def save_peak(self, centroids, z, peakmz, pks=None):
        """
        Save a peak to the MatchedCollection after its charge has been predicted.
        First, it optimizes the shift of the peak, then checks the matches and adds it to the MatchedCollection.
        :param centroids: Centroid data, m/z in first column, intensity in second
        :param z: Predicted charge
        :param peakmz: Peak m/z value
        :param pks: MatchedCollection peaks object
        :return: Matched indexes relative to the original centroid data
        """
        peaks = optimize_shift2(self.config, centroids, z, peakmz, tol=self.config.matchtol, maxshift=self.config.maxshift)
        if len(peaks) > 0:
            all_matched_indexes = []
            for m in peaks:
                for ind in m.matchedindexes:
                    if ind not in all_matched_indexes:
                        all_matched_indexes.append(ind)
                if pks is not None:
                    pks.add_peak(m)
                    #pks.add_pk_to_masses(m, 10)
                else:
                    self.pks.add_peak(m)
                    #self.pks.add_pk_to_masses(m, 10)
            return all_matched_indexes, peaks
        else:
            return [],[]

    def get_matches(self, centroids, z, peakmz, pks=None):
        """
        Get the matches for a peak
        :param centroids: Centroid data, m/z in first column, intensity in second
        :param z: Predicted charge
        :param peakmz: Peak m/z value
        :param pks: MatchedCollection peaks object
        :return: Indexes of matched peaks from the centroid data
        """
        if len(centroids) < self.config.minpeaks:
            return []
        if z == 0 or z > self.maxz:
            return []
        matchedindexes, peaks = self.save_peak(centroids, z, peakmz, pks=pks)
        return matchedindexes, peaks

    def batch_process_spectrum(self, data, window=None, threshold=0.001, centroided=False):
        """
        Process a spectrum and identify the peaks. It first identifies peak cluster, then predicts the charge,
        then checks the peaks. If all is good, it adds them to the MatchedCollection as a MatchedPeak object.

        :param data: Spectrum data, m/z in first column, intensity in second
        :param window: Window for peak selection
        :param threshold: Threshold for peak selection
        :param centroided: Whether the data is already centroided. If not, it will centroid it.
        :return: MatchedCollection of peaks
        """
        starttime = time.perf_counter()
        if window is None:
            window = self.config.peakwindow

        # TODO: Need a way to test for whether data is centroided already
        if centroided:
            centroids = data
        else:
            centroids = deepcopy(get_all_centroids(data, window=5, threshold=threshold * 0.1))

        if self.use_wrapper:
            self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.config)
        else:
            kwindow = window
            threshold = threshold
            for i in range(self.config.knockdown_rounds):

                if i <= 5:
                    self.config.css_thresh = 0.85
                else:
                    self.config.css_thresh = 0.75

                if i > 0:
                    kwindow = kwindow * 0.5
                self.config.current_KD_round = i
                peaks = fastpeakdetect(centroids, window=kwindow, threshold=threshold)
                # print("Knockdown:", i, "Peaks:", len(peaks))
                #print("Knockdown:", i, "Peaks:", peaks)
                if len(peaks) == 0:
                    break

                emats, peaks, centlist, indexes = encode_phase_all(centroids, peaks, lowmz=self.config.mzwindow[0],
                                                                   highmz=self.config.mzwindow[1], phaseres=self.phaseres)


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
                    matchedindexes, peaks = self.get_matches(centlist[j], z, p[0], pks=self.pks)


                    if len(matchedindexes) > 0:
                        ngood += 1
                        # Find matches
                        indval = indexes[j]
                        matchindvals = indval[matchedindexes]
                        # Knock them down
                        knockdown.extend(matchindvals)
                        #centroids = self.perform_modelling_kd(centroids, indval, peaks)
                self.config.acceptedclusters += ngood
                #print("NGood:", ngood)
                if len(knockdown) == 0:
                    continue
                knockdown = np.array(knockdown)
                centroids = np.delete(centroids, knockdown, axis=0)

                if len(centroids) < 3:
                    break
                #centroids = centroids[centroids[:, 1] > 0]
            # print("Time:", time.perf_counter() - starttime)
        return self.pks

    def perform_modelling_kd(self, centroids, indval, peaks):
        centroids = deepcopy(centroids)
        min_index = indval[0]
        max_ratio = 2

        full_kds = []
        partial_kd_lists = []

        for peak in peaks:
            partial_kds = []
            for i in range(len(peak.matchedindexes)):
                intensity_ratio = centroids[peak.matchedindexes[i], 1] / peak.isodist[peak.isomatches[i], 1]
                if intensity_ratio < max_ratio and peak.matchedindexes[i] not in full_kds:
                    full_kds.append(peak.matchedindexes[i])
                else:
                    partial_kds.append([peak.matchedindexes[i], intensity_ratio])
            partial_kd_lists.append(partial_kds)


        agg_partial_kds = []
        for i in range(len(partial_kd_lists[0])):
            ratios = [partial_kd_lists[0][i][1]]
            for j in range(len(peaks)):
                if j != 0:
                    for k in range(len(partial_kd_lists[j])):
                        if partial_kd_lists[j][k][0] == partial_kd_lists[0][i][0]:
                            ratios.append(partial_kd_lists[j][k][1])
            if len(ratios) == len(peaks):
                agg_partial_kds.append([partial_kd_lists[0][i][0], np.mean(np.array(ratios))])

        for i in range(len(agg_partial_kds)):
            print(peaks[0].centroids[agg_partial_kds[i][0], 0], agg_partial_kds[i][1])
            print(centroids[agg_partial_kds[i][0] + min_index, 0], centroids[agg_partial_kds[i][0] + min_index, 1])
            centroids[agg_partial_kds[i][0] + min_index, 1] *= (1 / agg_partial_kds[i][1])

        for i in range(len(full_kds)):
            centroids[full_kds[i] + min_index, 1] = 0

        return centroids



    def process_file(self, file, scans=None):
        starttime = time.perf_counter()
        self.config.filepath = file
        # Get importer and check it
        reader = ud.get_importer(file)
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

        #self.pks.save_pks()
        return reader

    def export_peaks(self, type="prosightlite", filename=None, reader=None, max_precursors=None):
        if filename is None:
            filename = "peaks.csv"

        if reader is None:
            reader = self.reader

        if type == "prosightlite":
            self.pks.export_prosightlite(filename)
        elif type == "msalign":
            self.pks.export_msalign(reader, filename, max_precursors=max_precursors)
        elif type == "pkl":
            self.pks.save_pks()
        else:
            raise ValueError("Unknown Export Type", type)


class IsoDecConfig:
    def __init__(self):
        """
        Configuration class for IsoDec Engine. Holds the key parameters for the data processing and deconvolution.
        """
        self.filepath = ""
        self.batch_size = 32
        self.test_batch_size = 2048
        self.peakwindow = 40
        self.acceptedclusters = 0
        self.matchtol = 0.005
        self.minpeaks = 3
        self.minmatchper = 0.67
        self.css_thresh = 0.80
        self.maxshift = 3  # This will get overwritten for smaller z, where it's dangerous to have more than 1 or 2
        self.mzwindow = [-1.5, 2.5]
        self.plusoneintwindow = [0.1, 0.6]
        self.knockdown_rounds = 10
        self.current_KD_round = 0
        self.activescan = -1
        self.activescanrt = -1
        self.activescanorder = -1

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


if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecEngine(phaseres=8)
    topdirectory = "C:\\Data\\IsoNN\\training"

    dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]
    eng.create_merged_dataloader(dirs, "phase83", noise_percent=0.0, batchsize=32, double_percent=0.4)
    # eng.train_model(epochs=3)
    eng.train_model(epochs=10, lossfn="crossentropy", forcenew=True)
    #eng.train_model(epochs=3, lossfn="focal", forcenew=False)

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

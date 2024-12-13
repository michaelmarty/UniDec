import time
import numpy as np
from itertools import chain
import os
import torch
from torch.utils.data import DataLoader

from unidec.IsoDec.models import example, PhaseModel
from unidec.IsoDec.datatools import fastpeakdetect, get_all_centroids, fastnearest, check_spacings, remove_noise_cdata
from unidec.IsoDec.match import optimize_shift2, IsoDecConfig, MatchedCollection
from unidec.IsoDec.encoding import data_dirs, encode_noise, encode_phase_all, small_data_dirs, \
    encode_double, encode_harmonic, extract_centroids

from copy import deepcopy
import pickle as pkl

from unidec.IsoDec.c_interface import IsoDecWrapper
from unidec.IsoDec.plots import cplot
import platform
from unidec.IsoDec.altdecon import thrash_predict

from unidec.UniDecImporter.ImporterFactory import *


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

# TODO: Inherit this from IsoDecRuntime
class IsoDecEngine:
    """
    Main class for IsoDec Engine
    """

    def __init__(self, phaseres=8, verbose=False, use_wrapper=False):
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

        self.config.phaseres = phaseres
        self.phasemodel = PhaseModel()
        self.maxz = 50
        self.phasemodel.dims = [self.maxz, self.config.phaseres]
        if self.config.phaseres == 8:
            self.phasemodel.modelid = 1
        elif self.config.phaseres == 4:
            self.phasemodel.modelid = 0
        else:
            self.phasemodel.modelid = 2

        self.use_wrapper = use_wrapper
        if platform.system() == "Linux":
            self.use_wrapper = False

        if self.use_wrapper:
            self.wrapper = IsoDecWrapper()
        else:
            self.wrapper = None

        self.reader = None
        self.predmode = 0

    def drop_ones(self, percentage=0.8):
        """
        Drop 80% of the training data with charge 1
        :param percentage: Percentage of data to drop
        :return: None
        """
        print("Dropping charge 1 data:", percentage)
        z = self.training_data[1]
        emat = self.training_data[0]
        centroids = self.training_data[2]
        keep = []
        for i in range(len(z)):
            if z[i] != 1 or np.random.rand() < percentage:
                keep.append(i)
        self.training_data = [emat[keep], z[keep], [centroids[i] for i in keep]]
        print("New Length:", len(self.training_data[0]))

    def add_harmonics(self, harmonic_percent=0.4):
        """
        Add harmonics at 2x charge to the training and test data
        :param harmonic_percent: Percent of total data to add as harmonics
        :return: None
        """
        ltraining = len(self.training_data[0])
        ltest = len(self.test_data[0])
        nharm = int(ltraining * harmonic_percent)
        print(f"Adding {nharm} harmonic samples to training data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(nharm):
            index = np.random.randint(0, ltraining - 1)
            centroid = self.training_data[2][index]
            z = self.training_data[1][index]
            emat, centroid2 = encode_harmonic(centroid, z, phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid2)

        self.training_data[0] = np.concatenate((self.training_data[0], emats), axis=0)
        self.training_data[1] = np.concatenate((self.training_data[1], zs), axis=0)
        self.training_data[2] = self.training_data[2] + tempcentroids

        nharmtest = int(ltest * harmonic_percent)
        print(f"Adding {nharmtest} harmonic samples to test data")
        emats = []
        zs = []
        tempcentroids = []
        for i in range(nharmtest):
            index = np.random.randint(0, ltest - 1)
            centroid = self.test_data[2][index]
            z = self.test_data[1][index]
            emat, centroid2 = encode_harmonic(centroid, z, phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid2)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + tempcentroids

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
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.config.phaseres)
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
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.config.phaseres)
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
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.config.phaseres)
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
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid)

        self.test_data[0] = np.concatenate((self.test_data[0], emats), axis=0)
        self.test_data[1] = np.concatenate((self.test_data[1], zs), axis=0)
        self.test_data[2] = self.test_data[2] + tempcentroids

    def load_training_data(self, training_path, test_path=None, noise_percent=0.0, double_percent=0.4,
                           harmonic_percent=0.0, onedrop_percent=0.0):
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

        if onedrop_percent > 0:
            self.drop_ones(percentage=onedrop_percent)

        if harmonic_percent > 0:
            self.add_harmonics(harmonic_percent)

        if double_percent > 0:
            self.add_doubles(double_percent)

        if noise_percent > 0:
            self.add_noise(noise_percent)

        print("Loaded:", len(self.training_data[0]), "Training Samples")

    def create_training_dataloader(self, training_path, test_path=None, noise_percent=0, batchsize=None,
                                   double_percent=0.4, harmonic_percent=0, one_drop_percent=0):
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
                                double_percent=double_percent, harmonic_percent=harmonic_percent,
                                onedrop_percent=one_drop_percent)

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.config.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.config.test_batch_size, shuffle=False,
                                          pin_memory=False)

    def create_merged_dataloader(self, dirs, training_path, noise_percent=0.0, batchsize=None, double_percent=0.4,
                                 harmonic_percent=0.0, onedrop_percent=0.0):
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
            self.load_training_data(training_path, noise_percent=noise_percent, double_percent=double_percent,
                                    harmonic_percent=harmonic_percent, onedrop_percent=onedrop_percent)
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
        return [int(z), 0]

    def thrash_predictor(self, centroids):
        return [thrash_predict(centroids), 1]

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
        pk = optimize_shift2(self.config, centroids, z, peakmz)
        if pk is not None:
            if pks is not None:
                pk.rt = self.config.activescanrt
                pk.scan = self.config.activescan
                pks.add_peak(pk)
                pks.add_pk_to_masses(pk, 10)
            else:
                self.pks.add_peak(pk)
            return pk.matchedindexes
        else:
            return []

    def get_matches_multiple_z(self, centroids, zs, peakmz, pks=None):
        if len(centroids) < self.config.minpeaks:
            return []
        if zs[0] == 0 or zs[0] > self.maxz:
            return []
        pk1 = optimize_shift2(self.config, centroids, zs[0], peakmz)
        pk2 = optimize_shift2(self.config, centroids, zs[1], peakmz)
        if pk1 is not None and pk2 is not None:
            # Retain the peak with the highest score
            pk1_maxscore = np.amax(pk1.acceptedshifts[:, 1])
            pk2_maxscore = np.amax(pk2.acceptedshifts[:, 1])
            if self.config.verbose:
                print("Pk1 score:", pk1_maxscore, "Pk2 score:", pk2_maxscore)
            if pk1_maxscore > pk2_maxscore:
                if pks is not None:
                    pks.add_peak(pk1)
                else:
                    self.pks.add_peak(pk1)
                return pk1.matchedindexes
            else:
                if pks is not None:
                    pks.add_peak(pk2)
                else:
                    self.pks.add_peak(pk2)
                return pk2.matchedindexes
        elif pk1 is not None and pk2 is None:
            if pks is not None:
                pks.add_peak(pk1)
            else:
                self.pks.add_peak(pk1)
            return pk1.matchedindexes
        elif pk1 is None and pk2 is not None:
            if pks is not None:
                pks.add_peak(pk2)
            else:
                self.pks.add_peak(pk2)
            return pk2.matchedindexes
        else:
            return []

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
        if self.config.verbose:
            print("Processing spectrum with prediction mode:", self.predmode)
        starttime = time.perf_counter()
        if window is None:
            window = self.config.peakwindow
        if threshold is None:
            threshold = self.config.peakthresh
        if self.config.css_thresh < 0.6:
            self.config.adjust_css = False

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

        if self.use_wrapper:
            self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.config)
        else:
            kwindow = window
            threshold = threshold
            for i in range(self.config.knockdown_rounds):
                # Adjust settings based on round
                if i >= 5 and self.config.adjust_css:
                    self.config.css_thresh = self.config.css_thresh * 0.90
                    if self.config.css_thresh < 0.6:
                        self.config.css_thresh = 0.6
                if self.config.verbose:
                    print("Spectrum length: ", len(centroids))
                if i > 0:
                    kwindow = kwindow * 0.5
                    if kwindow < 1:
                        kwindow = 1
                    threshold = threshold * 0.5
                    if threshold < 0.000001:
                        threshold = 0.000001
                self.config.current_KD_round = i

                # Pick peaks
                peaks = fastpeakdetect(centroids, window=int(kwindow), threshold=threshold)
                # print("Knockdown:", i, "Peaks:", len(peaks))
                if self.config.verbose:
                    print("\n\nKnockdown:", i, "NPeaks:", len(peaks), "Peaks:", peaks[:, 0], kwindow)
                if len(peaks) == 0:
                    break

                if self.predmode == 0:
                    # Encode phase of all
                    emats, peaks, centlist, indexes = encode_phase_all(centroids, peaks, lowmz=self.config.mzwindow[0],
                                                                       highmz=self.config.mzwindow[1],
                                                                       phaseres=self.config.phaseres,
                                                                       minpeaks=2, datathresh=self.config.datathreshold)


                    emats = [torch.as_tensor(e, dtype=torch.float32) for e in emats]
                    # emats = torch.as_tensor(emats, dtype=torch.float32).to(self.phasemodel.device)
                    data_loader = DataLoader(emats, batch_size=2048, shuffle=False, pin_memory=True)

                    # Predict Charge
                    preds = self.phasemodel.batch_predict(data_loader)
                elif self.predmode == 1:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindow[
                                                                                                0],
                                                                                            highmz=self.config.mzwindow[
                                                                                                1],
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.phase_predictor(c) for c in encodingcentroids]
                elif self.predmode == 2:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindow[
                                                                                                0],
                                                                                            highmz=self.config.mzwindow[
                                                                                                1],
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.thrash_predictor(c) for c in encodingcentroids]
                elif self.predmode == 3:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindow[
                                                                                                0],
                                                                                            highmz=self.config.mzwindow[
                                                                                                1],
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.phase_predictor(c) for c in encodingcentroids]
                    preds2 = [self.thrash_predictor(c) for c in encodingcentroids]
                    preds = [[preds[i][0], preds2[i][0]] for i in range(len(preds))]
                else:
                    raise ValueError("Unknown mode", self.predmode)

                knockdown = []
                ngood = 0


                # print(peaks, len(peaks))
                # Loop through all peaks to check if they are good
                for j, p in enumerate(peaks):
                    z = preds[j]
                    kindex = fastnearest(centroids[:, 0], p[0])

                    if self.config.verbose:
                        print("Peak:", p, z)

                    if kindex in knockdown:
                        continue
                    if z[0] == 0:
                        knockdown.append(kindex)
                        continue

                    # Get the centroids around the peak
                    if z[1] != 0 and z[1] != z[0]:
                        matchedindexes = self.get_matches_multiple_z(centlist[j], z, p[0], pks=self.pks)
                    else:
                        matchedindexes = self.get_matches(centlist[j], z[0], p[0], pks=self.pks)

                    if len(matchedindexes) > 0:
                        ngood += 1
                        # Find matches
                        indval = indexes[j]
                        matchindvals = indval[matchedindexes]
                        # Knock them down
                        knockdown.extend(matchindvals)
                    else:
                        knockdown.append(kindex)
                if len(knockdown) == 0:
                    continue

                knockdown = np.array(knockdown)
                centroids = np.delete(centroids, knockdown, axis=0)

                if len(centroids) < self.config.minpeaks:
                    break

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
                    spectrum = reader.get_single_scan(s)
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

    def export_peaks(self, type="prosightlite", filename="output", reader=None, act_type="HCD", max_precursors=None):
        if filename is None:
            filename = "peaks.csv"

        if reader is None:
            reader = self.reader

        if type == "prosightlite":
            self.pks.export_prosightlite(filename)
        elif type == "msalign":
            self.pks.export_msalign(self.config, reader, filename, act_type=act_type, max_precursors=max_precursors)
        elif type == "pkl":
            self.pks.save_pks()
        else:
            raise ValueError("Unknown Export Type", type)


if __name__ == "__main__":
    starttime = time.perf_counter()
    eng = IsoDecEngine(phaseres=8)
    topdirectory = "C:\\Data\\IsoNN\\training"

    dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]
    eng.create_merged_dataloader(dirs, "phase83", noise_percent=0.0, batchsize=32, double_percent=0.4,
                                 harmonic_percent=0.1, onedrop_percent=0.8)
    # eng.train_model(epochs=3)
    eng.train_model(epochs=10, lossfn="crossentropy", forcenew=True)
    # eng.train_model(epochs=3, lossfn="focal", forcenew=False)

    # eng.create_merged_dataloader([os.path.join(topdirectory, small_data_dirs[2])], "phase82", noise_percent=0.2,
    #                             batchsize=32, double_percent=0.2)
    # eng.save_bad_data()
    exit()
    import matplotlib.pyplot as plt
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

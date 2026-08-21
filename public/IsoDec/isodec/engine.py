import time
import numpy as np
from itertools import chain
import os
import sys
from pathlib import Path

import torch
from torch.utils.data import DataLoader

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "isodec"

from .encoding import encode_synthetic
from .models import example, PhaseModel
from .datatools import fastpeakdetect, get_all_centroids, fastnearest, check_spacings, remove_noise_cdata
from .match import optimize_shift2, MatchedCollection
from .config import IsoDecConfig
from .encoding import data_dirs, encode_noise, encode_phase_all, small_data_dirs, \
    encode_double, encode_harmonic, extract_centroids

from copy import deepcopy
import pickle as pkl

from .c_interface import IsoDecWrapper
from .plots import *
import platform
from .altdecon import thrash_predict

from .io import ImporterFactory
from ._version import __version__


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


# Note, this is primarily a tool for training, testing, and development.
# We recommend most people use IsoDecRuntime for routine use.


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
        self.version = __version__

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
        self.ematkey="emat"
        self.centroidkey="centroids"

    def _get_augmentation_data(self, dataset):
        """Return the mutable data collection selected for augmentation."""
        if dataset == "training":
            data = self.training_data
        elif dataset == "test":
            data = self.test_data
        else:
            raise ValueError("dataset must be 'training' or 'test'")

        if data is None:
            raise ValueError(f"{dataset.capitalize()} data has not been loaded")
        return data

    def drop_ones(self, percentage=0.8, dataset="training"):
        """
        Drop charge-1 samples from the selected data collection.
        :param percentage: Percentage of data to drop
        :param dataset: Data collection to modify: ``"training"`` or ``"test"``
        :return: None
        """
        data = self._get_augmentation_data(dataset)
        print(f"Dropping charge 1 data from {dataset} data:", percentage)
        z = data[1]
        emat = data[0]
        centroids = data[2]
        keep = []
        for i in range(len(z)):
            if z[i] != 1 or np.random.rand() < percentage:
                keep.append(i)
        data[:] = [emat[keep], z[keep], [centroids[i] for i in keep]]
        print("New Length:", len(data[0]))

    def add_harmonics(self, harmonic_percent=0.4, dataset="training"):
        """
        Add harmonics at 2x charge to the selected data collection.
        :param harmonic_percent: Percent of total data to add as harmonics
        :param dataset: Data collection to modify: ``"training"`` or ``"test"``
        :return: None
        """
        data = self._get_augmentation_data(dataset)
        data_length = len(data[0])
        nharm = int(data_length * harmonic_percent)
        print(f"Adding {nharm} harmonic samples to {dataset} data")
        if nharm == 0:
            return

        emats = []
        zs = []
        tempcentroids = []
        for _ in range(nharm):
            index = np.random.randint(data_length)
            centroid = data[2][index]
            z = data[1][index]
            emat, centroid2 = encode_harmonic(centroid, z, phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid2)

        data[0] = np.concatenate((data[0], emats), axis=0)
        data[1] = np.concatenate((data[1], zs), axis=0)
        data[2] = data[2] + tempcentroids

    def add_noise(self, noise_percent, dataset="training"):
        """
        Add noise to the selected data collection.
        :param noise_percent: Percent of total data to add as noise
        :param dataset: Data collection to modify: ``"training"`` or ``"test"``
        :return: None
        """
        data = self._get_augmentation_data(dataset)
        data_length = len(data[0])
        nnoise = int(data_length * noise_percent)
        print(f"Adding {nnoise} noise samples to {dataset} data")
        if nnoise == 0:
            return

        emats = []
        zs = []
        tempcentroids = []
        for _ in range(nnoise):
            index = np.random.randint(data_length)
            centroid = data[2][index]
            emat, centroid, z = encode_noise(centroid[0, 0], np.amax(centroid[:, 1]), phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(0)
            tempcentroids.append(centroid)

        data[0] = np.concatenate((data[0], emats), axis=0)
        data[1] = np.concatenate((data[1], zs), axis=0)
        data[2] = data[2] + tempcentroids

    def add_doubles(self, double_percent, dataset="training"):
        """
        Add double peaks to the selected data collection.
        :param double_percent: Percent of total data to add as double peaks
        :param dataset: Data collection to modify: ``"training"`` or ``"test"``
        :return: None
        """
        data = self._get_augmentation_data(dataset)
        data_length = len(data[0])
        ndouble = int(data_length * double_percent)
        print(f"Adding {ndouble} double samples in {dataset} data")
        if ndouble == 0:
            return

        emats = []
        zs = []
        tempcentroids = []
        for _ in range(ndouble):
            index = np.random.randint(data_length)
            centroid1 = data[2][index]
            index2 = np.random.randint(data_length)
            centroid2 = data[2][index2]
            z = data[1][index]
            emat, centroid = encode_double(centroid1, centroid2, phaseres=self.config.phaseres)
            emats.append(emat)
            zs.append(z)
            tempcentroids.append(centroid)

        data[0] = np.concatenate((data[0], emats), axis=0)
        data[1] = np.concatenate((data[1], zs), axis=0)
        data[2] = data[2] + tempcentroids

    def equilize_data(self, dataset="training"):
        """
        Equalize charge-state representation in the selected data collection.
        :param dataset: Data collection to modify: ``"training"`` or ``"test"``
        :return: None
        """
        data = self._get_augmentation_data(dataset)
        data_length = len(data[0])
        if data_length == 0:
            return

        zdata = data[1]
        zcounts = np.bincount(zdata)[1:]
        thebar = np.median(zcounts)
        uz = np.arange(1, len(zcounts) + 1)
        print(f"Equilizing {dataset} data. The Bar:", thebar)

        if len(uz) != len(zcounts):
            raise ValueError("Unique Zs and Z Counts are not the same length")

        if thebar <= 1:
            thebar = 10

        emats = []
        zs = []
        tempcentroids = []
        for i in range(1, 50):
            if i not in uz:
                num = 0
            else:
                num = zcounts[np.where(uz == i)[0][0]]

            if num < thebar:
                n = int(thebar - num)
                for _ in range(n):
                    index = np.random.randint(data_length)
                    centroid1 = data[2][index]
                    z = data[1][index]
                    emat, centroid = encode_synthetic(centroid1, z, i, phaseres=self.config.phaseres)
                    emats.append(emat)
                    zs.append(i)
                    tempcentroids.append(centroid)

        if not emats:
            return

        data[0] = np.concatenate((data[0], emats), axis=0)
        data[1] = np.concatenate((data[1], zs), axis=0)
        data[2] = data[2] + tempcentroids

    def load_data(self, data_path, dataset="training", noise_percent=0.0, double_percent=0.4,
                  harmonic_percent=0.0, onedrop_percent=0.0, equilize=False,
                  ematkey="emat", centroidkey="centroids"):
        """
        Load and augment one training or test data collection from a file.
        :param data_path: Path to an NPZ file or its file tag
        :param dataset: Data collection to load: ``"training"`` or ``"test"``
        :param noise_percent: The percent of noise samples to add
        :param double_percent: The percent of double-peak samples to add
        :param harmonic_percent: The percent of harmonic samples to add
        :param onedrop_percent: The percent of charge 1 data to drop
        :param equilize: Whether to equilize the data
        :param ematkey: The key for the encoded matrix in the npz file
        :param centroidkey: The key for the centroids in the npz file
        :return: None
        """
        if dataset not in ("training", "test"):
            raise ValueError("dataset must be 'training' or 'test'")

        data_path = os.fspath(data_path)
        ext = ".npz"
        if not data_path.lower().endswith(ext):
            directory, tag = os.path.split(data_path)
            data_path = os.path.join(directory, f"{dataset}_data_{tag}{ext}")

        with np.load(data_path, allow_pickle=True) as loaded:
            if centroidkey in ("centroids", "profiles"):
                data = [loaded[ematkey], loaded["z"], list(loaded[centroidkey])]
            elif centroidkey == "both":
                centroid_ematkey = ematkey.replace("emat_p_", "emat_c_").replace("emat_b_", "emat_c_")
                profile_ematkey = centroid_ematkey.replace("emat_c_", "emat_p_")
                centroid_data = [loaded[centroid_ematkey], loaded["z"], list(loaded["centroids"])]
                profile_data = [loaded[profile_ematkey], loaded["z"], list(loaded["profiles"])]
                data = [
                    np.concatenate((centroid_data[0], profile_data[0]), axis=0),
                    np.concatenate((centroid_data[1], profile_data[1]), axis=0),
                    centroid_data[2] + profile_data[2],
                ]
            else:
                raise ValueError("centroidkey must be 'centroids', 'profiles', or 'both'")

        if dataset == "training":
            self.training_data = data
        else:
            self.test_data = data

        if onedrop_percent > 0:
            self.drop_ones(percentage=onedrop_percent, dataset=dataset)

        if harmonic_percent > 0:
            self.add_harmonics(harmonic_percent, dataset=dataset)

        if double_percent > 0:
            self.add_doubles(double_percent, dataset=dataset)

        if noise_percent > 0:
            self.add_noise(noise_percent, dataset=dataset)

        if equilize:
            self.equilize_data(dataset=dataset)

        print(f"Loaded: {len(data[0])} {dataset.capitalize()} Samples")

    def create_training_dataloader(self, training_path, test_path=None, noise_percent=0, batchsize=None,
                                   double_percent=0.4, harmonic_percent=0, one_drop_percent=0, equalize=False,
                                   ematkey="emat", centroidkey="centroids"):
        """
        Create the training and test dataloaders from a single file path
        :param training_path: Path to the training data file or name of the file tag
        :param test_path: Optional path to the test data file or name of the file tag. Default is same as training
        :param noise_percent: Percent of noise to add to the training and test data
        :param batchsize: Batch size for training
        :param double_percent: Percent of double peaks to add to the training and test data
        :param harmonic_percent: Percent of harmonic peaks to add to the training and test data
        :param one_drop_percent: Percent of charge 1 data to drop
        :param equalize: Whether to equilize the data
        :return:
        """
        if batchsize is not None:
            self.config.batch_size = batchsize

        self.load_data(training_path, dataset="training", noise_percent=noise_percent,
                       double_percent=double_percent, harmonic_percent=harmonic_percent,
                       ematkey=ematkey, centroidkey=centroidkey, onedrop_percent=one_drop_percent,
                       equilize=equalize)

        if test_path is None:
            test_path = os.fspath(training_path)
            if test_path.lower().endswith(".npz"):
                directory, filename = os.path.split(test_path)
                test_path = os.path.join(directory, filename.replace("training", "test", 1))

        self.load_data(test_path, dataset="test", noise_percent=noise_percent,
                       double_percent=double_percent, harmonic_percent=harmonic_percent,
                       ematkey=ematkey, centroidkey=centroidkey, onedrop_percent=one_drop_percent,
                       equilize=equalize)

        self.training_data = IsoDecDataset(self.training_data[0], self.training_data[1])
        self.test_data = IsoDecDataset(self.test_data[0], self.test_data[1])

        self.train_dataloader = DataLoader(self.training_data, batch_size=self.config.batch_size, shuffle=True,
                                           pin_memory=True)
        self.test_dataloader = DataLoader(self.test_data, batch_size=self.config.test_batch_size, shuffle=False,
                                          pin_memory=False)

    def create_merged_dataloader(self, dirs, training_path, noise_percent=0.0, batchsize=None, double_percent=0.4,
                                 harmonic_percent=0.0, onedrop_percent=0.0, equilize=False, ematkey="emat",
                                 centroidkey="centroids", dataset="both"):
        """
        Create a merged dataloader from multiple directories. Looks for common file names and merges them together
        :param dirs: Directories to look in
        :param training_path: File name or tag, fed to load_data
        :param noise_percent: Percent of noise to add to the training and test data
        :param batchsize: Batch size for training
        :param double_percent: Percent of double peaks to add to the training and test data
        :param harmonic_percent: Percent of harmonic peaks to add to the training and test data
        :param onedrop_percent: Percent of charge 1 data to drop
        :param equilize: Whether to equilize the data
        :param dataset: Data to import: ``"training"``, ``"test"``, or ``"both"``
        :return:
        """
        if dataset not in ("training", "test", "both"):
            raise ValueError("dataset must be 'training', 'test', or 'both'")

        if batchsize is not None:
            self.config.batch_size = batchsize

        training_data = []
        test_data = []
        load_training = dataset in ("training", "both")
        load_test = dataset in ("test", "both")

        self.ematkey = ematkey
        self.centroidkey = centroidkey
        self.training_data = None
        self.test_data = None
        self.training_centroids = []
        self.test_centroids = []
        self.train_dataloader = None
        self.test_dataloader = None

        for d in dirs:
            os.chdir(d)
            print(d)
            if load_training:
                self.load_data(training_path, dataset="training", noise_percent=noise_percent,
                               double_percent=double_percent, harmonic_percent=harmonic_percent,
                               onedrop_percent=onedrop_percent, ematkey=ematkey,
                               centroidkey=centroidkey, equilize=equilize)
                training_data.append(self.training_data)

            if load_test:
                test_path = os.fspath(training_path)
                if test_path.lower().endswith(".npz"):
                    directory, filename = os.path.split(test_path)
                    test_path = os.path.join(directory, filename.replace("training", "test", 1))
                self.load_data(test_path, dataset="test", noise_percent=noise_percent,
                               double_percent=double_percent, harmonic_percent=harmonic_percent,
                               onedrop_percent=onedrop_percent, ematkey=ematkey,
                               centroidkey=centroidkey, equilize=equilize)
                test_data.append(self.test_data)

        if load_training:
            merged_training = [np.concatenate([data[0] for data in training_data], axis=0),
                               np.concatenate([data[1] for data in training_data], axis=0)]
            self.training_centroids = list(chain(*[data[2] for data in training_data]))
            self.training_data = IsoDecDataset(merged_training[0], merged_training[1])
            self.train_dataloader = DataLoader(
                self.training_data, batch_size=self.config.batch_size, shuffle=True, pin_memory=True
            )
            print(f"Training Data Length: {len(self.training_data)}")

        if load_test:
            merged_test = [np.concatenate([data[0] for data in test_data], axis=0),
                           np.concatenate([data[1] for data in test_data], axis=0)]
            self.test_centroids = list(chain(*[data[2] for data in test_data]))
            self.test_data = IsoDecDataset(merged_test[0], merged_test[1])
            self.test_dataloader = DataLoader(
                self.test_data, batch_size=self.config.test_batch_size, shuffle=False, pin_memory=True
            )
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

        accuracy_score=0

        if self.train_dataloader is None or self.test_dataloader is None:
            raise ValueError("DataLoaders not created. Run create_training_dataloader first.")
        self.phasemodel.get_class_weights(self.train_dataloader)
        for t in range(epochs):
            print(f"Epoch {t + 1}\n-------------------------------")
            self.phasemodel.train_model(self.train_dataloader, lossfn=lossfn, forcenew=forcenew)
            _, accuracy_score = self.phasemodel.evaluate_model(self.test_dataloader)

        if self.ematkey != "emat":
            output = self.ematkey.replace("emat_", "")
            savepath = f"phase_model_{output}.pth"
        else:
            savepath = None

        if save:
            self.phasemodel.save_model(savepath)
        print("Done! Time:", time.perf_counter() - starttime)
        return accuracy_score

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
            pk.matchedcentroids = centroids[pk.matchedindexes]
            if pks is not None:
                pk.rt = self.config.activescanrt
                pk.scan = self.config.activescan
                pks.add_peak(pk)
                pks.add_pk_to_masses(pk, self.config)
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
        if pk1 is not None:
            pk1.matchedcentroids = centroids[pk1.matchedindexes]
        pk2 = optimize_shift2(self.config, centroids, zs[1], peakmz)
        if pk2 is not None:
            pk2.matchedcentroids = centroids[pk2.matchedindexes]
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

    def get_matches_zloop(self, centroids, predvec, peakmz, pks=None):
        if len(centroids) < self.config.minpeaks:
            return []
        # Order the predvec by descending
        order = np.argsort(predvec[0])[::-1]
        # for i in order:
        #     print("Z:", i, "Prediction:", predvec[0][i])
        matchedindexes = []
        css_scores = []
        peaks = []
        for i in order:
            if self.config.verbose:
                print("Z:", i, "Prediction:", predvec[0][i])
            if i == 0:
                continue
            # if i == 1:
            #     #Charge 1 must be within 90% of the max score
            #     if predvec[0][i] < 0.5 * predvec[0][order[0]]:
            #         continue
            # #Check if the prediction score is within a threshold of the max score
            # elif i != 1 and predvec[0][i] < self.config.zscore_threshold * predvec[0][order[0]]:
            #     break

            pk = optimize_shift2(self.config, centroids, i, peakmz)
            if pk is not None:
                # Add unique matched indexes to the list
                matchedindexes.append(pk.matchedindexes)
                # Get the max score
                maxscore = np.amax(pk.acceptedshifts[:, 1])
                css_scores.append(maxscore)
                peaks.append(pk)
                # if pks is not None:
                #     pks.add_peak(pk)
                # else:
                #     self.pks.add_peak(pk)

        if len(peaks) > 0:
            # Find the peak with the highest score
            maxindex = np.argmax(css_scores)
            if pks is not None:
                pks.add_peak(peaks[maxindex])
            else:
                self.pks.add_peak(peaks[maxindex])
            return matchedindexes[maxindex]
        else:
            return []

    def batch_process_spectrum(self, data, type=None, window=None, threshold=None, centroided=False, refresh=False):
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
            self.pks = self.wrapper.process_spectrum(centroids, self.pks, self.config, type)
        else:
            kwindow = window
            threshold = threshold
            for i in range(self.config.knockdown_rounds):
                # Adjust settings based on round
                if i >= 5:
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
                    emats, peaks, centlist, indexes = encode_phase_all(centroids, peaks, lowmz=self.config.mzwindowlb,
                                                                       highmz=self.config.mzwindowub,
                                                                       phaseres=self.config.phaseres,
                                                                       minpeaks=2, datathresh=self.config.datathreshold)

                    emats = [torch.as_tensor(e, dtype=torch.float32) for e in emats]
                    # emats = torch.as_tensor(emats, dtype=torch.float32).to(self.phasemodel.device)
                    data_loader = DataLoader(emats, batch_size=2048, shuffle=False, pin_memory=True)

                    # Predict Charge
                    preds = self.phasemodel.batch_predict(data_loader)
                elif self.predmode == 1:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindowlb,
                                                                                            highmz=self.config.mzwindowub,
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.phase_predictor(c) for c in encodingcentroids]
                elif self.predmode == 2:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindowlb,
                                                                                            highmz=self.config.mzwindowub,
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.thrash_predictor(c) for c in encodingcentroids]
                elif self.predmode == 3:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindowlb,
                                                                                            highmz=self.config.mzwindowub,
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.phase_predictor(c) for c in encodingcentroids]
                    preds2 = [self.thrash_predictor(c) for c in encodingcentroids]
                    preds = [[preds[i][0], preds2[i][0]] for i in range(len(preds))]

                elif self.predmode == 4:
                    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(centroids, peaks,
                                                                                            lowmz=self.config.mzwindowlb,
                                                                                            highmz=self.config.mzwindowub,
                                                                                            minpeaks=2,
                                                                                            datathresh=self.config.datathreshold)
                    peaks = goodpeaks
                    centlist = outcentroids
                    preds = [self.phasemodel.predict_returnvec(c) for c in encodingcentroids]


                else:
                    raise ValueError("Unknown mode", self.predmode)

                knockdown = []
                ngood = 0

                # print(peaks, len(peaks))
                # Loop through all peaks to check if they are good
                for j, p in enumerate(peaks):
                    z = preds[j]
                    if self.predmode == 4:
                        z = 0
                    kindex = fastnearest(centroids[:, 0], p[0])

                    matchedindexes = []

                    if self.config.verbose:
                        print("Peak:", p, z)

                    if kindex in knockdown:
                        continue

                    if self.predmode != 4:
                        if z[0] == 0:
                            knockdown.append(kindex)
                            continue
                        # Get the centroids around the peak
                        if z[1] != 0 and z[1] != z[0]:
                            matchedindexes = self.get_matches_multiple_z(centlist[j], z, p[0], pks=self.pks)
                        else:
                            matchedindexes = self.get_matches(centlist[j], z[0], p[0], pks=self.pks)
                    elif self.predmode == 4:
                        matchedindexes = self.get_matches_zloop(centlist[j], preds[j], p[0], pks=self.pks)

                    if len(matchedindexes) > 0:
                        ngood += 1
                        # Find matches
                        indval = indexes[j]
                        matchindvals = indval[matchedindexes]
                        self.pks.peaks[-1].matchedindexes = np.array(matchindvals)
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

    def export_peaks(self, type="prosightlite", filename="output", reader=None, act_type="HCD", max_precursors=1):
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

    topdirectory = "C:\\Data\\IsoNN\\training"
    # topdirectory = "Z:\\Group Share\\JGP\\TrainingData"

    dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]

    # eng = IsoDecEngine(phaseres=8)
    # # eng.create_merged_dataloader(dirs, "phase83", noise_percent=0.0, batchsize=32, double_percent=0,
    # #                              harmonic_percent=0, onedrop_percent=0, equilize=True)
    # eng.create_merged_dataloader(dirs, "phase83", noise_percent=0.0, batchsize=32, double_percent=0,
    #                              harmonic_percent=0, onedrop_percent=0.9, equilize=True)
    # # eng.train_model(epochs=3)
    # eng.train_model(epochs=10, lossfn="crossentropy", forcenew=True)
    # eng.train_model(epochs=3, lossfn="focal", forcenew=False)

    # eng.create_merged_dataloader([os.path.join(topdirectory, small_data_dirs[2])], "phase82", noise_percent=0.2,
    #                             batchsize=32, double_percent=0.2)
    # eng.save_bad_data()

    eng = IsoDecEngine(phaseres=1)
    # eng.create_merged_dataloader(dirs, "phase83", noise_percent=0.0, batchsize=32, double_percent=0,
    #                              harmonic_percent=0, onedrop_percent=0, equilize=True)
    eng.create_merged_dataloader(dirs, "phase1", noise_percent=0.0, batchsize=32, double_percent=0,
                                 harmonic_percent=0, onedrop_percent=0.9, equilize=True)
    # eng.train_model(epochs=3)
    eng.train_model(epochs=3, lossfn="crossentropy", forcenew=True)

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

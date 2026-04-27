import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenc import fft_pep_mass_to_dist
import random


class IsoGenMassEngine(IsoGenEngineBase):
    def __init__(self, isolen=128):
        super().__init__()
        self.veclen = 5
        if isolen > 128:
            modelid = 1
        else:
            modelid = 0
        self.model = IsoGenModelBase(isolen=isolen, savename="isogenmass_model_", vectorlen=self.veclen, modelid=modelid)
        self.isolen = isolen

        self.massranges = np.array(
            [[10, 1200], [400, 12000], [5000, 60000], [8000, 120000], [80000, 1100000]]
        )
        self.lengths = np.array([8, 32, 64, 128, 1024])

        self.models = []
        for l in self.lengths:
            if l > 128:
                modelid = 1
            else:
                modelid = 0
            model = IsoGenModelBase(isolen=l, savename="isogenmass_model_", vectorlen=self.veclen, modelid=modelid)
            self.models.append(model)

    def train(self, masses=None, dists=None, n=100000, epochs=10, length=128, forcenew=False):
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
        except:
            print("Length not found in list")
            return
        massrange = self.massranges[modelindex]
        model = self.models[modelindex]
        print("Generating Training Data for isolen:" + str(length))
        if masses is None or dists is None:
            input, target = gen_training_data(n, isolen=length, massrange=massrange)
            testinput, testtarget = gen_training_data(int(n * 0.1), isolen=length, massrange=massrange)
        else:
            if dists.shape[1] < length:
                print("Warning: Isolen does not match training data. Changing to", dists.shape[1])
                self.isolen = dists.shape[1]

            elif dists.shape[1] > length:
                print("Warning: Isolen does not match training data. Truncating to", self.isolen)
                dists = dists[:, :self.isolen]

            # Shuffle once and split
            n = len(masses)
            indices = np.random.permutation(n)
            split = int(n * 0.9)

            train_idx = indices[:split]
            test_idx = indices[split:]

            input = masses[train_idx]
            target = dists[train_idx]
            testinput = masses[test_idx]
            testtarget = dists[test_idx]
        print("Created Data:", length, n)
        trd, ted = self.create_data_loaders([input, target], [testinput, testtarget])
        model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

    def transfer_train(self, trainfile, epochs=10, length=128):
        tdata = np.load(trainfile)
        dists = tdata["dists"]
        masses = tdata["masses"]
        trd, ted = self.setup_data(dists, masses)

        indices = np.where(self.lengths == length)
        if len(indices) > 0:
            model = self.models[indices[0][0]]
            model.run_training(trd, ted, epochs=epochs)
        else:
            print("No model with the requested isolen (" + str(length) + ") exists")
            raise ValueError("No model with selected isolen exists.")



    def train_all(self, n=100000, epochs=10, forcenew=False, ignore=[]):
        for i, l in enumerate(self.lengths):
            if l in ignore:
                continue
            self.train(n, epochs, l)

    def get_model_index(self, mass):
        for i in range(len(self.massranges)):
            if self.massranges[i][0] <= mass <= self.massranges[i][1]:
                return i

        print("Mass out of range", mass)
        raise ValueError("Mass out of range. Must be less than 1 MDa.")

    def inputs_to_vectors(self, inputs):
        return np.array([mass_to_vector(m) for m in inputs])

    def predict(self, mass, isolen=None):
        import unidec.IsoDec.IsoGen.isogen_tools as ig
        vec = ig.mass_to_vector(mass)
        if isolen is None:
            index = self.get_model_index(mass)
            model = self.models[index]
            return model.predict(vec)
        else:
            indices = np.where(self.lengths == isolen)
            if len(indices) > 0:
                model = self.models[indices[0][0]]
                return model.predict(vec)
            else:
                print("No model with the requested isolen (" + str(isolen) + ") exists")
                raise ValueError("No model with selected isolen exists.")

def gen_training_data(n, isolen=128, massrange=[100, 1000000], log=False):
    if log:
        rlog = np.random.uniform(np.log10(massrange[0]), np.log10(massrange[1]), n)
        randmasses = np.power(10, rlog)
    else:
        randmasses = np.random.uniform(massrange[0], massrange[1], n)
    # dists = [isomike(mass=m, length=isolen) for m in randmasses]
    dists = [fft_pep_mass_to_dist(m, isolen) for m in randmasses]
    return randmasses, np.array(dists)

if __name__ == "__main__":

    eng = IsoGenMassEngine()
    #os.chdir("Z:\\Group Share\\JGP\\PeptideTraining")
    #masses, dists = gen_training_data(1000000, isolen=1024, massrange=[80000, 1100000])
    #np.savez_compressed("massdata_"+str(len(masses))+".npz", masses=masses, dists=dists)

    if False:
        n = 60000
        # eng.train_all(n, 20, forcenew=False, ignore=[])
        eng.train(n, 10, 1024, forcenew=False)
        # eng.train_multiple(["massdata_100000.npz"], epochs=40, forcenew=False, inputname="masses")
    # exit()

    if False:
        eng.train_all(n=1000000, epochs=10, forcenew=False)

    if True:
        topdir = r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training"
        os.chdir(topdir)

        import fnmatch

        def match_files(directory, string, exclude=None):
            files = []
            for file in os.listdir(directory):
                if fnmatch.fnmatch(file, string):
                    if exclude is None or exclude not in file:
                        files.append(file)
            return np.array(files)

        def correct_dist_lengths(dists, length):
            corrected_dists = []
            for d in dists:
                if len(d) == length:
                    corrected_dists.append(d)
                if len(d) > length:
                    corrected = [d[i] for i in range(length)]
                    corrected_dists.append(corrected)
                if len(d) < length:
                    diff = length - len(d)
                    corrected = np.append(d, np.zeros(diff))
                    corrected_dists.append(corrected)

            return corrected_dists

        polyXdata = np.load(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training\poly_X_peptides_min_6_max_25.npz")
        polyXdists = polyXdata['dists']
        polyXmasses = polyXdata['masses']

        smallpepdata = np.load(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training\all_peptides_min_1_max_5.npz")
        smallpepdists = smallpepdata['dists']
        smallpepmasses = smallpepdata['masses']

        trainfiles = match_files(topdir, "*1100.npz")

        dists = []
        masses = []
        for file in trainfiles:
            data = np.load(file)
            dists.extend(data['dists'])
            masses.extend(data['masses'])

        data_indices = [[] for massrange in eng.massranges]

        for i in range(len(masses)):
            for j in range(len(eng.massranges)):
                if eng.massranges[j][0] <= masses[i] <= eng.massranges[j][1]:
                    data_indices[j].append(i)

        for i in range(len(eng.lengths)):
            curr_masses = np.array([masses[j] for j in data_indices[i]])
            # curr_masses = np.append(curr_masses, smallpepmasses)
            # curr_masses = np.append(curr_masses, polyXmasses)

            curr_dists = np.array([dists[j] for j in data_indices[i]])
            # curr_dists = np.concatenate((curr_dists, smallpepdists))
            # curr_dists = np.concatenate((curr_dists, polyXdists))

            curr_dists = np.array(correct_dist_lengths(curr_dists, eng.lengths[i]))

            eng.isolen = eng.lengths[i]
            eng.train(masses=curr_masses, dists=curr_dists, length=eng.lengths[i], epochs=10, forcenew=True)



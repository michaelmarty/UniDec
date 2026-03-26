import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenc import fft_pep_mass_to_dist
from unidec.modules.isotopetools import isomike


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

    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
        except:
            print("Length not found in list")
            return
        massrange = self.massranges[modelindex]
        model = self.models[modelindex]
        print("Generating Training Data for isolen:" + str(length))
        input, target = gen_training_data(n, isolen=length, massrange=massrange)
        testinput, testtarget = gen_training_data(
            int(n * 0.1), isolen=length, massrange=massrange
        )
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

    if True:
        eng.train_all(n=1000000, epochs=10, forcenew=False)


from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenc import fft_rna_mass_to_dist
from unidec.IsoDec.IsoGen.rna_compositional_model import *
import matplotlib.pyplot as plt
import multiprocessing
from unidec.modules.isotopetools import isojim_rna
import time
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

def process_sequence(formula_isolen_tuple, length):
    return isojim_rna(formula_isolen_tuple[0], length=formula_isolen_tuple[1])

class IsoGenRNAveragineEngine(IsoGenEngineBase):
    def __init__(self, isolen=64):
        super().__init__()
        self.veclen = 5
        self.model = IsoGenModelBase(isolen=isolen, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=0)
        self.lengths = np.array([32, 64, 128])
        self.massranges = np.array(
            [[200, 20000], [1000, 80000], [8000, 160000]])
        self.isolen = isolen
        self.modelindex = 0
        try:
            self.modelindex = np.argwhere(self.lengths == isolen)[0][0]
        except:
            print("Length not found in list")
            return
        self.models = []
        for l in self.lengths:
            model = IsoGenModelBase(isolen=l, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=0)
            self.models.append(model)


    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        print("Training RNA Averagine Model")
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
            model = self.models[modelindex]
            print("Generating Training Data")
            massrange = self.massranges[modelindex]
            print("Mass Range:", massrange)

            input, target = gen_training_data(n, isolen=length, massrange=massrange)
            testinput, testtarget = gen_training_data(int(n * 0.1), isolen=length, massrange=massrange)
            print("Created Data:", length, n)
            trd, ted = self.create_data_loaders([input, target], [testinput, testtarget])

            model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

        except:
            print("Length not found in list")
            return


    def train_all(self, n=100000, epochs=10, forcenew=False, ignore=[]):
        for i, l in enumerate(self.lengths):
            if l in ignore:
                continue
            self.train(n, epochs, length=l, forcenew=forcenew)
            print("Trained RNA Averagine Model for length", l)

    def inputs_to_vectors(self, inputs):
        return np.array([mass_to_vector(m) for m in inputs])


    def predict(self, mass):
        vec = mass_to_vector(mass)
        for i in range(len(vec)):
            print("Py vec", i, ":", vec[i])
        isolen = rnamass_to_isolen(mass)
        model = None
        for m in self.models:
            if m.isolen == isolen:
                model = m
                break
        if model == None:
            print("Unsupported mass.")
            return None
        else:
            return model.predict(vec)

def rna_mass_to_dist(input):
    dist = fft_rna_mass_to_dist(input[0], input[1])
    return dist

def gen_training_data(n, isolen=128, massrange=[100, 1000000]):
    masses = np.random.uniform(massrange[0], massrange[1], n)
    inputs = [(mass, isolen) for mass in masses]
    print("Calculating Isotope Dists...")
    dists = []
    for i in range(len(masses)):
        dist = fft_rna_mass_to_dist(masses[i], isolen)
        dists.append(dist)
    return np.array(masses), np.array(dists)

def rna_seq_to_mass(seq):
    # Get the molecular weight of the sequence
    mw = 0
    for base in seq:
        mw += mm.Formula(rna_dict[base]).mass
    return mw

def rna_seq_to_formula(seq):
    #Get the molecular formula of the sequence
    formula = np.zeros(5)
    for base in seq:
        formula += rna_dict_sep[base]
    return formula

def mass_to_rnaveragine_dist(mass_len_tuple):
    dist = fft_rna_mass_to_dist(mass_len_tuple[0], isolen=128)
    return dist



if __name__ == "__main__":
    eng = IsoGenRNAveragineEngine()
    eng.train_all(1000000, forcenew=True)









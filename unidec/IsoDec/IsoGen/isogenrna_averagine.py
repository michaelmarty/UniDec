from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogen_base import *
from unidec.IsoDec.IsoGen.isogenmass import IsoGenMassEngine
from unidec.IsoDec.IsoGen.rna_compositional_model import *
import matplotlib.pyplot as plt
import multiprocessing
from unidec.modules.isotopetools import isojim_rna
import time
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

def process_sequence(formula_isolen_tuple, length):
    return isojim_rna(formula_isolen_tuple[0], length=formula_isolen_tuple[1])

class IsoGenRNAveragineEngine(IsoGenEngineBase):
    def __init__(self, isolen=128):
        super().__init__()
        self.veclen = 5
        if isolen > 128:
            modelid = 1
        else:
            modelid = 0
        self.model = IsoGenModelBase(isolen=isolen, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=modelid)
        self.lengths = np.array([8, 32, 64])
        self.massranges = np.array(
            [[10, 2400], [800, 20000], [1000, 80000]])
        self.isolen = isolen
        self.modelindex = -1
        try:
            self.modelindex = np.argwhere(self.lengths == isolen)[0][0]
        except:
            print("Length not found in list")
            return
        self.models = []
        for l in self.lengths:
            if l > 128:
                modelid = 1
            else:
                modelid = 0
            model = IsoGenModelBase(isolen=l, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=modelid)
            self.models.append(model)


    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        print("Training RNA Averagine Model")
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
        except:
            print("Length not found in list")
            return
        print("Generating Training Data")
        massrange = self.massranges[modelindex]
        print("Mass Range:", massrange)

        input, target = gen_training_data(n, isolen=length, massrange=massrange)
        testinput, testtarget = gen_training_data(int(n * 0.1), isolen=length, massrange=massrange)
        print("Created Data:", length, n)
        trd, ted = self.create_data_loaders([input, target], [testinput, testtarget])

        self.model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

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
        return self.model.predict(vec)

def gen_training_data(n, isolen=128, massrange=[100, 1000000], rand_seqs=False):


    masses = []

    if rand_seqs:
        #determine the min and max lengths
        #Divide the min mass by the max nucleotide mass to get the min length
        min_length = np.floor(massrange[0] / 305.2)
        #Divide the max mass by the min nucleotide mass to get the max length
        max_length = np.ceil(massrange[1] / 345.2)

        #Generate n random lengths between min and max length
        lengths = np.random.randint(min_length, max_length, n)
        #Generate n random sequences of the lengths
        sequences = []

        for l in lengths:
            #Order doesn't matter so we can just do this
            num_A = np.random.randint(0, l)
            num_C = np.random.randint(0, l - num_A)
            num_G = np.random.randint(0, l - num_A - num_C)
            num_U = l - num_A - num_C - num_G
            #Generate a random sequence of length l
            seq = ""
            seq = seq + "A" * num_A
            seq = seq + "C" * num_C
            seq = seq + "G" * num_G
            seq = seq + "U" * num_U

            mass = 0
            mass += num_A * mm.Formula(rna_dict["A"]).mass
            mass += num_C * mm.Formula(rna_dict["C"]).mass
            mass += num_G * mm.Formula(rna_dict["G"]).mass
            mass += num_U * mm.Formula(rna_dict["U"]).mass

            masses.append(mass)
            sequences.append(seq)

        print("Determining formulas...")
        with multiprocessing.Pool() as pool:
            formulas = pool.map(rna_seq_to_formula, sequences)


        input_data= [(formula, isolen) for formula in formulas]

        print("Calculating Isotope Dists...")
        with multiprocessing.Pool() as pool:
            dists = pool.map(process_sequence, input_data)

    else:
        masses = np.random.uniform(massrange[0], massrange[1], n)
        inputs = [(mass, isolen) for mass in masses]
        print("Calculating Isotope Dists...")
        with multiprocessing.Pool() as pool:
            dists = pool.map(mass_to_rnaveragine_dist, inputs)

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
    num_rnaveragines = int(mass_len_tuple[0] / rnaveragine_mass)
    isolist = rnaveragine_comp_numerical * num_rnaveragines
    dist = isojim_rna(isolist, length=mass_len_tuple[1])
    return dist



if __name__ == "__main__":
    # ls = rna_mass_to_isolist(1285.760279992)



    #mass of GUAC
    # eng = IsoGenMassEngine()
    # m_td = eng.predict(1285.760279992)
    # rna_dist = isojim_rna(ls)
    #
    # # mass to dist approach
    # import matplotlib.pyplot as plt
    # import matplotlib as mpl
    #
    # mpl.use('WxAgg')
    #
    # true_mass = rna_seq_to_mass("GUAC")






    mpl.use("WxAgg")
    eng = IsoGenRNAveragineEngine(isolen=8)
    if True:
        # train the model
        eng.train(n=1000000, epochs=10, length=eng.isolen, forcenew=True)

    exit()
    #test the speed
    if True:
        n = 10000
        masses = np.random.uniform(eng.massranges[eng.modelindex][0], eng.massranges[eng.modelindex][1], n)

        #predict a distribution for each mass
        pred_dists = []
        avg_dists = []
        chis = []
        for m in masses:
            pred_dists.append(eng.predict(m))
            avg_dists.append(mass_to_rnaveragine_dist(mass_len_tuple=(m, eng.isolen)))
            chi = calculate_pearson_chisquared(pred_dists[-1], avg_dists[-1])
            chis.append(chi)
            if True:
                print("Chi squared:", chi)
                plt.plot(pred_dists[-1], label="AI")
                plt.plot(avg_dists[-1], label="RNAveragine")
                plt.legend()
                plt.title(f"Mass: {m:.2f} Da")
                plt.show()

        print("Average Chi Squared:", np.mean(chis))









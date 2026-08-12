import math

import numpy as np
from numpy import fft as fftpack
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.IsoDec.IsoGen.isogen_tools import *
import os
from unidec.IsoDec.IsoGen.isogenpep import IsoGenPepEngine
from collections import Counter
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import calculate_cosinesimilarity
import shutil
import random
import pickle
from matplotlib.pyplot import rcParams


mass_diff_c = 1.0033
massavgine = 111.1254
avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])

isotopes = np.array([[[12., 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
                     [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
                     [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
                     [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
                     [[31.9720707, 95.02], [32.9714585, 0.75], [33.9678668, 4.21], [35.9760809, 0.02]]])


isolength = 1024
buffer = np.zeros(isolength)

h = np.array([1, 0.000115, 0, 0])
c = np.array([1,0.010816,0,0])
#c = np.array([1,0.104,0,0]) #Weird A, G, and L have this isotope ratio
n = np.array([1,0.008,0,0])
o = np.array([1, 0.0004, 0.002, 0])
s = np.array([1, 0.0079, 0.044, 0])
n_heavy = np.array([0,1,0,0])
c_heavy = np.array([0,1,0,0])

h = np.append(h, buffer)
c = np.append(c, buffer)
n = np.append(n, buffer)
o = np.append(o, buffer)
s = np.append(s, buffer)
n_heavy = np.append(n_heavy, buffer)
c_heavy = np.append(c_heavy, buffer)

dt = np.dtype(np.complex128)
hft = fftpack.rfft(h).astype(dt)
cft = fftpack.rfft(c).astype(dt)
nft = fftpack.rfft(n).astype(dt)
oft = fftpack.rfft(o).astype(dt)
sft = fftpack.rfft(s).astype(dt)
c_heavy_ft = fftpack.rfft(c_heavy).astype(dt)
n_heavy_ft = fftpack.rfft(n_heavy).astype(dt)

pre_ft = np.array([cft, hft, nft, oft, sft])

#Now do the set-up by AA instead of atom.
aa_dists = {
    'A': np.array([1.0, 0.320975, 0.037256, 0.00206, 8.1e-05, 3e-06, 0.0, 0.0]),
    'C': np.array([1.0, 0.055875, 0.047244, 0.002236, 0.000129, 5e-06, 0.0, 0.0]),
    'D': np.array([1.0, 0.061775, 0.007538, 0.000388, 2.1e-05, 1e-06, 0.0, 0.0]),
    'E': np.array([1.0, 0.083005, 0.008958, 0.000555, 3e-05, 1e-06, 0.0, 0.0]),
    'F': np.array([1.0, 0.126435, 0.0092, 0.000495, 2e-05, 1e-06, 0.0, 0.0]),
    'G': np.array([1.0, 0.216745, 0.014641, 0.000529, 2.5e-05, 0.0, 0.0, 0.0]),
    'H': np.array([1.0, 0.103205, 0.006723, 0.000331, 1.2e-05, 0.0, 0.0, 0.0]),
    'I': np.array([1.0, 0.087665, 0.005303, 0.000244, 7e-06, 0.0, 0.0, 0.0]),
    'K': np.array([1.0, 0.09578, 0.006016, 0.000287, 9e-06, 0.0, 0.0, 0.0]),
    'L': np.array([1.0, 0.633665, 0.170286, 0.025341, 0.002311, 0.000138, 6e-06, 0.0]),
    'M': np.array([1.0, 0.082335, 0.048904, 0.003496, 0.000196, 8e-06, 0.0, 0.0]),
    'N': np.array([1.0, 0.06949, 0.006012, 0.000308, 1.2e-05, 0.0, 0.0, 0.0]),
    'P': np.array([1.0, 0.074205, 0.004299, 0.000186, 5e-06, 0.0, 0.0, 0.0]),
    'Q': np.array([1.0, 0.08272, 0.006935, 0.000387, 1.6e-05, 1e-06, 0.0, 0.0]),
    'R': np.array([1.0, 0.11178, 0.007612, 0.00039, 1.4e-05, 0.0, 0.0, 0.0]),
    'S': np.array([1.0, 0.074605, 0.006328, 0.000336, 1.4e-05, 0.0, 0.0, 0.0]),
    'T': np.array([1.0, 0.061605, 0.005527, 0.000264, 1e-05, 0.0, 0.0, 0.0]),
    'V': np.array([1.0, 0.074435, 0.004316, 0.000187, 5e-06, 0.0, 0.0, 0.0]),
    'W': np.array([1.0, 0.16055, 0.013895, 0.000859, 4e-05, 1e-06, 0.0, 0.0]),
    'Y': np.array([1.0, 0.126835, 0.011251, 0.000752, 3.8e-05, 2e-06, 0.0, 0.0])
}

# A = [1.0, 0.047975, 0.002865, 0.000102, 2e-06, 0.0, 0.0, 0.0]
# C = [1.0, 0.055875, 0.047244, 0.002236, 0.000129, 5e-06, 0.0, 0.0]
# D = [1.0, 0.061775, 0.007538, 0.000388, 2.1e-05, 1e-06, 0.0, 0.0]
# E = [1.0, 0.083005, 0.008958, 0.000555, 3e-05, 1e-06, 0.0, 0.0]
# F = [1.0, 0.126435, 0.0092, 0.000495, 2e-05, 1e-06, 0.0, 0.0]
# G = [1.0, 0.034745, 0.002403, 7e-05, 1e-06, 0.0, 0.0, 0.0]
# H = [1.0, 0.103205, 0.006723, 0.000331, 1.2e-05, 0.0, 0.0, 0.0]
# I = [1.0, 0.087665, 0.005303, 0.000244, 7e-06, 0.0, 0.0, 0.0]
# K = [1.0, 0.09578, 0.006016, 0.000287, 9e-06, 0.0, 0.0, 0.0]
# L = [1.0, 0.087665, 0.005303, 0.000244, 7e-06, 0.0, 0.0, 0.0]
# M = [1.0, 0.082335, 0.048904, 0.003496, 0.000196, 8e-06, 0.0, 0.0]
# N = [1.0, 0.06949, 0.006012, 0.000308, 1.2e-05, 0.0, 0.0, 0.0]
# P = [1.0, 0.074205, 0.004299, 0.000186, 5e-06, 0.0, 0.0, 0.0]
# Q = [1.0, 0.08272, 0.006935, 0.000387, 1.6e-05, 1e-06, 0.0, 0.0]
# R = [1.0, 0.11178, 0.007612, 0.00039, 1.4e-05, 0.0, 0.0, 0.0]
# S = [1.0, 0.074605, 0.006328, 0.000336, 1.4e-05, 0.0, 0.0, 0.0]
# T = [1.0, 0.061605, 0.005527, 0.000264, 1e-05, 0.0, 0.0, 0.0]
# V = [1.0, 0.074435, 0.004316, 0.000187, 5e-06, 0.0, 0.0, 0.0]
# W = [1.0, 0.16055, 0.013895, 0.000859, 4e-05, 1e-06, 0.0, 0.0]
# Y = [1.0, 0.126835, 0.011251, 0.000752, 3.8e-05, 2e-06, 0.0, 0.0]


K_heavy = [1.0, 0.01755, 0.00209, 3.4e-05, 0.0, 0.0, 0.0, 0.0]

def initialize_aa_preft(isolength=64):
    buffer = np.zeros(isolength)
    dt = np.dtype(np.complex128)

    aa_preft = {}

    for key,value in aa_dists.items():
        v = np.append(value, buffer)
        aa_preft[key] = fftpack.rfft(v).astype(dt)

    return aa_preft

def initialize_aa_preft_heavyK(isolength=64):
    buffer = np.zeros(isolength)
    dt = np.dtype(np.complex128)

    aa_preft = {}

    for key,value in aa_dists.items():
        if key == 'K':
            v = np.append(K_heavy, buffer)
            aa_preft[key] = fftpack.rfft(v).astype(dt)
        else:
            v = np.append(value, buffer)
            aa_preft[key] = fftpack.rfft(v).astype(dt)

    return aa_preft

def pepseq_to_isolist(seq):
    isolist = np.array([0,2,0,1,0])
    for aa in seq:
        if aa in aa_elementalcomp_dict:
            isolist += aa_elementalcomp_dict[aa]
        else:
            print("Unexpected amino acid found: ", str(aa))
    return isolist


def isojim_local(isolist, length=16):
    """Thanks to Jim Prell for Sketching this Code"""
    '''
    Takes as in input a list of numbers of elements in the order above [C, H, N, O, S]
    Calculates the aggregated isotope distribution based on the input list.
    '''
    num_atoms = len(isolist)
    allft = 0

    for i in range(num_atoms):
        if i == 0:
            allft = pre_ft[i] ** isolist[i]
        else:
            allft = allft * (pre_ft[i] ** isolist[i])

    allift = np.abs(fftpack.irfft(allft))

    allift = allift / np.amax(allift)
    return allift[:length]


def isojim_aa(sequence, aa_preft, length=16):
    aa_counts = Counter(sequence)

    aas_processed = 0

    allft = 0

    for aa, count in aa_counts.items():
        if aas_processed == 0:
            allft = aa_preft[aa] ** count
        else:
            allft = allft * (aa_preft[aa] ** count)

        aas_processed += 1

    allift = np.abs(fftpack.irfft(allft))

    allift = allift / np.amax(allift)
    return allift[:length]


def add_noise_to_dist(dist, noise_level=0.2):
    noisy_dist = [0 for _ in dist]
    noise_multipliers = np.random.normal(loc = 1, scale=noise_level, size=len(dist))
    for i in range(len(dist)):
        noisy_dist[i] = dist[i] * noise_multipliers[i]
    return noisy_dist


def get_aa_isodist(aa):
    isolist = aa_elementalcomp_dict[aa]
    dist = isojim_local(isolist, 16)
    print("AA: ", aa)
    print("Dist: [" + str(round(dist[8], 6)) + ", " + str(round(dist[9], 6)) + ", " + str(round(dist[10], 6)) + ", " + str(round(dist[11], 6)) + ", " +
          str(round(dist[12], 6)) + ", " + str(round(dist[13], 6)) + ", " + str(round(dist[14], 6)) + ", " + str(round(dist[15], 6)) + "]")


def get_tryptic_peptides(min_length=7, max_length=50):
    data = np.load(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training\human_protein_seqs.npz")
    seqs = data['seqs']

    tryptic_peptides = []

    for seq in seqs:
        if len(seq) > 300:
            continue
        curr_pep = ''
        for aa in seq:
            curr_pep = curr_pep + aa
            if aa == 'K' or aa == 'R':
                if len(curr_pep) >= min_length and len(curr_pep) <= max_length:
                    tryptic_peptides.append(curr_pep)
                    curr_pep = ''

    print("Generated ", len(tryptic_peptides), " tryptic peptide sequences.")
    return tryptic_peptides


def select_random_seqs(all_seqs, n=1000):
    inds = random.sample(range(len(all_seqs)), n)
    selected = [all_seqs[i] for i in inds]
    return selected


def assess_model_performance(all_seqs, aa_preft, isogenpep_engine):
    css_vals = []
    for seq in all_seqs:
        true_dist = isojim_aa(seq, aa_preft)
        nn_dist = isogenpep_engine.predict(seq)
        css = calculate_cosinesimilarity(nn_dist, true_dist)
        css_vals.append(css)

    return np.array(css_vals)


def assess_model_performance_inputdists(seqs, dists, engine):
    css_vals = []
    for i in range(len(dists)):
        nn_dist = engine.predict(seqs[i])
        true_dist = dists[i]
        css = calculate_cosinesimilarity(nn_dist, true_dist)
        css_vals.append(css)
    return np.array(css_vals)

class Histogram:
    def __init__(self, min, max, step):
        self.min = min
        self.max = max
        self.step = step
        self.bin_count = math.ceil((max - min)/ step)
        self.bins = [0 for i in range(self.bin_count+1)]
        self.bin_mins = self.populate_bin_mins()
        self.bin_fracs = None


    def populate_bin_mins(self):
        bin_mins = [0 for i in range(self.bin_count + 1)]
        bin_mins[0] = self.min
        for i in range(1, self.bin_count+1):
            bin_mins[i] = bin_mins[i-1] + self.step
        return bin_mins

    def add_value(self, value):
        bin = math.floor((value - self.min) / self.step)
        self.bins[bin] += 1

    def add_list(self, list):
        for value in list:
            self.add_value(value)

    def add_existing_hist(self, hist):
        if hist.min != self.min or hist.max != self.max or hist.step != self.step:
            print("Unable to merge histograms with different parameters.")
            return

        for i, val in enumerate(hist.bins):
            self.bins[i] += val

    def populate_bin_fracs(self, normalize=True):
        total = sum(self.bins)
        self.bin_fracs = [[] for i in range(self.bin_count + 1)]
        for i,v in enumerate(self.bins):
            self.bin_fracs[i] = self.bins[i] / total

        if normalize:
            self.bin_fracs = np.array(self.bin_fracs) / max(self.bin_fracs)



def reset_model():
    original_pth_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\OriginalModels\isogenpep_model_16.pth"
    original_bin_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\OriginalModels\isogenpep_model_16.pth"

    existing_pth_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\isogenpep_model_16.pth"
    existing_bin_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\isogenpep_model_16.bin"

    shutil.copy(original_bin_path, existing_bin_path)
    shutil.copy(original_pth_path, existing_pth_path)

    return

def delete_model():
    existing_pth_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\isogenpep_model_16.pth"
    existing_bin_path = r"C:\Python\UniDecDev\unidec\IsoDec\IsoGen\isogenpep_model_16.bin"

    os.remove(existing_pth_path)
    os.remove(existing_bin_path)

def transfer_train_and_test_inputdists(training_seqs, training_dists, test_seqs, test_dists, engine=None, forcenew=False):
    if forcenew:
        delete_model()

    if engine is None:
        engine = IsoGenPepEngine(isolen=16)

    training_seqs = np.array(training_seqs)
    training_dists= np.array(training_dists)

    engine.transfer_train(training_seqs, training_dists, epochs=10, length=16, forcenew=forcenew)

    transfertrained_css_vals = assess_model_performance_inputdists(test_seqs, test_dists, engine)

    #reset_model()
    return transfertrained_css_vals

def transfer_train_and_test(all_seqs, aa_preft, n_training=1000, forcenew=False, input_eng=None):
    if forcenew:
        delete_model()
        eng = IsoGenPepEngine(isolen=16)
    else:
        if input_eng is None:
            eng = IsoGenPepEngine(isolen=16)
        else:
            eng = input_eng

    training_seqs = select_random_seqs(all_seqs, n=n_training)
    training_dists = []

    for seq in training_seqs:
        dist = isojim_aa(seq, aa_preft)
        noisy_dist = add_noise_to_dist(dist, noise_level=0.2)
        training_dists.append(noisy_dist)

    training_seqs = np.array(training_seqs)
    training_dists= np.array(training_dists)

    print("Transfer-training model.")
    eng.transfer_train(training_seqs, training_dists, epochs=10, length=16, forcenew=forcenew)

    print("Assessing transfer-trained model.")
    transfertrained_css_vals = assess_model_performance(all_seqs, aa_preft, eng)

    print("Average trained css: ", np.average(transfertrained_css_vals))

    # plt.hist(untrained_css_vals, bins=20, alpha=0.6, label="Untrained")
    # plt.hist(transfertrained_css_vals, bins=20, alpha=0.6, label="Trained")
    # plt.legend()
    # plt.show()

    #reset_model()
    return transfertrained_css_vals


if __name__ == "__main__":
    # get_aa_isodist('K')
    # exit()

    os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover")

    # Set backend to Agg
    mpl.use('WxAgg')
    aa_preft = initialize_aa_preft(1024)
    aa_preft_heavyk = initialize_aa_preft_heavyK(1024)

    if False:
        seq = "GDGKLVVIWPPR"
        light_isodist = isojim_aa(seq, aa_preft)
        heavy_isodist = isojim_aa(seq, aa_preft_heavyk)

        plt.plot(light_isodist, label="Light")
        plt.plot(heavy_isodist, label="Heavy")
        plt.legend()
        plt.show()
        exit()

    if False:
        rcParams['ps.useafm'] = True
        rcParams['ps.fonttype'] = 42
        rcParams['pdf.fonttype'] = 42
        rcParams['font.family'] = 'Arial'

        tick_size = 16
        label_size = 20
        linewidth = 3

        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']

        save_figs = True

        with open("trained.txt", "rb") as myFile:
            trained = pickle.load(myFile)

        with open("random_init_trained.txt", "rb") as myFile:
            random_init = pickle.load(myFile)


        ns = []
        trained_avg_css = []
        trained_css_std_dev = []

        random_init_ns = []
        random_init_avg_css = []
        random_init_css_std_dev = []

        for k,v in trained.items():
            ns.append(k)
            all_css_vals = []
            for it in v:
                all_css_vals.extend(it)
            trained_avg_css.append(np.average(all_css_vals))
            trained_css_std_dev.append(np.std(all_css_vals))

        for k,v in random_init.items():
            random_init_ns.append(k)
            all_css_vals = []
            for it in v:
                all_css_vals.extend(it)
            random_init_avg_css.append(np.average(all_css_vals))
            random_init_css_std_dev.append(np.std(all_css_vals))

        ns.append(0)
        trained_avg_css.append(0.7566234)
        trained_css_std_dev.append(0)

        n_ordered_indices = np.argsort(ns)

        ns = [ns[i] for i in n_ordered_indices]
        trained_avg_css = [trained_avg_css[i] for i in n_ordered_indices]
        trained_css_std_dev = [trained_css_std_dev[i] for i in n_ordered_indices]
        trained_css_95ci = 1.96 * (np.array(trained_css_std_dev) / np.sqrt(10))

        random_init_n_ordered_indices = np.argsort(random_init_ns)
        random_init_ns = [random_init_ns[i] for i in random_init_n_ordered_indices]
        random_init_avg_css = [random_init_avg_css[i] for i in random_init_n_ordered_indices]
        random_init_css_std_dev = [random_init_css_std_dev[i] for i in random_init_n_ordered_indices]
        random_init_95ci = 1.96 * (np.array(random_init_css_std_dev) / np.sqrt(10))

        plt.errorbar(ns, trained_avg_css, yerr=trained_css_95ci, capsize=2, marker="o", color='b', label="Trained Model")
        plt.errorbar(random_init_ns, random_init_avg_css, yerr=random_init_95ci, capsize=2, marker='o', color='r', label="New Model")
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.hlines(1, 0, 5000, colors='k')
        plt.tick_params(length=5, width=2)
        plt.xlim(0, 5100)
        plt.ylim(bottom=0.5)
        plt.xticks([0,1000,2000,3000,4000,5000], fontsize=tick_size)
        plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=tick_size)
        plt.xlabel("Training Peptides", fontsize=label_size)
        plt.ylabel("Average CSS", fontsize=label_size)
        plt.legend(ncol=2, frameon=True, edgecolor="black", fontsize=16, loc='lower right')
        plt.tight_layout()
        if save_figs:
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\traincount_vs_css.pdf", dpi=300, transparent=True)
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\traincount_vs_css.png", dpi=300, transparent=True)
        plt.show()
        plt.show()

        exit()

    if False:
        reset_model()

        tryptic_peps = get_tryptic_peptides(min_length=7, max_length=10)

        eng = IsoGenPepEngine(isolen=16)
        print("Assessing untrained model.")
        untrained_css_vals = assess_model_performance(tryptic_peps, aa_preft, eng)
        print("Average untrained css: ", np.average(untrained_css_vals))

        n_values = [50, 100, 250, 500, 1000, 2500, 5000]
        #n_values = [100, 500]
        iterations = 10

        ns = []
        trained_css_values = {}
        random_init_trained_css_values = {}

        for n in n_values:
            for i in range(iterations):
                trained = transfer_train_and_test(tryptic_peps, aa_preft, n, forcenew=False)
                random_init_trained = transfer_train_and_test(tryptic_peps, aa_preft, n, forcenew=True)

                if n in trained_css_values:
                    trained_css_values[n].append(trained)
                else:
                    trained_css_values[n] = [trained]

                if n in random_init_trained_css_values:
                    random_init_trained_css_values[n].append(random_init_trained)
                else:
                    random_init_trained_css_values[n] = [random_init_trained]


        with open("trained.txt", "wb") as myFile:
            pickle.dump(trained_css_values, myFile)

        with open("random_init_trained.txt", "wb") as myFile:
            pickle.dump(random_init_trained_css_values, myFile)
        exit()

    if True:
        reset_model()
        data = np.load("transfer_learning_dists_8762.npz")
        dists = data['dists']
        seqs = data['seqs']

        max_len = 50

        filtered = [
            (seq, dist)
            for seq, dist in zip(seqs, dists)
            if len(seq) < max_len
        ]

        seqs = [seq for seq, dist in filtered]
        dists = [dist for seq, dist in filtered]

        iterations = 1
        untrained_avg_css_vals = []
        trained_avg_css_vals = []

        untrained_histogram = Histogram(0,1, 0.01)
        trained_histogram = Histogram(0,1,0.01)

        for i in range(iterations):
            trainper = 0.9
            ntrain = int(round(trainper*len(dists)))
            ntest = len(dists) - ntrain

            untrained_hist = Histogram(0, 1, 0.01)
            trained_hist = Histogram(0, 1, 0.01)

            train_indices= random.sample(range(0,len(dists)), ntrain)
            all_indices = range(0,len(dists))

            training_dists = [dists[i] for i in train_indices]
            training_seqs = [seqs[i] for i in train_indices]

            test_dists = [dists[i] for i in all_indices if i not in train_indices]
            test_seqs = [seqs[i] for i in all_indices if i not in train_indices]

            eng = IsoGenPepEngine(isolen=16)

            untrained_css_vals = assess_model_performance_inputdists(seqs, dists, eng)
            untrained_hist.add_list(untrained_css_vals)
            untrained_histogram.add_existing_hist(untrained_hist)

            untrained_avg_css_vals.append(np.average(untrained_css_vals))

            eng = None

            trained_css_vals = transfer_train_and_test_inputdists(training_seqs, training_dists, test_seqs, test_dists, engine=eng, forcenew=False)
            trained_hist.add_list(trained_css_vals)
            trained_histogram.add_existing_hist(trained_hist)
            trained_avg_css_vals.append(np.average(trained_css_vals))

        #reset_model()
        avg_diff = np.array(trained_avg_css_vals) - np.array(untrained_avg_css_vals)
        diff_std = np.std(avg_diff)
        diff_ci = 1.96 * (diff_std / np.sqrt(iterations))
        print("Average difference: ", np.average(avg_diff))
        print("95% CI for Diff: ", diff_ci)
        print("Average untrained: ", np.average(untrained_avg_css_vals))
        print("Average trained: ", np.average(trained_avg_css_vals))

        untrained_histogram.populate_bin_fracs(normalize=False)
        trained_histogram.populate_bin_fracs(normalize=False)

        # Set backend to Agg
        mpl.use('WxAgg')
        rcParams['ps.useafm'] = True
        rcParams['ps.fonttype'] = 42
        rcParams['pdf.fonttype'] = 42
        rcParams['font.family'] = 'Arial'

        tick_size = 16
        label_size = 20
        linewidth = 3

        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']

        plt.plot(untrained_histogram.bin_mins, untrained_histogram.bin_fracs, linewidth=2, color='r', label="Untrained")
        plt.plot(trained_histogram.bin_mins, trained_histogram.bin_fracs, linewidth=2, color='b', label="Trained")
        plt.xlim(0.7,1.01)
        plt.xlabel("CSS", fontsize=label_size)
        plt.ylabel("Fraction", fontsize=label_size)
        plt.xticks([0.7,0.8,0.9,1], labels=["0.7", "0.8", "0.9", "1"], fontsize=tick_size)
        plt.yticks([0,0.3,0.6], labels=["0", "0.3", "0.6"], fontsize=tick_size)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.tick_params(length=5, width=2)
        plt.legend(ncol=2, frameon=True, edgecolor="black", fontsize=16, loc='upper center')
        plt.tight_layout()
        # plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\IsoGen_Figure7.pdf", dpi=300, transparent=True)
        # plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\IsoGen_Figure7.png", dpi=300, transparent=True)
        plt.show()


    if False:
        # Set backend to Agg
        mpl.use('WxAgg')

        rcParams['ps.useafm'] = True
        rcParams['ps.fonttype'] = 42
        rcParams['pdf.fonttype'] = 42
        rcParams['font.family'] = 'Arial'

        tick_size = 16
        label_size = 20
        linewidth = 3

        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']

        reset_model()

        tryptic_peps = get_tryptic_peptides(10,50)
        pep_num = 10

        normal_dist = fft_pep_seq_to_dist(tryptic_peps[pep_num])
        perturbed_dist = isojim_aa(tryptic_peps[pep_num], aa_preft)

        eng = IsoGenPepEngine(isolen=16)

        pre_transfer_learning = eng.predict(tryptic_peps[pep_num])

        transfer_train_and_test(tryptic_peps, aa_preft, n_training=500, forcenew=False, input_eng=eng)

        post_transfer_learning = eng.predict(tryptic_peps[pep_num])


        plt.plot(normal_dist, color=CB_color_cycle[0], linewidth=2, label="Natural")
        plt.plot(perturbed_dist, color=CB_color_cycle[1], linewidth=2, label="Perturbed")
        plt.plot(pre_transfer_learning, color=CB_color_cycle[2], linewidth=2, linestyle='--', label="NN-Pre")
        plt.plot(post_transfer_learning, color=CB_color_cycle[3], linewidth=2, linestyle='--', label="NN-Post")
        plt.title(tryptic_peps[pep_num], fontsize=24)
        plt.xlabel("Isotopologue #", fontsize=label_size)
        plt.ylabel("Abundance", fontsize=label_size)
        plt.xlim(0,8)
        plt.ylim(0,1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.tick_params(length=5, width=2)
        plt.xticks(fontsize=tick_size)
        plt.yticks(labels=["0", "0.5", "1"], ticks=[0,0.5,1], fontsize=tick_size)
        plt.legend(ncol=2, frameon=True, edgecolor="black", fontsize=14)
        plt.tight_layout()
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pre_post_transferlearning.pdf",dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pre_post_transferlearning.png",dpi=300, transparent=True)
        plt.show()

        reset_model()








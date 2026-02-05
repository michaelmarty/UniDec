import matplotlib.pyplot as plt
import numpy as np
import gc
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.IsoDec.IsoGen.isogenmass import IsoGenMassEngine
from unidec.IsoDec.IsoGen.isogenpep import *
from unidec.IsoDec.IsoGen.isogen_tools import *
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *
import fnmatch
import time
import random
from matplotlib.pyplot import rcParams
from Scripts.JGP.IsoGen_Analysis.BRAIN_wrapper import *
import pickle
import math


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def remove_outliers(data, multiplier=2):
    """
    Remove outliers from a list using the IQR method.

    Args:
        data (list or iterable of float): The input list of numbers.
        multiplier (float): The IQR multiplier (default 1.5). Larger = less strict.

    Returns:
        list: A list with outliers removed.
    """
    if len(data) < 4:
        # Too few data points to determine outliers reliably
        return data

    sorted_data = sorted(data)
    q1 = sorted_data[len(sorted_data) // 4]
    q3 = sorted_data[(len(sorted_data) * 3) // 4]
    iqr = q3 - q1

    lower_bound = q1 - multiplier * iqr
    upper_bound = q3 + multiplier * iqr

    return [x for x in data if lower_bound <= x <= upper_bound]

def get_s_count_from_pepmass(mass):
    formula, minmassint, isolist = pep_makemass(mass)
    return isolist[4]

def get_aacount_scount_array():
    aa_counts = np.arange(5,300,1)

    #Store the AA counts where S count increases
    transition_aa_counts = []
    s_counts = []
    previous_s_count = 0

    for aa_count in aa_counts:
        mass = aa_count * massavgine
        formula, minmassint, isolist = pep_makemass(mass)
        if isolist[4] > previous_s_count:
            transition_aa_counts.append(aa_count)
        previous_s_count = isolist[4]

    return transition_aa_counts

def get_key_from_mass(mass, bin_width):
    key = int(math.floor(mass / bin_width) * bin_width)
    return key


# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

save_figs = False

eng = IsoGenPepEngine()
eng_avg = IsoGenMassEngine()

tick_size = 16
label_size = 20

if False:
    masses = [1000, 5000, 10000, 50000]
    isolen = 128

    for mass in masses:
        py_fft_dist = pepmass_to_dist(mass, isolen=128)
        py_nn_dist = eng_avg.predict(mass)
        c_fft_dist = fft_pep_mass_to_dist(mass, isolen=128)
        c_nn_dist = nn_pep_mass_to_dist(mass)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.legend()
        plt.show()

if False:
    seqs = ["ADAKCGS", "AACCGKHIYHHHCCHG"]

    for seq in seqs:
        py_fft_dist = peptide_to_dist(seq, 64)
        py_nn_dist = eng.predict(seq)
        c_fft_dist = fft_pep_seq_to_dist(seq)
        c_nn_dist = nn_pep_seq_to_dist(seq)

        plt.plot(py_fft_dist, color='r', label="Py-FFT")
        plt.plot(py_nn_dist, color='orange', label='Py-NN')
        plt.plot(c_fft_dist, color='b', linestyle=':', label="C-FFT")
        plt.plot(c_nn_dist, color='purple', linestyle='--', label="C-NN")
        plt.legend()
        plt.show()

if True:
    topdir = r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Assessment"
    os.chdir(topdir)
    assessment_files = match_files(topdir, "*1000.npz")

    seqs = []
    dists = []
    vecs = []
    masses = []


    for file in assessment_files:
        data = np.load(file)
        seqs.extend(data["seqs"])
        dists.extend(data["dists"])
        vecs.extend(data["vecs"])
        masses.extend(data["masses"])

    # with open("masses_seqs.txt", "w") as f:
    #     for i in range(len(masses)):
    #         f.write(f"{masses[i]:.4f} {seqs[i]}\n")

    all_aa_counts = np.zeros(len(masses))
    s_counts = np.zeros(len(masses))

    fft_seq_times = {}
    fft_mass_times = {}
    nn_seq_times = {}
    nn_mass_times = {}
    brain_mass_times = {}
    brain_seq_times = {}

    fft_mass_css_vals = {}
    nn_seq_css_vals = {}
    nn_mass_css_vals = {}
    brain_mass_css_vals = {}
    brain_seq_css_vals = {}
    py_mass_css_vals = {}


    #These are the values separated into a mass dictionary (100 Da / bin, and the key is the bottom)
    fft_seq_times_bymass = {}
    fft_mass_times_bymass = {}
    nn_seq_times_bymass = {}
    nn_mass_times_bymass = {}
    brain_mass_times_bymass = {}
    brain_seq_times_bymass = {}

    fft_mass_css_vals_bymass = {}
    nn_seq_css_vals_bymass = {}
    nn_mass_css_vals_bymass = {}
    brain_mass_css_vals_bymass = {}
    brain_seq_css_vals_bymass = {}


    shuffled_indices = np.arange(0, len(masses))
    random.shuffle(shuffled_indices)


    num_processed = 0
    mass_bin_width = 100
    curr_length = 50

    for i in shuffled_indices:
        num_processed += 1
        percent_processed = round((num_processed / len(masses)) * 100, 2)
        if num_processed % 10000 == 0:
            print("Processed ", percent_processed, "% of sequences")

        s_counts[i] = get_s_count_from_pepmass(masses[i])

        aa_count = peptide_to_aacount(seqs[i])


        all_aa_counts[i] = aa_count

        pep_isolen = pepmass_to_isolen(masses[i])

        mass_key = get_key_from_mass(masses[i], mass_bin_width)

        timestart = time.perf_counter()
        fft_seq_dist = fft_pep_seq_to_dist(seqs[i])
        fft_seq_time = round(time.perf_counter() - timestart, 6)

        if aa_count in fft_seq_times:
            fft_seq_times[aa_count].append(fft_seq_time)
        else:
            fft_seq_times[aa_count] = [fft_seq_time]

        if mass_key in fft_seq_times_bymass:
            fft_seq_times_bymass[mass_key].append(fft_seq_time)
        else:
            fft_seq_times_bymass[mass_key] = [fft_seq_time]


        timestart = time.perf_counter()
        fft_mass_dist = fft_pep_mass_to_dist(masses[i])
        fft_mass_time = round(time.perf_counter() - timestart, 6)
        fft_mass_css = calculate_cosinesimilarity(fft_mass_dist, fft_seq_dist)

        if aa_count in fft_mass_times:
            fft_mass_times[aa_count].append(fft_mass_time)
            fft_mass_css_vals[aa_count].append(fft_mass_css)
        else:
            fft_mass_times[aa_count] = [fft_mass_time]
            fft_mass_css_vals[aa_count] = [fft_mass_css]

        if mass_key in fft_mass_times_bymass:
            fft_mass_times_bymass[mass_key].append(fft_mass_time)
            fft_mass_css_vals_bymass[mass_key].append(fft_mass_css)
        else:
            fft_mass_times_bymass[mass_key] = [fft_mass_time]
            fft_mass_css_vals_bymass[mass_key] = [fft_mass_css]


        timestart = time.perf_counter()
        nn_seq_dist = nn_pep_seq_to_dist(seqs[i])
        nn_seq_time = round(time.perf_counter() - timestart, 6)
        nn_seq_css = calculate_cosinesimilarity(nn_seq_dist, fft_seq_dist)

        if aa_count in nn_seq_times:
            nn_seq_times[aa_count].append(nn_seq_time)
            nn_seq_css_vals[aa_count].append(nn_seq_css)
        else:
            nn_seq_times[aa_count] = [nn_seq_time]
            nn_seq_css_vals[aa_count] = [nn_seq_css]

        if mass_key in nn_seq_times_bymass:
            nn_seq_times_bymass[mass_key].append(nn_seq_time)
            nn_seq_css_vals_bymass[mass_key].append(nn_seq_css)
        else:
            nn_seq_times_bymass[mass_key] = [nn_seq_time]
            nn_seq_css_vals_bymass[mass_key] = [nn_seq_css]

        timestart = time.perf_counter()
        nn_mass_dist = nn_pep_mass_to_dist(masses[i])
        nn_mass_time = round(time.perf_counter() - timestart, 6)
        nn_mass_css = calculate_cosinesimilarity(nn_mass_dist, fft_seq_dist)

        if aa_count in nn_mass_times:
            nn_mass_times[aa_count].append(nn_mass_time)
            nn_mass_css_vals[aa_count].append(nn_mass_css)
        else:
            nn_mass_times[aa_count] = [nn_mass_time]
            nn_mass_css_vals[aa_count] = [nn_mass_css]

        if mass_key in nn_mass_times_bymass:
            nn_mass_times_bymass[mass_key].append(nn_mass_time)
            nn_mass_css_vals_bymass[mass_key].append(nn_mass_css)
        else:
            nn_mass_times_bymass[mass_key] = [nn_mass_time]
            nn_mass_css_vals_bymass[mass_key] = [nn_mass_css]

        timestart = time.perf_counter()
        brain_mass_dist = brain_pepmass_to_dist(masses[i])
        brain_mass_time = round(time.perf_counter() - timestart, 6)
        brain_mass_css = calculate_cosinesimilarity(brain_mass_dist, fft_seq_dist)

        if aa_count in brain_mass_times:
            brain_mass_times[aa_count].append(brain_mass_time)
            brain_mass_css_vals[aa_count].append(brain_mass_css)
        else:
            brain_mass_times[aa_count] = [brain_mass_time]
            brain_mass_css_vals[aa_count] = [brain_mass_css]

        if mass_key in brain_mass_times_bymass:
            brain_mass_times_bymass[mass_key].append(brain_mass_time)
            brain_mass_css_vals_bymass[mass_key].append(brain_mass_css)
        else:
            brain_mass_times_bymass[mass_key] = [brain_mass_time]
            brain_mass_css_vals_bymass[mass_key] = [brain_mass_css]


        timestart = time.perf_counter()
        brain_seq_dist = brain_pepseq_to_dist(seqs[i])
        brain_seq_time = round(time.perf_counter() - timestart, 6)
        brain_seq_css = calculate_cosinesimilarity(brain_seq_dist, fft_seq_dist)

        if aa_count in brain_seq_times:
            brain_seq_times[aa_count].append(brain_seq_time)
            brain_seq_css_vals[aa_count].append(brain_seq_css)
        else:
            brain_seq_times[aa_count] = [brain_seq_time]
            brain_seq_css_vals[aa_count] = [brain_seq_css]

        if mass_key in brain_seq_times_bymass:
            brain_seq_times_bymass[mass_key].append(brain_seq_time)
            brain_seq_css_vals_bymass[mass_key].append(brain_seq_css)
        else:
            brain_seq_times_bymass[mass_key] = [brain_seq_time]
            brain_seq_css_vals_bymass[mass_key] = [brain_seq_css]


        # py_isojim = pepmass_to_dist(masses[i])
        # py_isojim_css = calculate_cosinesimilarity(py_isojim, fft_seq_dist)
        #
        # if aa_count in py_mass_css_vals:
        #     py_mass_css_vals[aa_count].append(py_isojim_css)
        # else:
        #     py_mass_css_vals[aa_count] = [py_isojim_css]
        #
        # py_nn_dist = eng_avg.predict(masses[i])
        # #Ensure that the length of py_nn_dist is 128
        # py_nn_dist = np.array(list(py_nn_dist[:128]) + [0] * max(0, 128 - len(py_nn_dist)))
        #
        #
        # if aa_count == curr_length:
        #     plt.plot(fft_seq_dist, label="FFT-Seq")
        #     plt.plot(fft_mass_dist, label="FFT-Mass")
        #     plt.plot(py_nn_dist, label="Py-NN-Mass")
        #     plt.plot(nn_seq_dist, label="NN-Seq")
        #     plt.plot(nn_mass_dist, label="NN-Mass")
        #     plt.title("Length: " + str(curr_length) + " Mass: " + str(masses[i]))
        #     plt.legend()
        #     plt.show()
        #
        #
        #     plt.plot(fft_seq_dist-fft_mass_dist, label="FFT-Mass Diff")
        #     plt.plot(fft_seq_dist-nn_mass_dist, label="NN-Mass Diff")
        #     plt.plot(fft_seq_dist-py_nn_dist, label="Py-NN-Mass Diff")
        #     plt.title("Length: " + str(curr_length) + " Mass: " + str(masses[i]))
        #     plt.legend()
        #     plt.show()
        #
        #     curr_length += 20


    with open("fft_seq_times.txt", "wb") as myFile:
        pickle.dump(fft_seq_times, myFile)
    for k,v in fft_seq_times.items():
        fft_seq_times[k] = np.average(remove_outliers(v))

    with open("fft_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(fft_seq_times_bymass, myFile)
    for k,v in fft_seq_times_bymass.items():
        fft_seq_times_bymass[k] = np.average(remove_outliers(v))


    with open("fft_mass_times.txt", "wb") as myFile:
        pickle.dump(fft_mass_times, myFile)
    for k,v in fft_mass_times.items():
        fft_mass_times[k] = np.average(remove_outliers(v))

    with open("fft_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(fft_mass_times_bymass, myFile)
    for k,v in fft_mass_times_bymass.items():
        fft_mass_times_bymass[k] = np.average(remove_outliers(v))


    with open("nn_seq_times.txt", "wb") as myFile:
        pickle.dump(nn_seq_times, myFile)
    for k,v in nn_seq_times.items():
        nn_seq_times[k] = np.average(remove_outliers(v))

    with open("nn_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(nn_seq_times_bymass, myFile)
    for k,v in nn_seq_times_bymass.items():
        nn_seq_times_bymass[k] = np.average(remove_outliers(v))


    with open("nn_mass_times.txt", "wb") as myFile:
        pickle.dump(nn_mass_times, myFile)
    for k,v in nn_mass_times.items():
        nn_mass_times[k] = np.average(remove_outliers(v))

    with open("nn_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(nn_mass_times_bymass, myFile)
    for k,v in nn_mass_times_bymass.items():
        nn_mass_times_bymass[k] = np.average(remove_outliers(v))


    with open("fft_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(fft_mass_css_vals, myFile)
    for k,v in fft_mass_css_vals.items():
        fft_mass_css_vals[k] = np.average(v)

    with open("fft_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(fft_mass_css_vals_bymass, myFile)
    for k,v in fft_mass_css_vals_bymass.items():
        fft_mass_css_vals_bymass[k] = np.average(v)


    with open("nn_seq_css_vals.txt", "wb") as myFile:
        pickle.dump(nn_seq_css_vals, myFile)
    for k,v in nn_seq_css_vals.items():
        nn_seq_css_vals[k] = np.average(v)

    with open("nn_seq_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(nn_seq_css_vals_bymass, myFile)
    for k,v in nn_seq_css_vals_bymass.items():
        nn_seq_css_vals_bymass[k] = np.average(v)


    with open("nn_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(nn_mass_css_vals, myFile)
    for k,v in nn_mass_css_vals.items():
        nn_mass_css_vals[k] = np.average(v)

    with open("nn_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(nn_mass_css_vals_bymass, myFile)
    for k,v in nn_mass_css_vals_bymass.items():
        nn_mass_css_vals_bymass[k] = np.average(v)


    with open("brain_mass_times.txt", "wb") as myFile:
        pickle.dump(brain_mass_times, myFile)
    for k,v in brain_mass_times.items():
        brain_mass_times[k] = np.average(remove_outliers(v))

    with open("brain_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(brain_mass_times_bymass, myFile)
    for k,v in brain_mass_times_bymass.items():
        brain_mass_times_bymass[k] = np.average(remove_outliers(v))


    with open("brain_seq_times.txt", "wb") as myFile:
        pickle.dump(brain_seq_times, myFile)
    for k,v in brain_seq_times.items():
        brain_seq_times[k] = np.average(remove_outliers(v))

    with open("brain_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(brain_seq_times_bymass, myFile)
    for k,v in brain_seq_times_bymass.items():
        brain_seq_times_bymass[k] = np.average(remove_outliers(v))


    with open("brain_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(brain_mass_css_vals, myFile)
    for k,v in brain_mass_css_vals.items():
        brain_mass_css_vals[k] = np.average(v)

    with open("brain_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(brain_mass_css_vals_bymass, myFile)
    for k,v in brain_mass_css_vals_bymass.items():
        brain_mass_css_vals_bymass[k] = np.average(v)


    with open("brain_seq_css_vals.txt", "wb") as myFile:
        pickle.dump(brain_seq_css_vals, myFile)
    for k,v in brain_seq_css_vals.items():
        brain_seq_css_vals[k] = np.average(v)

    with open("brain_seq_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(brain_seq_css_vals_bymass, myFile)
    for k,v in brain_seq_css_vals_bymass.items():
        brain_seq_css_vals_bymass[k] = np.average(v)

    for k,v in py_mass_css_vals.items():
        py_mass_css_vals[k] = np.average(v)


    aa_counts = list(fft_seq_times.keys())
    fft_seq_avgtimes = list(fft_seq_times.values())
    fft_mass_avgtimes = list(fft_mass_times.values())
    nn_seq_avgtimes = list(nn_seq_times.values())
    nn_mass_avgtimes = list(nn_mass_times.values())
    brain_mass_avgtimes = list(brain_mass_times.values())
    brain_seq_avgtimes = list(brain_seq_times.values())

    fft_mass_avgcss = list(fft_mass_css_vals.values())
    nn_mass_avgcss = list(nn_mass_css_vals.values())
    nn_seq_avgcss = list(nn_seq_css_vals.values())
    brain_mass_avgcss = list(brain_mass_css_vals.values())
    brain_seq_avgcss = list(brain_seq_css_vals.values())
    py_mass_avgcss = list(py_mass_css_vals.values())

    masses = list(fft_seq_times_bymass.keys())
    fft_seq_avgtimes_bymass = list(fft_seq_times_bymass.values())
    fft_mass_avgtimes_bymass = list(fft_mass_times_bymass.values())
    nn_seq_avgtimes_bymass = list(nn_seq_times_bymass.values())
    nn_mass_avgtimes_bymass = list(nn_mass_times_bymass.values())
    brain_mass_avgtimes_bymass = list(brain_mass_times_bymass.values())
    brain_seq_avgtimes_bymass = list(brain_seq_times_bymass.values())

    fft_mass_avgcss_bymass = list(fft_mass_css_vals_bymass.values())
    nn_mass_avgcss_bymass = list(nn_mass_css_vals_bymass.values())
    nn_seq_avgcss_bymass = list(nn_seq_css_vals_bymass.values())
    brain_mass_avgcss_bymass = list(brain_mass_css_vals_bymass.values())
    brain_seq_avgcss_bymass = list(brain_seq_css_vals_bymass.values())


    plt.scatter(aa_counts, fft_seq_avgtimes, label="FFT-Seq")
    plt.scatter(aa_counts, fft_mass_avgtimes, label="FFT-Mass")
    plt.scatter(aa_counts, nn_seq_avgtimes, label="NN-Seq")
    plt.scatter(aa_counts, nn_mass_avgtimes, label="NN-Mass")
    plt.scatter(aa_counts, brain_mass_avgtimes, label="BRAIN-Mass")
    plt.scatter(aa_counts, brain_seq_avgtimes, label="BRAIN-Seq")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlim(0,1000)
    plt.ylim(bottom=0)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel("AA Count", fontsize=label_size)
    plt.ylabel("Average Time (sec)", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.png", dpi=300, transparent=True)
    plt.show()


    aa_count_s_transitions = get_aacount_scount_array()

    plt.scatter(aa_counts, fft_mass_avgcss, label="FFT-Mass")
    plt.scatter(aa_counts, nn_seq_avgcss, label="NN-Seq")
    plt.scatter(aa_counts, nn_mass_avgcss, label="NN-Mass")
    plt.scatter(aa_counts, brain_mass_avgcss, label="BRAIN-Mass")
    plt.scatter(aa_counts, brain_seq_avgcss, label="BRAIN-Seq")
    # for aa_count in aa_count_s_transitions:
    #     plt.vlines(aa_count, 0.995, 1, color='r', linestyle='-')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 1000)
    plt.ylim(top=1)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel("AA Count",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_css.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_css.png", dpi=300, transparent=True)
    plt.show()


    plt.scatter(masses, fft_seq_avgtimes_bymass, label="FFT-Seq")
    plt.scatter(masses, fft_mass_avgtimes_bymass, label="FFT-Mass")
    plt.scatter(masses, nn_seq_avgtimes_bymass, label="NN-Seq")
    plt.scatter(masses, nn_mass_avgtimes_bymass, label="NN-Mass")
    plt.scatter(masses, brain_mass_avgtimes_bymass, label="BRAIN-Mass")
    plt.scatter(masses, brain_seq_avgtimes_bymass, label="BRAIN-Seq")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylim(bottom=0)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel("Mass (Da)", fontsize=label_size)
    plt.ylabel("Average Time (sec)", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.png", dpi=300, transparent=True)
    plt.show()


    aa_count_s_transitions = get_aacount_scount_array()

    plt.scatter(masses, fft_mass_avgcss_bymass, label="FFT-Mass")
    plt.scatter(masses, nn_seq_avgcss_bymass, label="NN-Seq")
    plt.scatter(masses, nn_mass_avgcss_bymass, label="NN-Mass")
    plt.scatter(masses, brain_mass_avgcss_bymass, label="BRAIN-Mass")
    plt.scatter(masses, brain_seq_avgcss_bymass, label="BRAIN-Seq")

    for aa_count in aa_count_s_transitions:
        plt.vlines(aa_count, 0.995, 1, color='r', linestyle='-')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.ylim(top=1)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel("Mass (Da)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_css.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_css.png", dpi=300, transparent=True)
    plt.show()




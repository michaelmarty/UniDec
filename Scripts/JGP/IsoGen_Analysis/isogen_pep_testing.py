import collections

import matplotlib.pyplot as plt
import molmass
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
import pickle
import math
import pyopenms

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

def get_mass_s_count_array():
    masses = np.arange(0, 100000, 10)
    transition_masses = []
    previous_s_count = 0
    for mass in masses:
        s_count = get_s_count_from_pepmass(mass)
        if s_count > previous_s_count:
            transition_masses.append(mass)
            previous_s_count = s_count
    return transition_masses

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

def molmass_pep_seq_to_dist(sequence):
    seq = pyopenms.AASequence(sequence)
    formula = seq.getFormula().toString()
    dist = get_dist_from_formula(formula)
    dist = np.array(dist) / max(dist)
    return dist

def molmass_pep_mass_to_dist(mass):
    formula, minmassint, intnum = pep_makemass(mass)
    dist = get_dist_from_formula(formula)
    dist = np.array(dist) / max(dist)
    return dist

def get_key_from_mass(mass, bin_width):
    key = int(math.floor(mass / bin_width) * bin_width)
    return key

def order_dict_by_aacount(dict):
    o_dict = collections.OrderedDict(sorted(dict.items()))
    return o_dict

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

include_molmass = False

if True:
    topdir = r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Assessment"
    os.chdir(topdir)
    #assessment_files = match_files(topdir, "*1000.npz")
    assessment_files = ["human_protein_seqs.npz", "mouse_protein_seqs.npz", "ecoli_protein_seqs.npz", "yeast_protein_seqs.npz"]
    #assessment_files = ["training_random_human_proteins_5_min_5_max_1000.npz"]

    seqs = []
    dists = []
    vecs = []
    masses = []

    new_avgine_comp = np.array([5.117756, 7.982229, 1.431965, 1.487561, 0.040586])
    new_avgine_mass = 114.60070342105581

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
    molmass_mass_css_vals = {}
    molmass_seq_css_vals = {}


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

    molmass_mass_css_vals_bymass = {}
    molmass_seq_css_vals_bymass = {}

    nn_mass_human_css_vals_bymass = {}

    fft_newavgine_css_vals_bymass = {}


    shuffled_indices = np.arange(0, len(masses))
    random.shuffle(shuffled_indices)


    num_processed = 0
    mass_bin_width = 1000
    curr_length = 50

    eng_mass = IsoGenMassEngine()

    for i in shuffled_indices:
        num_processed += 1
        percent_processed = round((num_processed / len(masses)) * 100, 2)
        if num_processed % 500 == 0:
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
        brain_mass_dist = brain_pep_mass_to_dist(masses[i])
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
        brain_seq_dist = brain_pep_seq_to_dist(seqs[i])
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

        if include_molmass:
            molmass_mass_dist = molmass_pep_mass_to_dist(masses[i])
            molmass_mass_css_val = calculate_cosinesimilarity(molmass_mass_dist, fft_seq_dist)
            if aa_count in molmass_mass_css_vals:
                molmass_mass_css_vals[aa_count].append(molmass_mass_css_val)
            else:
                molmass_mass_css_vals[aa_count] = [molmass_mass_css_val]

            if mass_key in molmass_mass_css_vals_bymass:
                molmass_mass_css_vals_bymass[mass_key].append(molmass_mass_css_val)
            else:
                molmass_mass_css_vals_bymass[mass_key] = [molmass_mass_css_val]

            molmass_seq_dist = molmass_pep_seq_to_dist(seqs[i])
            molmass_seq_css_val = calculate_cosinesimilarity(molmass_seq_dist, fft_seq_dist)
            if aa_count in molmass_seq_css_vals:
                molmass_seq_css_vals[aa_count].append(molmass_seq_css_val)
            else:
                molmass_seq_css_vals[aa_count] = [molmass_seq_css_val]

            if mass_key in molmass_seq_css_vals_bymass:
                molmass_seq_css_vals_bymass[mass_key].append(molmass_seq_css_val)
            else:
                molmass_seq_css_vals_bymass[mass_key] = [molmass_seq_css_val]


        nn_human_mass_dist = eng_mass.predict(masses[i])
        nn_human_mass_css = calculate_cosinesimilarity(nn_human_mass_dist, fft_seq_dist)
        if mass_key in nn_mass_human_css_vals_bymass:
            nn_mass_human_css_vals_bymass[mass_key].append(nn_human_mass_css)
        else:
            nn_mass_human_css_vals_bymass[mass_key] = [nn_human_mass_css]


        fft_newavgine_mass_dist = pepmass_to_dist_avgine_input(masses[i], new_avgine_comp, new_avgine_mass)
        fft_newavgine_mass_css = calculate_cosinesimilarity(fft_newavgine_mass_dist, fft_seq_dist)
        if mass_key in fft_newavgine_css_vals_bymass:
            fft_newavgine_css_vals_bymass[mass_key].append(fft_newavgine_mass_css)
        else:
            fft_newavgine_css_vals_bymass[mass_key] = [fft_newavgine_mass_css]

        # if aa_count == curr_length:
        #     plt.plot(fft_seq_dist, label="FFT-Seq")
        #     plt.plot(fft_mass_dist, label="FFT-Mass")
        #     plt.plot(nn_seq_dist, label="NN-Seq")
        #     plt.plot(nn_mass_dist, label="NN-Mass")
        #     plt.plot(molmass_mass_dist, label='molmass-Mass')
        #     plt.plot(molmass_seq_dist, label='molmass-Seq')
        #     plt.title("Length: " + str(curr_length) + " Mass: " + str(masses[i]))
        #     plt.legend()
        #     plt.show()
        #
        #
        #     curr_length += 20

    fft_seq_times = order_dict_by_aacount(fft_seq_times)
    with open("fft_seq_times.txt", "wb") as myFile:
        pickle.dump(fft_seq_times, myFile)
    for k,v in fft_seq_times.items():
        fft_seq_times[k] = np.average(remove_outliers(v))

    fft_seq_times_bymass = order_dict_by_aacount(fft_seq_times_bymass)
    with open("fft_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(fft_seq_times_bymass, myFile)
    for k,v in fft_seq_times_bymass.items():
        fft_seq_times_bymass[k] = np.average(remove_outliers(v))

    fft_mass_times = order_dict_by_aacount(fft_mass_times)
    with open("fft_mass_times.txt", "wb") as myFile:
        pickle.dump(fft_mass_times, myFile)
    for k,v in fft_mass_times.items():
        fft_mass_times[k] = np.average(remove_outliers(v))

    fft_mass_times_bymass = order_dict_by_aacount(fft_mass_times_bymass)
    with open("fft_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(fft_mass_times_bymass, myFile)
    for k,v in fft_mass_times_bymass.items():
        fft_mass_times_bymass[k] = np.average(remove_outliers(v))


    nn_seq_times = order_dict_by_aacount(nn_seq_times)
    with open("nn_seq_times.txt", "wb") as myFile:
        pickle.dump(nn_seq_times, myFile)
    for k,v in nn_seq_times.items():
        nn_seq_times[k] = np.average(remove_outliers(v))

    nn_seq_times_bymass = order_dict_by_aacount(nn_seq_times_bymass)
    with open("nn_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(nn_seq_times_bymass, myFile)
    for k,v in nn_seq_times_bymass.items():
        nn_seq_times_bymass[k] = np.average(remove_outliers(v))

    nn_mass_times = order_dict_by_aacount(nn_mass_times)
    with open("nn_mass_times.txt", "wb") as myFile:
        pickle.dump(nn_mass_times, myFile)
    for k,v in nn_mass_times.items():
        nn_mass_times[k] = np.average(remove_outliers(v))

    nn_mass_times_bymass = order_dict_by_aacount(nn_mass_times_bymass)
    with open("nn_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(nn_mass_times_bymass, myFile)
    for k,v in nn_mass_times_bymass.items():
        nn_mass_times_bymass[k] = np.average(remove_outliers(v))

    fft_mass_css_vals = order_dict_by_aacount(fft_mass_css_vals)
    with open("fft_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(fft_mass_css_vals, myFile)
    for k,v in fft_mass_css_vals.items():
        fft_mass_css_vals[k] = np.average(v)

    fft_mass_css_vals_bymass = order_dict_by_aacount(fft_mass_css_vals_bymass)
    with open("fft_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(fft_mass_css_vals_bymass, myFile)
    for k,v in fft_mass_css_vals_bymass.items():
        fft_mass_css_vals_bymass[k] = np.average(v)


    nn_seq_css_vals = order_dict_by_aacount(nn_seq_css_vals)
    with open("nn_seq_css_vals.txt", "wb") as myFile:
        pickle.dump(nn_seq_css_vals, myFile)
    for k,v in nn_seq_css_vals.items():
        nn_seq_css_vals[k] = np.average(v)

    nn_seq_css_vals_bymass = order_dict_by_aacount(nn_seq_css_vals_bymass)
    with open("nn_seq_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(nn_seq_css_vals_bymass, myFile)
    for k,v in nn_seq_css_vals_bymass.items():
        nn_seq_css_vals_bymass[k] = np.average(v)

    nn_mass_css_vals = order_dict_by_aacount(nn_mass_css_vals)
    with open("nn_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(nn_mass_css_vals, myFile)
    for k,v in nn_mass_css_vals.items():
        nn_mass_css_vals[k] = np.average(v)

    nn_mass_css_vals_bymass = order_dict_by_aacount(nn_mass_css_vals_bymass)
    with open("nn_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(nn_mass_css_vals_bymass, myFile)
    for k,v in nn_mass_css_vals_bymass.items():
        nn_mass_css_vals_bymass[k] = np.average(v)

    brain_mass_times = order_dict_by_aacount(brain_mass_times)
    with open("brain_mass_times.txt", "wb") as myFile:
        pickle.dump(brain_mass_times, myFile)
    for k,v in brain_mass_times.items():
        brain_mass_times[k] = np.average(remove_outliers(v))

    brain_mass_times_bymass = order_dict_by_aacount(brain_mass_times_bymass)
    with open("brain_mass_times_bymass.txt", "wb") as myFile:
        pickle.dump(brain_mass_times_bymass, myFile)
    for k,v in brain_mass_times_bymass.items():
        brain_mass_times_bymass[k] = np.average(remove_outliers(v))

    brain_seq_times = order_dict_by_aacount(brain_seq_times)
    with open("brain_seq_times.txt", "wb") as myFile:
        pickle.dump(brain_seq_times, myFile)
    for k,v in brain_seq_times.items():
        brain_seq_times[k] = np.average(remove_outliers(v))

    brain_seq_times_bymass = order_dict_by_aacount(brain_seq_times_bymass)
    with open("brain_seq_times_bymass.txt", "wb") as myFile:
        pickle.dump(brain_seq_times_bymass, myFile)
    for k,v in brain_seq_times_bymass.items():
        brain_seq_times_bymass[k] = np.average(remove_outliers(v))

    brain_mass_css_vals = order_dict_by_aacount(brain_mass_css_vals)
    with open("brain_mass_css_vals.txt", "wb") as myFile:
        pickle.dump(brain_mass_css_vals, myFile)
    for k,v in brain_mass_css_vals.items():
        brain_mass_css_vals[k] = np.average(v)

    brain_mass_css_vals_bymass = order_dict_by_aacount(brain_mass_css_vals_bymass)
    with open("brain_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(brain_mass_css_vals_bymass, myFile)
    for k,v in brain_mass_css_vals_bymass.items():
        brain_mass_css_vals_bymass[k] = np.average(v)

    brain_seq_css_vals = order_dict_by_aacount(brain_seq_css_vals)
    with open("brain_seq_css_vals.txt", "wb") as myFile:
        pickle.dump(brain_seq_css_vals, myFile)
    for k,v in brain_seq_css_vals.items():
        brain_seq_css_vals[k] = np.average(v)

    brain_seq_css_vals_bymass = order_dict_by_aacount(brain_seq_css_vals_bymass)
    with open("brain_seq_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(brain_seq_css_vals_bymass, myFile)
    for k,v in brain_seq_css_vals_bymass.items():
        brain_seq_css_vals_bymass[k] = np.average(v)

    if include_molmass:
        with open("molmass_mass_css_vals_bymass.txt", "wb") as myFile:
            pickle.dump(molmass_mass_css_vals_bymass, myFile)
        for k, v in molmass_mass_css_vals_bymass.items():
            molmass_mass_css_vals_bymass[k] = np.average(v)

        with open("molmass_mass_css_vals.txt", "wb") as myFile:
            pickle.dump(molmass_mass_css_vals, myFile)
        for k, v in molmass_mass_css_vals.items():
            molmass_mass_css_vals[k] = np.average(v)

        with open("molmass_seq_css_vals_bymass.txt", "wb") as myFile:
            pickle.dump(molmass_seq_css_vals_bymass, myFile)
        for k, v in molmass_seq_css_vals_bymass.items():
            molmass_seq_css_vals_bymass[k] = np.average(v)

        with open("molmass_seq_css_vals.txt", "wb") as myFile:
            pickle.dump(molmass_seq_css_vals, myFile)
        for k, v in molmass_seq_css_vals.items():
            molmass_seq_css_vals[k] = np.average(v)

    nn_mass_human_css_vals_bymass = order_dict_by_aacount(nn_mass_human_css_vals_bymass)
    with open("nn_human_mass_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(nn_mass_human_css_vals_bymass, myFile)
    for k,v in nn_mass_human_css_vals_bymass.items():
        nn_mass_human_css_vals_bymass[k] = np.average(v)

    fft_newavgine_css_vals_bymass = order_dict_by_aacount(fft_newavgine_css_vals_bymass)
    with open("fft_newavgine_css_vals_bymass.txt", "wb") as myFile:
        pickle.dump(fft_newavgine_css_vals_bymass, myFile)
    for k, v in fft_newavgine_css_vals_bymass.items():
        fft_newavgine_css_vals_bymass[k] = np.average(v)


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
    if include_molmass:
        molmass_mass_avgcss = list(molmass_mass_css_vals.values())
        molmass_seq_avgcss = list(molmass_seq_css_vals.values())
        molmass_mass_avgcss_bymass = list(molmass_mass_css_vals_bymass.values())
        molmass_seq_avgcss_bymass = list(molmass_seq_css_vals_bymass.values())

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
    nn_human_mass_avgcss_bymass = list(nn_mass_human_css_vals_bymass.values())


    if False:
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

    if False:
        plt.scatter(aa_counts, fft_mass_avgcss, label="FFT-Mass")
        plt.scatter(aa_counts, nn_seq_avgcss, label="NN-Seq")
        plt.scatter(aa_counts, nn_mass_avgcss, label="NN-Mass")
        plt.scatter(aa_counts, brain_mass_avgcss, label="BRAIN-Mass")
        plt.scatter(aa_counts, brain_seq_avgcss, label="BRAIN-Seq")
        if include_molmass:
            plt.scatter(aa_counts, molmass_mass_avgcss, label="molmass-Mass")
            plt.scatter(aa_counts, molmass_seq_avgcss, label="molmass-Seq")
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

    if False:
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
    mass_s_transitions = get_mass_s_count_array()

    plt.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass")
    #plt.plot(masses, nn_seq_avgcss_bymass, label="NN-Seq")
    #plt.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass")
    #plt.plot(masses, brain_mass_avgcss_bymass, label="BRAIN-Mass")
    #plt.plot(masses, brain_seq_avgcss_bymass, label="BRAIN-Seq")
    plt.plot(masses, nn_human_mass_avgcss_bymass, label="NN-Human-Mass")
    for mass in mass_s_transitions:
        plt.vlines(mass, 0.975, 1, color='r', linestyle='-')
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




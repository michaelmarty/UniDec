import pickle
from matplotlib.pyplot import rcParams
import matplotlib as mpl
import os
import numpy as np
from unidec.IsoDec.IsoGen.isogen_tools import pep_makemass, massavgine
import matplotlib.pyplot as plt
import collections

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

def get_mass_s_count_array():
    masses = np.arange(0, 100000, 10)
    transition_masses = []
    previous_s_count = 0
    for mass in masses:
        s_count = get_s_count_from_pepmass(mass)
        if s_count > previous_s_count:
            transition_masses.append(mass)
            previous_s_count = s_count
    transition_masses = np.array(transition_masses)/ 1e3
    return transition_masses

def order_dict_by_aacount(dict):
    o_dict = collections.OrderedDict(sorted(dict.items()))
    return o_dict

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



save_figs = True

include_molmass = False

topdir = r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Assessment"
os.chdir(topdir)

print("Loading data...")

#region Load Data
with open("fft_seq_times.txt", "rb") as myFile:
    fft_seq_times = pickle.load(myFile)
    fft_seq_times = order_dict_by_aacount(fft_seq_times)

with open("fft_mass_times.txt", "rb") as myFile:
    fft_mass_times = pickle.load(myFile)
    fft_mass_times = order_dict_by_aacount(fft_mass_times)

with open("nn_seq_times.txt", "rb") as myFile:
    nn_seq_times = pickle.load(myFile)
    nn_seq_times = order_dict_by_aacount(nn_seq_times)

with open("nn_mass_times.txt", "rb") as myFile:
    nn_mass_times = pickle.load(myFile)
    nn_mass_times = order_dict_by_aacount(nn_mass_times)

with open("brain_seq_times.txt", "rb") as myFile:
    brain_seq_times = pickle.load(myFile)
    brain_seq_times = order_dict_by_aacount(brain_seq_times)

with open("brain_mass_times.txt", "rb") as myFile:
    brain_mass_times = pickle.load(myFile)
    brain_mass_times = order_dict_by_aacount(brain_mass_times)

with open("fft_mass_css_vals.txt", "rb") as myFile:
    fft_mass_css_vals = pickle.load(myFile)
    fft_mass_css_vals = order_dict_by_aacount(fft_mass_css_vals)

with open("nn_seq_css_vals.txt", "rb") as myFile:
    nn_seq_css_vals = pickle.load(myFile)
    nn_seq_css_vals = order_dict_by_aacount(nn_seq_css_vals)

with open("nn_mass_css_vals.txt", "rb") as myFile:
    nn_mass_css_vals = pickle.load(myFile)
    nn_mass_css_vals = order_dict_by_aacount(nn_mass_css_vals)

with open("brain_seq_css_vals.txt", "rb") as myFile:
    brain_seq_css_vals = pickle.load(myFile)
    brain_seq_css_vals = order_dict_by_aacount(brain_seq_css_vals)

with open("brain_mass_css_vals.txt", "rb") as myFile:
    brain_mass_css_vals = pickle.load(myFile)
    brain_mass_css_vals = order_dict_by_aacount(brain_mass_css_vals)


with open("fft_seq_times_bymass.txt", "rb") as myFile:
    fft_seq_times_bymass = pickle.load(myFile)
    fft_seq_times_bymass = order_dict_by_aacount(fft_seq_times_bymass)

with open("fft_mass_times_bymass.txt", "rb") as myFile:
    fft_mass_times_bymass = pickle.load(myFile)
    fft_mass_times_bymass = order_dict_by_aacount(fft_mass_times_bymass)

with open("nn_seq_times_bymass.txt", "rb") as myFile:
    nn_seq_times_bymass = pickle.load(myFile)
    nn_seq_times_bymass = order_dict_by_aacount(nn_seq_times_bymass)

with open("nn_mass_times_bymass.txt", "rb") as myFile:
    nn_mass_times_bymass = pickle.load(myFile)
    nn_mass_times_bymass = order_dict_by_aacount(nn_mass_times_bymass)

with open("brain_seq_times_bymass.txt", "rb") as myFile:
    brain_seq_times_bymass = pickle.load(myFile)
    brain_seq_times_bymass = order_dict_by_aacount(brain_seq_times_bymass)

with open("brain_mass_times_bymass.txt", "rb") as myFile:
    brain_mass_times_bymass = pickle.load(myFile)
    brain_mass_times_bymass = order_dict_by_aacount(brain_mass_times_bymass)

with open("fft_mass_css_vals_bymass.txt", "rb") as myFile:
    fft_mass_css_vals_bymass = pickle.load(myFile)
    fft_mass_css_vals_bymass = order_dict_by_aacount(fft_mass_css_vals_bymass)

with open("nn_seq_css_vals_bymass.txt", "rb") as myFile:
    nn_seq_css_vals_bymass = pickle.load(myFile)
    nn_seq_css_vals_bymass = order_dict_by_aacount(nn_seq_css_vals_bymass)

with open("nn_mass_css_vals_bymass.txt", "rb") as myFile:
    nn_mass_css_vals_bymass = pickle.load(myFile)
    nn_mass_css_vals_bymass = order_dict_by_aacount(nn_mass_css_vals_bymass)

with open("brain_seq_css_vals_bymass.txt", "rb") as myFile:
    brain_seq_css_vals_bymass = pickle.load(myFile)
    brain_seq_css_vals_bymass = order_dict_by_aacount(brain_seq_css_vals_bymass)

with open("brain_mass_css_vals_bymass.txt", "rb") as myFile:
    brain_mass_css_vals_bymass = pickle.load(myFile)
    brain_mass_css_vals_bymass = order_dict_by_aacount(brain_mass_css_vals_bymass)

with open("nn_human_mass_css_vals_bymass.txt", "rb") as myFile:
    nn_human_mass_css_vals_bymass = pickle.load(myFile)
    nn_human_mass_css_vals_bymass = order_dict_by_aacount(nn_human_mass_css_vals_bymass)

with open("fft_newavgine_css_vals_bymass.txt", "rb") as myFile:
    fft_newavgine_css_vals_bymass = pickle.load(myFile)
    fft_newavgine_css_vals_bymass = order_dict_by_aacount(fft_newavgine_css_vals_bymass)
#endregion Load Data

print("Averaging and removing outliers...")

#region Average and Remove Outliers
for k,v in fft_seq_times.items():
    fft_seq_times[k] = np.average(remove_outliers(v))

for k,v in fft_mass_times.items():
    fft_mass_times[k] = np.average(remove_outliers(v))

for k,v in nn_seq_times.items():
    nn_seq_times[k] = np.average(remove_outliers(v))

for k,v in nn_mass_times.items():
    nn_mass_times[k] = np.average(remove_outliers(v))

for k,v in fft_mass_css_vals.items():
    fft_mass_css_vals[k] = np.average(v)

for k,v in nn_seq_css_vals.items():
    nn_seq_css_vals[k] = np.average(v)

for k,v in nn_mass_css_vals.items():
    nn_mass_css_vals[k] = np.average(v)


for k,v in fft_seq_times_bymass.items():
    fft_seq_times_bymass[k] = np.average(remove_outliers(v))

for k,v in fft_mass_times_bymass.items():
    fft_mass_times_bymass[k] = np.average(remove_outliers(v))

for k,v in nn_seq_times_bymass.items():
    nn_seq_times_bymass[k] = np.average(remove_outliers(v))

for k,v in nn_mass_times_bymass.items():
    nn_mass_times_bymass[k] = np.average(remove_outliers(v))

for k,v in fft_mass_css_vals_bymass.items():
    fft_mass_css_vals_bymass[k] = np.average(v)

for k,v in nn_seq_css_vals_bymass.items():
    nn_seq_css_vals_bymass[k] = np.average(v)

for k,v in nn_mass_css_vals_bymass.items():
    nn_mass_css_vals_bymass[k] = np.average(v)

for k,v in brain_mass_times.items():
    brain_mass_times[k] = np.average(remove_outliers(v))

for k,v in brain_seq_times.items():
    brain_seq_times[k] = np.average(remove_outliers(v))

for k,v in brain_mass_times_bymass.items():
    brain_mass_times_bymass[k] = np.average(remove_outliers(v))

for k,v in brain_seq_times_bymass.items():
    brain_seq_times_bymass[k] = np.average(remove_outliers(v))

for k,v in brain_mass_css_vals.items():
    brain_mass_css_vals[k] = np.average(remove_outliers(v))

for k,v in brain_seq_css_vals.items():
    brain_seq_css_vals[k] = np.average(remove_outliers(v))

keys_to_remove = []
for k,v in brain_mass_css_vals_bymass.items():
    brain_mass_css_vals_bymass[k] = np.average(v)

for k,v in brain_seq_css_vals_bymass.items():
    brain_seq_css_vals_bymass[k] = np.average(v)

for k,v in nn_human_mass_css_vals_bymass.items():
    nn_human_mass_css_vals_bymass[k] = np.average(v)

for k,v in fft_newavgine_css_vals_bymass.items():
    fft_newavgine_css_vals_bymass[k] = np.average(v)


if include_molmass:
    with open("molmass_seq_css_vals_bymass.txt", "rb") as myFile:
        molmass_seq_css_vals_bymass = pickle.load(myFile)
        molmass_seq_css_vals_bymass = order_dict_by_aacount(molmass_seq_css_vals_bymass)

    with open("molmass_mass_css_vals_bymass.txt", "rb") as myFile:
        molmass_mass_css_vals_bymass = pickle.load(myFile)
        molmass_mass_css_vals_bymass = order_dict_by_aacount(brain_mass_css_vals_bymass)

    with open("molmass_seq_css_vals.txt", "rb") as myFile:
        molmass_seq_css_vals = pickle.load(myFile)
        molmass_seq_css_vals = order_dict_by_aacount(molmass_seq_css_vals)

    with open("molmass_mass_css_vals.txt", "rb") as myFile:
        molmass_mass_css_vals = pickle.load(myFile)
        molmass_mass_css_vals = order_dict_by_aacount(molmass_mass_css_vals)

    for k, v in molmass_mass_css_vals.items():
        molmass_mass_css_vals[k] = np.average(remove_outliers(v))

    for k, v in molmass_seq_css_vals.items():
        molmass_seq_css_vals[k] = np.average(remove_outliers(v))

    for k, v in molmass_mass_css_vals_bymass.items():
        molmass_mass_css_vals_bymass[k] = np.average(v)

    for k, v in molmass_seq_css_vals_bymass.items():
        molmass_seq_css_vals_bymass[k] = np.average(v)

    molmass_mass_avgcss = list(brain_mass_css_vals.values())
    molmass_seq_avgcss = list(brain_seq_css_vals.values())
    molmass_mass_avgcss_bymass = list(brain_mass_css_vals_bymass.values())
    molmass_seq_avgcss_bymass = list(brain_seq_css_vals_bymass.values())

aa_counts = list(fft_seq_times.keys())

masses = np.array(list(fft_seq_times_bymass.keys()))/1e3

fft_seq_avgtimes = list(fft_seq_times.values())
fft_seq_avgtimes = [x * 1000000.0 for x in fft_seq_avgtimes]

fft_seq_avgtimes_bymass = list(fft_seq_times_bymass.values())
fft_seq_avgtimes_bymass = [x * 1000000.0 for x in fft_seq_avgtimes_bymass]


fft_mass_avgtimes = list(fft_mass_times.values())
fft_mass_avgtimes = [x * 1000000.0 for x in fft_mass_avgtimes]

fft_mass_avgtimes_bymass = list(fft_mass_times_bymass.values())
fft_mass_avgtimes_bymass = [x * 1000000.0 for x in fft_mass_avgtimes_bymass]


nn_seq_avgtimes = list(nn_seq_times.values())
nn_seq_avgtimes = [x * 1000000.0 for x in nn_seq_avgtimes]

nn_seq_avgtimes_bymass = list(nn_seq_times_bymass.values())
nn_seq_avgtimes_bymass = [x * 1000000.0 for x in nn_seq_avgtimes_bymass]


nn_mass_avgtimes = list(nn_mass_times.values())
nn_mass_avgtimes = [x * 1000000.0 for x in nn_mass_avgtimes]

nn_mass_avgtimes_bymass = list(nn_mass_times_bymass.values())
nn_mass_avgtimes_bymass = [x * 1000000.0 for x in nn_mass_avgtimes_bymass]


brain_mass_avgtimes = list(brain_mass_times.values())
brain_mass_avgtimes = [x * 1000000.0 for x in brain_mass_avgtimes]

brain_mass_avgtimes_bymass = list(brain_mass_times_bymass.values())
brain_mass_avgtimes_bymass = [x * 1000000.0 for x in brain_mass_avgtimes_bymass]


brain_seq_avgtimes = list(brain_seq_times.values())
brain_seq_avgtimes = [x * 1000000.0 for x in brain_seq_avgtimes]

brain_seq_avgtimes_bymass = list(brain_seq_times_bymass.values())
brain_seq_avgtimes_bymass = [x * 1000000.0 for x in brain_seq_avgtimes_bymass]


fft_mass_avgcss = list(fft_mass_css_vals.values())
nn_mass_avgcss = list(nn_mass_css_vals.values())
nn_seq_avgcss = list(nn_seq_css_vals.values())
brain_mass_avgcss = list(brain_mass_css_vals.values())
brain_seq_avgcss = list(brain_seq_css_vals.values())


fft_mass_avgcss_bymass = list(fft_mass_css_vals_bymass.values())
nn_mass_avgcss_bymass = list(nn_mass_css_vals_bymass.values())
nn_seq_avgcss_bymass = list(nn_seq_css_vals_bymass.values())
brain_mass_avgcss_bymass = list(brain_mass_css_vals_bymass.values())
brain_seq_avgcss_bymass = list(brain_seq_css_vals_bymass.values())
nn_human_mass_avgcss_bymass = list(nn_human_mass_css_vals_bymass.values())
fft_newavgine_avgcss_bymass = list(fft_newavgine_css_vals_bymass.values())
#endregion Average and Remove Outliers

transition_masses = get_mass_s_count_array()

# Plot by AA count:
if False:
    plt.plot(aa_counts, fft_seq_avgtimes, label="FFT-Seq", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(aa_counts, fft_mass_avgtimes, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
    plt.plot(aa_counts, nn_seq_avgtimes, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(aa_counts, nn_mass_avgtimes, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
    plt.plot(aa_counts, brain_mass_avgtimes, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
    plt.plot(aa_counts, brain_seq_avgtimes, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
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
    plt.ylabel("Average Time\n(microseconds)", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\aa_count_speed.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\aa_count_speed.png", dpi=300, transparent=True)
    plt.show()


    aa_count_s_transitions = get_aacount_scount_array()

    plt.plot(aa_counts, fft_mass_avgcss, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
    plt.plot(aa_counts, nn_seq_avgcss, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(aa_counts, nn_mass_avgcss, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
    plt.plot(aa_counts, brain_mass_avgcss, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
    plt.plot(aa_counts, brain_seq_avgcss, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])

    for aa_count in aa_count_s_transitions:
        plt.vlines(aa_count, 0.99, 1, color='r', linestyle='--')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 1000)
    plt.ylim(top=1.00001)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel("AA Count",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\aa_count_css.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\aa_count_css.png", dpi=300, transparent=True)
    plt.show()

#Plot by mass
if False:
    plt.plot(masses, fft_seq_avgtimes_bymass, label="FFT-Seq", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(masses, fft_mass_avgtimes_bymass, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
    plt.plot(masses, nn_seq_avgtimes_bymass, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, nn_mass_avgtimes_bymass, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
    plt.plot(masses, brain_mass_avgtimes_bymass, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
    plt.plot(masses, brain_seq_avgtimes_bymass, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlim(0, 100)
    plt.ylim(0, 400)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xticks(fontsize=tick_size)
    plt.yticks([0,100,200,300,400], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)", fontsize=label_size)
    plt.ylabel("Average Time\n(microseconds)", fontsize=label_size)
    plt.legend(ncol=3, frameon=True, edgecolor="black", fontsize=12)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pep_mass_v_speed.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pep_mass_v_speed.png", dpi=300, transparent=True)
    plt.show()

    mass_scount_transitions = np.array(get_mass_scount_array(list(fft_seq_times_bymass.keys())))/1e3

    plt.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
    plt.plot(masses, nn_seq_avgcss_bymass, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
    plt.plot(masses, brain_mass_avgcss_bymass, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
    plt.plot(masses, brain_seq_avgcss_bymass, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])

    # for mass in mass_scount_transitions:
    #     plt.vlines(mass, 0.95, 1, color='r', linestyle='--')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(bottom=0.98, top=1.001)
    plt.xticks(fontsize=tick_size)
    plt.yticks(ticks=[0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pep_mass_v_css.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\pep_mass_v_css.png", dpi=300, transparent=True)
    plt.show()

#Plot separately the mass and sequence methods by mass
if False:
    #Mass methods avg times
    plt.plot(masses, fft_mass_avgtimes_bymass, label="FFT", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(masses, nn_mass_avgtimes_bymass, label="NN", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, brain_mass_avgtimes_bymass, label="BRAIN", linewidth=linewidth, color=CB_color_cycle[5])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(0, 400)
    plt.xticks([0,25,50,75,100], fontsize=tick_size)
    plt.yticks([0,100,200,300,400], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)", fontsize=label_size)
    plt.ylabel("Average Time\n(microseconds)", fontsize=label_size)
    #plt.legend(ncol=3, frameon=True, edgecolor="black", fontsize=16, loc='upper center')
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_speed_massmethods.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_speed_massmethods.png", dpi=300, transparent=True)
    plt.show()

    #Seq methods avg times
    plt.plot(aa_counts, fft_seq_avgtimes, label="FFT-Seq", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(aa_counts, nn_seq_avgtimes, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(aa_counts, brain_seq_avgtimes, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 1000)
    plt.ylim(0, 400)
    plt.xticks([0,250,500,750,1000], fontsize=tick_size)
    plt.yticks([0,100,200,300, 400], fontsize=tick_size)
    plt.xlabel("Amino Acids", fontsize=label_size)
    plt.ylabel("Average Time\n(microseconds)", fontsize=label_size)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_aacount_v_speed_seqmethods.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_aacount_v_speed_seqmethods.png", dpi=300, transparent=True)
    plt.show()


    #Mass methods avgcss
    plt.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, brain_mass_avgcss_bymass, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[5])
    # for mass in transition_masses:
    #     plt.vlines(mass, 0.95, 1, linestyle='--', color='r')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(bottom=0.975, top=1.001)
    plt.xticks(ticks=[0, 25, 50, 75, 100], fontsize=tick_size)
    plt.yticks(ticks=[0.975, 0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods.png", dpi=300, transparent=True)
    plt.show()

    #seq methods avg css
    plt.plot(aa_counts, nn_seq_avgcss, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(aa_counts, brain_seq_avgcss, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 1000)
    plt.ylim(bottom=0.975, top=1.001)
    plt.xticks(ticks=[0, 250, 500, 750, 1000], fontsize=tick_size)
    plt.yticks(ticks=[0.975, 0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Amino Acids",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_aacount_v_css_seqmethods.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_aacount_v_css_seqmethods.png", dpi=300, transparent=True)
    plt.show()


#plot mass and sequence methods separately by mass, with molmass included, only CSS
if False:
    #Mass methods avgcss
    plt.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, brain_mass_avgcss_bymass, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[5])
    plt.plot(masses, molmass_mass_avgcss_bymass, label="molmass-Mass", linewidth=linewidth, color=CB_color_cycle[1])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(bottom=0.98, top=1.001)
    plt.xticks(fontsize=tick_size)
    plt.yticks(ticks=[0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods_molmass.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods_molmass.png", dpi=300, transparent=True)
    plt.show()

    #seq methods avg css
    plt.plot(masses, nn_seq_avgcss_bymass, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, brain_seq_avgcss_bymass, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
    plt.plot(masses, molmass_seq_avgcss_bymass, label="molmass-Seq", linewidth=linewidth, color=CB_color_cycle[1])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(bottom=0.98, top=1.001)
    plt.xticks(fontsize=tick_size)
    plt.yticks(ticks=[0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend()
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_seqmethods_molmass.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_seqmethods_molmass.png", dpi=300, transparent=True)
    plt.show()

#plot mass methods v avgcss including the specifically-trained NN-Mass
if True:
    #Mass methods avgcss
    plt.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass (Averagine)", linewidth=linewidth, color=CB_color_cycle[0])
    plt.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass (Averagine)", linewidth=linewidth, color=CB_color_cycle[2])
    plt.plot(masses, nn_human_mass_avgcss_bymass, label="NN-Mass (Select Proteomes)", linewidth=linewidth, color=CB_color_cycle[1])
    plt.plot(masses, fft_newavgine_avgcss_bymass, label="FFT-Mass (Select Proteomes)", linewidth=linewidth, color=CB_color_cycle[3])
    # for mass in transition_masses:
    #     plt.vlines(mass, 0.95, 1, linestyle='--', color='r')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.xticks()
    plt.tick_params(length=5, width=2)
    plt.xlim(0, 100)
    plt.ylim(bottom=0.975, top=1.005)
    plt.xticks(ticks=[0, 25, 50, 75, 100], fontsize=tick_size)
    plt.yticks(ticks=[0.975, 0.980, 0.985, 0.990, 0.995, 1.000], fontsize=tick_size)
    plt.xlabel("Protein Mass (kDa)",fontsize=label_size)
    plt.ylabel("CSS", fontsize=label_size)
    plt.legend(ncol=2, frameon=True, loc='upper center', edgecolor="black", fontsize=10)
    plt.tight_layout()
    if save_figs:
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods_spec_nnmass.pdf", dpi=300, transparent=True)
        plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\protein_mass_v_css_massmethods_spec_nnmass.png", dpi=300, transparent=True)
    plt.show()

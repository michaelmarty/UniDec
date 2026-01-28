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



save_figs = False

topdir = r"C:\Users\Admin\Desktop\martylab\IsoGen\IntactProtein\Assessment"
os.chdir(topdir)

print("Loading data...")

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


print("Removing outliers...")
for k,v in fft_seq_times.items():
    fft_seq_times[k] = np.average(remove_outliers(v))

for k,v in fft_mass_times.items():
    fft_mass_times[k] = np.average(remove_outliers(v))

for k,v in nn_seq_times.items():
    nn_seq_times[k] = np.average(remove_outliers(v))

for k,v in nn_mass_times.items():
    nn_mass_times[k] = np.average(remove_outliers(v))

for k,v in fft_mass_css_vals.items():
    fft_mass_css_vals[k] = np.average(remove_outliers(v))

for k,v in nn_seq_css_vals.items():
    nn_seq_css_vals[k] = np.average(remove_outliers(v))

for k,v in nn_mass_css_vals.items():
    nn_mass_css_vals[k] = np.average(remove_outliers(v))

for k,v in brain_mass_times.items():
    brain_mass_times[k] = np.average(remove_outliers(v))

for k,v in brain_seq_times.items():
    brain_seq_times[k] = np.average(remove_outliers(v))

for k,v in brain_mass_css_vals.items():
    brain_mass_css_vals[k] = np.average(remove_outliers(v))

for k,v in brain_seq_css_vals.items():
    brain_seq_css_vals[k] = np.average(remove_outliers(v))


aa_counts = list(fft_seq_times.keys())
fft_seq_avgtimes = list(fft_seq_times.values())
fft_seq_avgtimes = [x * 1000000.0 for x in fft_seq_avgtimes]

fft_mass_avgtimes = list(fft_mass_times.values())
fft_mass_avgtimes = [x * 1000000.0 for x in fft_mass_avgtimes]

nn_seq_avgtimes = list(nn_seq_times.values())
nn_seq_avgtimes = [x * 1000000.0 for x in nn_seq_avgtimes]

nn_mass_avgtimes = list(nn_mass_times.values())
nn_mass_avgtimes = [x * 1000000.0 for x in nn_mass_avgtimes]

brain_mass_avgtimes = list(brain_mass_times.values())
brain_mass_avgtimes = [x * 1000000.0 for x in brain_mass_avgtimes]

brain_seq_avgtimes = list(brain_seq_times.values())
brain_seq_avgtimes = [x * 1000000.0 for x in brain_seq_avgtimes]

fft_mass_avgcss = list(fft_mass_css_vals.values())
nn_mass_avgcss = list(nn_mass_css_vals.values())
nn_seq_avgcss = list(nn_seq_css_vals.values())
brain_mass_avgcss = list(brain_mass_css_vals.values())
brain_seq_avgcss = list(brain_seq_css_vals.values())


plt.plot(aa_counts, fft_seq_avgtimes, label="FFT-Seq", linewidth=linewidth, color=CB_color_cycle[0])
plt.plot(aa_counts, fft_mass_avgtimes, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
plt.plot(aa_counts, nn_seq_avgtimes, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
plt.plot(aa_counts, nn_mass_avgtimes, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
plt.plot(aa_counts, brain_mass_avgtimes, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
plt.plot(aa_counts, brain_seq_avgtimes, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlim(0,300)
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
    plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.pdf", dpi=300, transparent=True)
    plt.savefig(r"C:\Users\Admin\Desktop\martylab\IsoGen\Figures\aa_count_speed.png", dpi=300, transparent=True)
plt.show()


aa_count_s_transitions = get_aacount_scount_array()

plt.plot(aa_counts, fft_mass_avgcss, label="FFT-Mass", linewidth=linewidth, color=CB_color_cycle[1])
plt.plot(aa_counts, nn_seq_avgcss, label="NN-Seq", linewidth=linewidth, color=CB_color_cycle[2])
plt.plot(aa_counts, nn_mass_avgcss, label="NN-Mass", linewidth=linewidth, color=CB_color_cycle[3])
plt.plot(aa_counts, brain_mass_avgcss, label="BRAIN-Mass", linewidth=linewidth, color=CB_color_cycle[4])
plt.plot(aa_counts, brain_seq_avgcss, label="BRAIN-Seq", linewidth=linewidth, color=CB_color_cycle[5])

# for aa_count in aa_count_s_transitions:
#     plt.vlines(aa_count, 0.99, 1, color='r', linestyle='-')

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.xticks()
plt.tick_params(length=5, width=2)
plt.xlim(0, 300)
plt.ylim(top=1.001)
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
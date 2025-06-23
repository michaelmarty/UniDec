import numpy as np
import gc
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.IsoDec.IsoGen.isogenmass import IsoGenMassEngine
from unidec.IsoDec.IsoGen.isogenpep import *
from unidec.IsoDec.IsoGen.isogen_tools import *
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *
import fnmatch
import time


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenPepEngine()
eng_avg = IsoGenMassEngine()

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
    topdir = r"C:\Users\Admin\Desktop\martylab\IsoGen\IntactProtein\Assessment"
    os.chdir(topdir)
    assessment_files = match_files(topdir, "*0.npz")
    assessment_files = ["assessment_random_human_proteins_1000_min_5_max_200_mods0.npz"]

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


    fft_seq_times = {}
    fft_mass_times = {}
    nn_seq_times = {}
    nn_mass_times = {}

    fft_mass_css_vals = {}
    nn_seq_css_vals = {}
    nn_mass_css_vals = {}


    for i in range(len(masses)):
        if i % 10000 == 0:
            print("Processing sequence", i, "of", len(masses))

        aa_count = peptide_to_aacount(seqs[i])
        pep_isolen = pepmass_to_isolen(masses[i])
        rna_isolen = rnamass_to_isolen(masses[i])

        timestart = time.perf_counter()
        fft_seq_dist = fft_pep_seq_to_dist(seqs[i])
        fft_seq_time = round(time.perf_counter() - timestart, 6)

        if aa_count in fft_seq_times:
            fft_seq_times[aa_count].append(fft_seq_time)
        else:
            fft_seq_times[aa_count] = [fft_seq_time]


        timestart = time.perf_counter()
        fft_mass_dist = fft_pep_mass_to_dist(masses[i])
        fft_mass_time = round(time.perf_counter() - timestart, 6)
        fft_mass_css = calculate_pearson_chisquared(fft_mass_dist, fft_seq_dist)

        if aa_count in fft_mass_times:
            fft_mass_times[aa_count].append(fft_mass_time)
            fft_mass_css_vals[aa_count].append(fft_mass_css)
        else:
            fft_mass_times[aa_count] = [fft_mass_time]
            fft_mass_css_vals[aa_count] = [fft_mass_css]


        timestart = time.perf_counter()
        nn_seq_dist = nn_pep_seq_to_dist(seqs[i])
        nn_seq_time = round(time.perf_counter() - timestart, 6)
        nn_seq_css = calculate_pearson_chisquared(nn_seq_dist, fft_seq_dist)

        if aa_count in nn_seq_times:
            nn_seq_times[aa_count].append(nn_seq_time)
            nn_seq_css_vals[aa_count].append(nn_seq_css)
        else:
            nn_seq_times[aa_count] = [nn_seq_time]
            nn_seq_css_vals[aa_count] = [nn_seq_css]


        timestart = time.perf_counter()
        nn_mass_dist = nn_pep_mass_to_dist(masses[i])
        nn_mass_time = round(time.perf_counter() - timestart, 6)
        nn_mass_css = calculate_pearson_chisquared(nn_mass_dist, fft_seq_dist)

        if aa_count in nn_mass_times:
            nn_mass_times[aa_count].append(nn_mass_time)
            nn_mass_css_vals[aa_count].append(nn_mass_css)
        else:
            nn_mass_times[aa_count] = [nn_mass_time]
            nn_mass_css_vals[aa_count] = [nn_mass_css]

        # if i % 100000 == 0:
        #     plt.plot(fft_seq_dist, label="FFT-Seq")
        #     plt.plot(fft_mass_dist, label="FFT-Mass")
        #     plt.plot(nn_seq_dist, label="NN-Seq")
        #     plt.plot(nn_mass_dist, label="NN-Mass")
        #     plt.legend()
        #     plt.show()




    for k,v in fft_seq_times.items():
        fft_seq_times[k] = np.average(v)

    for k,v in fft_mass_times.items():
        fft_mass_times[k] = np.average(v)

    for k,v in nn_seq_times.items():
        nn_seq_times[k] = np.average(v)

    for k,v in nn_mass_times.items():
        nn_mass_times[k] = np.average(v)

    for k,v in fft_mass_css_vals.items():
        fft_mass_css_vals[k] = np.average(v)

    for k,v in nn_seq_css_vals.items():
        nn_seq_css_vals[k] = np.average(v)

    for k,v in nn_mass_css_vals.items():
        nn_mass_css_vals[k] = np.average(v)

    aa_counts = list(fft_seq_times.keys())
    fft_seq_avgtimes = list(fft_seq_times.values())
    fft_mass_avgtimes = list(fft_mass_times.values())
    nn_seq_avgtimes = list(nn_seq_times.values())
    nn_mass_avgtimes = list(nn_mass_times.values())
    fft_mass_avgcss = list(fft_mass_css_vals.values())
    nn_mass_avgcss = list(nn_mass_css_vals.values())
    nn_seq_avgcss = list(nn_seq_css_vals.values())

    plt.scatter(aa_counts, fft_seq_avgtimes, label="FFT-Seq")
    plt.scatter(aa_counts, fft_mass_avgtimes, label="FFT-Mass")
    plt.scatter(aa_counts, nn_seq_avgtimes, label="NN-Seq")
    plt.scatter(aa_counts, nn_mass_avgtimes, label="NN-Mass")
    plt.legend()
    plt.show()

    plt.scatter(aa_counts, fft_mass_avgcss, label="FFT-Mass")
    plt.scatter(aa_counts, nn_seq_avgcss, label="NN-Seq")
    plt.scatter(aa_counts, nn_mass_avgcss, label="NN-Mass")
    plt.legend()
    plt.show()




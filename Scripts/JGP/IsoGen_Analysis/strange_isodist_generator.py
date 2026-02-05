import numpy as np
from numpy import fft as fftpack
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.IsoDec.IsoGen.isogen_tools import *
import os

from unidec.IsoDec.IsoGen.isogenmass import IsoGenMassEngine

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

#This is the natural isotope ratio
#c = np.array([1, 0.01082, 0, 0])

#This is the ratio assuming complete isotopic/atomic scrambling
#c = np.array([1,0.0225,0,0])

#Let's tone down the scrambling to be somewhat more realistic..
c = np.array([1,0.013,0,0])

c_heavy = np.array([0,1,0,0])

#This is the natural isotope ratio
#n = np.array([1, 0.0037, 0, 0])

#This is the ratio assuming complete isotopic/atomic scrambling
#n = np.array([1,0.106,0,0])

#Let's tone down the scrambling to be somewhat more realistic..
n = np.array([1,0.008,0,0])

n_heavy = np.array([0,1,0,0])
o = np.array([1, 0.0004, 0.002, 0])
s = np.array([1, 0.0079, 0.044, 0])

h = np.append(h, buffer)
c = np.append(c, buffer)
c_heavy = np.append(c_heavy, buffer)
n = np.append(n, buffer)
n_heavy = np.append(n_heavy, buffer)
o = np.append(o, buffer)
s = np.append(s, buffer)

dt = np.dtype(np.complex128)
hft = fftpack.rfft(h).astype(dt)
cft = fftpack.rfft(c).astype(dt)
nft = fftpack.rfft(n).astype(dt)
oft = fftpack.rfft(o).astype(dt)
sft = fftpack.rfft(s).astype(dt)

pre_ft = np.array([cft, hft, nft, oft, sft])

def pepseq_to_isolist(seq):
    isolist = np.array([0,2,0,1,0])
    for aa in seq:
        if aa in aa_elementalcomp_dict:
            isolist += aa_elementalcomp_dict[aa]
        else:
            print("Unexpected amino acid found: ", str(aa))
    return isolist


def isojim_scrambled(isolist, length):
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


def calc_modified_isotope_ratio():
    """"
    Back-of-the-envelope calculation of a modified isotope ratio in
    """
    dmem_aa_comp_dict = {
        "G":0.4,
        "R":0.398,
        "C":0.201,
        "Q":4.0,
        "H":0.2,
        "I":0.802,
        "L":0.802,
        "K":0.798,
        "M":0.201,
        "F":0.4,
        "S":0.4,
        "T":0.798,
        "W":0.078,
        "Y":0.398,
        "V":0.803
    }

    sum_conc = 0
    sum_isolist = np.array([0.0,0.0,0.0,0.0,0.0])

    #We can sum the concentrations of each AA and basically treat them as moles since they're in the same volume
    for k,v in dmem_aa_comp_dict.items():
        sum_conc += v

    #Add in the dextrose
    sum_conc += 25

    for k,v in dmem_aa_comp_dict.items():
        sum_isolist += (v/sum_conc) * aa_elementalcomp_dict[k]


    dextrose = np.array([6,12,0,6,0])
    #Add in the dextrose (again)
    sum_isolist += (25/sum_conc) * dextrose

    #These are fixed isotopically, so must be subtracted from the sum isolist and treated differently
    #This will eventually allow for the calculation of a modified N and C isotope ratio in the SILAC-labelled proteins
    R_heavyN_contribution = (dmem_aa_comp_dict["R"] / sum_conc) * 4
    R_heavyC_contribution = (dmem_aa_comp_dict["R"] / sum_conc) * 6

    c_total = sum_isolist[0]
    h_total = sum_isolist[1]
    n_total = sum_isolist[2]
    o_total = sum_isolist[3]
    s_total = sum_isolist[4]

    sum_isolist[0] -= R_heavyC_contribution
    sum_isolist[2] -= R_heavyN_contribution

    C_isoratio = ((sum_isolist[0]/c_total) * c) + ((R_heavyC_contribution/c_total) * c_heavy)
    C_isoratio = C_isoratio / np.max(C_isoratio)
    print("Modified Carbon isotope ratio: ", str(C_isoratio))

    N_isoratio = ((sum_isolist[2] / n_total) * n) + ((R_heavyN_contribution / n_total) * n_heavy)
    N_isoratio = N_isoratio / np.max(N_isoratio)
    print("Modified Nitrogen isotope ratio: ", str(N_isoratio))


def add_noise_to_dist(dist, noise_level=0.2):
    noisy_dist = [0 for _ in dist]
    noise_multipliers = np.random.normal(loc = 1, scale=noise_level, size=len(dist))
    for i in range(len(dist)):
        noisy_dist[i] = dist[i] * noise_multipliers[i]
    return noisy_dist


if __name__ == "__main__":
    # calc_modified_isotope_ratio()
    # exit()

    if False:
        pepseq = "AAAVVVKRAAAC"
        isolist = pepseq_to_isolist(pepseq)

        standard_dist = fft_pep_seq_to_dist(pepseq)
        altered_dist = isojim(isolist, length=128)

        plt.plot(standard_dist, label="Standard")
        plt.plot(altered_dist, label="SILAC")
        plt.xlim(0)
        plt.ylim(0,1)
        plt.legend()
        plt.show()

    if False:
        pepmass = 4000

        formula,minmassint,isolist = pep_makemass(pepmass)

        standard_dist = fft_pep_mass_to_dist(pepmass)
        altered_dist = isojim_scrambled(isolist, 128)
        noisy_dist = add_noise_to_dist(altered_dist)

        for i in range(len(standard_dist)):
            plt.vlines(i+0.05,0,standard_dist[i], color='r')
        for i in range(len(altered_dist)):
            plt.vlines(i+0.1,0,altered_dist[i], color='b')
        for i in range(len(noisy_dist)):
            plt.vlines(i+0.15, 0, noisy_dist[i], color='g')

        plt.xlim(0)
        plt.ylim(0)
        plt.legend()
        plt.show()

    if True:
        os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Assessment")

        min_mass = 500
        max_mass = 3000
        n = 50000

        pepmasses = np.random.uniform(min_mass, max_mass, n)
        dists = [0 for _ in pepmasses]

        for i in range(len(pepmasses)):
            formula,minmassint,isolist = pep_makemass(pepmasses[i])
            dist = isojim_scrambled(isolist, length=128)
            dists[i] = dist

        dists = np.array(dists)
        np.savez_compressed(str(n) + "_scrambled_assessment_dists_min_" + str(min_mass) + "_max_" + str(max_mass) + ".npz", dists=dists, masses=pepmasses)

    if False:
        os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling")
        from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *
        eng = IsoGenMassEngine(isolen=32)

        trainfile = r"Training\1000_scrambled_training_dists_min_500_max_3000.npz"

        assessmentfile = r"Assessment\50000_scrambled_assessment_dists_min_500_max_3000.npz"
        assessmentdata = np.load(assessmentfile)
        assessment_masses = assessmentdata["masses"]
        assessment_dists = assessmentdata["dists"]

        untrained_css_vals = []
        for i in range(len(assessment_masses)):
            nn_dist = eng.predict(assessment_masses[i], isolen=32)
            css = calculate_cosinesimilarity(nn_dist, assessment_dists[i])
            untrained_css_vals.append(css)

        untrained_css_vals = np.array(untrained_css_vals)
        print("Average Untrained CSS:" + str(np.average(untrained_css_vals)))

        eng.transfer_train(trainfile, epochs=10, length=32)

        trained_css_vals = []

        for i in range(len(assessment_masses)):
            nn_dist = eng.predict(assessment_masses[i], isolen=32)
            css = calculate_cosinesimilarity(nn_dist, assessment_dists[i])
            trained_css_vals.append(css)

        trained_css_vals = np.array(trained_css_vals)
        print("Average transfer-trained CSS:" + str(np.average(trained_css_vals)))







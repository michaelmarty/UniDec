import matplotlib.pyplot as plt

from unidec.IsoDec.IsoGen.isogenpep import *
from unidec.IsoDec.IsoGen.isogen_tools import *
import time
import os
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

# Set backend to Agg
mpl.use('WxAgg')

unimod_obo_file = "C:\\Users\\Admin\\Desktop\\martylab\\IsoGen\\unimod.obo.txt"



eng = IsoGenPepEngine()

data = np.load("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\Training\\training_random_ecoli_proteins_10000_min_50_max_200_mods2.npz")

seqs = data["seqs"]
masses = data["masses"]
lengths = [len(seq) for seq in seqs]
dists = data["dists"]

# plt.hist(masses, bins=100)
# plt.show()

chis = []
nn_times = []
ft_times = []

for i, seq in enumerate(seqs):
    ai_tstart = time.perf_counter()
    pred_dist = eng.predict(seq)
    ai_tend = time.perf_counter()
    nn_times.append(ai_tend - ai_tstart)

    ft_tstart = time.perf_counter()
    true_dist = peptide_to_dist(seq)
    ft_tend = time.perf_counter()
    ft_times.append(ft_tend - ft_tstart)

    #Normalize the two distributions to peak at 1
    pred_dist = pred_dist / np.max(pred_dist)
    true_dist = true_dist / np.max(true_dist)

    mass = peptide_to_mass(seq)
    averagine_dist = mass_to_dist(mass)

    if i % 1000 == 0 and False:
        plt.plot(pred_dist * 100, label="AI", color="b")
        plt.plot(true_dist * 100, color="k", label="True")
        plt.legend()
        plt.title(str(len(seq)))
        plt.xlabel("Isotope Number")
        plt.ylabel("%")
        plt.xlim(0, 20)
        plt.show()

    # if len(seq) > 30 and len(seq) < 50:
    #     plt.plot(pred_dist * 100, label="AI", color="b")
    #     plt.plot(true_dist * 100, color="k", label="True")
    #     plt.legend()
    #     plt.title(str(len(seq)))
    #     plt.xlabel("Isotope Number")
    #     plt.ylabel("%")
    #     plt.xlim(0, 20)
    #     plt.show()

    chi = calculate_pearson_chisquared(pred_dist, true_dist)
    chis.append(chi)

chis_lengths = {}
nn_times_lengths = {}
ft_times_lengths = {}

for i in range(len(chis)):
    if lengths[i] in chis_lengths:
        chis_lengths[lengths[i]].append(chis[i])
        nn_times_lengths[lengths[i]].append(nn_times[i])
        ft_times_lengths[lengths[i]].append(ft_times[i])
    else:
        chis_lengths[lengths[i]] = [chis[i]]
        nn_times_lengths[lengths[i]] = [nn_times[i]]
        ft_times_lengths[lengths[i]] = [ft_times[i]]

meanchis_lengths = {}
mean_nntimes_lengths = {}
mean_fttimes_lengths = {}

for length in chis_lengths.keys():
    meanchis_lengths[length] = np.mean(np.array(chis_lengths[length]))
    mean_nntimes_lengths[length] = np.mean(np.array(nn_times_lengths[length]))
    mean_fttimes_lengths[length] = np.mean(np.array(ft_times_lengths[length]))

plt.scatter(meanchis_lengths.keys(), meanchis_lengths.values())
plt.xlabel("Peptide Length")
plt.ylabel("Mean Chi")
#plt.savefig("Z:\\Group Share\\JGP\\IsoGen\\Figures\\isogenpep_chidist.png", dpi=300, transparent=True)
plt.show()


plt.scatter(mean_nntimes_lengths.keys(), mean_nntimes_lengths.values(), label="NN")
plt.scatter(mean_fttimes_lengths.keys(), mean_fttimes_lengths.values(), label="FT")
plt.xlabel("Peptide Length")
plt.ylabel("Average Time (seconds)")
plt.legend()
#plt.savefig("Z:\\Group Share\\JGP\\IsoGen\\Figures\\isogenpep_timecomp.png", dpi=300, transparent=True)
plt.show()



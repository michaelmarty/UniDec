from unidec.IsoDec.IsoGen.isogenmass import *
import matplotlib.pyplot as plt
from unidec.modules.isotopetools import *
import numpy as np
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenMassEngine()

data = np.load("Z:\\Group Share\\JGP\\PeptideTraining\\IntactProtein\\training_random_proteins_100_min_1_max_300.npz")

seqs = data["seqs"]
masses = data["masses"]
lengths = [len(seq) for seq in seqs]
dists = ["dists"]


nn_chis = []
nn_times = []
ft_chis = []
ft_times = []

for i, mass in enumerate(masses):
    ai_tstart = time.perf_counter()
    pred_dist = eng.predict(mass)
    ai_tend = time.perf_counter()
    nn_times.append(ai_tend - ai_tstart)

    ft_tstart = time.perf_counter()
    ft_dist = pepmass_to_dist(mass)
    ft_tend = time.perf_counter()
    ft_times.append(ft_tend - ft_tstart)


    #Normalize the two distributions to peak at 1
    pred_dist = pred_dist / np.max(pred_dist)
    ft_dist = ft_dist / np.max(ft_dist)
    true_dist = peptide_to_dist(seqs[i])
    true_dist = true_dist / np.max(true_dist)



    nn_chi = calculate_pearson_chisquared(pred_dist, true_dist)
    nn_chis.append(nn_chi)

    ft_chi = calculate_pearson_chisquared(ft_dist, true_dist)
    ft_chis.append(ft_chi)

nn_chis_lengths = {}
nn_times_lengths = {}
ft_chis_lengths = {}
ft_times_lengths = {}

for i in range(len(nn_chis)):
    if lengths[i] in nn_chis_lengths:
        nn_chis_lengths[lengths[i]].append(nn_chis[i])
        nn_times_lengths[lengths[i]].append(nn_times[i])
        ft_chis_lengths[lengths[i]].append(ft_chis[i])
        ft_times_lengths[lengths[i]].append(ft_times[i])
    else:
        nn_chis_lengths[lengths[i]] = [nn_chis[i]]
        nn_times_lengths[lengths[i]] = [nn_times[i]]
        ft_chis_lengths[lengths[i]] = [ft_chis[i]]
        ft_times_lengths[lengths[i]] = [ft_times[i]]

nn_meanchis_lengths = {}
mean_nntimes_lengths = {}
ft_meanchis_lengths = {}
mean_fttimes_lengths = {}

for length in nn_chis_lengths.keys():
    nn_meanchis_lengths[length] = np.mean(np.array(nn_chis_lengths[length]))
    mean_nntimes_lengths[length] = np.mean(np.array(nn_times_lengths[length]))
    ft_meanchis_lengths[length] = np.mean(np.array(ft_chis_lengths[length]))
    mean_fttimes_lengths[length] = np.mean(np.array(ft_times_lengths[length]))

plt.scatter(nn_meanchis_lengths.keys(), nn_meanchis_lengths.values(), label="NN")
plt.scatter(ft_meanchis_lengths.keys(), ft_meanchis_lengths.values(), label="FT")
plt.xlabel("Peptide Length")
plt.ylabel("Mean Chi")
plt.legend()
plt.show()


plt.scatter(mean_nntimes_lengths.keys(), mean_nntimes_lengths.values(), label="NN")
plt.scatter(mean_fttimes_lengths.keys(), mean_fttimes_lengths.values(), label="FT")
plt.xlabel("Peptide Length")
plt.ylabel("Average Time (seconds)")
plt.legend()
plt.show()

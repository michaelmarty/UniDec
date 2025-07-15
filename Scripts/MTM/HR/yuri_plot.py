import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import unidec.IsoDec.match as match

dir = "C:\\Data\\Yuri\\"
os.chdir(dir)

datafile = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3.csv"

data = np.loadtxt(datafile, delimiter=",", skiprows=1)

mz = data[:, 0]
intensity = data[:, 1]

# df = pd.read_csv("shared_peaks.tsv", sep="\t")
df = pd.read_csv("unmatched_decon.tsv", sep="\t")

mc = match.df_to_matchedcollection(df, ion=None, checkmultiplemonoisos=True)

for pk in mc.peaks:
    if pk.mz > 0:
        plt.axvline(pk.mz, color='red', linestyle='--', label=f'Peak {pk.z} {pk.mz:.4f}')

plt.plot(mz, intensity, label='Data', color='blue')
plt.xlabel('m/z')
plt.show()
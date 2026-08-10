import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

import unidec.IsoDec.compare_tools
import unidec.IsoDec.match as match
from unidec.IsoDec.plots import cplot

dir = "C:\\Data\\Yuri\\"
os.chdir(dir)

datafile = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3.csv"
datafile = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms.csv"

data = np.loadtxt(datafile, delimiter=",", skiprows=1)

cdatafile = os.path.splitext(datafile)[0] + "_1_centroided.txt"
cdata = np.loadtxt(cdatafile)

mz = data[:, 0]
intensity = data[:, 1]

# df = pd.read_csv("shared_peaks.tsv", sep="\t")
df = pd.read_csv("unmatched_decon.tsv", sep="\t")

mc = unidec.IsoDec.compare_tools.df_to_matchedcollection(df, ion=None, checkmultiplemonoisos=True)

# for pk in mc.peaks:
#     if pk.mz > 0:
#         plt.axvline(pk.mz, color='red', linestyle='--', label=f'Peak {pk.z} {pk.mz:.4f}')

plt.plot(mz, intensity, label='Data', color='blue')
cplot(cdata, color="red")
plt.xlabel('m/z')
plt.show()
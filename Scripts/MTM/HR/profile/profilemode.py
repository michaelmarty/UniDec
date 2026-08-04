import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import plot_pks
import time
import matplotlib as mpl

os.chdir("C:\\Data\\Yuri\\")
# file = "C:\\Data\\Yuri\\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3.csv"
# data = np.genfromtxt(file, delimiter=',')

file = "etd_spectrum_profile.txt"
data = np.loadtxt(file)
# data = ud.datachop(data, 1415, 1445)

eng = IsoDecRuntime(verbose=False)
# eng.config.mzwindow[1] = 4.05
print(eng.phase_predictor(data))
stime = time.perf_counter()
pks = eng.batch_process_spectrum(data)
print(pks, "Time:", time.perf_counter() - stime)

# Write data to file
output_file = "processed.txt"
np.savetxt(output_file, data)

isodat = pks.to_merged_isodists()
print("Isodist:", isodat)

# export isodist
isodist_file = "isodist.txt"
np.savetxt(isodist_file, isodat)


plot_pks(pks, data)
plt.show()

# plt.plot(data[:, 0], data[:, 1])
# plt.show()


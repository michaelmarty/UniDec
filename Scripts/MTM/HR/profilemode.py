import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import plot_pks

file = "C:\\Data\\Yuri\\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3.csv"

data = np.genfromtxt(file, delimiter=',')
data = ud.datachop(data, 1415, 1422)

eng = IsoDecEngine(verbose=True)
print(eng.phase_predictor(data))
pks = eng.batch_process_spectrum(data, centroided=False)
print(pks)

plot_pks(pks, data)
plt.show()

# plt.plot(data[:, 0], data[:, 1])
# plt.show()


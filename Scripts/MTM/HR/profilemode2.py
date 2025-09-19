import numpy as np
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.plots import cplot
import unidec.tools as ud
import time

test1 = "C:\\Data\\Volker noisy protein spectra\\Protein+cov-binder_#212#141.txt"

data = np.loadtxt(test1)
data = ud.datachop(data, 800, 2000)  # Chop data to a reasonable range
# data = ud.datachop(data, 890,    900)  # Chop data to a reasonable range
data = ud.datacompsub(data, 10)

eng = IsoDecRuntime()
eng.config.verbose = False
eng.config.matchtol = 50
eng.batch_process_spectrum(data, centroided=False)
print("N Peaks Found:", len(eng.pks.peaks))
st = time.perf_counter()
eng.pks.fill_fwhm(data)
simdata, sumsim, huskleft = eng.pks.gen_peaks(data)
print("Set FWHM for Peaks", time.perf_counter() - st)
#
# fwhmdat = [[p.mz, p.avgsigma] for p in eng.pks.peaks]
# fwhmdat = np.array(fwhmdat)

# plt.plot(data[:, 0], data[:, 1]/np.amax(data[:,1]), label='Processed Data')
# # plt.plot(fwhmdat[:, 0], fwhmdat[:, 1], 'r.', label='FWHM')
# plt.show()
#
#
#
# exit()
dmax = np.amax(data[:, 1])
plt.plot(data[:, 0], data[:, 1], color="k")
cplot(eng.pks.to_merged_isodists())
# for i, d in enumerate(simdata):
#     plt.plot(data[:, 0], d, label=f'Peak {i+1}')
plt.plot(data[:, 0], sumsim, 'r--', label='Simulated Data')
plt.plot(data[:, 0], huskleft, 'b:', label='Leftover Data')
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.show()

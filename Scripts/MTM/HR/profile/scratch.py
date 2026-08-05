import numpy as np
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.plots import cplot
import unidec.tools as ud
import time
from unidec.modules.fitting import *
from unidec.modules.fwhmtools import *

# def calc_sigma(data):
#     ymax = np.amax(data[:, 1])
#     diffs = np.diff(data[:,0])
#     ys = data[:-1,1]
#     ys2 = data[1:,1]
#     ys = (ys + ys2) / 2.0  # Average the y values
#     sum = np.sum(ys * diffs)
#     return sum/(np.sqrt(2 * np.pi) * ymax)

test1 = "C:\\Data\\Volker noisy protein spectra\\Protein+cov-binder_#212#141.txt"

data = np.loadtxt(test1)
data = ud.datachop(data, 800, 2000)  # Chop data to a reasonable range
data = ud.datachop(data, 1265,    1266)  # Chop data to a reasonable range
# data = ud.datacompsub(data, 10)
data[:,1] = data[:,1] - np.amin(data[:,1])  # Subtract baseline
data[:, 1] = data[:, 1] /np.amax(data[:,1]) # Normalize data

fit, fitdat = isolated_peak_fit(data[:,0], data[:,1], 0)

print("Fitted Peak:", fit)
print("Fit Sigma:", fit[0][0] / 2.35482)
sigma = calc_sigma(data)
print("Sigma:", sigma)
m, s, a = calc_gauss(data)
print("Calculated Gaussian:", m, s, a)


f2 = ndis_std(data[:, 0], data[np.argmax(data[:, 1]), 0], sigma, 1)
f3 = ndis_std(data[:, 0], m, s, a)

offset = 0
plt.plot(data[:, 0], data[:, 1], label='Processed Data')
plt.plot(data[:, 0], fitdat-1*offset, 'r-', label='Fitted Peak')
plt.plot(data[:, 0], f2-2*offset, 'g-', label='Fitted Peak')
plt.plot(data[:, 0], f3-3*offset, 'b-', label='Calculated Gaussian')
plt.show()

exit()

plt.plot(data[:, 0], data[:, 1])
cplot(eng.pks.to_merged_isodists())
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.show()

import numpy as np
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.plots import cplot
import unidec.tools as ud
import time
from unidec.modules.fitting import *

# def calc_sigma(data):
#     ymax = np.amax(data[:, 1])
#     diffs = np.diff(data[:,0])
#     ys = data[:-1,1]
#     ys2 = data[1:,1]
#     ys = (ys + ys2) / 2.0  # Average the y values
#     sum = np.sum(ys * diffs)
#     return sum/(np.sqrt(2 * np.pi) * ymax)

from numba import njit

@njit(fastmath=True)
def calc_sigma(data):
    s = 0
    x1, y1 = data[0]
    for d in data[1:]:
        x2, y2 = d
        s += (y1 + y2) * (x2 - x1) / 2.0
        x1, y1 = x2, y2
    return s/(np.sqrt(2 * np.pi) * np.amax(data[:, 1]))

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

f2 = ndis_std(data[:, 0], data[np.argmax(data[:, 1]), 0], sigma, 1)

plt.plot(data[:, 0], data[:, 1], label='Processed Data')
plt.plot(data[:, 0], fitdat-1, 'r-', label='Fitted Peak')
plt.plot(data[:, 0], f2-2, 'g-', label='Fitted Peak')
plt.show()



exit()

plt.plot(data[:, 0], data[:, 1])
cplot(eng.pks.to_merged_isodists())
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.show()

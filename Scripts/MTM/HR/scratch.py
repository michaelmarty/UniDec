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

@njit(fastmath=True)
def calc_ab(data, c):
    b1 = data[:,1] > np.amax(data[:,1]) * 0.25  # Remove points with very low intensity
    data = data[b1]  # Remove zero intensity points
    ys = data[:, 1]**2
    m11 = np.sum(ys)
    m12 = np.sum(data[:,0] * ys)
    m22 = np.sum(data[:,0]**2 * ys)

    M = np.array([[m11, m12], [m12, m22]])

    logdat = np.log(data[:, 1])

    v1 = np.sum(ys * logdat) - c * m22
    v2 = np.sum(data[:,0] * ys * logdat) - c * np.sum(data[:,0] **3 * ys)

    v = np.array([v1, v2])

    ab = np.linalg.solve(M, v)
    return ab

@njit(fastmath=True)
def calc_gauss(data):
    """
    Calculate the Gaussian fit parameters for the data.
    https://arxiv.org/pdf/1907.07241
    :param data: Data array with m/z and intensity
    :return: Gaussian fit parameters (amplitude, mean, sigma)
    """
    # mean = np.sum(data[:, 0] * data[:, 1]) / np.sum(data[:, 1])
    sigma = calc_sigma(data)
    # amplitude = np.amax(data[:, 1])
    c = -1 / (2 * sigma ** 2)
    a, b = calc_ab(data, c)
    mean = -b / (2 * c)
    amplitude = np.exp(a - (b ** 2 / (4*c)))
    return mean, sigma, amplitude

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

from pyopenms import *
import os
import Scripts.JGP.xic_xtractor as xic_xtractor
import matplotlib.pyplot as plt
import time
import matplotlib as mpl
from scipy.optimize import curve_fit
import numpy as np

def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def multi_gaussian(x, *params):
    # Expecting params in the order: A1, mu1, sigma1, A2, mu2, sigma2, ...
    # Reshape params to separate each Gaussian
    n = len(params) // 3
    gaussians = np.zeros_like(x)
    for i in range(n):
        A = params[3 * i]
        mu = params[3 * i + 1]
        sigma = params[3 * i + 2]
        gaussians += A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    return gaussians
# Set backend to Agg
mpl.use('WxAgg')


file = "Z:\\Group Share\\JGP\\CS MS1 Runs\\mzML\\2024_09_04_AS_28_CS_05_ESIneg.mzML"

exp = MSExperiment()
MzMLFile().load(file, exp)
print("Creating indexed file")
ind_file = xic_xtractor.IndexedFile(exp)

ppm_tol = 40

original_mz = 1399.9648
isotopes_to_check = 400
width = original_mz * ppm_tol / 1e6



mzs = np.zeros(2*isotopes_to_check+1)
mzs[isotopes_to_check] = original_mz
for i in range(1, isotopes_to_check+1):
    mzs[isotopes_to_check + i] = original_mz + i*1.00335
    mzs[isotopes_to_check - i] = original_mz - i*1.00335

#make an xic for each of the mzs and plot them
timestart= time.perf_counter()
for mz in mzs:
    xic = ind_file.extract_xic(mz, width, 0, 20)
    plt.plot(xic[:,0]/60, xic[:,2], label=str(mz))
print("Time taken: ", time.perf_counter() - timestart)
plt.legend()
plt.show()
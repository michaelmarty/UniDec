import numpy as np
import time
# from scipy import fftpack
from numpy import fft as fftpack
from numba import njit
from copy import deepcopy
from unidec.IsoDec.IsoGen.isogenwrapper import IsoGenWrapper

# import numba as nb
# import pyteomics.mass as ms

mass_diff_c = 1.0033
massavgine = 111.1254
avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])
isotopes = np.array([[[12., 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
                     [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
                     [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
                     [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
                     [[31.9720707, 95.02], [32.9714585, 0.75], [33.9678668, 4.21], [35.9760809, 0.02]]])

isoparams = [1.00840852e+00, 1.25318718e-03, 2.37226341e+00, 8.19178000e-04, -4.37741951e-01, 6.64992972e-04,
             9.94230511e-01, 4.64975237e-01, 1.00529041e-02, 5.81240792e-01]


@njit(fastmath=True)
def isotopemid(mass: float):
    a = -4.37741951e-01  # isoparams[4]
    b = 6.64992972e-04  # isoparams[5]
    c = 9.94230511e-01  # isoparams[6]
    return a + b * pow(mass, c)


@njit(fastmath=True)
def isotopesig(mass: float):
    a = 4.64975237e-01  # isoparams[7]
    b = 1.00529041e-02  # isoparams[8]
    c = 5.81240792e-01  # isoparams[9]
    return a + b * pow(mass, c)


@njit(fastmath=True)
def isotopealpha(mass: float):
    a = 1.00840852e+00  # isoparams[0]
    b = 1.25318718e-03  # isoparams[1]
    return a * np.exp(-mass * b)


@njit(fastmath=True)
def isotopebeta(mass: float):
    a = 2.37226341e+00  # isoparams[2]
    b = 8.19178000e-04  # isoparams[3]
    return a * np.exp(-mass * b)


@njit(fastmath=True)
def isomike(mass: float, length=128) -> np.ndarray:
    mid = isotopemid(mass)
    sig = isotopesig(mass)
    alpha = isotopealpha(mass)
    amp = (1.0 - alpha) / (sig * 2.50662827)
    beta = isotopebeta(mass)
    maxval = 0
    isoindex = np.arange(0, length)
    isotopeval = np.zeros(length)
    for k in range(length):
        e = alpha * np.exp(-isoindex[k] * beta)
        g = amp * np.exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)))
        temp = e + g
        if temp > maxval:
            maxval = temp
        isotopeval[k] = temp
    return isotopeval / maxval


@njit(fastmath=True)
def makemassmike(testmass: float) -> float:
    num = testmass / massavgine * avgine
    intnum = np.array([int(round(n)) for n in num])
    x = intnum * isotopes[:, 0, 0]
    minmassint = np.sum(x)

    return minmassint


#@njit(fastmath=True)
def fast_calc_averagine_isotope_dist(mass, charge=1, adductmass=1.007276467):
    # Predict Isotopic Intensities
    intensities = isogenmass(mass)
    # Calculate masses for these
    masses = np.arange(0, len(intensities)) * mass_diff_c + mass

    # Load Into Array
    dist = np.zeros((len(masses), 2))
    dist[:, 0] = masses
    dist[:, 1] = intensities

    # Filter Low Intensities
    b1 = intensities > np.amax(intensities) * 0.01
    dist = dist[b1]

    # Convert to m/z
    dist[:, 0] = (dist[:, 0] + float(charge) * adductmass) / float(charge)

    return dist

#@njit(fastmath=True)
def isogenmass(mass):
    eng = IsoGenWrapper()
    return eng.gen_isodist(mass)


#@njit(fastmath=True)
def fast_calc_averagine_isotope_dist_dualoutput(mass, charge=1, adductmass=1.007276467, isotopethresh: float = 0.01):
    # Predict Isotopic Intensities
    intensities = isogenmass(mass)
    # Calculate masses for these
    masses = np.arange(0, len(intensities)) * mass_diff_c + mass

    # Load Into Array
    dist = np.zeros((len(masses), 2))
    dist[:, 0] = masses
    dist[:, 1] = intensities

    # Filter Low Intensities
    b1 = intensities > np.amax(intensities) * isotopethresh

    dist = dist[b1]

    # Convert to m/z
    massdist = dist.copy()
    dist[:, 0] = (dist[:, 0] + float(charge) * adductmass) / float(charge)

    return dist, massdist


# @njit(fastmath=True)
def makemass(testmass):
    num = testmass / massavgine * avgine
    intnum = np.array([int(round(n)) for n in num])
    x = intnum * isotopes[:, 0, 0]
    minmassint = np.sum(x)
    formula = ""
    if intnum[0] != 0:
        formula = formula + "C" + str(intnum[0])
    if intnum[1] != 0:
        formula = formula + "H" + str(intnum[1])
    if intnum[2] != 0:
        formula = formula + "N" + str(intnum[2])
    if intnum[3] != 0:
        formula = formula + "O" + str(intnum[3])
    if intnum[4] != 0:
        formula = formula + "S" + str(intnum[4])
    return formula, minmassint, intnum


isolength = 128
buffer = np.zeros(isolength)
h = np.array([1, 0.00015, 0, 0])
c = np.array([1, 0.011, 0, 0])
n = np.array([1, 0.0037, 0, 0])
o = np.array([1, 0.0004, 0.002, 0])
s = np.array([1, 0.0079, 0.044, 0])
h = np.append(h, buffer)
c = np.append(c, buffer)
n = np.append(n, buffer)
o = np.append(o, buffer)
s = np.append(s, buffer)

dt = np.dtype(np.complex128)
hft = fftpack.rfft(h).astype(dt)
cft = fftpack.rfft(c).astype(dt)
nft = fftpack.rfft(n).astype(dt)
oft = fftpack.rfft(o).astype(dt)
sft = fftpack.rfft(s).astype(dt)


# @njit(fastmath=True)
def isojim(isolist, length=isolength):
    """Thanks to Jim Prell for Sketching this Code"""
    numc = isolist[0]
    numh = isolist[1]
    numn = isolist[2]
    numo = isolist[3]
    nums = isolist[4]

    allft = cft ** numc * hft ** numh * nft ** numn * oft ** numo * sft ** nums

    # with nb.objmode(allift='float64[:]'):
    #    allift = fftpack.irfft(allft)
    allift = np.abs(fftpack.irfft(allft))
    # allift = np.abs(allift)
    allift = allift / np.amax(allift)
    return allift[:isolength]  # .astype(nb.float64)


# @njit(fastmath=True)
def predict_charge(mass):
    """
    Give predicted native charge state of species of defined mass
    https://www.pnas.org/content/107/5/2007.long
    :param mass: Mass in Da
    :return: Float, predicted average native charge state
    """
    m = 0.0467
    s = 0.533
    nativez = m * (mass ** s)
    return nativez


# @njit(fastmath=True)
def calc_averagine_isotope_dist(mass, mono=False, charge=None, adductmass=1.007276467, crop=False, fast=True,
                                length=isolength, **kwargs):
    if fast:
        minmassint = makemassmike(mass)
        intensities = isomike(mass)
    else:
        _, minmassint, isolist = makemass(mass)
        # print(isolist)
        intensities = isojim(isolist)

    if mono:
        minmassint = mass
    masses = np.arange(0, len(intensities)) + minmassint
    dist = np.zeros((len(masses), 2))
    dist[:, 0] = masses
    dist[:, 1] = intensities
    if not mono:
        dist = correct_avg(dist, mass)

    if charge == "Auto":
        z = predict_charge(mass)
    elif charge is not None:
        z = float(charge)
    else:
        z = 1.0

    if adductmass is None:
        adductmass = 0.0

    if z is not None and z != 0:
        dist[:, 0] = (dist[:, 0] + z * adductmass) / z

    if crop:
        b1 = dist[:, 1] > np.amax(dist[:, 1]) * 0.01
        dist = dist[b1]

    return dist


@njit(fastmath=True)
def correct_avg(dist, mass):
    """Note, switched from weighted average to peak correction"""
    avg = dist[np.argmax(dist[:, 1]), 0]  # np.average(dist[:, 0], weights=dist[:, 1])
    dist[:, 0] = dist[:, 0] - avg + mass
    return dist


def get_apex_mono_diff(mass):
    dist = calc_averagine_isotope_dist(mass)
    apex = dist[np.argmax(dist[:, 1]), 0]
    mono = dist[0, 0]
    diff = apex - mono
    return diff


def predict_apex_mono_diff(mass):
    m = 0.00063139
    b = -0.53143767
    fit = m * mass + b
    if fit < 0:
        fit = 0
    return round(fit)


if __name__ == "__main__":

    m = 1000
    print(isogenmass(m))
    exit()

    m = 1000
    formula, minmassint, isolist = makemass(m)
    x = isojim(isolist)[:10]
    y = isomike(m)[:10]

    print(x, y)

    exit()
    x = 10 ** np.arange(2, 6, 0.1)
    y = [get_apex_mono_diff(m) for m in x]
    fit = np.polyfit(x, y, 1)
    print(fit)
    fitdat = [predict_apex_mono_diff(m) for m in x]
    import matplotlib.pyplot as plt

    plt.plot(x, y)
    plt.plot(x, fitdat)
    plt.show()

    exit()
    mval = 1000000
    startime = time.perf_counter()
    dist = calc_averagine_isotope_dist(mval)
    endtime = time.perf_counter()
    print(endtime - startime)
    plt.plot(dist[:, 0], dist[:, 1])
    plt.show()

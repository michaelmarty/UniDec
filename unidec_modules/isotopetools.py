import numpy as np
import time
from scipy import fftpack
import matplotlib.pyplot as plt

# import pyteomics.mass as ms

massavgine = 111.1254
avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])
isotopes = np.array([[[12., 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
                     [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
                     [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
                     [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
                     [[31.9720707, 95.02], [32.9714585, 0.75], [33.9678668, 4.21], [35.9760809, 0.02]]])


def makemass(testmass):
    num = testmass / massavgine * avgine
    intnum = [int(round(n)) for n in num]
    minmassint = (intnum * isotopes[:, 0, 0]).sum()
    formula = ""
    if intnum[0] is not 0:
        formula = formula + "C" + str(intnum[0])
    if intnum[1] is not 0:
        formula = formula + "H" + str(intnum[1])
    if intnum[2] is not 0:
        formula = formula + "N" + str(intnum[2])
    if intnum[3] is not 0:
        formula = formula + "O" + str(intnum[3])
    if intnum[4] is not 0:
        formula = formula + "S" + str(intnum[4])
    return formula, minmassint, intnum


def isojim(isolist):
    '''Thanks to Jim Prell for Sketching this Code'''
    numc = isolist[0]
    numh = isolist[1]
    numn = isolist[2]
    numo = isolist[3]
    nums = isolist[4]

    buffer = np.zeros(700)
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
    hft = fftpack.fft(h).astype(dt)
    cft = fftpack.fft(c).astype(dt)
    nft = fftpack.fft(n).astype(dt)
    oft = fftpack.fft(o).astype(dt)
    sft = fftpack.fft(s).astype(dt)

    allft = cft ** numc * hft ** numh * nft ** numn * oft ** numo * sft ** nums
    allift = np.abs(fftpack.ifft(allft))
    allift = allift / np.amax(allift)
    return allift


def calc_averagine_isotope_dist(mass, mono=False):
    formula, minmassint, isolist = makemass(mass)
    # print(isolist)
    intensities = isojim(isolist)
    masses = np.arange(0, len(intensities)) + minmassint
    dist = np.transpose([masses, intensities])
    if not mono:
        dist = correct_avg(dist, mass)
    # print(ms.calculate_mass(formula=formula, average=True))
    # print(np.average(dist[:, 0], weights=dist[:, 1]))
    # print(minmassint)
    return np.array(dist)


def correct_avg(dist, mass):
    '''Note, switched from weighted average to peak correction'''
    avg = dist[np.argmax(dist[:, 1]), 0]  # np.average(dist[:, 0], weights=dist[:, 1])
    dist[:, 0] = dist[:, 0] - avg + mass
    return dist

if __name__ == "__main__":
    mval = 1000000
    startime = time.perf_counter()
    dist = calc_averagine_isotope_dist(mval)
    endtime = time.perf_counter()
    print(endtime - startime)
    plt.plot(dist[:, 0], dist[:, 1])
    plt.show()

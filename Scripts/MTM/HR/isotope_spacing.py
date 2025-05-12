import numpy as np
from numpy import fft as fftpack
from unidec.modules.isotopetools import makemass
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.use("wxAgg")

def merge_isotope_spacing(x1, x2, y1, y2, isolength):
    diffs = np.add.outer(x1, x2)
    probs = np.outer(y1, y2)

    x = []
    y = []
    for i in range(isolength):
        b1 = diffs < i + 0.5
        b2 = diffs > i - 0.5
        b1 = np.logical_and(b1, b2)
        intensity = 0
        weighted_average = i
        if np.any(b1):
            intensity = np.sum(probs[b1])
        if intensity!=0:
            weighted_average = np.sum(diffs[b1] * probs[b1]) / intensity
        x.append(weighted_average)
        y.append(intensity)

    return np.array(x), np.array(y)

def threshold(array, thresh):
    b1 = array<np.amax(array) * thresh
    array[b1] = 0
    return array

def merge_all(deltas, dists, isolength, flen=7):
    x1 = deltas[0]
    y1 = dists[0]

    for i in range(flen - 1):
        x2 = deltas[i + 1]
        y2 = dists[i + 1]
        # y2 = threshold(y2, 0.00001)
        # y1 = threshold(y1, 0.00001)
        x1, y1 = merge_isotope_spacing(x1, x2, y1, y2, isolength)
    return x1, y1

def calc_mzs(mass, diffs, pre_ft, isolength):
    formula, minmassint, intnum = makemass(mass)
    intnum = np.concatenate((intnum, [intnum[3], intnum[4]]))
    print(intnum)

    dists = []
    deltas = []
    for i in range(len(intnum)):
        allft = pre_ft[i] ** intnum[i]
        allift = np.abs(fftpack.irfft(allft))
        allift = allift / np.amax(allift)
        dists.append(allift[:isolength])
        diffs_temp = diffs[i] * np.arange(isolength)
        deltas.append(diffs_temp)

    x1, y1 = merge_all(deltas, dists, isolength)

    x2 = np.diff(x1)
    y2 = y1[:-1]

    weighted_avg = np.sum(x2 * y2) / np.sum(y2)
    return weighted_avg


avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])
isotopes = np.array([[[12., 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
                     [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
                     [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
                     [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
                     [[31.9720707, 95.02], [32.9714585, 0.75], [33.9678668, 4.21], [35.9760809, 0.02]]])

isolength = 128
buffer = np.zeros(isolength)
h = np.array([1, 0.00015, 0, 0])
c = np.array([1, 0.011, 0, 0])
n = np.array([1, 0.0037, 0, 0])
o = np.array([1, 0, 0.002, 0])
o17 = np.array([1, 0.0004, 0, 0])
s = np.array([1, 0, 0.044, 0])
s33 = np.array([1, 0.0079, 0, 0])


h = np.append(h, buffer)
c = np.append(c, buffer)
n = np.append(n, buffer)
o = np.append(o, buffer)
s = np.append(s, buffer)
o17 = np.append(o17, buffer)
s33 = np.append(s33, buffer)


dt = np.dtype(np.complex128)
hft = fftpack.rfft(h).astype(dt)
cft = fftpack.rfft(c).astype(dt)
nft = fftpack.rfft(n).astype(dt)
oft = fftpack.rfft(o).astype(dt)
sft = fftpack.rfft(s).astype(dt)
o17ft = fftpack.rfft(o17).astype(dt)
s33ft = fftpack.rfft(s33).astype(dt)

pre_ft = np.array([cft, hft, nft, oft, sft, o17ft, s33ft])

dh = 2.0141018-1.0078250
dc = 13.0033548-12.0000000
dn = 15.0001089-14.0030740
do = 17.9991605-15.9949146
ds = 33.9678668-31.9720707
ds33 = 32.9714585-31.9720707
do17 = 16.9991315-15.9949146

diffs = [dc, dh, dn, do, ds, do17, ds33]

a = calc_mzs(100, diffs, pre_ft, isolength)

x = np.arange(2, 5, 0.1)
masses = 10**x
ints = [calc_mzs(m, diffs, pre_ft, isolength) for m in masses]

plt.plot(masses, ints)
plt.show()

# plt.plot(np.diff(x1))
# plt.plot(y1)
# plt.show()




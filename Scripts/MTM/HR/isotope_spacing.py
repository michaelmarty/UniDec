import numpy as np
from numpy import fft as fftpack
from unidec.modules.isotopetools import makemass

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

formula, minmassint, intnum = makemass(1000)
print(formula, minmassint, intnum)
intnum = np.concatenate((intnum, [intnum[3], intnum[4]]))
print(intnum)

dists = []
for i in range(len(intnum)):
    allft = pre_ft[i] ** intnum[i]
    allift = np.abs(fftpack.irfft(allft))
    allift = allift / np.amax(allift)
    dists.append(allift[:isolength])




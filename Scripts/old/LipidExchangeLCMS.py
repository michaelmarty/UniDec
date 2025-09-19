__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
import pymzml
import os

ps = [760., 765.]
dps = [788., 794.]
pc = [804., 809.]
dpc = [832., 839.]

file = "MTM150515_PSPC_50.mzML"
dir = "C:\\Users\\Michael Marty\\Box Sync\\Teaching\\525 Fall 2016\\Data"

ranges = [ps, dps, pc, dpc]

mids = [(r[0] + r[1]) / 2. for r in ranges]
spans = [(r[1] - r[0]) / 2. for r in ranges]

path = os.path.join(dir, file)
msrun = pymzml.run.Reader(path)
times = []
cints = []
plt.figure()
for spectrum in msrun:
    try:
        scantime = float(spectrum["scan start time"])
    except:
        scantime = -1
        
    if scantime > 0 and scantime < 30:
        times.append(scantime)
        # impdat=np.transpose([spectrum.mz,spectrum.i])
        mz = np.array(spectrum.mz)
        intensity = np.array(spectrum.i)

        ints = []
        for i in range(0, len(ranges)):
            mid = mids[i]
            span = spans[i]
            mztest = np.abs(mz - mid)
            btest = mztest < span
            inttot = np.sum(intensity[btest])
            ints.append(inttot)
            plt.plot(mz[btest], intensity[btest])
        cints.append(ints)

cints = np.transpose(cints)
times = np.array(times)
plt.show()
maxes = np.amax(cints, axis=1)
psr = maxes[0] / maxes[1]
pcr = maxes[2] / maxes[3]
# print psr,pcr
tot = psr + pcr
psp = psr / tot * 100
pcp = pcr / tot * 100
print("POPS: %.2f  POPC: %.2f " % (psp, pcp))

maxes = np.sum(cints, axis=1)
psr = maxes[0] / maxes[1]
pcr = maxes[2] / maxes[3]
# print psr,pcr
tot = psr + pcr
psp = psr / tot * 100
pcp = pcr / tot * 100
print("POPS: %.2f  POPC: %.2f " % (psp, pcp))

plt.figure()
for i in range(0, len(ranges)):
    plt.plot(times, cints[i])
plt.show()

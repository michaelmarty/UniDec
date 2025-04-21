import numpy as np

import matplotlib.pyplot as plt
import os
import matplotlib
'''
minmzdiff = []
for i, m in enumerate(mz):
    diff = mz - m
    diff = diff[diff > 0]
    if len(diff)>0:
        minmzdiff.append(np.amin(diff))
        if i % 10000 == 0:
            print(i)

np.savetxt("negdiffs.txt", minmzdiff)
print("Mean Minimum Diff:", np.mean(minmzdiff))
print("Min Min Diff:", np.amin(minmzdiff))
'''
matplotlib.use('WxAgg')

libcsvpath = "Z:\Group Share\MTM\\neglib.csv"
os.chdir(os.path.dirname(libcsvpath))

posdat = np.loadtxt("posdiffs.txt")
negdat = np.loadtxt("negdiffs.txt")

posdat = posdat[posdat<0.2]
negdat = negdat[negdat<0.2]

plt.subplot(211)
plt.title("Positive Mode")
plt.hist(posdat, bins=100)
plt.xlabel('Minimum Difference in Mass (Da)')
plt.subplot(212)
plt.title("Negative Mode")
plt.hist(negdat, bins=100)
plt.xlabel('Minimum Difference in Mass (Da)')
plt.tight_layout()
plt.show()

print(posdat)

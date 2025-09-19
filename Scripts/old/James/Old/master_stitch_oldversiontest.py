import numpy as np
import matplotlib.pyplot as plt
import os
import UniDec.unidec_modules.unidectools as ud

dir = "Z:\\Group Share\\James Keener\\LXMS Paper\\Stitch\\60_320V_method5"
os.chdir(dir)

files = ud.match_files(dir, "*_lipav.txt")
print(files)

colors = ['purple', 'navy', 'forestgreen']

nrange = [2, 175]

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
for i, f in enumerate(files):
    if "mix" in f:
        color = colors[int(i / 2)]
        l = np.loadtxt(f)

        b1 = l[0] <= nrange[1]
        b2 = l[1] >= nrange[0]
        b3 = np.logical_and(b1, b2)

        ax2.scatter(l[0, b3], l[1, b3], c=color)
        ax2.errorbar(l[0, b3], l[1, b3], yerr=l[2, b3], c=color, capsize=3, capthick=1, elinewidth=1)
        ax2.set_xlim(nrange)

        f2 = f[:-10] + "_datavg.txt"
        print(f2)

        d = np.loadtxt(f2)

        b1 = d[:, 0] <= np.amax(l[0])
        b2 = d[:, 0] >= np.amin(l[0])
        b3 = np.logical_and(b1, b2)
        d[:, 1] -= np.amin(d[b3, 1])
        d[:, 1] /= np.amax(d[b3, 1])

        ax1.plot(d[b3, 0], d[b3, 1], color=color)

        ax1.set_xlim(nrange)

popc = 760.076
popg = 749.02

ax2.hlines((0.5 * popc + 0.5 * popg), nrange[0], nrange[1])
ax2.hlines(popg, nrange[0], nrange[1])
ax2.hlines(popc, nrange[0], nrange[1])
ax2.hlines((0.25 * popc + 0.75 * popg), nrange[0], nrange[1])
ax2.hlines(753.6)

plt.show()
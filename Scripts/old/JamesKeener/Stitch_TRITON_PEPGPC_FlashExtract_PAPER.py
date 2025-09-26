import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.tools as ud

dir = "C:\\Users\\Lab.DESKTOP-HJRL65I\\University of Arizona\\martylab - Documents\\Papers\\Flash Extraction MP Nanodisc\\Data\\TRITON_POPE_POPC_POPG"
os.chdir(dir)

files = ud.match_files(dir, "*_lipav.txt")
print(files)

colors = ['k', 'navy', 'r']

nrange = [-1, 15]

plt.figure(figsize=(9.5,11.5))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
#ax3 = plt.subplot(313)
for i, f in enumerate(files):
    if "mix" in f:
        color = colors[int(i / 2)]
        l = np.loadtxt(f)

        b1 = l[0] <= nrange[1]
        b2 = l[1] >= nrange[0]
        b3 = np.logical_and(b1, b2)

        ax2.scatter(l[0, b3], l[1, b3], c="yellow", marker="o", edgecolor=color, zorder=200, s=100)
        ax2.errorbar(l[0, b3], l[1, b3], yerr=l[2, b3], c=color, capsize=3, capthick=1, elinewidth=1)
        ax2.set_xlim(nrange)
        ax2.set_ylim(717.996, 760.076)

        #ax3.plot(l[0,b3][1:], filt.gaussian_filter1d(np.diff(l[1,b3]),10))
       # ax3.hlines(0, nrange[0], nrange[1], linestyle="--")

        f2 = f[:-10] + "_datavg.txt"
        print(f2)

        d = np.loadtxt(f2)

        b1 = d[:, 0] <= np.amax(l[0])+0.5
        b2 = d[:, 0] >= np.amin(l[0])-1.5
        b3 = np.logical_and(b1, b2)
        d[:, 1] -= np.amin(d[b3, 1])
        d[:, 1] /= np.amax(d[b3, 1])

        ax1.plot(d[b3, 0], d[b3, 1]*100, color=color)

        ax1.set_xlim(nrange)
        ax1.set_ylim(0,100)

popc = 760.076
popg = 748.999
ax1.get_xaxis().set_visible(False)
ax2.hlines((742.357), nrange[0], nrange[1], linestyle="--", color="black")
ax2.hlines((748.397), nrange[0], nrange[1], linestyle="--", color="grey")
ax2.hlines((744.018), nrange[0], nrange[1], linestyle="--", color="green")
#ax2.hlines(popg, nrange[0], nrange[1])
#ax2.hlines(popc, nrange[0], nrange[1])
#ax2.hlines((0.25 * popc + 0.75 * popg), nrange[0], nrange[1])
plt.tight_layout()
plt.savefig("finalstitch_method5.pdf")
plt.show()
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

rcParams['font.size'] = 16

nfile = "z:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\example spectrum\\EX_n.txt"
pfile = "z:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\example spectrum\\EX_P.txt"

ndat = np.loadtxt(nfile)
pdat = np.loadtxt(pfile)

pdat[:, 1] /= np.amax(pdat[:int(len(pdat) / 2), 1]) / 100
ndat[:, 1] /= np.amax(ndat[:, 1]) / 100

nboo1 = ndat[:, 0] < 670
nboo2 = ndat[:, 0] > 666
nboo3 = np.logical_and(nboo1, nboo2)

ndat[nboo3, 1] /= np.amax(ndat[:int(len(ndat) / 2), 1]) / 100

pmax2 = np.amax(pdat[int(len(pdat) / 2):, 1])
nmax2 = np.amax(ndat[int(len(ndat) / 2):, 1])
ratio = 0.2

ndat[int(len(ndat) / 2):, 1] *= ratio * pmax2 / nmax2

plt.figure(figsize=(3, 4))
plt.ylim(0, 100)
plt.plot(ndat[:, 0], ndat[:, 1], color="g", label="PG (-)")
plt.plot(pdat[:, 0], pdat[:, 1], color="k", label="PC (+)")
plt.xlabel("m/z")
# plt.ylabel("%")


ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.get_yaxis().set_ticks([0, 50, 100])
ax.get_yaxis().set_ticklabels(["0", '%', "100"])

# plt.legend(loc=9)
plt.tight_layout()

if True:
    dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
    path = os.path.join(dir, "Example Spectrum_2_10.png")
    plt.savefig(path, transparent=True)

plt.show()

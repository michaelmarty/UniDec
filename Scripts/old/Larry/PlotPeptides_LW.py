import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
# rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\MassSpecData\\Melittin\\DMPG\\MassDefectExtracts\\"
os.chdir(path)

files = "MelittinDMPG_Average.txt"

ydata = np.loadtxt(files)
yerr = np.loadtxt("MelittinDMPG_errors.txt")

for i in range(0, len(ydata[0])):
    factor = np.sum(ydata[:, i]) / 100
    ydata[:, i] /= factor
    yerr[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("MelittinDMPG_Xvalues.txt")

labels = ["0", "1", "2", "3", "4", "5", "6", "7"]
colors = ["k", "b", "g", "r", "c", "m", "y", "#808080"]
markers = ["o", ">", "^", "s", "x", "+", "*", "D"]

cmap = matplotlib.cm.get_cmap('inferno')
colors = cmap(np.arange(0, len(labels)) / (len(labels)-1))


ticks = np.arange(0, 21, 5)

for i in range(0, len(ydata)):
    y = ydata[i]
    ye = yerr[i]
    plt.errorbar(xdata, y, yerr=ye, marker=markers[i], label=labels[i], color=colors[i])
# plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
# plt.errorbar(x, y-0.4, yerr=yerr, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(0, 20)
plt.xticks(ticks)

plt.xlabel("Peptide:Nanodisc")
plt.ylabel("Relative Intensity (%)")
plt.legend()
plt.show()

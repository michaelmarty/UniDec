import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\Figures\\Paper\\GramVoltageRamps\\"
os.chdir(path)

files = "GA_DMPC_RAMP_Averages.txt"

ydata = np.loadtxt(files)
yerr = np.loadtxt("GA_DMPC_RAMP_Errors.txt")

for i in range(0, len(ydata[0])):
    factor = np.sum(ydata[:, i]) / 100
    ydata[:, i] /= factor
    yerr[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("GramDMPC_Xvalues.txt")

bool=xdata<150

labels = ["0", "2", "4", "6", "7"]
#colors = ["k", "#808080", "b", "r", "c", "m", "y", "g"]
markers = ["o", ">", "^", "s", "x", "+", "*", "D"]

#cmap = matplotlib.cm.get_cmap('tab10')
#colors = cmap(np.arange(0, len(labels)) / (len(labels)-1))

colors = [[0.267004, 0.004874, 0.329415, 1.],
          [0.275191, 0.194905, 0.496005, 1.],
          [0.212395, 0.359683, 0.55171, 1.],
          [0.153364, 0.497, 0.557724, 1.],
          [0.122312, 0.633153, 0.530398, 1.],
          [0.288921, 0.758394, 0.428426, 1.],
          [0.626579, 0.854645, 0.223353, 1.],
          [0.993248, 0.906157, 0.143936, 1.]]
ticks = np.arange(0, 101, 20)

for i in range(0, len(ydata)):
    y = ydata[i]
    ye = yerr[i]
    plt.errorbar(xdata[bool], y[bool], yerr=ye[bool], marker=markers[i], label=labels[i], color=colors[i])
#plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
#plt.errorbar(x, y-0.4, yerr=yerr, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(0, 101)
ax.set_ylim(0,100)
plt.xticks(ticks)

ticks = np.arange(0, 101, 50)
#plt.yticks(ticks)
plt.yticks(ticks, labels=["0","%","100"])

plt.xlabel("Voltage (V)")
#plt.ylabel("%")
#plt.legend()
plt.show()

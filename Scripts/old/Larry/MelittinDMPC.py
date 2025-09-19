import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import UniDec.unidec_modules.unidectools as ud

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\Figures\\Paper\\MelittinTitrations\\"
os.chdir(path)

files = "MelittinDMPC_Averages2.txt"

ydata = np.loadtxt(files)
yerr = np.loadtxt("MelittinDMPC_Errors2.txt")

for i in range(0, len(ydata[0])):
    factor = np.sum(ydata[:, i]) / 100
    ydata[:, i] /= factor
    yerr[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("MelittinDMPC_Xvalues.txt")

bool=xdata<20

labels = ["0", "1", "2", "3", "4", "5", "6", "7"]
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
ticks = np.arange(0, 18.5, 5)

for i in range(0, len(ydata)):
    y = ydata[i]
    ye = yerr[i]
    plt.errorbar(xdata[bool], y[bool], yerr=ye[bool], marker=markers[i], label=labels[i], color=colors[i])
# plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
# plt.errorbar(x, y-0.4, yerr=yerr, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(0, 18.5)
ax.set_ylim(0,100)
plt.xticks(ticks)

ticks = np.arange(0, 101, 50)
plt.yticks(ticks)
plt.yticks(ticks, labels=["0","%","100"])
plt.xlabel("Melittin/Nanodisc")
#plt.ylabel("Relative Intensity (%)")
#plt.legend()
plt.show()

#plt.figure()
#ydata2=np.array([y/np.amax(y)*100. for y in ydata.transpose()]).transpose()
#plt.contourf(xdata, sdat, ydata)
#norm=cm.colors.Normalize(vmax=100, vmin=0)
#plt.imshow(ydata, aspect='auto', norm=norm, interpolation='nearest', origin="lower", cmap="ocean_r") #viridis, ocean_r, plasma, gnuplot, jet
#cbar=plt.colorbar(use_gridspec=True, ticks=[0, 50, 100])
#cbar.ax.get_yaxis().set_tick_params(direction='out')
#cbar.ax.set_yticklabels(["0", "%", "100"])
#plt.ylabel("# of Peptides / Nanodisc")
#plt.xlabel("Peptide Ratio")
#xticklabels = np.arange(0, 21, 5)
#plt.xticks(np.arange(0,len(xdata)),xdata)
#plt.show()


x = xdata[bool]
for i in range(0, len(x)):
    y = ydata[:,i]
    n = np.array(labels[:len(y)]).astype(float)
    print(n, i , len(x))
    print(y)
    ye = yerr[:, i]
    fits, fitdat = ud.binomial_fit(n,y)
    print(fits, fitdat)
    err=np.abs(fitdat-y)
    color=colors[i]
    plt.errorbar(n,y, yerr=ye, color=color)
    plt.plot(n,fitdat, color=color, linestyle="--")

    plt.show()
exit()
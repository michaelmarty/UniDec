import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib
import UniDec.unidec_modules.unidectools as ud

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\Figures\\Paper\\Data for figures\\Gramicidin Titrations\\DMPC\\"
os.chdir(path)

files = "DMPC_Gram_Averages.txt"

ydata = np.loadtxt(files)
yerr = np.loadtxt("DMPC_Gram_errors.txt")

# for i in range(0, len(ydata[0])):
#   factor = np.sum(ydata[:, i]) / 100
#  ydata[:, i] /= factor
# yerr[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("DMPC_Gram_Xvalues.txt")
labels = ["0", "2", "4", "6", "8", "10", "12", "14"]  # ,"8"]
colors = ["k", "#808080", "b", "r", "c", "m", "y", "g"]
markers = ["o", ">", "^", "s", "x", "+", "*", "D"]

cmap = matplotlib.cm.get_cmap('viridis')
colors = cmap(np.arange(0, len(labels)) / (len(labels) - 1))
colors = [[0.267004, 0.004874, 0.329415, 1.],
          [0.275191, 0.194905, 0.496005, 1.],
          [0.212395, 0.359683, 0.55171, 1.],
          [0.153364, 0.497, 0.557724, 1.],
          [0.122312, 0.633153, 0.530398, 1.],
          [0.288921, 0.758394, 0.428426, 1.],
          [0.626579, 0.854645, 0.223353, 1.],
          [0.993248, 0.906157, 0.143936, 1.]]
'''colors=[[1.46200e-03, 4.66000e-04, 1.38660e-02, 1.00000e+00],
 [1.55850e-01, 4.45590e-02, 3.25338e-01, 1.00000e+00],
 [3.97674e-01, 8.32570e-02, 4.33183e-01, 1.00000e+00],
 [6.21685e-01, 1.64184e-01, 3.88781e-01, 1.00000e+00],
 [8.32299e-01, 2.83913e-01, 2.57383e-01, 1.00000e+00],
 [9.61293e-01, 4.88716e-01, 8.42890e-02, 1.00000e+00],
 [9.81173e-01, 7.59135e-01, 1.56863e-01, 1.00000e+00],
 [.95, .95, 6.44924e-01,1.00000e+00]]'''

ticks = np.arange(0, 12.5, 5)
boo1 = xdata < 25

x = xdata[boo1]
n = np.array(labels[:len(ydata[0])]).astype(float)
print(n)
for i in range(0, len(x)):
    y = ydata[:,i]
    ye = yerr[:, i]
    print(y)
    fits, fitdat = ud.binomial_fit(n,y)

    err=np.abs(fitdat-y)
    color=colors[i]
    plt.errorbar(n,y, yerr=ye, color=color)
    plt.plot(n,fitdat, color=color, linestyle="--")

    plt.show()
exit()


for i in range(0, len(ydata)):
    y = ydata[i]
    ye = yerr[i]
    plt.errorbar(xdata[boo1], y[boo1], yerr=ye[boo1], marker=markers[i], label=labels[i], color=colors[i])
# plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
# plt.errorbar(x, y-0.4, yerr=yerr, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(0, 12.5)
ax.set_ylim(0, 1)
plt.xticks(ticks)

ticks = np.arange(0, 1.01, 0.5)
plt.yticks(ticks)

plt.xlabel("Peptide Ratio")
plt.ylabel("Relative Intensity (%)")
plt.legend()
plt.show()

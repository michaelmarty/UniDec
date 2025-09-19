import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib
from scipy import stats
from scipy.optimize import curve_fit

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\Figures\\Paper\\LL37Titrations\\"
os.chdir(path)

files = "LL37_DMPG_Titration.txt"

ydata = np.loadtxt(files)
yerr = np.loadtxt("LL37_DMPG_Titration_errors.txt")

ydata = np.vstack((ydata, ydata[0]*0, ydata[0]*0))
yerr = np.vstack((yerr, ydata[0]*0, ydata[0]*0))
print(ydata)

for i in range(0, len(ydata[0])):
    factor = np.sum(ydata[:, i]) / 100
    ydata[:, i] /= factor
    yerr[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("LL37_DMPG_XVALS.txt")
labels = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]  # ,"8"]
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
          [0.993248, 0.906157, 0.143936, 1.],
          [0.993248, 0.906157, 0.143936, 1.],
          [0.993248, 0.906157, 0.143936, 1.]]
'''colors=[[1.46200e-03, 4.66000e-04, 1.38660e-02, 1.00000e+00],
 [1.55850e-01, 4.45590e-02, 3.25338e-01, 1.00000e+00],
 [3.97674e-01, 8.32570e-02, 4.33183e-01, 1.00000e+00],
 [6.21685e-01, 1.64184e-01, 3.88781e-01, 1.00000e+00],
 [8.32299e-01, 2.83913e-01, 2.57383e-01, 1.00000e+00],
 [9.61293e-01, 4.88716e-01, 8.42890e-02, 1.00000e+00],
 [9.81173e-01, 7.59135e-01, 1.56863e-01, 1.00000e+00],
 [.95, .95, 6.44924e-01,1.00000e+00]]'''

ticks = np.arange(0, 18.5, 5)
boo1 = xdata < 25
# for i in range(0, len(ydata)):
#    y = ydata[i]
#    ye = yerr[i]
#    plt.errorbar( xdata[boo1], y[boo1]+i*100,yerr=ye[boo1]*0, marker=markers[i], label=labels[i], color=colors[i])

print(ydata)

x = xdata[boo1]


def binomial2(x, n, p, A, n2, p2, A2):
    mono = stats.binom.pmf(x, n, p) * A
    di = stats.binom.pmf(x[::2], n2, p2) * A2
    mono[::2] += di
    return mono


def binomial_fit(xvals, yvals):
    xvals = np.array(xvals)
    Aguess = np.amax(yvals)
    Aguess = 80
    nguess = np.amax(xvals)
    pguess = float(xvals[np.argmax(yvals)]) / (1.5*float(nguess))
    guess = [nguess, pguess, Aguess, nguess, pguess, Aguess]
    # print(guess)
    fits = curve_fit(binomial2, xvals, yvals, p0=guess, maxfev=1000000)[0]
    fitdat = binomial2(xvals, fits[0], fits[1], fits[2], fits[3], fits[4], fits[5])
    return fits, fitdat


for i in range(0, len(x)):
    y = ydata[:, i]
    n = np.array(labels[:len(y)]).astype(float)
    # print(n, i , len(x))
    # print(y)
    ye = yerr[:, i]
    fits, fitdat = binomial_fit(n, y)
    print(fits)
    err = np.abs(fitdat - y)
    color = colors[i]
    plt.errorbar(n, y, yerr=ye, color=color)
    plt.plot(n, fitdat, color=color, linestyle="--")
plt.show()
exit()

x = np.arange(0, len(ydata))
print(x)
for i in range(0, len(ydata[0])):
    y = ydata[:, i]
    ye = yerr[:, i]
    plt.bar(x, y, yerr=ye, bottom=-100 * i, width=0.4, color=colors, label=xdata[i])
    bars = plt.bar(x + 0.5, y, yerr=ye, bottom=-100 * i, width=0.4, color=colors, label=xdata[i], edgecolor="w")
    for bar in bars:
        bar.set_hatch("/////")
    plt.text(7.5, 20 - 100 * i, str(xdata[i]))
plt.text(5.5, 80, "Peptide/Nanodisc")
# plt.errorbar(x, y-0.2, yerr=yerr, marker="o")
# plt.errorbar(x, y-0.4, yerr=yerr, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ax.set_xlim(0, 18.5)
ax.set_ylim(-500, 100)
# plt.xticks(ticks)

ticks = np.arange(0, 101, 50)
plt.yticks(ticks, labels=["0", "%", "100"])

# plt.xlabel("Peptide Ratio")
# plt.ylabel("Relative Intensity (%)")
# plt.legend()
plt.show()

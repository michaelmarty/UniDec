import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib

# rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

path = "Z:\\Group Share\\Larry Walker\\Figures\\Paper\\Data for figures\\Gramicidin Titrations\\"
os.chdir(path)

file1 = "DMPG_Averages.txt"
file2 = "DMPC_Averages.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DMPG_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DMPC_Errors.txt")

for i in range(0, len(ydata1[0])):
    factor = np.sum(ydata1[:, i]) / 100
    ydata1[:, i] /= factor
    yerr1[:, i] /= factor

for i in range(0, len(ydata2[0])):
    factor = np.sum(ydata2[:,i]) / 100
    ydata2[:, i] /= factor
    yerr2[:, i] /= factor
# factor = np.amax(ydata) / 100
# ydata /= factor

xdata = np.loadtxt("DMPG_XValues.txt")
labels = ["0", "2", "4", "6", "8", "10", "12", "14"]
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

ticks = np.arange(0, 18.5, 5)
boo1 = xdata < 25
# for i in range(0, len(ydata)):
 #   y = ydata[i]
  #  ye = yerr1[i]
   # plt.errorbar( xdata[boo1], y[boo1]+i*100,yerr1=ye[boo1]*0, marker=markers[i], label=labels[i], color=colors[i])

print(ydata1)
print(ydata2)

x = np.arange(0,len(ydata1))*2
wid=0.8
print(x)
for i in range(0, len(ydata1[0])):
    y = ydata1[:, i]
    ye = yerr1[:, i]
    plt.bar(x, y, yerr=ye, bottom=-100*i, width = wid, color=colors, label=xdata[i])
    #bars=plt.bar(x, y, yerr=ye, bottom=-100 * i, width=wid, color=colors, label=xdata[i], edgecolor="w")
    #for bar in bars:
    #    bar.set_hatch("/////")
    plt.text(12.5, 20-100*i,str(xdata[i]))

x2 = np.arange(0,len(ydata2))*2
print(x2)
for i in range(0, len(ydata2[0])):
    y2 = ydata2[:, i]
    ye2 = yerr2[:, i]
    plt.bar(x2+1, y2, yerr=ye2, bottom=-100*i, width = wid, color=colors, label=xdata[i])
    bars2=plt.bar(x2+1, y2, yerr=ye2, bottom=-100 * i, width=wid, color=colors, label=xdata[i], edgecolor="w")
    for bar in bars2:
        bar.set_hatch("/////")
    #plt.text(7.5, 20-100*i,str(xdata[i]))

plt.text(10.5, 80,"Peptide/Nanodisc")
# plt.errorbar(x, y-0.2, yerr1=yerr1, marker="o")
# plt.errorbar(x, y-0.4, yerr1=yerr1, marker="o")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ax.set_xlim(0, 18.5)
ax.set_ylim(-500, 100)
# plt.xticks(ticks)

ticks = np.arange(0, 101, 50)
plt.yticks(ticks, labels=["0","%","100"])

# plt.xlabel("Peptide Ratio")
#plt.ylabel("Relative Intensity (%)")
#plt.legend()
plt.show()

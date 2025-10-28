import numpy as np
import matplotlib.pyplot as plt
import os

path = "Z:\Group Share\Deseree\Overlay Data"

files = ["T0T1T2 Overlay.asc", "T1 60.asc", "T2 no void.asc","G1 60.asc", "R1 60.asc"]    # "Mixedbelt.asc"]
row = [1, 1, 1, 1, 1]
corr = [0,0,0,0,0]
plt.figure(figsize=(9, 9))
for i, f in enumerate(files):
    fname = os.path.join(path, f)
    data = np.loadtxt(fname, skiprows=3)
    ydata = data[:, row[i]]
    ydata = ydata - np.amin(ydata)
    #ydata=ydata-i
    xdata = data[:, 0]-corr[i]

    b1 = xdata > 8
    b2 = xdata < 20
    b3 = np.logical_and(b1, b2)
    ydata = ydata / np.amax(ydata[b3])
    labels = ['T0','T1',' T2', 'G1','R1']
    colors = ['b', 'k', 'red', 'orange','green']

    plt.plot(xdata[b3], ydata[b3], label=labels[i], color=colors[i])

plt.legend(loc='upper right', title="MSP")
plt.xlabel("Retention Time (mL)")
plt.ylabel("Normalized A280 (mAU)")
plt.show()

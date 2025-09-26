import numpy as np
import matplotlib.pyplot as plt
import os

path = "C:\\Users\\Marty Lab\\Box Sync\\LabShared\\Lab Data\\Nanodisc Spectra"

files = ["N1.asc"]  # ,"D1T3 POPC.asc"]
row = [1]
corr = [0]
plt.figure(figsize=(9, 9))
for i, f in enumerate(files):
    fname = os.path.join(path, f)
    data = np.loadtxt(fname, skiprows=3)
    ydata = data[:, row[i]]
    ydata = ydata - np.amin(ydata)
    ydata = ydata / np.amax(ydata)
    xdata = data[:, 0] - corr[i]

    b1 = xdata > 5
    b2 = xdata < 20
    b3 = np.logical_and(b1, b2)


    plt.plot(xdata[b3], ydata[b3])

plt.xlabel("Retention Time (mL)")
plt.ylabel("Normalized Absorbance")
plt.savefig("C:\\Users\Marty Lab\Box Sync\LabShared\Papers\James\Aubin\N1212.pdf")
plt.show()



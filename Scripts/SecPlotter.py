import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42



#path = "Z:\Group Share\Deseree"
path = "E:"

#files = ["G1 R1 Overlay.asc", "T0T1T2 Overlay.asc", "P1B.asc", "PB.asc", "PAX.asc", "P1D.asc"]    # "Mixedbelt.asc"]
files = ["JEK_20210401_AMTB_PCPGPE_R2.csv"]
column = [1, 1, 1, 1, 1, 1]
corr = [0, 0, 0, 0, 0, 0]

plt.figure(figsize=(6, 4))

for i, f in enumerate(files):
    # Create the file path
    fname = os.path.join(path, f)
    # Load the data
    data = np.loadtxt(fname, skiprows=1, delimiter=",")

    # Extract the x data
    xdata = data[:, 0]
    # Extract the y data
    ydata = data[:, column[i]]

    # Setting plotting range
    b1 = xdata > 5
    b2 = xdata < 20
    b3 = np.logical_and(b1, b2)

    #Slice it
    xdata = xdata[b3]
    ydata = ydata[b3]

    # Set the minimum to 0 and normalize to a max of 1
    ydata = ydata - np.amin(ydata)
    #ydata = ydata /np.amax(ydata)

    # Set the Labels and Colors
    #labels = ['AmtB', 'AqpZ', '1', '1-2', '2', '0-2']
    colors = ['k', 'k', 'y', 'orange', 'r', 'purple']

    plt.plot(xdata, ydata, color=colors[i])

#plt.legend(loc='upper right')#, title="Extra T's")
plt.xlabel("Retention Time (mL)")
plt.ylabel("A280 (mAU)")
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.ylim(0,400)
plt.xlim(5,20)
plt.show()

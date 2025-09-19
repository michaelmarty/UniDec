import numpy as np
import matplotlib.pyplot as plt
import os

path= "Z:\wresager\SEC"

files = ["rhoddm.asc","rhoddm_bleached.asc"]#,"D1T3 POPC.asc"]
row=[1,1,1,1]
corr=[0,0,0,1.52]
labels=["DDM", "DDM Bleached"]
plt.figure(figsize=(2,2))
for i,f in enumerate(files):
    fname = os.path.join(path,f)
    data=np.loadtxt(fname,skiprows=3)
    ydata = data[:,row[i]]
    ydata= ydata - np.amin(ydata)
    ydata= ydata /np.amax(ydata)
    xdata=data[:,0]-corr[i]

    b1=xdata>5
    b2=xdata<20
    b3=np.logical_and(b1,b2)

    plt.plot(xdata[b3],ydata[b3], label=labels[i])

plt.xlabel("Retention Time (min)")
plt.ylabel("Normalized A280")
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os

path= "E:\\"

files = ["6132017rhodopsin chaps.asc", "6132017rhodopsin chaps bleached.asc", "6152017rhodopsin c8e4.asc", "6152017rhodopsin c8e4 bleached.asc"]#,"D1T3 POPC.asc"]
row=[1,1,1,1]
corr=[0,0,0,0]
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

    plt.plot(xdata[b3],ydata[b3])

plt.xlabel("Retention Time (mL)")
plt.ylabel("Normalized A280")
plt.show()

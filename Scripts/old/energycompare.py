import numpy as np
import matplotlib.pyplot as plt
import os

topdir="C:\\Data\\New"
sample=["POPC_D1T1m","POPC_D1T0m","POPC_D1T0-1m","POPC_D1T0-2m","POPC_D1T1-2m","POPC_D1T2m"]
types=["HCD","ISTRAP"]
energies=[np.arange(0,250,20),np.arange(0,300,25)]
colors=["r","b"]

plt.figure()
for j,s in enumerate(sample):
    for i,type in enumerate(types):
        print(i, type)
        energy=energies[i]
        newdir = os.path.join(topdir, s + "_" + types[i])
        extracts=np.loadtxt(os.path.join(newdir,"extracts.txt"))

        #grid=np.loadtxt(os.path.join(newdir,"ExtractFull2D_Da.txt"))

        b1=extracts<149400
        if j==0:
            label=type
        else:
            label=None
        plt.plot(energy[b1],extracts[b1]/1000, marker='o', color=colors[i],label=label)


plt.xlabel("Volts")
plt.ylabel("50% Threshold Center of Mass (kDa)")
plt.ylim([110.000,150.000])
plt.legend()
plt.savefig("C:\\Users\Michael Marty\Box Sync\LabShared\MS_Data\dissociation.pdf")
plt.show()
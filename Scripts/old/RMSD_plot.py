__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
#rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
rcParams['font.size']=8


residues=[78,812,445]

topdir="C:\\VMD\\AmtB"

plt.figure(figsize=(3.13,2))
for i,r in enumerate(residues):
    file=os.path.join(topdir,"rmsdres"+str(r)+".dat")
    data=np.loadtxt(file,skiprows=2)
    plt.plot(data[:,0]/10.,data[:,1])
plt.xlabel("Time (ns)",fontsize=8)
plt.ylabel("RMSD (${\AA}$)",fontsize=8)

plt.gca().set_xlim(20,100)
plt.tight_layout()
plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\RMSD.png")
plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\RMSD.pdf")
plt.show()

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
from scipy import interpolate
from scipy.interpolate import interp1d
from matplotlib import rcParams


# Data for plotting
path = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\FIGURES"
os.chdir(path)

file1 = "x values for figure 2.txt"
file2 = "POPC SOPE y values.txt"
file3 = "POPC SOPE error.txt"
file4 ="POPC SOPE POPS y values.txt"
file5 ="POPC SOPE POPS error.txt"

t = np.loadtxt(file1)
y1 = np.loadtxt(file2)
e = np.loadtxt(file3)
y2 = np.loadtxt(file4)
e2 = np.loadtxt(file5)



#rcParams['font.size']=10
fig, ax = plt.subplots()


#f2=interp1d(t,y1,kind='cubic')
#f2=interp1d(t,y2,kind='cubic')
xnew = np.linspace(1, 36,num=360,endpoint=True)
ax.set(ylim=(746,762))
ax.plot(t, y1)
ax.plot(t, y2)
plt.plot(t,y1, 'o')
plt.errorbar(t,y1,yerr=e,fmt='none', capsize=3, capthick=1, elinewidth=1)

plt.plot(t,y2, 'o')
plt.errorbar(t,y2,yerr=e2,fmt='none', capsize=3, capthick=1, elinewidth=1)


ax.set_ylabel('Average Lipid Mass', size=16)
ax.set_xlabel('Collision Voltage', size=16)
ax.set_title('GC 60-200V of AmtB 50%PG Nanodisc', size=20)
#ax.grid()
plt.show()
#fig.savefig("test.png")
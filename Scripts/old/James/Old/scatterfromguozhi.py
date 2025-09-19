import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
from scipy import interpolate
from scipy.interpolate import interp1d
from matplotlib import rcParams


# Data for plotting
path = "Z:\\Group Share\\James Keener\\PYTHON SCATTER 1"
os.chdir(path)

file1 = "xaxislipids.txt"
file2 = "yaxismass.txt"
file3 = "error.txt"

t = np.loadtxt(file1)
y1 = np.loadtxt(file2)
e = np.loadtxt(file3)



#rcParams['font.size']=10
fig, ax = plt.subplots()


f2=interp1d(t,y1,kind='cubic')
xnew = np.linspace(1, 36,num=360,endpoint=True)
ax.set(ylim=(749,763))
#ax.plot(t, y)
plt.plot(t,y1, 'o')
plt.errorbar(t,y1,yerr=e,fmt='none', capsize=3, capthick=1, elinewidth=1)



ax.set_ylabel('Average Lipid Mass', size=16)
ax.set_xlabel('Number of Lipids', size=16)
ax.set_title('GC 60-200V of AmtB 50%PG Nanodisc', size=20)
#ax.grid()
plt.show()
#fig.savefig("test.png")
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(t, a, b, c):
    #return a * np.log(b * (d*t)) + c
    return a * (1-np.exp(-b * t))+ c

path = "p:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\Empty Exchange\\Size at high conc"
os.chdir(path)

# A and E
file1 = "Time.txt"
file2 = "D1.txt"
file3 = "E3.txt"

t = np.loadtxt(file1)/60.
y1 = np.loadtxt(file2)*100
y1err1 = np.loadtxt("D1_Error.txt")*100
y2 = np.loadtxt(file3)*100
y2err2 = np.loadtxt("E3_Error.txt")*100

p0=[np.amax(y1), 1/np.mean(t), np.amin(y1)]

rcParams['font.size']=16

tfit = np.arange(0, np.amax(t), 0.1)

final1=34.91
ferr1=6.15

popt, pcov = curve_fit(func, t, y1, p0=p0)
plt.figure()

plt.errorbar(t, y1, color="r", marker="o", linestyle="", label="MSP1D1 (9.8 nm)", yerr=y1err1, capsize=3)
plt.fill_between(t, t*0+ final1-ferr1, t*0+ final1+ferr1, color="r", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), color="r", linestyle="--", label=None)


final2=38.28
ferr2=3.49

popt, pcov = curve_fit(func, t, y2, p0=p0)

plt.errorbar(t, y2, color="navy", marker="o", linestyle="", label="MSP1E3D1 (13.1 nm)", yerr=y2err2, capsize=3)
plt.fill_between(t, t*0+ final2-ferr2, t*0+ final2+ferr2, color="navy", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), color='navy', linestyle='--', label=None)

plt.xlabel('Time (hr)', size=20)
plt.ylabel('% POPC', size=20)

ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.legend(loc=4)
plt.tight_layout()
plt.title('[Lipid] (1657 \u00B5M)')
#dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
#file="size"
#plt.savefig(os.path.join(dir, file+".pdf"), transparent=True, dpi=1200)
#plt.savefig(os.path.join(dir, file+".png"), transparent=True, dpi=1200)

plt.show()
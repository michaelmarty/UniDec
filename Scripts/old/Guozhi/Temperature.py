import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(t, a, b, c):
    #return a * np.log(b * (d*t)) + c
    return a * (1-np.exp(-b * t))+ c

path = "z:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\Empty Exchange\\Temperature"
os.chdir(path)


file1 = "Time1.txt"
file3 = "4 deg C.txt"

file2 = "Time2.txt"
file4 = "roomT.txt"

t1 = np.loadtxt(file1)/60.
t2 = np.loadtxt(file2)/60.

y1 = np.loadtxt(file3)*100
y1err1 = np.loadtxt("4 deg C_Errors.txt")*100
y2 = np.loadtxt(file4)*100
y2err2 = np.loadtxt("roomT_Errors.txt")*100

p0=[np.amax(y1), 1/np.mean(t1), np.amin(y1)]

rcParams['font.size']=16



final1=43.91
ferr1=0.38


popt, pcov = curve_fit(func, t1, y1, p0=p0)
print(popt)
tfit1 = np.arange(0, np.amax(t1), 0.1)

plt.figure()
plt.errorbar(t1, y1, color="navy", marker="o", linestyle="", label="4 \u00B0C", yerr=y1err1, capsize=3)
#plt.hlines(final1, 0, np.amax(t1), color="navy", linestyle="-")
#plt.hlines(final1-ferr1, 0, np.amax(t1), color="navy", linestyle="--")
#plt.hlines(final1+ferr1, 0, np.amax(t1), color="navy", linestyle="--")
plt.fill_between(tfit1, tfit1*0+ final1-ferr1, tfit1*0+ final1+ferr1, color="navy", alpha=0.2)
plt.plot(tfit1, func(tfit1, *popt), color='navy', linestyle="--", label=None)
plt.legend()


final2=36.43
ferr2=2.10

popt, pcov = curve_fit(func, t2, y2, p0=p0)
print(popt)
plt.errorbar(t2, y2, color="b", marker="o", linestyle="", label="22 \u00B0C", yerr=y2err2, capsize=3)
plt.fill_between(tfit1, tfit1*0+ final2-ferr2, tfit1*0+ final2+ferr2, color="b", alpha=0.2)
#plt.hlines(final2, 0, np.amax(t1), color="b", linestyle="-")
#plt.hlines(final2-ferr2, 0, np.amax(t1), color="b", linestyle="--")
#plt.hlines(final2+ferr2, 0, np.amax(t1), color="b", linestyle="--")


tfit2 = np.arange(0, np.amax(t2), 0.1)
plt.plot(tfit2, func(tfit2, *popt), 'b--', label=None)

plt.xlabel('Time (hr)', size=20)
plt.ylabel('% POPC', size=20)

ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.legend(loc=4)
plt.tight_layout()

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
file="temp"
plt.savefig(os.path.join(dir, file+".pdf"), transparent=True, dpi=1200)
plt.savefig(os.path.join(dir, file+".png"), transparent=True, dpi=1200)

plt.show()
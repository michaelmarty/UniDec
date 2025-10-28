import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(t, a, b, c):
    #return a * np.log(b * (d*t)) + c
    return a * (1-np.exp(-b * t))+ c

path = "p:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\Empty Exchange\\Concentration\\E3"
os.chdir(path)

file1 = "TIME.txt"
file2 = "Ori.txt"
file3 = "Dil.txt"
#file4 = ".txt"

t = np.loadtxt(file1)/60.

tfit = np.arange(0, np.amax(t), 0.1)

y1 = np.loadtxt(file2)*100
y1err1 = np.loadtxt("Ori_Errors.txt")*100
y2 = np.loadtxt(file3)*100
y2err2 = np.loadtxt("Dil_Errors.txt")*100
#y3 = np.loadtxt(file4)*100
#y3err3 = np.loadtxt("very low_Error.txt")*100

p0=[np.amax(y1), 1/np.mean(t), np.amin(y1)]

rcParams['font.size']=16

final1=38.28
ferr1=3.49

popt, pcov = curve_fit(func, t, y1, p0=p0)
plt.figure()
plt.errorbar(t, y1, color="navy", marker="o", linestyle="", label="6.4 \u00B5M", yerr=y1err1, capsize=3)
#plt.hlines(final1, 0, np.amax(t), color="navy", linestyle="-")
plt.fill_between(t, t*0+ final1-ferr1, t*0+ final1+ferr1, color="navy", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), color="navy", linestyle="--", label=None)

final2=42.93
ferr2=4.17


popt, pcov = curve_fit(func, t, y2, p0=p0)
plt.errorbar(t, y2, color="b", marker="o", linestyle="", label="2.1 \u00B5M", yerr=y2err2, capsize=3)
#plt.hlines(final2, 0, np.amax(t), color="b", linestyle="-")
plt.fill_between(t, t*0+ final2-ferr2, t*0+ final2+ferr2, color="b", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), 'b--', label=None)


#final3=
#ferr3=


#popt, pcov = curve_fit(func, t, y3, p0=p0)
#plt.errorbar(t, y3, color="skyblue", marker="o", linestyle="", label="1.1 \u00B5M", yerr=y3err3, capsize=3)
#plt.fill_between(t, t*0+ final3-ferr3, t*0+ final3+ferr3, color="skyblue", alpha=0.2)
#plt.hlines(final3, 0, np.amax(t), color="skyblue", linestyle="-")
#plt.plot(tfit, func(tfit, *popt), color="skyblue", linestyle="--", label=None)


plt.xlabel('Time (hr)', size=20)
plt.ylabel('% POPC', size=20)

ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.legend(loc=4, title="[MSP1E3D1]", fontsize=14)
plt.tight_layout()

#dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
#file="conc"
#plt.savefig(os.path.join(dir, file+".pdf"), transparent=True, dpi=1200)
#plt.savefig(os.path.join(dir, file+".png"), transparent=True, dpi=1200)

plt.show()
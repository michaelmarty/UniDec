import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit, newton
from matplotlib import rcParams

from copy import deepcopy
#import pandas as pd
#import statsmodels.api as sm



def func(t, a, b, c):
    #return a * np.log(b * (d*t)) + c
    y = a * (1-np.exp(-b * t))+ c
    return y

def get_hl(t, y1):
    p0 = [np.amax(y1), 1 / np.mean(t), np.amin(y1)]
    popt, pcov = curve_fit(func, t, y1, p0=p0)
    k = popt[1]
    hl = np.log(2)/k
    return hl

def bootstrap(t, ydata, leaveout=1, n=100):
    l = len(t)
    cos=[]
    for i in range(n):
        x = deepcopy(t)
        y = deepcopy(ydata)
        choice = np.random.choice(l, l-leaveout, replace=True)
        #drop = np.random.randint(0, l)
        #newx = np.delete(x, drop)
        #newy = np.delete(y, drop)
        newx = x[choice]
        newy = y[choice]
        try:
            co = get_hl(newx, newy)
            cos.append(co)
        except:
            pass
    #plt.hist(cos)
    #plt.show()
    print(np.mean(cos),"+/-", np.std(cos))


#halflife = np.log(2) / 0.079404
#halflife1= np.log(2) / (0.079404+0.01168162)
#halflife2= np.log(2) / (0.079404-0.01168162)
#print('halflife', halflife, halflife1, halflife2)

t = np.array([4, 49, 214, 395, 902, 2622])/60.
ydata = np.array ([0.095368, 0.06988, 0.097227, 0.121098043, 0.182946392, 0.344652,
                   0.039601, 0.10238, 0.153579, 0.180558186, 0.225052, 0.347536,
                   0.049713, 0.084317, 0.133108686, 0.146389472, 0.23493, 0.36619])
#([0.064421576, 0.083878092, 0.159556773, 0.183348, 0.24558688, 0.319126423,   # batch a
                 # 0.029828785, 0.073093514, 0.147610869, 0.167533354, 0.236034583, 0.352859551,  # batch b
                  #0.018291219, 0.039799601, 0.082531323, 0.136846317, 0.217191265, 0.299744443])  # batch c

tdata = np.ravel([t,t,t])


ay = ydata[:6]
by = ydata[6:12]
cy = ydata[12:]

yav = np.mean([ay, by, cy], axis=0)*100
yerr1 = np.std([ay, by, cy], axis=0)*100

hls = [get_hl(t, ay), get_hl(t, by), get_hl(t, cy)]
print(np.mean(hls), np.std(hls))

print(t, yav)
p0=[np.amax(yav), 1/np.mean(t), np.amin(yav)]
popt, pcov = curve_fit(func, tdata, ydata*100, p0=p0)

print('parameter', popt)
p_sigma = np.sqrt(np.diag(pcov))
print('sqrt error', p_sigma)
print("Global Fit", np.log(2)/popt[1] )
bootstrap(tdata, ydata)


rcParams['font.size']=16

tfit = np.arange(0, np.amax(t), 0.1)

final1=36.43
ferr1=2.1


plt.figure()

#plt.plot(t, ay*100)
#plt.plot(t, by*100)
#plt.plot(t, cy*100)

plt.errorbar(t, yav, color="navy", marker="o", linestyle="", label="MSPD1 (6.4 \u00B5M)", yerr=yerr1, capsize=3)
plt.fill_between(t, t*0+ final1-ferr1, t*0+ final1+ferr1, color="navy", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), color="navy", linestyle="--", label=None)






plt.xlabel('Time (hr)', size=20)
plt.ylabel('% POPC', size=20)

ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.legend(loc=4)
plt.tight_layout()

#dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
#file="size"
#plt.savefig(os.path.join(dir, file+".pdf"), transparent=True, dpi=1200)
#plt.savefig(os.path.join(dir, file+".png"), transparent=True, dpi=1200)

plt.show()




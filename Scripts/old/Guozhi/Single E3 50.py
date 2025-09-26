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

def func2(t, b):
    #return a * np.log(b * (d*t)) + c
    y = -26 * (1-np.exp(-b * t))+ 55
    return y

def get_hl2(t, y1):
    #p0 = [np.amax(y1), 1 / np.mean(t), np.amin(y1)]
    p0 = [1/np.mean(t)]
    popt, pcov = curve_fit(func2, t, y1, p0=p0)
    k = popt[0]
    hl = np.log(2)/k
    return hl

def get_hl(t, y1):
    p0 = [np.amax(y1), 1 / np.mean(t), np.amin(y1)]
    #p0 = [1/np.mean(t)]
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

t = np.array([0, 4, 49, 214, 395, 1521, 2643, 3983])/60.
ydata = np.array([0.557742264, 0.469541209, 0.459261728, 0.424115368, 0.397571634, 0.374643974, 0.358313764, 0.290636321,	 # batch a
                  0.545783287, 0.419092582, 0.424415491, 0.392409428, 0.35904848, 0.346535932, 0.344401301, 0.295479672,  # batch b
                  0.554598383, 0.423801979, 0.431666285, 0.427627786, 0.391957897, 0.350120913, 0.343860789, 0.317144315])  # batch c

tdata = np.ravel([t,t,t])


ay = ydata[:8]
by = ydata[8:16]
cy = ydata[16:]

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
print("Global Fit", get_hl(tdata, ydata*100))
bootstrap(tdata, ydata*100)


rcParams['font.size']=16

tfit = np.arange(0, np.amax(t), 0.1)

final1=28.89
ferr1=0.7



plt.figure()

#plt.plot(t, ay*100)
#plt.plot(t, by*100)
#plt.plot(t, cy*100)

plt.errorbar(t, yav, color="navy", marker="o", linestyle="", label="MSP1E3D1 (2.4 \u00B5M)", yerr=yerr1, capsize=3)
plt.fill_between(t, t*0+ final1-ferr1, t*0+ final1+ferr1, color="navy", alpha=0.2)
plt.plot(tfit, func(tfit, *popt), color="navy", linestyle="--", label=None)



plt.xlabel('Time (hr)', size=20)
plt.ylabel('% POPC', size=20)

ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.legend(loc=1)
plt.tight_layout()

#dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\LXMS and ENDMS"
#file="size"
#plt.savefig(os.path.join(dir, file+".pdf"), transparent=True, dpi=1200)
#plt.savefig(os.path.join(dir, file+".png"), transparent=True, dpi=1200)

plt.show()




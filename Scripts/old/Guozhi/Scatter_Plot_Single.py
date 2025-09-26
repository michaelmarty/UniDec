

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(x, a, b, c,d):
    return a * np.log(b * (d*x)) + c

# Data for plotting
path = "Z:\\Group Share\\Guozhi Zhang\\SynaptXS\\4 year presentation figure"
os.chdir(path)

file1 = "x_time_mins.txt"
file2 = "y_percent PC.txt"

x = np.loadtxt(file1)
y1 = np.loadtxt(file2)



#rcParams['font.size']=10
fig, ax = plt.subplots()

#final1=100
#plt.hlines(final1, 0, np.amax(t), color="navy", linestyle="--")
ax.scatter(x, y1, marker='o', label="Mixed nanodiscs", color="black")
#ax.plot(x, y)

rcParams['font.size']=16
popt, pcov = curve_fit(func, x, y1)
#plt.figure()

plt.plot(x, func(x, *popt), 'r--', label="Fitted Curve")

ax.set_ylabel('% PC', size=16)
ax.set_xlabel('Time (mins)', size=16)
ax.set_title('Exchange Plot', size=20)
#ax.grid()
plt.legend()
plt.show()
#fig.savefig("test.png")
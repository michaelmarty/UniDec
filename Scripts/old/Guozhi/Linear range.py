

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(x, a, b):
    return (a*x) + b

# Data for plotting
path = "Z:\\Group Share\\Guozhi Zhang\\SynaptXS\\4 year presentation figure"
os.chdir(path)

file1 = "x_log_conc.txt"
file2 = "y_popc_log_intensity.txt"
file3 = "y_dmpg_log_intensity.txt"
file4 = "y_popg_log_intensity.txt"
x = np.loadtxt(file1)
y1 = np.loadtxt(file2)
y2 = np.loadtxt(file3)
y3 = np.loadtxt(file4)

rcParams['font.size']=20
fig, ax = plt.subplots()


ax.scatter(x, y1, marker='o', label="POPC, positive mode", color="black")
popt, pcov = curve_fit(func, x, y1)
plt.plot(x, func(x, *popt), color="black", linestyle="--")#, label="PC Fitted Curve")


ax.scatter(x, y2, marker='^', label="DMPG, negative mode", color="g")
popt, pcov = curve_fit(func, x, y2)
plt.plot(x, func(x, *popt), color="g", linestyle="--")#, label="PG Fitted Curve")

ax.scatter(x, y3, marker='^', label="POPG, negative mode", color="r")
popt, pcov = curve_fit(func, x, y3)
plt.plot(x, func(x, *popt), color="r", linestyle="--")

ax.set_ylabel("Log(intensity)")
ax.set_xlabel("Log(concentration)")
#ax.set_title('Response Curve')

rcParams['font.size']=16
plt.legend()
plt.show()


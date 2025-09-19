

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

file1 = "x_ther0.txt"
file2 = "y_actual_pc.txt"
file3 = "y_actual_pg.txt"
x = np.loadtxt(file1)
y1 = np.loadtxt(file2)
y2 = np.loadtxt(file3)

rcParams['font.size']=20
fig, ax = plt.subplots()


#ax.scatter(x, y1, marker='o', label="PC, positive mode", color="black")
#popt, pcov = curve_fit(func, x, y1)
#plt.plot(x, func(x, *popt), color="black", linestyle="--")#, label="PC Fitted Curve")


ax.scatter(x, y2, marker='^', label="PG, negative mode", color="g")
popt, pcov = curve_fit(func, x, y2)
plt.plot(x, func(x, *popt), color="g", linestyle="--")#, label="PG Fitted Curve")

ax.set_ylabel('Predicted Response')
ax.set_xlabel('Theoretical Response')
#ax.set_title('Response Curve')

rcParams['font.size']=16
plt.legend()
plt.show()


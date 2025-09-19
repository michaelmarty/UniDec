import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams

def func(t, a, b, c, d):
    return a * np.log(b * (d*t)) + c
    #return a * np.exp(-b * x) + c

path = "P:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\Empty Exchange\\D1"
os.chdir(path)

file1 = "D1 Room Temp Time.txt"
file2 = "D1 Room Temp y B1.txt"
t = np.loadtxt(file1)
y = np.loadtxt(file2)

rcParams['font.size']=16
popt, pcov = curve_fit(func, t, y)
plt.figure()


plt.plot(t, y, 'ko', label="D1 Room T Exchange")
plt.plot(t, func(t, *popt), 'r--', label="Fitted Curve")
plt.xlabel('Time (m)', size=20)
plt.ylabel('%PC', size=20)
plt.legend()
plt.show()
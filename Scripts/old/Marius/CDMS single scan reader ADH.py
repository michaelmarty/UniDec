import numpy
import numpy as np
import h5py
import os
import scipy
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import ms_deisotope
from ms_deisotope import plot



dir = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\Data"
os.chdir(dir)


data = np.loadtxt(fname="ADH single scan.txt", delimiter="\t", skiprows=(7))
#load data and skip first 7 rows

x = (data[:, 0])
#pull x data for column 1
y = (data[:, 1])
#pull y fata for column 2
y = (0.25*y)/70.79778994140625
#convert to S/N with Gmax
# Setting plotting range
b1 = x > 5000
b2 = x < 6200
b3 = np.logical_and(b1, b2)
#slice
x = x[b3]
y = y[b3]
#set graph limit
#S/N = [0,16]
#set m/z limit
#m/z = [10000,12500]
#Plot
rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
plt.figure()
#plt.scatter(x,y, color = "black")
plt.vlines(x,y, ymax= 0, colors = "black")
ax=plt.gca()
ax.set(ylim=(0,6))
ax.set(xlim=(5000,6200))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.xlim(10000,12500)
#plt.ylim(10000,12500)
plt.xlabel("m/z")
plt.ylabel("S/N")
plt.show()
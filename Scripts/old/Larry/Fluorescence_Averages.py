import easygui
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import os
from matplotlib import rcParams

easygui.msgbox("This Script will make a Normalized X,Y Scatter Plot", "Fluorescence", "Continue")
file=easygui.fileopenbox(default="Z:/Group Share/Larry Walker/Fluorescence/")
Time_steps=easygui.integerbox(msg="Enter time step for each addition separated by a comma")

base=os.path.basename(file)
filename=os.path.splitext(base)
print(base)
print(filename[0])

Time = ()
Intensity = ()

with open(file) as f:
    n = len(list(f))
print("n=", n)

Time = np.append(Time, np.loadtxt(file, skiprows=4, dtype=str, delimiter='\t', max_rows=n-5, usecols=(0)))
print("Time=", len(Time))
Intensity = np.append(Intensity, np.loadtxt(file, skiprows=4, delimiter='\t', max_rows=n-5, usecols=(1)))

Time = pd.to_numeric(Time)
Time[:] = [x / 60 for x in Time]
for i in Intensity:
    Intensity = (Intensity - Intensity.min()) / (Intensity.max() - Intensity.min())

fig, ax=plt.subplots()
ax.plot(Time, Intensity)

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
fig.suptitle(filename[0], fontsize=30)
ax.set_xlabel("Time (Mins)", fontsize=20)
ax.set_ylabel("Normalized Intensity", fontsize=20)
ax.set_xlim(0, Time[-1])
ax.set_xticks(np.arange(0, Time[-1], step=10))
ax.set_ylim(0, 1.01)
ax.set_yticks(np.arange(0, 1.01, step=0.25))
ax.set_yticklabels(np.arange(0, 101, step=25))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()
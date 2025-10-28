import easygui
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
import os
from matplotlib import rcParams

#easygui.msgbox("This Script will make a Normalized X,Y Scatter Plot", "Fluorescence", "Continue")
file1=easygui.fileopenbox(default="D:/Group Share/Fluorescence/")
file2=easygui.fileopenbox(default="D:/Group Share/Fluorescence/")
file3=easygui.fileopenbox(default="D:/Group Share/Fluorescence/")
file4=easygui.fileopenbox(default="D:/Group Share/Fluorescence/")

base1=os.path.basename(file1)
filename1=os.path.splitext(base1)
base2=os.path.basename(file2)
filename2=os.path.splitext(base2)
base3=os.path.basename(file3)
filename3=os.path.splitext(base3)
base4=os.path.basename(file4)
filename4=os.path.splitext(base4)
#print(base1, base2, base3, base4)
#print(filename1[0], filename2[0], filename3[0], filename4[0])


Time1 = ()
Intensity1 = ()
Time2 = ()
Intensity2 = ()
Time3 = ()
Intensity3 = ()
Time4 = ()
Intensity4 = ()

with open(file1) as f1:
    n = len(list(f1))

with open(file2) as f2:
    m = len(list(f2))

with open(file3) as f3:
    o = len(list(f3))

with open(file4) as f4:
    p = len(list(f4))

print(n,m,o,p)

Time1 = np.append(Time1, np.loadtxt(file1, skiprows=4, dtype=str, delimiter='\t', max_rows=n-5, usecols=(0)))
Intensity1 = np.append(Intensity1, np.loadtxt(file1, skiprows=4, delimiter='\t', max_rows=n-5, usecols=(1)))

Time2 = np.append(Time2, np.loadtxt(file2, skiprows=4, dtype=str, delimiter='\t', max_rows=m-5, usecols=(0)))
Intensity2 = np.append(Intensity2, np.loadtxt(file2, skiprows=4, delimiter='\t', max_rows=m-5, usecols=(1)))

Time3 = np.append(Time3, np.loadtxt(file3, skiprows=4, dtype=str, delimiter='\t', max_rows=o-5, usecols=(0)))
Intensity3 = np.append(Intensity3, np.loadtxt(file3, skiprows=4, delimiter='\t', max_rows=o-5, usecols=(1)))

Time4 = np.append(Time4, np.loadtxt(file4, skiprows=4, dtype=str, delimiter='\t', max_rows=p-5, usecols=(0)))
Intensity4 = np.append(Intensity4, np.loadtxt(file4, skiprows=4, delimiter='\t', max_rows=p-5, usecols=(1)))

Time1 = pd.to_numeric(Time1)
Time1[:] = [x / 60 for x in Time1]
for i in Intensity1:
    Intensity1 = (Intensity1 - Intensity1.min()) / (Intensity1.max() - Intensity1.min())

Time2 = pd.to_numeric(Time2)
Time2[:] = [x / 60 for x in Time2]
for i in Intensity2:
    Intensity2 = (Intensity2 - Intensity2.min()) / (Intensity2.max() - Intensity2.min())

Time3 = pd.to_numeric(Time3)
Time3[:] = [x / 60 for x in Time3]
for i in Intensity3:
    Intensity3 = (Intensity3 - Intensity3.min()) / (Intensity3.max() - Intensity3.min())

Time4 = pd.to_numeric(Time4)
Time4[:] = [x / 60 for x in Time4]
for i in Intensity4:
    Intensity4 = (Intensity4 - Intensity4.min()) / (Intensity4.max() - Intensity4.min())

fig, ax = plt.subplots(1,1)
ax.plot(Time1, Intensity1, label=filename1[0])#, color='grey')
ax.plot(Time2, Intensity2, label=filename2[0])#, color='blue')
ax.plot(Time3, Intensity3, label=filename3[0])#, color='pink')
ax.plot(Time4, Intensity4, label=filename4[0])#, color='black')

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
#fig.suptitle(filename1[0], fontsize=30)
ax.set_xlabel("Time (Mins)", fontsize=20)
ax.set_ylabel("Normalized Intensity", fontsize=20)
#ax.set_xlim(0, Time1[-1])
ax.set_xticks(np.arange(0, Time1[-1], step=10))
ax.set_ylim(0, 1.01)
ax.set_yticks(np.arange(0, 1.01, step=0.25))
ax.set_yticklabels(np.arange(0, 101, step=25))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend()
plt.show()
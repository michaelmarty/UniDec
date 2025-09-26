import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\\Group Share\\Larry Walker\\MassSpecData\\All Peptides\\Melittin\\"
os.chdir(path)

file1 = "DMPG_3to1_Melittin.txt"
file2 = "DMPG_9to1_Melittin.txt"
file3 = "DMPG_18to1_Melittin.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DMPG_3to1_Melittin_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DMPG_9to1_Melittin_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("DMPG_18to1_Melittin_Errors.txt")

labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8']#, '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-width, ydata1, width, color='c', label='3:1')#, yerr=yerr1)
rects2 = ax.bar(x, ydata2, width, color='m', label='9:1')#, yerr=yerr2)
rects3 = ax.bar(x+width, ydata3, width, color='y', label='18:1')#, yerr=yerr3)

rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
#rcParams['lines.linewidth'] =100
#rcParams['axes.linewidth'] =100
rcParams['font.size']=16
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('%', size=24)
ax.set_xlabel('Number of Peptides Incorporated', size=24)
ax.set_yticklabels(np.arange(0, 75, step=(17.5)))
ax.set_title('Melittin')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=0.7)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,0.75, step=0.175), size=20)
fig.tight_layout()

plt.show()
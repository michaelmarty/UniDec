import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "F:\\Group Share\\MassSpecData\\All Peptides\\LL37\\"
os.chdir(path)

file1 = "DPPG_3to1_LL37.txt"
file2 = "DPPG_6to1_LL37.txt"
file3 = "DPPG_9to1_LL37.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DPPG_3to1_LL37_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DPPG_6to1_LL37_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("DPPG_9to1_LL37_Errors.txt")

labels = ['0', '1', '2', '3', '4', '5', '6', '7']#, '8', '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-width, ydata1, width, color='lightblue', label='3:1', yerr=yerr1, capsize=3.0)
rects2 = ax.bar(x, ydata2, width, color='slategray', label='6:1', yerr=yerr2, capsize=3.0)
rects3 = ax.bar(x+width, ydata3, width,color='steelblue', label='9:1', yerr=yerr3, capsize=3.0)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('%', size=24)
ax.set_xlabel('LL-37 Peptides Incorporated', size=24)
ax.set_yticklabels(np.arange(0, 101, step=(25)))
ax.set_title('DPPG')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0, 1.01, step=0.25), size=20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout()

plt.show()
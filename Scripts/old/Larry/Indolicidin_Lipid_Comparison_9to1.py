import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "F:\\Group Share\\Figures\\Lipid_Tail_comparison\\Indolicidin\\Text Files\\"
os.chdir(path)

file1 = "DLPG_9to1_Indolicidin.txt"
file2 = "DMPG_9to1_Indolicidin.txt"
file3 = "DPPG_9to1_Indolicidin.txt"
file4 = "POPG_9to1_Indolicidin.txt"
file5 = "SOPG_9to1_Indolicidin.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DLPG_9to1_Indolicidin_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DMPG_9to1_Indolicidin_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("DPPG_9to1_Indolicidin_Errors.txt")
ydata4 = np.loadtxt(file4)
yerr4 = np.loadtxt("POPG_9to1_Indolicidin_Errors.txt")
ydata5 = np.loadtxt(file5)
yerr5 = np.loadtxt("SOPG_9to1_Indolicidin_Errors.txt")

labels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']#, '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.15  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-2*width, ydata1, width, color='lightblue', label='DLPG')#, yerr=yerr1, capsize=3.0)
rects2 = ax.bar(x-width, ydata2, width, color='slategray', label='DMPG')#, yerr=yerr2, capsize=3.0)
rects3 = ax.bar(x, ydata3, width, color='steelblue', label='DPPG')#, yerr=yerr3, capsize=3.0)
rects4 = ax.bar(x+width, ydata4, width, color='navy', label='POPG')#, yerr=yerr4, capsize=3.0)
rects5 = ax.bar(x+2*width, ydata5, width, color='aqua', label='SOPG')#, yerr=yerr5, capsize=3.0)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('%', size=24)
ax.set_xlabel('Peptides Incorporated', size=24)
ax.set_yticklabels(np.arange(0, 101, step=(25)))
ax.set_title('9:1 Indolicidin:Nanodisc')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,1.01, step=0.25), size=20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#fig.tight_layout()

plt.show()
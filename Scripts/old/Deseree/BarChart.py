import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\\Group Share\\FPOP DATA\\"
os.chdir(path)

file1 = "CombinedData.txt"
#file2 = "BarChartData2.txt"
#file3 = "DMPG_12to1_LactoferricinB_DTT.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("CombinedError.txt")
#ydata2 = np.loadtxt(file2)
#yerr2 = np.loadtxt("BarChartErrors2.txt")
#ydata3 = np.loadtxt(file3)
#yerr3 = np.loadtxt("DMPG_12to1_LactoferricinB_DTT_Errors.txt")

labels = ['Free', 'Nanodiscs']#, '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14']

x = [1, 1.6]  # the label locations
width = 0.5 # the width of the bars
#plt.figure(figsize=(2,4))

fig, ax = plt.subplots(figsize=(3,4))
rects1 = ax.bar(x, ydata1, width, color='c', yerr=yerr1, capsize=3)
#rects2 = ax.bar(x+width, 16.71, width, color='m', label='9:1', yerr=yerr2)
#rects3 = ax.bar(x+width, ydata3, width, color='y', label='12:1', yerr=yerr3)

rcParams['font.size']=10
rcParams['ps.fonttype'] =42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('% Oxidation', size=24)
#ax.set_xlabel('Number of Peptides Incorporated', size=24)
ax.set_yticklabels(np.arange(0, 50, step=(15)))
#ax.set_title('DMPG')
ax.set_xticks(x)
ax.set_xticklabels(labels)
#ax.legend()
plt.xticks(size=20)
#plt.xlim(0,len(labels))
plt.ylim(top=50)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0, 50, step=15), size=20)
fig.tight_layout()
ax.spines['top'].set_visible(False) #Removes box around plot
ax.spines['right'].set_visible(False)

plt.show()
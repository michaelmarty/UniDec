import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

#path = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\Peptide\\LL37\\hamilton Bacterial"
#os.chdir(path)

#file1 = "DPPG 3to1 yvalues .txt"
#file2 = "DPPG SOPE POCL 3to1 yvalues.txt"
#file3 = "DPPG_9to1_LL37.txt"

ydataA = 6986.42
yerrA = 3754.511
ydataB = 18926.6747
yerrB = 6019.495
ydataC = 14573.54
yerrC = 7677.48

labels = ['BSA']# '6', '7', '8', '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-width, ydataA, width, color='lightblue', label='A', yerr=yerrA, capsize=3.0)
rects2 = ax.bar(x, ydataB, width, color='slategray', label='B', yerr=yerrB, capsize=3.0)
rects3 = ax.bar(x+width, ydataC, width,color='steelblue', label='C', yerr=yerrC, capsize=3.0)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('S/N', size=24)
ax.set_xlabel('Needle type', size=24)
#ax.set_yticklabels(np.arange(0, 30001, step=(25)))
ax.set_title('S/N ratio for BSA +15')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=30000.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0, 30000.1, step=10000), size=20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout()

plt.show()
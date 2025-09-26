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

ydataA = 15.94
yerrA = 8.23
ydataB = 32.05
yerrB = 6.07
ydataC =9.16
yerrC = 6.81
#ydataC = 14573.545
#yerrC = 7677.48

labels = ['Flowrate'] #, 'control']# '6', '7', '8', '9', '10', '11', '12', '13', '14']

yvals1 = [11.51, 12.33, 11.67, 28.27]
yvals2 = [25.97, 29.59, 32.38, 40.27]
yvals3 = [12.00, 4.44, 17.42, 2.79]

x = np.arange(len(labels))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-width, ydataA, width, color='steelblue', label='PFDCS', yerr=yerrA, capsize=3.0)
rects2 = ax.bar(x, ydataB, width, color='lightblue', label='PEG', yerr=yerrB, capsize=3.0)
rects3 = ax.bar(x+width, ydataC, width,color='slategrey', label='Control', yerr=yerrC, capsize=3.0)
ax.plot((x-width)* np.ones(len(yvals1)), yvals1, marker="o", linestyle=" ", color="black")
ax.plot((x)* np.ones(len(yvals2)), yvals2, marker="o", linestyle=" ", color="black")
ax.plot((x+width)* np.ones(len(yvals3)), yvals3, marker="o", linestyle=" ", color="black")
rcParams['font.size']=12
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('nL/min', size=24)
ax.set_xlabel('Needle type', size=24)
#ax.set_yticklabels(np.arange(0, 30001, step=(25)))
ax.set_title('BSA Coated vs control')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=45)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0, 45.1, step=5), size=20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout()

plt.show()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib



# data to plot
n_groups = 4
path = "Z:\\Group Share\\Larry Walker\\Figures\\Seattle\\BioSpinComp\\"
os.chdir(path)

files = "GA8to1_DMPG_PRE_4to10.txt"
files1 = "GA8to1_DMPG_POST_4to10.txt"
pre = np.loadtxt(files)
post = np.loadtxt(files1)
#yerrDMPC = np.loadtxt("GA6to1BarDMPCErrors.txt")
#yerrDMPG = np.loadtxt("GA6to1BarDMPGErrors.txt")


#rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 12

# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.35
opacity = 0.8

rects1 = plt.bar(index, pre, bar_width,
                 alpha=opacity,
                 color='midnightblue',
                 label='PRE')

rects2 = plt.bar(index + bar_width, post, bar_width,
                 alpha=opacity,
                 color='r',
                 label='POST')

plt.xlabel('Gramicidin Incorporated')
plt.ylabel('Relative Abundance')
#plt.title('Scores by person')
plt.xticks(index + bar_width/2, ('4','6','8','10'))
plt.legend()
plt.ylim(0,0.55)
plt.tight_layout()
plt.show()

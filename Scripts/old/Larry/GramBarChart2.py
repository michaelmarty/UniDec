import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib
import matplotlib



# data to plot
n_groups = 4
#means_DMPC = (.30, 0.55, 0.12, 0.03)
#means_DMPG = (0.57, 0.41, 0.02, 0.0)
path = "Z:\\Group Share\\Julia\\Imidazole analogues\\190818\\"
os.chdir(path)


files = "averagecontrol.txt"
files1 = "averageimidazole.txt"
files2 = "average2iso.txt"
means_control = np.loadtxt(files)
means_Imidazole = np.loadtxt(files1)
means_2Iso = np.loadtxt(files2)
yerrControl = np.loadtxt("controlstddev.txt")
yerrImidazole = np.loadtxt("imidazolestddev.txt")
yerr2Iso = np.loadtxt("2isostddev.txt")


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
pos = list(range(len('averagecontrol')))

#rects1 = plt.bar(index, means_control, bar_width, yerr=yerrControl, alpha=opacity, color="k")

#rects2 = plt.bar(index + bar_width, means_Imidazole, bar_width, yerr=yerrImidazole,
                # alpha=opacity,
                 #color='r',
                 #label='Imidazole')
#rects3 = plt.bar(index+bar_width, means_2Iso, bar_width, yerr=yerr2Iso,
                 #alpha=opacity,
                 #color='b',
                 #label="2-Isopropyl")
figs, ax = plt.subplots()
rects1 =plt.bar([pos], "average control", bar_width, label='control')
rects2 = plt.bar([p+ bar_width for p in pos], "averageimidazole", bar_width, label = 'imidazole')
rects3 = ax.bar([p+bar_width*2 for p in pos], "average2iso", bar_width, label= "2iso")




plt.xlabel('Gramicidin Incorporated')
plt.ylabel('Relative Abundance')
#plt.title('Scores by person')
plt.xticks(index + bar_width/2, ('0', '2', '4', '6'))
plt.legend()
plt.ylim(0,0.7)
plt.tight_layout()
plt.show()

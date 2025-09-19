import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\\Group Share\\Deseree\\New FPOP\\Indolicidin\\"
os.chdir(path)

file1 = "Indo 3 to 1 Average Ox.txt"


ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("Indo 3 to 1 St Dev.txt")


labels = ['DMPG Con', 'DMPG Ox', 'DMPC Con', 'DMPC Ox']#, '8', '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars



plt.bar(labels, ydata1)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.

plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=0.5)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,0.51, step=0.10), size=20)

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.show()

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\\Group Share\\Deseree\\New FPOP\\"
os.chdir(path)

file1 = "LL37 9to1 PG PC ecoli.txt"


ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("LL37 9to1 PG PC ecoli Errors.txt")


labels = ['DMPG Con', 'DMPG Ox', 'DMPC Con', 'DMPC Ox', 'e. Coli Con', 'e.Coli Ox']

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars



plt.bar(labels, ydata1, yerr=yerr1, capsize=3)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.

plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=0.1)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,0.11, step=0.05), size=20)

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.show()

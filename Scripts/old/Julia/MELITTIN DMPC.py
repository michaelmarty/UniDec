import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\Julia/3rdyearBCP\MELITTIN DMPC"
os.chdir(path)

file1 = "AVERAGE_3.TXT"
file2 = "AVERAGE_6.txt"
file3 = "AVERAGE_9.txt"
file4 = "average_12.txt"
file5 = "average_18.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("Error_3.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("error_6.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("error_9.txt")
ydata4 = np.loadtxt(file4)
yerr4 = np.loadtxt('error_12.txt')
ydata5 = np.loadtxt(file5)
yerr5 = np.loadtxt("error_18.txt")

labels = ['0', '1', '2', '3', '4'] #'7']'8', '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.15  # the width of the bars


fig, ax = plt.subplots()
rects1 = ax.bar(x, ydata1, width, color='lightblue', label='3:1', yerr=yerr1, capsize=3.0)
rects2 = ax.bar(x+0.15, ydata2, width, color='c', label='6:1', yerr=yerr2, capsize=3.0)
rects3 = ax.bar(x+0.3, ydata3, width, color='seagreen', label='9:1', yerr=yerr3, capsize=3.0)
rects4 = ax.bar(x+0.45, ydata4, width, color='royalblue', label='12:1', yerr=yerr4, capsize=3.0)
rects5 = ax.bar(x+0.6, ydata5, width, color='darkslateblue', label='18:1', yerr=yerr5, capsize=3.0)



rcParams['font.size']=28
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_ylabel('%', size=24)
ax.set_xlabel('Melittin Peptides Incorporated', size=24)
ax.set_yticklabels(np.arange(0, 101, step=(25)))
#ax.set_title('DMPG')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=20)
plt.xlim(-0.5,len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,1.01, step=0.5), size=20)
ax.get_yaxis().set_ticklabels([0, "%", 100], size=30)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout()

plt.show()
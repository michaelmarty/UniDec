import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "Z:\Larry Walker\Figures\Imidazole_Titration"
os.chdir(path)

file1 = "0mM.txt"
file2 = "6mM.txt"
file3 = "12mM.txt"
file4 = "24mM.txt"
file5 = "36mM.txt"
file6 = "48mM.txt"

ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("0mM_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("6mM_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("12mM_Errors.txt")
ydata4 = np.loadtxt(file4)
yerr4 = np.loadtxt("24mM_Errors.txt")
ydata5 = np.loadtxt(file5)
yerr5 = np.loadtxt("36mM_Errors.txt")
ydata6 = np.loadtxt(file6)
yerr6 = np.loadtxt("48mM_Errors.txt")


labels = ('0', '1', '2', '3', '4', '5', '6') #, '7', '8', '9', '10', '11', '12', '13', '14']


x = np.arange(len(labels))  # the label locations
width = 0.12  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-2.5*width, ydata1, width, color='lightskyblue', label='0 mM', yerr=yerr1, capsize=3.0)
rects2 = ax.bar(x-1.5*width, ydata2, width, color='c', label='6 mM', yerr=yerr2, capsize=3.0)
rects3 = ax.bar(x-0.5*width, ydata3, width, color='seagreen', label='12 mM', yerr=yerr3, capsize=3.0)
rects4 = ax.bar(x+0.5*width, ydata4, width, color='teal', label='24 mM', yerr=yerr4, capsize=3.0)
rects5 = ax.bar(x+1.5*width, ydata5, width, color='darkslateblue', label='36 mM', yerr=yerr5 ,capsize=3.0)
rects6 = ax.bar(x+2.5*width, ydata6, width, color='indigo', label='48 mM', yerr=yerr6, capsize=3.0)

rcParams['font.size']=28
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
# Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_ylabel('%', size=24)
ax.set_xlabel('Melittin Peptides Incorporated', size=30)
#ax.set_ylabels["0", "%", "100"]
#ax.set_yticklabels(np.arrange(0, 101, step=(50)))
#ax.set_title('Imidazole Titration', size=35)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(size=30)
plt.xlim(-0.5, len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0, 1.01, step=0.5), size=20)
#plt.yticklabels(0, "%", 100)
ax.get_yaxis().set_ticklabels([0, "%", 100], size=30)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#fig.tight_layout()

plt.show()
plt.savefig("plotter.pdf", transparent="True")
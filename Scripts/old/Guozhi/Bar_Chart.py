import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "P:\\Group Share\\Guozhi Zhang\\Lipid exchange new\\PYTHON GRAPH\\same msp conc titration"
os.chdir(path)

file1 = "Empty ND.txt"
file2 = "Total Lipid Ratio.txt"
file3 = "AmtB ND.txt"
#file4 = "DMPG_3to1_LL37.txt"
#file5 = "DMPG_3to1_Melittin.txt"
#file6 = "DMPG_3to1_Gramicidin.txt" # load file

ydata1 = np.loadtxt(file1)
#yerr1 = np.loadtxt("DMPG_3to1_Indolicidin_Errors.txt")
ydata2 = np.loadtxt(file2)
#yerr2 = np.loadtxt("DMPG_3to1_TemporinA_Errors.txt")
ydata3 = np.loadtxt(file3)
#yerr3 = np.loadtxt("DMPG_3to1_TemporinL_Errors.txt")
#ydata4 = np.loadtxt(file4)
#yerr4 = np.loadtxt("DMPG_3to1_LL37_Errors.txt")
#ydata5 = np.loadtxt(file5)
#yerr5 = np.loadtxt("DMPG_3to1_Melittin_Errors.txt")
#ydata6 = np.loadtxt(file6)
#yerr6 = np.loadtxt("DMPG_3to1_Gramicidin_Errors.txt") # load ave and stdev

labels = ['1', '2', '3', '4', '5', '6']#$, '7', '8', '9' , '10', '11', '12', '13', '14'], x labels

x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars, might need change depends on

fig, ax = plt.subplots()
rects1 = ax.bar(x+width, ydata1, width, color='gray', label='Empty ND')#, yerr=yerr1)
rects2 = ax.bar(x+2*width, ydata2, width, color='b', label='Total Lipid Ratio')#, yerr=yerr2)
rects3 = ax.bar(x+3*width, ydata3, width, color='y', label='AmtB ND')#, yerr=yerr3)
#rects4 = ax.bar(x+4*width, ydata4, width, color='r', label='LL-37', yerr=yerr4)
#rects5 = ax.bar(x+5*width, ydata5, width, color='b', label='Melittin', yerr=yerr5)
#rects6 = ax.bar(x+6*width, ydata6, width, color='gray', label='Gramicidin', yerr=yerr6)

rcParams['font.size']=16
#rcParams['legend.loc']='lower center'
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('% PC', size=20)
ax.set_xlabel('Titration Set', size=20)
ax.set_title('Crude Titration')
ax.set_xticks(x+2*width)
ax.set_xticklabels(labels)
ax.set_yticklabels(np.arange(0, 101, step=(25)))#change the scale of y
ax.legend()
plt.xticks(size=20)
plt.xlim(0,len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,1.1, step=0.25), size=20)
fig.tight_layout()


plt.show()
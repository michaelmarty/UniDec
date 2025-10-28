import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "F:\\Group Share\\Figures\\Peptide_Paper2\\All Peptides\\3to1\\DMPG\\"
os.chdir(path)

file1 = "DMPG_3to1_TemporinA.txt"
file2 = "DMPG_3to1_Bactenecin_DTT.txt"
file3 = "DMPG_3to1_Bactenecin.txt"
file4 = "DMPG_3to1_TemporinL.txt"
file5 = "DMPG_3to1_LactoferricinB_DTT.txt"
file6 = "DMPG_3to1_Magainin2.txt"
file7= "DMPG_3to1_KR12.txt"
file8 = "DMPG_3to1_LactoferricinB.txt"
file9 = "DMPG_3to1_Indolicidin.txt"

file10 = "DMPG_3to1_LL37.txt"
file11 = "DMPG_3to1_Melittin.txt"
file12 = "DMPG_3to1_Gramicidin.txt"

ydata9 = np.loadtxt(file9)
yerr9 = np.loadtxt("DMPG_3to1_Indolicidin_Errors.txt")
ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DMPG_3to1_TemporinA_Errors.txt")
ydata4 = np.loadtxt(file4)
yerr4 = np.loadtxt("DMPG_3to1_TemporinL_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("DMPG_3to1_Bactenecin_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DMPG_3to1_Bactenecin_DTT_Errors.txt")
ydata7 = np.loadtxt(file7)
yerr7 = np.loadtxt("DMPG_3to1_KR12_Errors.txt")
ydata8 = np.loadtxt(file8)
yerr8 = np.loadtxt("DMPG_3to1_LactoferricinB_Errors.txt")
ydata5 = np.loadtxt(file5)
yerr5 = np.loadtxt("DMPG_3to1_LactoferricinB_DTT_Errors.txt")
ydata6 = np.loadtxt(file6)
yerr6 = np.loadtxt("DMPG_3to1_Magainin2_Errors.txt")
ydata10 = np.loadtxt(file10)
yerr10 = np.loadtxt("DMPG_3to1_LL37_Errors.txt")
ydata11 = np.loadtxt(file11)
yerr11 = np.loadtxt("DMPG_3to1_Melittin_Errors.txt")
ydata12 = np.loadtxt(file12)
yerr12 = np.loadtxt("DMPG_3to1_Gramicidin_Errors.txt")

labels = ['0', '1', '2', '3', '4', '5', '6']#, '7', '8', '9', '10', '11', '12', '13', '14']

x = np.arange(len(labels))  # the label locations
width = 0.09  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x-(5.5*width), ydata1, width, label='TempA', yerr=yerr1)
rects2 = ax.bar(x-(4.5*width), ydata2, width, label='Bactenecin (Red)', yerr=yerr2)
rects3 = ax.bar(x-(3.5*width), ydata3, width, label='Bactenecin (Ox)', yerr=yerr3)
rects4 = ax.bar(x-(2.5*width), ydata4, width, label='TempL', yerr=yerr4)
rects5 = ax.bar(x-(1.5*width), ydata5, width, label='LacB (Red)', yerr=yerr5)
rects6 = ax.bar(x-(0.5*width), ydata6, width, label='Mag2', yerr=yerr6)
rects7 = ax.bar(x+(0.5*width), ydata7, width, label='KR12', yerr=yerr7)
rects8 = ax.bar(x+(1.5*width), ydata8, width, label='LacB (Ox)', yerr=yerr8)
rects9 = ax.bar(x+(2.5*width), ydata9, width, label='Indolicidin', yerr=yerr9)
#rects10 = ax.bar(x+(3.5*width), ydata10, width, label='LL37', yerr=yerr10)
#rects11 = ax.bar(x+(4.5*width), ydata11, width, label='Melittin', yerr=yerr11)
#rects12 = ax.bar(x+(5.5*width), ydata12, width, label='Gram', yerr=yerr12)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
#rcParams['legend.loc']='lower center'
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('%', size=24)
ax.set_xlabel('Number of Peptides Incorporated', size=24)
ax.set_title('3:1 Peptide:DMPG Nanodisc')
ax.set_xticks(x-2*width)
ax.set_xticklabels(labels)
ax.set_yticklabels(np.arange(0, 101, step=(25)))
ax.legend()
plt.xticks(size=20)
#plt.xlim(0,len(labels))
plt.ylim(top=1.0)
plt.ylim(bottom=0.0)
plt.yticks(np.arange(0,1.1, step=0.25), size=20)
fig.tight_layout()

plt.show()
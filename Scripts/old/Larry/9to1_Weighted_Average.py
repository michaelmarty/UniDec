import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

path = "F:\\Group Share\\Figures\\Peptide_Paper2\\All Peptides\\9to1\\"
os.chdir(path)

file1 = "DMPG_9to1_TemporinA.txt"
file2 = "DMPG_9to1_Bactenecin_DTT.txt"
file3 = "DMPG_9to1_Bactenecin.txt"
file4 = "DMPG_9to1_TemporinL.txt"
file5 = "DMPG_9to1_LactoferricinB_DTT.txt"
file6 = "DMPG_9to1_Magainin2.txt"
file7 = "DMPG_9to1_KR12.txt"
file8 = "DMPG_9to1_LactoferricinB.txt"
file9 = "DMPG_9to1_Indolicidin.txt"
file10 = "DMPG_9to1_LL37.txt"
file11 = "DMPG_9to1_Melittin.txt"
file12 = "DMPG_9to1_Gramicidin.txt"

Pepnames=()
Pepnames=np.append(Pepnames, ['TemporinA', 'Bactenecin(Red)', 'Bactenecin(Ox)', 'LL37',  'Gramicidin', 'TemporinL',  'Melittin', 'KR12', 'Mag2', 'LacB(Ox)', 'Indolicidin'])

ydata9 = np.loadtxt(file9)
yerr9 = np.loadtxt("DMPG_9to1_Indolicidin_Errors.txt")
ydata1 = np.loadtxt(file1)
yerr1 = np.loadtxt("DMPG_9to1_TemporinA_Errors.txt")
ydata4 = np.loadtxt(file4)
yerr4 = np.loadtxt("DMPG_9to1_TemporinL_Errors.txt")
ydata3 = np.loadtxt(file3)
yerr3 = np.loadtxt("DMPG_9to1_Bactenecin_Errors.txt")
ydata2 = np.loadtxt(file2)
yerr2 = np.loadtxt("DMPG_9to1_Bactenecin_DTT_Errors.txt")
ydata7 = np.loadtxt(file7)
yerr7 = np.loadtxt("DMPG_9to1_KR12_Errors.txt")
ydata8 = np.loadtxt(file8)
yerr8 = np.loadtxt("DMPG_9to1_LactoferricinB_Errors.txt")
ydata5 = np.loadtxt(file5)
yerr5 = np.loadtxt("DMPG_9to1_LactoferricinB_DTT_Errors.txt")
ydata6 = np.loadtxt(file6)
yerr6 = np.loadtxt("DMPG_9to1_Magainin2_Errors.txt")
ydata10 = np.loadtxt(file10)
yerr10 = np.loadtxt("DMPG_9to1_LL37_Errors.txt")
ydata11 = np.loadtxt(file11)
yerr11 = np.loadtxt("DMPG_9to1_Melittin_Errors.txt")
ydata12 = np.loadtxt(file12)
yerr12 = np.loadtxt("DMPG_9to1_Gramicidin_Errors.txt")

Pepnum=(0,1,2,3,4,5,6,7,8,9)

ydata01=np.average(Pepnum, weights=ydata1)
ydata02=np.average(Pepnum, weights=ydata2)
ydata03=np.average(Pepnum, weights=ydata3)
ydata04=np.average(Pepnum, weights=ydata4)
ydata05=np.average(Pepnum, weights=ydata5)
ydata06=np.average(Pepnum, weights=ydata6)
ydata07=np.average(Pepnum, weights=ydata7)
ydata08=np.average(Pepnum, weights=ydata8)
ydata09=np.average(Pepnum, weights=ydata9)
ydata010=np.average(Pepnum, weights=ydata10)
ydata011=np.average(Pepnum, weights=ydata11)
ydata012=np.average(Pepnum, weights=ydata12)

yerr01 = np.sqrt((np.sum(ydata1*(Pepnum - ydata01)**2))/(((len(ydata1)-1)*np.sum(ydata1))/len(ydata1)))
yerr02 = np.sqrt((np.sum(ydata2*(Pepnum - ydata02)**2))/(((len(ydata2)-1)*np.sum(ydata2))/len(ydata2)))
yerr03 = np.sqrt((np.sum(ydata3*(Pepnum - ydata03)**2))/(((len(ydata3)-1)*np.sum(ydata3))/len(ydata3)))
yerr04 = np.sqrt((np.sum(ydata4*(Pepnum - ydata04)**2))/(((len(ydata4)-1)*np.sum(ydata4))/len(ydata4)))
yerr05 = np.sqrt((np.sum(ydata5*(Pepnum - ydata05)**2))/(((len(ydata5)-1)*np.sum(ydata5))/len(ydata5)))
yerr06 = np.sqrt((np.sum(ydata6*(Pepnum - ydata06)**2))/(((len(ydata6)-1)*np.sum(ydata6))/len(ydata6)))
yerr07 = np.sqrt((np.sum(ydata7*(Pepnum - ydata07)**2))/(((len(ydata7)-1)*np.sum(ydata7))/len(ydata7)))
yerr08 = np.sqrt((np.sum(ydata8*(Pepnum - ydata08)**2))/(((len(ydata8)-1)*np.sum(ydata8))/len(ydata8)))
yerr09 = np.sqrt((np.sum(ydata9*(Pepnum - ydata09)**2))/(((len(ydata9)-1)*np.sum(ydata9))/len(ydata9)))
yerr010 = np.sqrt((np.sum(ydata10*(Pepnum - ydata010)**2))/(((len(ydata10)-1)*np.sum(ydata10))/len(ydata10)))
yerr011 = np.sqrt((np.sum(ydata11*(Pepnum - ydata011)**2))/(((len(ydata11)-1)*np.sum(ydata11))/len(ydata11)))
yerr012 = np.sqrt((np.sum(ydata12*(Pepnum - ydata012)**2))/(((len(ydata12)-1)*np.sum(ydata12))/len(ydata12)))


mylist=()
mylist=np.append(mylist, [ydata01, ydata02, ydata03, ydata010,  yerr012, ydata04, ydata011, ydata07,  ydata06, ydata08, ydata09])
Errors=()
Errors=np.append(Errors, [yerr01, yerr02, yerr03, yerr010, yerr012, yerr04, yerr011, yerr07,  yerr06, yerr08, yerr09])

fig, ax = plt.subplots()
rects = ax.bar(Pepnames,mylist, yerr=Errors)


rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Average Peptides Incorporated', size=24)
#ax.set_xlabel('Number of Peptides Incorporated', size=24)
ax.set_title('9:1 Peptide:DMPG Nanodisc')
ax.set_xticks(Pepnames)
#ax.set_xticklabels(labels)
#ax.set_yticklabels(np.arange(0, 101, step=(25)))
#ax.legend()
plt.xticks(size=20)
#ax.set_xticklabels(Pepnames, rotation=45)
#plt.xlim(0,len(labels))
plt.ylim(top=6.5)
plt.ylim(bottom=0.0)
#plt.yticks(np.arange(0,1.1, step=0.25), size=20)
fig.tight_layout()

print(mylist)
#print(Errors)

plt.show()
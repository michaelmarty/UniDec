import numpy as np
import numpy
import os
import matplotlib.pyplot as plt

#set the directory
dir = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\MAMMALIAN\\POPC"
os.chdir(dir)
#load csv file with numpy
POPC = numpy.loadtxt(fname="POPC all.csv", delimiter = ",")
POPC_pos = POPC[0:3]
#print(POPC_SOPE_pos)
POPC_neg = POPC[3:6]
#print(POPC_SOPE_neg)
POPC_imid = POPC[6:9]
#print(POPC_SOPE_imid)
POPC_methyl = POPC[9:12]
#print(POPC_SOPE_methyl)
#x4 is an input for what to print
#calculate average and standard deviations for input x4
x4 = POPC_pos
print (x4.mean(axis=0))
print (x4.std(axis=0))

file1 = "x values for figure 2.txt"
methyl = x4.mean(axis=0)
methyl_deviation = x4.std(axis=0)
#file4 ="POPC SOPE POPS 2mimid y values.txt"
#file5 ="POPC SOPE POPS 2mimid error.txt"
#file6 ="POPC SOPE POPS SSM 2mimid y values.txt"
#file7 ="POPC SOPE POPS SSM 2mimid y errors.txt"
file8 = "POPC line.txt"
file9 = "POPC SOPE line.txt"
file10 = "POPC SOPE POPS line.txt"
file11 = "POPC SOPE POPS SSM line.txt"

t = np.loadtxt(file1)
y1 = methyl
e = methyl_deviation
#y2 = np.loadtxt(file4)
#e2 = np.loadtxt(file5)
#y3 = np.loadtxt(file6)
#e3 = np.loadtxt(file7)
y4 =np.loadtxt(file8)
y5 =np.loadtxt(file9)
y6 =np.loadtxt(file10)
y7 =np.loadtxt(file11)


#rcParams['font.size']=10
fig, ax = plt.subplots()


#f2=interp1d(t,y1,kind='cubic')
#f2=interp1d(t,y2,kind='cubic')
xnew = np.linspace(1, 36,num=360,endpoint=True)
ax.set(ylim=(750,762))
ax.set(xlim=(0,200))
ax.plot(t, y1)
#ax.plot(t, y2)
#ax.plot(t, y3)
ax.plot(t,y4,"--")
ax.plot(t,y5,"--")
ax.plot(t,y6,"--")
ax.plot(t,y7,"--")
plt.plot(t,y1, 'o')
plt.errorbar(t,y1,yerr=e,fmt='none', capsize=3, capthick=1, elinewidth=1)
#plt.plot(t,y2, 'o')
#plt.errorbar(t,y2,yerr=e2,fmt='none', capsize=3, capthick=1, elinewidth=1)
#plt.plot(t,y3, 'o')
#!plt.errorbar(t,y3,yerr=e3,fmt='none', capsize=3, capthick=1, elinewidth=1)

ax.set_ylabel('Average Lipid Mass', size=16)
ax.set_xlabel('Collision Voltage', size=16)
ax.set_title('2-methylimidazole', size=20)
#ax.grid()
#plt.show()
#fig.savefig("test.png")
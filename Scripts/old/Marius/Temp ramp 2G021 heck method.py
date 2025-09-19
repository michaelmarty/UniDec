import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

#set the directory
dir = "Z:\\Group Share\\Marius Kostelic\\Baker Lab AAVs"
os.chdir(dir)
#load csv file with numpy
#POPC_SOPE = numpy.loadtxt(fname="POPC SOPE all.csv", delimiter = ",")#
#POPC_SOPE_pos = POPC_SOPE [0:3]
#print(POPC_SOPE_pos)
#POPC_SOPE_neg = POPC_SOPE[3:6]
#print(POPC_SOPE_neg)
#POPC_SOPE_imid = POPC_SOPE[6:9]
#print(POPC_SOPE_imid)
#POPC_SOPE_methyl = POPC_SOPE[9:12]
#print(POPC_SOPE_methyl)
#x4 is an input for what to print
#calculate average and standard deviations for input x4
#x4 = POPC_SOPE_pos
#print (x4.mean(axis=0))
#print (x4.std(axis=0))

file1 = "2G021 x values.txt"
#methyl = x4.mean(axis=0)
#methyl_deviation = x4.std(axis=0)
file2 = "2G021 empty y values.txt"
file3 = "2G021 empty y error values.txt"
file4 ="2G021 partial y values.txt"
file5 ="2G021 partial y error values.txt"
file6 ="2G021 full y values.txt"
file7 ="2G021 full y error values.txt"
#file8 = "POPC line.txt"
#file9 = "POPC SOPE line.txt"
#file10 = "POPC SOPE POPS line.txt"
#file11 = "POPC SOPE POPS SSM line.txt"

t = np.loadtxt(file1)
y1 = np.loadtxt(file2)
e = np.loadtxt(file3)
y2 = np.loadtxt(file4)
e2 = np.loadtxt(file5)
y3 = np.loadtxt(file6)
e3 = np.loadtxt(file7)
#y4 =np.loadtxt(file8)
#5 =np.loadtxt(file9)
#y6 =np.loadtxt(file10)
#y7 =np.loadtxt(file11)

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

#rcParams['font.size']=10
fig, ax = plt.subplots()

#f2=interp1d(t,y1,kind='cubic')
#f2=interp1d(t,y2,kind='cubic')
xnew = np.linspace(1, 36,num=360,endpoint=True)
ax.set(ylim=(0,1))
ax.set(xlim=(10,80))
ax.legend()
ax.plot(t, y1, color='darkblue')
ax.plot(t, y2, color='grey')
ax.plot(t, y3, color='red')
#ax.plot(t,y4,"--"
#ax.plot(t,y5,"--")
#ax.plot(t,y6,"--")
#ax.plot(t,y7,"--")
plt.plot(t,y1, 'o',color='darkblue', markersize = 8)
plt.errorbar(t,y1,yerr=e,fmt='none',color='darkblue', capsize=3, capthick=1, elinewidth=1)
plt.plot(t,y2, 'D', color='grey', markersize = 8)
plt.errorbar(t,y2,yerr=e2,fmt='none',color='grey', capsize=3, capthick=1, elinewidth=1)
plt.plot(t,y3, 's', color='red', markersize = 8)
plt.errorbar(t,y3,yerr=e3,fmt='none', color='red', capsize=3, capthick=1, elinewidth=1)



ax.set_ylabel('Normalized Intensity', size=16)
ax.set_xlabel('Temperature', size=16)
ax.set_title('Empty, Partial, and Fully filled capsids', size=20)
fig.tight_layout()
#ax.grid()
plt.show()
#fig.savefig("test.png")
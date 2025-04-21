import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

directory = "Z:\\Group Share\\Annika\\LXMS, Fall 2024\\December Exchange"

os.chdir(directory)

file1 = "cholesterol1(-).txt"

data1 = np.loadtxt(file1)

file2 = "popc1(-).txt"

data2 = np.loadtxt(file2)

#file3 = "40uM_ToFMRM.txt"

#data3 = np.loadtxt(file3)

i1 = ud.integrate(data1, 5.3, 6.3)[0]
i2 = ud.integrate(data2, 5.3, 6.3)[0]
#i3 = ud.integrate(data3, 5.3, 6.3)[0]

#print(i1)
#print(i1/(np.sum([i1, i2])))
print(i1/i2)


#plt.subplot(211)
plt.plot(data1[:, 0], data1[:, 1], label="Cholesterol" , color = "royalblue", linewidth=3)
plt.plot(data2[:,0], data2[:,1], label="POPC", color="darkgreen", linewidth=3)
#plt.plot(data3[:,0], data3[:,1], label="40 uM Cholesterol ToF MRM", color="blue")
plt.legend()
plt.xlim(6,12.0)
#plt.ylim(0,6000000)
plt.xlabel("Retention Time (min)")
plt.ylabel("Intensity")

#plt.subplot(211)
#plt.plot(data3[:,0], data3[:,1], color="k")
#plt.xlabel("m/z")
#plt.ylabel("Intensity")

plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.tight_layout()
plt.savefig("2024_12_19_LX1(-).pdf")

plt.show()

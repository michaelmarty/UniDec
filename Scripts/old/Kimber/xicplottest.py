import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

directory = "Z:\\Group Share\\Kimber\\OzID - Current"

os.chdir(directory)

file1 = "spectest.txt"

data1 = np.loadtxt(file1)

#file2 = "spectestplz.txt"

#data2 = np.loadtxt(file2)

#file3 = "40uM_ToFMRM.txt"

#data3 = np.loadtxt(file3)

i1 = ud.integrate(data1, 5.3, 6.3)[0]
#i2 = ud.integrate(data2, 5.3, 6.3)[0]
#i3 = ud.integrate(data3, 5.3, 6.3)[0]

print(i1)
#print(i1/(np.sum([i1, i2])))
#print(i1/i2)


plt.subplot(212)
plt.plot(data1[:, 0], data1[:, 1], label="OzID of PE 16:1_18:1")
#plt.plot(data2[:,0], data2[:,1], label="40 uM Cholesterol MSMS", color="red")
#plt.plot(data3[:,0], data3[:,1], label="40 uM Cholesterol ToF MRM", color="blue")
plt.legend()
#plt.xlim(0,12.0)
#plt.ylim(0,6000000)
plt.xlabel("Retention Time (min)")
plt.ylabel("%")

#plt.subplot(211)
#plt.plot(data3[:,0], data3[:,1], color="k")
#plt.xlabel("m/z")
#plt.ylabel("Intensity")

#plt.tight_layout()
#plt.savefig("Ecoli_1D_HILIC.pdf")

plt.show()

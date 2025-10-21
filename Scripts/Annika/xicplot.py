import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

directory = "C:\\Users\\as244475\\OneDrive - The University of Texas at Austin\\Desktop"

os.chdir(directory)

file1 = "5min_DMPC_678.txt"  # Change this for each TIC

file2 = "5min_PE_16-0_18-1.txt"

data1 = np.loadtxt(file1)

data2 = np.loadtxt(file2)
data1[:, 1] = data1[:, 1] / np.max(data1[:, 1])
data2[:, 1] = data2[:, 1] / np.max(data2[:, 1])

# --- Check ---
print("Max height data1 =", np.max(data1[:, 1]))
print("Max height data2 =", np.max(data2[:, 1]))

#file3 = "40uM_ToFMRM.txt"

#data3 = np.loadtxt(file3)

#i1 = ud.integrate(data1, 5.3, 6.3)[0]
#i2 = ud.integrate(data2, 5.3, 6.3)[0]
#i3 = ud.integrate(data3, 5.3, 6.3)[0]

#print(i1)
#print(i1/(np.sum([i1, i2])))
#print(i1/i2)


#plt.subplot(211)
plt.plot(data1[:, 0], data1[:, 1], label="DMPC" , color = "royalblue", linewidth=3)
plt.plot(data2[:,0], data2[:,1], label="POPG", color="darkgreen", linewidth=3)
#plt.plot(data3[:,0], data3[:,1], label="40 uM Cholesterol ToF MRM", color="blue")
plt.legend()
plt.xlim(0,5.0)
#plt.ylim(0,6000000)
plt.xlabel("Retention Time (min)")
plt.ylabel("Intensity")

#plt.subplot(211)
#plt.plot(data3[:,0], data3[:,1], color="k")
#plt.xlabel("m/z")
#plt.ylabel("Intensity")

plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.tight_layout()
plt.savefig("Overlaid_C8_runs2.pdf")

plt.show()

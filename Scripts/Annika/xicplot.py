import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter

directory = "Z:\\Group Share\\Annika\\Stellar"

os.chdir(directory)

file1 = "AqpZ_Membrane_Extract.txt"
file2 = "HEK_Acidified_NegTrace.txt"  # Change this for each TIC

data1 = np.loadtxt(file1)

# data2 = np.loadtxt(file2)
#data1[:, 1] = data1[:, 1] / np.max(data1[:, 1])
#data2[:, 1] = data2[:, 1] / np.max(data2[:, 1])

# --- Check ---
#print("Max height data1 =", np.max(data1[:, 1]))
#print("Max height data2 =", np.max(data2[:, 1]))

#file3 = "40uM_ToFMRM.txt"

#data3 = np.loadtxt(file3)

#i1 = ud.integrate(data1, 5.3, 6.3)[0]
#i2 = ud.integrate(data2, 5.3, 6.3)[0]
#i3 = ud.integrate(data3, 5.3, 6.3)[0]

#print(i1)
#print(i1/(np.sum([i1, i2])))
#print(i1/i2)

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14
#plt.subplot(211)
fig, ax = plt.subplots(figsize=(11, 5))
plt.plot(data1[:, 0], data1[:, 1], label="ESI+" , color = "mediumvioletred", linewidth=2)
plt.plot(data2[:,0], data2[:,1], label="ESI-", color="indigo", linewidth=2)
#plt.plot(data3[:,0], data3[:,1], label="40 uM Cholesterol ToF MRM", color="blue")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_linewidth(2)  # Set bottom spine thickness
plt.gca().spines['left'].set_linewidth(2)  # Set left spine thickness
plt.tick_params(axis='both', width=2, length=5)  # Tick mark thickness

#plt.yticks(np.arange(0.0, 50001.0, 10000.0), fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend()
plt.xlim(0,15)
plt.ylim(0, None)
plt.xlabel("Retention Time (min)", fontsize=18)
plt.ylabel("Intensity", fontsize=18)

#plt.subplot(211)
#plt.plot(data3[:,0], data3[:,1], color="k")
# plt.xlabel("Retention Time (min)")
# plt.ylabel("Intensity")

plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.tight_layout()
#plt.savefig("Chol_D7Chol_XICs_Stellar.pdf")

plt.show()

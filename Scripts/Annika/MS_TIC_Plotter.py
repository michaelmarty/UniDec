import os
import numpy as np
import matplotlib.pyplot as plt

# Define the working directory
directory = "Z:\\Group Share\\Annika\\Stellar"
os.chdir(directory)

# Manually enter the TIC file name each time you run the script
file1 = "NEWEcoliChrom.txt"  # Change this for each TIC

#file2 = "(+)Ecoli_t-1_withFE.txt"
#file3 = "popc-postexch.txt"

#file4 = "D7pc-postexch.txt"

# Load the data from the text file
data = np.loadtxt(file1)
#data2 = np.loadtxt(file2)
#data3 = np.loadtxt(file3)
#data4 = np.loadtxt(file4)

# Assuming the first column is time and the second column is TIC intensity
time = data[:, 0]  # Retention time
intensity = data[:, 1]  # TIC intensity
# Plot the TIC
plt.figure(figsize=(8, 5))
plt.plot(data[:,0], data[:,1], label= "No Folch", color = "k", linewidth=2)
#plt.plot(data2[:,0], data2[:,1], label="With Folch", color="red", linewidth=2)
#plt.plot(data3[:,0], data3[:,1], label="POPC", color="darkgreen", linewidth=2)
#plt.plot(data4[:,0], data4[:,1], label="D7 PC", color="lime", linewidth=3)
plt.xlim(5, 12)
plt.ylim(bottom=0)
plt.ylim(0, 80000000)
plt.xlabel("Retention Time (min)")
plt.ylabel("Intensity")
#plt.title("E.col)")
plt.legend()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig("(-)EPL_Ctrl_Chrom.pdf")
plt.show()


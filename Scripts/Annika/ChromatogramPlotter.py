
import os
import numpy as np
import matplotlib.pyplot as plt



# Load the file (assumes comma-separated and correct headers)
directory = "Z:\\Group Share\\Annika\\Stellar\\Cholesterol"
os.chdir(directory)

file1 = "D7Chol_TIC.txt"  # Change this for each TIC

file2 = "Chol_TIC.txt"
#file3 = "popc-postexch.txt"

#file4 = "D7pc-postexch.txt"

# Load the data from the text file
data = np.genfromtxt(file1, delimiter='\t')

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(data[:,0], data[:,1], color = "k", linewidth=3)
#plt.xlim(5,8)
plt.xlim(2,9)
plt.ylim(0, 20000)
# You choose axis labels here
plt.xlabel("Retention Time (min)", fontsize=24)
plt.ylabel("Intensity", fontsize=24)
#plt.title("(+)Ecoli without Folch Extraction", fontsize=28)
#plt.legend()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.tight_layout()
#plt.savefig("Z:\\Group Share\\Annika\\Sample Prep Troubleshooting\\(+)Ecoli_t-1_noFE.pdf", bbox_inches='tight')


plt.show()

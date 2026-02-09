import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# Define the working directory
directory = "Z:\\Group Share\\Annika\\Figures"
os.chdir(directory)

# TIC files
file1 = "HEK_pos_TIC.txt"     #pos
file2 = "HEK_neg_TIC.txt"    #neg

# Load data
data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14


# Plot both on the same figure
plt.figure(figsize=(8, 5))
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_linewidth(2)  # Set bottom spine thickness
plt.gca().spines['left'].set_linewidth(2)  # Set left spine thickness
plt.tick_params(axis='both', width=2, length=5)  # Tick mark thickness


plt.plot(data1[:, 0], data1[:, 1],
         label="Positive", color="royalblue", linewidth=2)

plt.plot(data2[:, 0], data2[:, 1],
         label="Negative", color="fuchsia", linewidth=2)

# Formatting
plt.xlim(5, 12)
# plt.ylim(0, 8e7)
plt.xlabel("Retention Time (min)")
plt.ylabel("Intensity")
plt.legend()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

#plt.savefig("Ecoli_Overlay_Chrom.pdf")
plt.show()
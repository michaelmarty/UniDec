import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configure matplotlib settings
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.sans-serif'] = "Arial"

# Load data from the first text file
data1 = np.loadtxt("Z:\\Group Share\\Seun\\LX-MS\\Fluorometry data\\20241029 data\\mcherry emmision corrected edited.txt")
plt.plot(data1[:, 0], data1[:, 1], label="mCherry emission", color="tomato")

# Load data from the second text file
data2 = np.loadtxt("Z:\\Group Share\\Seun\\LX-MS\\Fluorometry data\\20241029 data\\mcherry blank corrected excitation edited.txt")  # Replace with the actual path to the second file
plt.plot(data2[:, 0], data2[:, 1], label="mCherry excitation", color="black")  # Change color as needed

# Set plot limits and labels
plt.gca().set_xlim(220, 800)
plt.ylim(0.4, )
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.title('Emission Spectrum')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure if needed
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20241029_mCherry Ex Em.png", dpi=1000)

plt.show()

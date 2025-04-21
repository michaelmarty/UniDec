from LCImporter import *
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.sans-serif'] = "Arial"


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-13 13-56-56\\20240613_Blank PMT 6.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-13 12-34-21\\20240613_Blank1_ 5 ul injection.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-13 14-31-42\\20240613_Blank PMT 12.D"

d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)
d3 = extract_data(p3, index = 0)

#plt.subplot(3, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="PMT 6", color="red")
plt.ylim(0, 50)
plt.xlim(0, 7)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.legend()

#plt.subplot(3, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="PMT 10", color="grey")
plt.ylim(0, 50)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.legend()

#plt.subplot(3, 1, 3)
plt.plot(d3[:,0], d3[:,1], label="PMT 12", color="black", linestyle="--")
plt.ylim(-5, 50)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240617_PMT optimization.png", dpi=1000)


plt.show()
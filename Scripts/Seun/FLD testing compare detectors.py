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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-12 14-05-11\\20240612_sfGFP_QCC_5 ul injection.D"



d1 = extract_data(p1, index = 3)
d2 = extract_data(p1, index = 0)





plt.subplot(2, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="MWD sfGFP", color="green")
plt.ylim(-50, 3000)
plt.xlim(0, 7)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="FLD sfGFP", color="green")
plt.ylim(-50, 600)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.legend()


# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240612_sfGFP detectors compared without column.png", dpi=1000)


plt.show()
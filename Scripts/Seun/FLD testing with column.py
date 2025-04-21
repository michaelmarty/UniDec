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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-12 16-42-18\\20240612_Blank1_ 5 ul injection.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-06-12 16-42-18\\20240612_sfGFP_QCC_5 ul injection.D"



d1 = extract_data(p1)
d2 = extract_data(p2)

#plt.subplot(3, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="Blank before", color="black")


#plt.subplot(3, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="sfGFP", color="green")


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
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240621_FLD sf GFP with column.png", dpi=1000)


plt.show()
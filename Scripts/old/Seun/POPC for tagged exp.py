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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-29 12-19-35\\20240529_(+)POPC_rep1_5 ul injection.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-29 12-19-35\\20240529_(+)POPC_rep2_5 ul injection.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-29 12-19-35\\20240529_(+)POPC_rep3_5 ul injection.D"

d1 = extract_data(p1, index = 4)
d2 = extract_data(p2, index = 4)
d3 = extract_data(p3, index = 4)


#plt.subplot(3, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="Rep 3", color="red")

#plt.subplot(3, 1, 2)
plt.plot(d2[:,0], d2[:,1] - 100, label="Rep 2", color="red" , linestyle="--")


#plt.subplot(3, 1, 3)
plt.plot(d3[:,0], d3[:,1] - 200, label="Rep 1", color="maroon")
plt.ylim(-250, 100)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('520 nm POPC All data stacked')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240529_520 nm POPC for tagged exp triplicates.png", dpi=1000)


plt.show()
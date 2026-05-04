from LCImporter import *
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 18
mpl.rcParams['font.sans-serif'] = "Arial"


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\Seun\\OF_20240309_control test 2024-05-15 15-10-48\\20240515_sfGFP 7 min M 5 ul injection.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-20 14-16-04\\20240520_sfGFP_QCC_5 ul injection.D"




d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)


#plt.subplot(2, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="Before", color="green")

#plt.subplot(2, 1, 2)
plt.plot(d2[:,0], d2[:,1] - 100, label="After", color="green" , linestyle="--")



# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('280 nm sf GFP All data stacked')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240520_280 nm sf GFP duplicates.png", dpi=1000)
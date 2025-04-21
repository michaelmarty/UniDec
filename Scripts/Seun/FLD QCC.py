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



p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2025-01-14 08-59-54\\20250114_Blank.D"
p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2025-01-14 08-59-54\\20250114_(+)sf GFP.D"

d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)


plt.plot(d1[:,0], d1[:,1], label="blank", color="black")




plt.plot(d2[:,0], d2[:,1], label="sf GFP", color="green")






plt.ylim(0.19, 5)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FLD (LU)')
plt.title('FLD QCC 2025-01-14 before run two columns stacked')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20250114_QCC b4 run 2 columns stacked.png", dpi=1000)


plt.show()
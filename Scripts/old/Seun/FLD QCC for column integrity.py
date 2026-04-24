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



p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-06-25 15-31-42\\20240625_Blank1.D"
p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-06-25 15-31-42\\20240625_diluted sf GFP.D"
p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-06-25 15-31-42\\20240625_Blank2.D"
p4 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-06-25 15-31-42\\20240625_sf GFP.D"

d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)
d3 = extract_data(p3, index = 0)
d4 = extract_data(p4, index = 0)


plt.plot(d1[:,0], d1[:,1], label="blank1", color="black")
plt.plot(d2[:,0], d2[:,1], label="diluted sf GFP", color="green")
plt.plot(d3[:,0], d3[:,1], label="blank2", color="black", linestyle="--")
plt.plot(d4[:,0], d4[:,1], label="sf GFP", color="green", linestyle="--")






plt.ylim(0, 16)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.title('FLD QCC')
plt.legend ()

# Adjust layout
plt.legend (loc="upper right")

# Save figure
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240626_QCC.png", dpi=1000)


plt.show()
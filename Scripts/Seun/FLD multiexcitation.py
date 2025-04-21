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


#p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-10-08 12-28-26\\20241008_Blank4 stack 2C.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240612_FLD testing 2024-12-12 15-18-10\\20241212_(+)D1_mCherry 16 uM.D"

#p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-10-08 12-28-26\\20241008_(-)sf GFP 500_510.D"

#p4 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240625_QCC 2024-10-08 12-28-26\\20241008_(+)mRuby3+(-)sf GFP 500_510.D"

#d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 2)
#d3 = extract_data(p3, index = 0)
#4 = extract_data(p4, index = 0)

#plt.plot(d1[:,0], d1[:,1], label="blank", color="black")



plt.plot(d2[:,0], d2[:,1], label="(+)MSP1D1_mCherry", color="indianred")




#plt.plot(d3[:,0], d3[:,1], label="(-)sf GFP", color="green")


#plt.plot(d4[:,0], d4[:,1], label="sf GFP + mRuby 3", color="orchid")


plt.ylim(0.14, 0.8)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.title('FLD with Ex at 490 nm and Em at 700 nm')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240106_FLD MSP1D1_mCherry490_700.png", dpi=1000)


plt.show()
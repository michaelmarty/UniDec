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


#p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20250107_Estevan LX 2025-01-07 11-43-06\\20250107_Blank.D"
#p2 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20241115_LX-MS IMAC 2024-11-15 18-18-30\\20241115_(+)POPC 7 µM.D"
#p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20241115_LX-MS IMAC 2024-12-05 16-23-30\\20241205_(-)DMPC 7 µM.D"
p4 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20250107_Estevan LX 2025-01-14 10-42-41\\20250114_Estevan LX 17 min.D"

#d1 = extract_data(p1, index = 0)
#d2 = extract_data(p2, index = 0)
#d3 = extract_data(p3, index = 0)
d4  = extract_data(p4, index = 0)
#plt.subplot(3, 1, 1)
#plt.plot(d1[:,0], d1[:,1], label="blank", color="black")

#plt.subplot(3, 1, 2)
#plt.plot(d2[:,0], d2[:,1], label="(+)POPC", color="red")

#plt.subplot(3, 1, 3)
#plt.plot(d3[:,0], d3[:,1], label="(-)DMPC", color="blue")

plt.plot(d4[:,0], d4[:,1], label="Estevan LX", color="green")

plt.ylim(0.2, 12)
plt.xlim(0, 32)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.title('LX-MS IMAC')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20250114_ Estevan LX IMAC 32 min.png", dpi=1000)


plt.show()
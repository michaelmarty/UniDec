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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-05 13-02-20\\20240405_sf GFP var_flwrate 5 ul.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-09 16-02-07\\20240409_sf GFP binding check 5 ul.D"



d1 = extract_data(p1, index = 2)
d2 = extract_data(p2, index = 2)


plt.subplot(2, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="sf GFP", color="green")
plt.ylim(-50, 300)
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')

plt.subplot(2, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="sf GFP", color="green")
plt.ylim(-50, 300)
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')

plt.tight_layout()

#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240410_490 nm 40 min method binding check compared.png", dpi=1000)


plt.show()
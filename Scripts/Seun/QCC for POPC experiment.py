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

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-20 17-13-27\\20240520_sfGFP_QCC_5 ul injection.D"




d1 = extract_data(p1, index = 2)
d2 = extract_data(p2, index = 2)


# Plotting in separate figures
#plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.subplot (2, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="sf GFP before", color="green")
plt.ylim(-50, 500)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
#plt.title('280 nm Blank for rep 3')
plt.legend (loc="upper right")


# Plotting in separate figures
#plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.subplot (2, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="sf GFP after", color="green")
plt.ylim(-50, 500)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
#plt.title('280 nm MSP D1T0 (-) rep 3')
plt.legend (loc="upper right")


plt.tight_layout()

plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240520_490 nm POPC experiment QCC.png", dpi=1000)

plt.show()
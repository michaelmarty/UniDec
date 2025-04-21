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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-22 15-01-51\\20240522_Blank1_ 5 ul injection.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-22 15-01-51\\20240522_4uM D1T0(+)_rep1_5 ul injection.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-22 15-01-51\\20240522_18uM D1T0(+)rep1_5 ul injection.D"

p4 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240520_7 minutes method 2024-05-22 15-01-51\\20240522_sfGFP_QCC_5 ul injection.D"


d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)
d3 = extract_data(p3, index = 0)
d4 = extract_data(p4, index = 2)


# Plotting in separate figures
plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.plot(d1[:, 0], d1[:, 1], label="blank", color="black")
plt.ylim(-50, 400)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('280 nm Blank')
plt.legend()
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240520_520 nm blank for POPC rep 2.png", dpi=1000)

# Plotting in separate figures
plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.plot(d2[:, 0], d2[:, 1], label="4 uM D1T0 (+)", color="blue")
plt.ylim(-50, 400)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('280 nm 4 uM D1T0 (+)')
plt.legend()
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240522_280 nm 4 uM D1T0 (+).png", dpi=1000)


# Plotting in separate figures
plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.plot(d3[:, 0], d3[:, 1], label="18 uM D1T0 (+)", color="blue")
plt.ylim(-50, 400)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('280 nm 18 uM D1T0 (+)')
plt.legend()
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240522_280 nm 18 uM D1T0 (+).png", dpi=1000)

# Plotting in separate figures
plt.figure(figsize=(8, 6))  # Create a new figure for the first plot
plt.plot(d4[:, 0], d4[:, 1], label="sf GFP", color="green")
plt.ylim(-50, 400)
plt.xlim(0, 7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('490 nm sf GFP')
plt.legend()
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240522_490 nm sf GFP.png", dpi=1000)


plt.show()
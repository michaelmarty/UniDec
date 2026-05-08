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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240701_SR7B nanodiscs 2024-07-11 15-37-02\\20240711_Blank_(-)POPC SR7B wc.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240701_SR7B nanodiscs 2024-07-11 15-37-02\\20240711_(-)POPC SR7B wc.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240701_SR7B nanodiscs 2024-07-11 15-37-02\\20240711_Blank_(+)POPC SR7B wc.D"

p4 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240701_SR7B nanodiscs 2024-07-11 15-37-02\\20240711_(+)POPC SR7B wc.D"


d1 = extract_data(p1, index = 3)
d2 = extract_data(p2, index = 3)
d3 = extract_data(p3, index = 3)
d4 = extract_data(p4, index = 3)

plt.plot(d1[:,0], d1[:,1], label="Blank for (-)POPC SR7B", color="black")


plt.plot(d2[:,0], d2[:,1], label="(-)POPC SR7B", color="blue")


plt.plot(d3[:,0], d3[:,1], label="Blank for (+)POPC SR7B", color="black", linestyle='--')

plt.plot(d4[:,0], d4[:,1], label="(+)POPC SR7B", color="red")

plt.ylim(-10, 40)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('Detection of SR7B nanodiscs at 520 nm without column')
plt.legend(loc="upper right")

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240712_SR7B nanodiscs 7 mins wc.png", dpi=1000)


plt.show()
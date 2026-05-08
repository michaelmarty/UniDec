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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\Blank.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\GFP.D"

d1 = extract_data(p1)
d2 = extract_data(p2)

plt.subplot(2, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="Blank", color="black")
plt.ylim(-50, 600)
plt.xlim(0, 30)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')

plt.subplot(2, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="GFP", color="green")
plt.ylim(-50, 600)
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')

plt.tight_layout()

#plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\BlankvsGFP.png", dpi=300)
#plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\BlankvsGFP.pdf")

plt.show()
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

p1 = "Z:\\Group Share\Seun\\LX-MS\LC\\LC Data\OF_vary vol inject 2024-03-02 09-42-57\\026-D1F-E1-Blue Food dye.D"

p2 = "Z:\\Group Share\Seun\\LX-MS\LC\\LC Data\OF_vary vol inject 2024-03-02 09-42-57\\026-D1F-E1-Blue Food dye.D"



d1 = extract_data(p1, index=2)
d2 = extract_data(p2, index=3)

ylims = [-20, 200]

plt.subplot(2, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="Blank", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')

plt.subplot(2, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="GFP", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')


plt.tight_layout()

plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\Bluefooddyeacross.png", dpi=300)
plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\Bluefooddyeacross.pdf")

plt.show()

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

p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-02 17-29-59\\20240402_blue food dye 0.2 ml_min 50 ul.D"
p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blank 0.2 ml_min 30 ul.D"
p3 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blue food dye 0.2 ml_min 20 ul.D"



d1 = extract_data(p1, index=4)
d2 = extract_data(p2, index=4)
d3 = extract_data(p3, index=4)







plt.subplot(3, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="Blue food dye 0.2 ml_min 50 ul", color='#0A2DA4')
plt.ylim(0, 100)
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(3, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="Blue food dye 0.2 ml_min 30 ul", color='#0A2DA4')
plt.ylim(0, 100)
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(3, 1, 3)
plt.plot(d3[:, 0], d3[:, 1], label="Blue food dye 0.2 ml_min 20 ul", color='#0A2DA4')
plt.ylim(0, 100)
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")



plt.tight_layout()

plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240404_bluefooddyecompared100.png", dpi=1000)
#plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\Blankvsbluefooddyecompare.pdf")

plt.show()

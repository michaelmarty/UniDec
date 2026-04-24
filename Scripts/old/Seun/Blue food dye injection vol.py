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

p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blank 0.2 ml_min 5 ul.D"
p2 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blue food dye 0.2 ml_min 5 ul.D"
p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blank 0.2 ml_min 10 ul.D"
p4 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blue food dye 0.2 ml_min 10 ul.D"
p5 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blank 0.2 ml_min 20 ul.D"
p6 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blue food dye 0.2 ml_min 20 ul.D"
p7 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blank 0.2 ml_min 30 ul.D"
p8 =  "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-04 10-24-05\\20240404_Blue food dye 0.2 ml_min 30 ul.D"

d1 = extract_data(p1, index=4)
d2 = extract_data(p2, index=4)
d3 = extract_data(p3, index=4)
d4 = extract_data(p4, index=4)
d5 = extract_data(p5, index=4)
d6 = extract_data(p6, index=4)
d7 = extract_data(p7, index=4)
d8 = extract_data(p8, index=4)


ylims = [0, 100]

plt.subplot(8, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="Blank 0.2 ml_min 5 ul", color="black")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="Blue food dye 0.2 ml_min 5 ul", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 3)
plt.plot(d3[:, 0], d3[:, 1], label="Blank 0.2 ml_min 10 ul", color="black")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 4)
plt.plot(d4[:, 0], d4[:, 1], label="Blue food dye 0.2 ml_min 10 ul", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 5)
plt.plot(d5[:, 0], d5[:, 1], label="Blank 0.2 ml_min 20 ul", color="black")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 6)
plt.plot(d6[:, 0], d6[:, 1], label="Blue food dye 0.2 ml_min 20 ul", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

plt.subplot(8, 1, 7)
plt.plot(d7[:, 0], d7[:, 1], label="Blank 0.2 ml_min 30 ul", color="black")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")


plt.subplot(8, 1, 8)
plt.plot(d8[:, 0], d8[:, 1], label="Blue food dye 0.2 ml_min 30 ul", color='#0A2DA4')
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.legend (loc="upper right")

#plt.tight_layout()

#plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\Blankvsbluefooddyecompare.png", dpi=1000)
#plt.savefig("C:\\Users\\admin\\University of Arizona\\martylab - Documents\\Presentations\\Seun\\Blankvsbluefooddyecompare.pdf")

plt.show()

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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-15 13-50-43\\20240411_rep 1 sf GFP 14 min method 5 ul.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-15 14-54-07\\20240415_Blank after sf GFP binding 5 ul.D"



d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)


plt.subplot(2, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="sf GFP", color="green")
plt.ylim(-50, 600)
plt.xlim(0, 15)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title ('280 nm sf GFP binding check')

plt.subplot(2, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="blank", color="black")
plt.ylim(-50, 600)
plt.xlim(0, 15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('280 nm blank after sf GFP binding check')

plt.tight_layout()

#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240415_280 nm sf GFP binding check.png", dpi=1000)


plt.show()
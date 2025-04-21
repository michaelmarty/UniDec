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


p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-11 15-09-20\\20240411_rep 1 sf GFP 14 min method 5 ul.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-11 15-09-20\\20240411_rep 2 sf GFP 14 min method 5 ul.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-04-11 15-09-20\\20240411_rep 3 sf GFP 14 min method 5 ul.D"

d1 = extract_data(p1, index = 2)
d2 = extract_data(p2, index = 2)
d3 = extract_data(p3, index = 2)


#plt.subplot(3, 1, 1)
plt.plot(d1[:,0], d1[:,1], label="Rep 1", color="green")

#plt.subplot(3, 1, 2)
plt.plot(d2[:,0], d2[:,1] - 100, label="Rep 2", color="green" , linestyle="--")


#plt.subplot(3, 1, 3)
plt.plot(d3[:,0], d3[:,1] - 200, label="Rep 3", color="lime")
plt.ylim(-250, 400)
plt.xlim(0, 15)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
plt.title('490 nm sf GFP 14 min method_All data stacked')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240411_490 nm 14 min method triplicates.png", dpi=1000)


plt.show()
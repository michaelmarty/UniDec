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

p1 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-02-28 15-20-17\\GFP.D"

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-21 17-21-01\\20240321_sf GFP 0.2ml_min 5 ul.D"

d1 = extract_data(p1, index=2)
d2 = extract_data(p2, index=2)


ylims = [-5, 1000]

plt.subplot(2, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="20 ul inject 0.1 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
#plt.title('20 ul injection 0.1 ml/min flowrate', loc='right')
plt.legend (loc="upper right")


plt.subplot(2, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="5 ul inject 0.2 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (mAU)')
#plt.title('5 ul injection 0.2 ml/min flow rate', loc='right')
plt.legend (loc="upper right")

plt.tight_layout()

plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\flowrate_0.1_0.2 ml_mincompared.png", dpi=1000)
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\flowrate_0.1_0.2 ml_mincompared.pdf")

plt.show()

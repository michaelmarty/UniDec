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

p2 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-26 13-01-30\\20240326_sf GFP 0.1 ml_min 5 ul.D"

p3 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-26 13-01-30\\20240326_sf GFP 0.2 ml_min 5 ul.D"

p4 = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_20240309_control test 2024-03-26 13-01-30\\20240326_sf GFP 0.3 ml_min 5 ul.D"

d1 = extract_data(p1, index=2)
d2 = extract_data(p2, index=2)
d3 = extract_data(p3, index=2)
d4 = extract_data(p4, index=2)

d = [d1, d2, d3, d4]
vols = [20, 5, 5, 5]
maxes = [np.amax(data) for data in d]

for data in d:
    #data[:, 1] = data[:, 1] - np.amin(data[:, 1])
    data[:, 1] = data[:, 1] / np.amax(data[:, 1]) * 100


#plt.plot(vols[1:], maxes[1:], 'o')
#plt.show()
#exit()

ylims = [-5, 100]

plt.subplot(4, 1, 1)
plt.plot(d1[:, 0], d1[:, 1], label="20 ul inject 0.1 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (%)')
#plt.title('20 ul injection 0.1 ml/min flowrate', loc='right')
plt.legend (loc="upper right")

plt.subplot(4, 1, 2)
plt.plot(d2[:, 0], d2[:, 1], label="5 ul inject 0.1 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (%)')
#plt.title('20 ul injection 0.5 ml/min flow rate', loc='right')
plt.legend (loc="upper right")


plt.subplot(4, 1, 3)
plt.plot(d3[:, 0], d3[:, 1], label="5 ul inject 0.2 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (%)')
#plt.title('10 ul injection 0.5 ml/min flowrate', loc='right')
plt.legend (loc="upper right")

plt.subplot(4, 1, 4)
plt.plot(d4[:, 0], d4[:, 1], label="5 ul inject 0.3 ml/min", color="green")
plt.ylim(ylims[0], ylims[1])
plt.xlim(0, 30)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('Abs (%)')
#plt.title('1 ul injection 0.5 ml/min flow rate', loc='right')
plt.legend (loc="upper right")

plt.tight_layout()

#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\Mar26integrity check.png", dpi=1000)
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\labelednanodiscscompare.pdf")
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\labelednanodiscscompare.svg")

plt.show()

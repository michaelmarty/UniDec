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


p1 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\\20240701_Blank.D"

p2 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\\20240701_(-)POPC SR7B.D"

p3 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\\20240701_(-)POPC SR7B.D"

p4 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\\20240701_(-)POPC SR7B.D"

p5 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\\20240701_(-)POPC SR7B.D"

#p6 = C:\\Users\Public\Documents\ChemStation\1\Data\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\20240701_(-)POPC SR7B.D"

#p7 = "C:\Users\Public\Documents\ChemStation\1\Data\OF_20240701_FLD optimization for SR7B labeled nanodiscs 2024-07-01 12-45-57\20240701_(-)POPC SR7B.D"

d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)
d3 = extract_data(p3, index = 1)
d4 = extract_data(p4, index = 2)
d5 = extract_data(p5, index = 3)
#d6 = extract_data(p6, index = 5)
#d7 = extract_data(p7, index = 6)



plt.plot(d1[:,0], d1[:,1], label="Blank 510 nm", color="black")

plt.plot(d2[:,0], d2[:,1], label=" FLD 510 nm", color="green")

plt.plot(d3[:,0], d3[:,1], label="FLD 490 nm", color="red")

plt.plot(d4[:,0], d4[:,1], label="FLD 530 nm", color="orchid")

plt.plot(d5[:,0], d5[:,1], label="FLD 540 nm", color="blue")

#plt.plot(d6[:,0], d6[:,1], label="Abs 280 nm", color="lime")

#plt.plot(d7[:,0], d7[:,1], label="Abs 620 nm", color="grey")

plt.ylim(0, 5)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.title('SR7B excitation scan')
plt.legend (loc="upper right")

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240701_FLD SR7B excitation scan.png", dpi=1000)


plt.show()
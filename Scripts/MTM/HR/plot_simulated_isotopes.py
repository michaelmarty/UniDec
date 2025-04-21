import numpy as np
import matplotlib.pyplot as plt
from unidec.modules.isotopetools import *
import os
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')

os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")

m1 = 750
z1 = 1

m2 = 3020
z2 = 4

dist1 = fast_calc_averagine_isotope_dist(m1, z1)

dist2 = fast_calc_averagine_isotope_dist(m2, z2)

mzaxis = np.arange(750, 760, 0.01)
int1 = np.zeros_like(mzaxis)
int2 = np.zeros_like(mzaxis)

sticks1 = np.zeros_like(mzaxis)
sticks2 = np.zeros_like(mzaxis)

sig=0.02

for i, d in enumerate(dist1):
    mz = d[0]
    inten = d[1]
    idx = np.argmin(np.abs(mzaxis - mz))
    sticks1[idx] = inten
    # Gaussian peak of width 0.05 for int1
    int1 += inten * np.exp(-0.5 * (mzaxis - mz) ** 2 / sig ** 2)

for i, d in enumerate(dist2):
    mz = d[0]
    inten = d[1]
    idx = np.argmin(np.abs(mzaxis - mz))
    sticks2[idx] = inten
    # Gaussian peak of width 0.05 for int2
    int2 += inten * np.exp(-0.5 * (mzaxis - mz) ** 2 / sig ** 2)

# Set font to Arial
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 16

plt.plot(mzaxis, sticks1, color="b")
plt.plot(mzaxis, sticks2, color="orange")
plt.plot(mzaxis, int1 * 0, color="k")
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Simplify y axis
plt.yticks([])
plt.ylabel("Intensity")
plt.xlabel("m/z")
plt.savefig("isotopes2.png", dpi=300, transparent=True)
plt.savefig("isotopes2.pdf")
#plt.legend()
plt.show()


plt.plot(mzaxis, int1, color="b")
plt.plot(mzaxis, int2, color="orange")
plt.plot(mzaxis, int1 * 0, color="k")
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Simplify y axis
plt.yticks([])
plt.ylabel("Intensity")
plt.xlabel("m/z")
plt.savefig("isotopes.png", dpi=300, transparent=True)
plt.savefig("isotopes.pdf")
#plt.legend()
plt.show()
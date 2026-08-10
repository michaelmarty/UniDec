import matplotlib.pyplot as plt
import numpy as np
import os

from rdkit.sping.colors import transparent

from unidec.IsoDec.encoding import *
from unidec.IsoDec.c_interface import example
from unidec.IsoDec.plots import *

outdir = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures"
os.chdir(outdir)

example = example[10:]

emat, example = encode_harmonic(example, 11)

cplot(example, color="k")
# shut off top and right axiws
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# set y axis to start at 0
plt.ylim(bottom=0)
# Simplify x ticks
plt.xticks(np.arange(568, 572, 1))

plt.savefig("harmonic.png", dpi=300, transparent=True)

plt.show()


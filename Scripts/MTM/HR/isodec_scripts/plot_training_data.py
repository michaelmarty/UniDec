import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from unidec.IsoDec.encoding import *

# Set backend to Agg
mpl.use('WxAgg')

path = "Z:\\Group Share\\JGP\\PXD019247\\20170105_L_MaD_ColC4_Ecoli20161108-12_01.pkl"

tr, te, zd = encode_phase_file(path)

zvals = range(1, 31)

for z in zvals:
    plt.subplot(10, 3, z)
    plt.ylabel("+" + str(z))
    plt.xticks([])
    plt.yticks([])
    # turn off axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    for t in tr:
        emat = t[0]
        centroids = t[1]
        zval = t[2]

        if zval == z:
            for c in centroids:
                plt.vlines(c[0], 0, c[1], color="k")
            break

os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")
plt.savefig("training_data.png", dpi=300, transparent=True)
plt.savefig("training_data.pdf")
plt.tight_layout()
plt.show()
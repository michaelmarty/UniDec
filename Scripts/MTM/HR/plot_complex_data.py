import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import unidec.tools as ud

from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import plot_pks

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoDecEngine(phaseres=8)

file2 = "C:\\Data\\IsoNN\\test2.txt"


os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")

if os.path.splitext(file2)[1] == ".txt":
    spectrum = np.loadtxt(file2, skiprows=0)
    spectrum = ud.datachop(spectrum, 953.5, 959)
    eng.use_wrapper = False
    eng.batch_process_spectrum(spectrum, centroided=True)
    plot_pks(eng.pks, centroids=spectrum, ccolor="k", plotmass=False)
    plt.tight_layout()
    plt.savefig("test2.png", dpi=300, transparent=True)
    plt.savefig("test2.pdf")

    plt.show()
    exit()


# Set backend to Agg
mpl.use('WxAgg')

file2 = "C:\\Data\\IsoNN\\test2.txt"


data = np.loadtxt(file2, skiprows=0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 16

data = ud.datachop(data, 953.5, 959)

for d in data:
    plt.vlines(d[0], 0, d[1], color="k")

plt.hlines(0, np.amin(data[:, 0]), np.amax(data[:, 0]), color="k")
plt.xlim
#plt.plot(data[:, 0], data[:, 1], color="k")
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
# Simplify y axis
plt.yticks([])
plt.ylabel("Intensity")
plt.xlabel("m/z")

plt.show()
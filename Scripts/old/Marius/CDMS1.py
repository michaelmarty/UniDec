import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.tools as ud

mpl.use("TkAgg")
from unidec import modules as TI

file = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\092292020\\20200929_MK_DPPG_pos_10xdilution_singleion4_pressure05_actual.raw"
TDI = TI.ThermoDataImporter(file)

pairs = None
max_scan = TDI.get_max_scans()
for s in range(1, max_scan):

    data = TDI.grab_scan_data(s)

    boo1 = data[:,1] > 2000

    data = data[boo1]

    if pairs is None:
        pairs=data
    else:
        pairs = np.concatenate((pairs, data), axis=0)

plt.subplot(121)
plt.scatter(pairs[:,0], pairs[:,1], marker="o", s=0.5)

slice = ud.datachop(pairs, 5750, 5950)
plt.subplot(122)
plt.hist(slice[:,1], bins=30)
plt.show()



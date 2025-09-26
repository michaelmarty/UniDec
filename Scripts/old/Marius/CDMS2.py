import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import unidec.tools as ud

mpl.use("TkAgg")
from unidec import modules as TI


def get_cdf(data):
    return np.cumsum(data)


def get_noise_level(data):
    std = np.std(data)
    return std * 3 + np.mean(data)


file =  "Z:\\Group Share\\Marius Kostelic\\CD-MS\\10192020 DPPC\\201019_mk_BSA_pos_singleion_longaquisition.raw"
#file = "C:\\Data\\CDMS\\LW_20200922_MK_BSA_pos_10x_dilution2_singleiontune2.RAW"
fname = os.path.splitext(file)[0]+".txt"
dir = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\10192020 DPPC"
os.chdir(dir)

if True:
    TDI = TI.ThermoDataImporter(file)

    pairs = None
    max_scan = TDI.get_max_scans()
    for s in range(1, max_scan):

        data = TDI.grab_scan_data(s)

        nl = get_noise_level(data[:, 1])
        boo1 = data[:, 1] > nl

        data = data[boo1]

        if pairs is None:
            pairs = data
        else:
            pairs = np.concatenate((pairs, data), axis=0)

    nl = get_noise_level(pairs[:, 1])
    boo1 = pairs[:, 1] > nl
    pairs = pairs[boo1]

    np.savetxt(fname, pairs)
else:
    pairs = np.loadtxt(fname)

print("Imported:", len(pairs))

plt.subplot(131)
plt.scatter(pairs[:, 0], pairs[:, 1])

bsa = 66430

slices = [3900, 4150, 4450, 5550, 5780, 6050]
z = [17, 16, 15, 24, 23, 22]
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)
fits = []
for s in slices:
    slice = ud.datachop(pairs, s, s + 100)
    # ax2.hist(slice[:,1], bins=30)
    cdf = (get_cdf(slice[:, 1])) / len(slice)
    x = np.arange(len(cdf)) / len(cdf)
    fit = np.polyfit(x, cdf, 1)
    min = np.amin(slice[:,1])
    avg = np.mean(slice[:,1])
    #fits.append(fit[0])
    fits.append(avg)
    print(s, fit[0], len(cdf), avg)
    ax3 = plt.plot(x, cdf)

fits = np.array(fits)
ax2.plot(z, fits, marker="o", linestyle=" ")
plt.show()

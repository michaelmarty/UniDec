import numpy as np
import h5py
import matplotlib.pyplot as plt
import unidec.tools as ud

file = "Z:\\Group Share\\JamesKeener Keener\\20190501_JEK_AmtB_PGPC_E3TO_B1_GC5_Large_HCDAIF_Pos_Ramp_0_200_20_testing.hdf5"
file = "Z:\Group Share\James Keener\September AmtB mix data\PREP 1\GC\GOOD DATA\\20190501_JEK_AmtB_PGPC_E3TO_B1_GC5_Large_HCDAIF_Pos_Ramp_0_200_20_testing.hdf5"

file = "C:\\Users\\Michael Marty\\Desktop\\20190501_JEK_AmtB_PGPC_E3TO_B1_GC5_Large_HCDAIF_Pos_Ramp_0_200_20_testing.hdf5"
"""
Methods:
3 = Center of Mass
4 = Local Max Position
5 = CoM 50%
6 = CoM 10%

"""


hdf = h5py.File(file)

agroup = hdf.get("ms_dataset")
mass_axis = get_dataset(agroup, "mass_axis")
mass_sum = get_dataset(agroup, "mass_sum")
hdf.close()

data = np.transpose([mass_axis, mass_sum])

indexes = np.arange(0, 50)
popc = 760.076
popg = 749.02

lmass = (popc + popg) / 2.
basemass = 126796
window = 200
method = 4

isodat = ud.datachop(data, basemass - window, basemass + window)
basemass = ud.data_extract(isodat, basemass, method, window=window)
basemass=126797.63
print(basemass)

lmasses = basemass + lmass * indexes

coms = []
comsor = []

colormap = plt.get_cmap('viridis')
colors = [colormap(float(k) / len(indexes)) for k in indexes]
plt.figure()
plt.subplot(121)
for i, l in enumerate(lmasses[1:]):
    isodat = ud.datachop(data, l - window, l + window)
    comsor.append(ud.data_extract(isodat, l, method, window=window))

    isodat[:, 0] = isodat[:, 0] - l
    # isodat2 = deepcopy(isodat)
    isodat[:, 0] /= i + 1
    isodat[:, 0] += lmass

    com = ud.data_extract(isodat, lmass, method, window=window)
    coms.append(com)
    plt.plot(isodat[:, 0], isodat[:, 1], color=colors[i])

plt.vlines(lmass, 0, 1)
plt.xlim(popg-10, popc+10)
plt.subplot(122)

coms = np.array(coms)
comsor = np.array(comsor)
fit=np.polyfit(indexes[1:], comsor,1)
cgrad = np.gradient(comsor)
print(fit)
lfit = (fit[0]-popg)/(popc-popg)
# coms2 = coms / indexes[1:] + lmass
'''
plt.plot(coms, marker="o")
plt.hlines(lmass, 0, len(indexes))
plt.hlines(popc, 0, len(indexes))
plt.hlines(popg, 0, len(indexes))'''
coms = (coms-popg) / (popc-popg)
cgrad = (cgrad-popg) / (popc-popg)
plt.plot(coms, marker="o")
plt.plot(cgrad, marker="o")
plt.hlines(0.5, 0, len(indexes))
plt.hlines(1, 0, len(indexes))
plt.hlines(0, 0, len(indexes))
plt.hlines(lfit, 0, len(indexes))
plt.show()

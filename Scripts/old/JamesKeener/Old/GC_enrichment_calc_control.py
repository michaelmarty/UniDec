import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec import metaunidec as me

# dir = "Z:\\Group Share\\JamesKeener Keener\\AmtB 50PG 50PC\\GC P1_P2_P3"
dir = "C:\\Data\\JamesKeener\\GC_cont"
os.chdir(dir)
files = ud.match_files(dir, "*.hdf5")

dir2 = "C:\\Data\\JamesKeener\\GC_mixed"
os.chdir(dir2)
files2 = ud.match_files(dir2, "*.hdf5")

popc = 760.076
popg = 749.02

lmass = (popc + popg) / 2.
basemass = 126796
window = 250
method = 5


def get_basemass(d, guess=126796):
    isodat = ud.datachop(d, guess - window, guess + window)
    basemass = ud.data_extract(isodat, guess, method, window=window)
    return basemass


def get_basemass2(d, guess=126796, lmass=754.5):
    indexes = np.arange(0, 40)
    lmasses = guess + lmass * indexes
    commass = []
    ints = []
    for i, l in enumerate(lmasses[1:]):
        isodat = ud.datachop(d, l - window, l + window)
        ints.append(np.amax(isodat[:, 1]))
        commass.append(ud.data_extract(isodat, l, method, window=window))

    fits = np.polyfit(indexes[1:], commass, 1, w=np.array(ints) ** 2)
    # print(fits)
    basemass = fits[1]
    lmass = fits[0]
    return basemass, lmass


def lr(d, basemass=126802, lmass=754.5):
    indexes = np.arange(2, 40)
    basemass, lmass = get_basemass2(d, lmass=lmass)
    basemass, lmass = get_basemass2(d, guess=basemass, lmass=lmass)
    print(basemass)

    lmasses = basemass + lmass * indexes

    commass = []
    ints = []
    for i, l in enumerate(lmasses[1:]):
        isodat = ud.datachop(d, l - window, l + window)
        ints.append(np.amax(isodat[:, 1]))
        commass.append(ud.data_extract(isodat, l, method, window=window))

    coms = (np.array(commass) - basemass) / indexes[1:]

    return coms, ints, indexes


def calc_e(file, lmass=754.5):
    eng = me.MetaUniDec()
    eng.open(file)
    eng.data.import_grids()

    data = ud.datachop(eng.data.massdat, 125000, 180000)
    # basemass = get_basemass(eng.data.massdat, guess=126805)
    basemass, lmass = get_basemass2(eng.data.massdat, guess=126805, lmass=lmass)
    # basemass = 126792+17

    #cav, ints = lr(data, basemass=basemass, lmass=lmass)
    # print(basemass)
    if True:
        dataset = []
        cv = []
        for s in eng.data.spectra:
            if s.var1 >= 60 and s.var1 <= 320:
                isodat = ud.datachop(s.massdat, 125000, 180000)
                dataset.append(isodat)
                cv.append(s.var1)

        coms = []
        ints = []
        for d in dataset:
            com, int, indexes = lr(d, lmass=lmass)
            coms.append(com)
            ints.append(int)

        cav = np.average(coms, axis=0, weights=np.array(ints) ** 2)

    return cav, basemass, indexes


cavs = []
bmasses = []
for f in files:
    cav, bm, indexes = calc_e(os.path.join(dir, f), lmass=popc)
    cavs.append(cav)
    bmasses.append(bm)

cav = np.average(cavs, axis=0)
stds = np.std(cavs, axis=0)

print(np.mean(bmasses))

plt.scatter(indexes[1:],cav, c='orange')
plt.errorbar(indexes[1:], cav, yerr=stds, c='black', capsize=3, capthick=1, elinewidth=1)

cavs2 = []
bmasses2 = []
for f in files2:
    cav, bm, indexes = calc_e(os.path.join(dir2, f), lmass=754.5)
    bmasses2.append(bm)
    cavs2.append(cav)

cav2 = np.average(cavs2, axis=0)
stds2 = np.std(cavs2, axis=0)
print(np.mean(bmasses2))

plt.ylim((751, 762))
#f2=interp1d(indexes[1:],cav2,kind='cubic')
#xnew = np.linspace(3, 39, num=1000,endpoint=True)
#plt.plot(indexes[1:],cav2,'o', xnew,f2(xnew),'-')
plt.scatter(indexes[1:], cav2, c='yellow')
plt.errorbar(indexes[1:], cav2, yerr=stds2, c='black', capsize=3, capthick=1, elinewidth=1)
plt.hlines(lmass, 0, len(cav))
#plt.hlines(popg, 0, len(cav))
plt.hlines(popc, 0, len(cav))
#plt.hlines((0.25 * popc + 0.75 * popg), 0, len(cav))

plt.show()
exit()

calc_e(os.path.join(dir, files[0]))
calc_e(os.path.join(dir2, files2[0]))

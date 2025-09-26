import os
import numpy as np
import matplotlib.pyplot as plt
import UniDec.unidec_modules.unidectools as ud
from UniDec import metaunidec as me

# dir = "Z:\\Group Share\\James Keener\\AmtB 50PG 50PC\\GC P1_P2_P3"


dir2 = "C:\\JEK data copy\\25PEPG"
os.chdir(dir2)
files2 = ud.match_files(dir2, "*.hdf5")

popg = 749.02
SOPE = 746.050
popc= 760.076
POPE= 717.996

lmass = 741.264
msp = 29982
amtb = 126796
basemass = amtb
window = 250
method = 5
nrange = [0, 20]
basemasscorrection = True
outheader = "GC"
print(files2)

bm1 = [basemass, basemass, basemass]
bm2 = [basemass, basemass, basemass]


# exit()


def get_basemass(d, guess=126796):
    isodat = ud.datachop(d, guess - window, guess + window)
    basemass = ud.data_extract(isodat, guess, method, window=window)
    return basemass


def get_basemass2(d, guess=126796, lmass=741.264, nrange=nrange):
    indexes = np.arange(nrange[0], nrange[1])
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
    if basemasscorrection:
        guess = basemass
    return guess, lmass


def lr(d, guess=126802, lmass=741.264, nrange=nrange):
    indexes = np.arange(nrange[0], nrange[1])
    basemass, lmass = get_basemass2(d, guess=guess, lmass=lmass, nrange=nrange)
    basemass, lmass = get_basemass2(d, guess=basemass, lmass=lmass, nrange=nrange)
    print(basemass, lmass)
    lmasses = basemass + lmass * indexes

    commass = []
    ints = []
    for i, l in enumerate(lmasses[1:]):
        isodat = ud.datachop(d, l - window, l + window)
        ints.append(np.amax(isodat[:, 1]))
        commass.append(ud.data_extract(isodat, l, method, window=window))

    coms = (np.array(commass) - basemass) / indexes[1:]

    return coms, ints, indexes


def calc_e(file, lmass=741.264, guess=126805, nrange=nrange):
    eng = me.MetaUniDec()
    eng.open(file)
    eng.data.import_grids()

    data = ud.datachop(eng.data.massdat, 0, 350000)

    if ud.isempty(data):
        eng.make_grids()
        eng.data.import_grids()
        data = ud.datachop(eng.data.massdat, guess-1000, 350000)

    # basemass = get_basemass(eng.data.massdat, guess=126805)
    basemass, lmass = get_basemass2(eng.data.massdat, guess=guess, lmass=lmass, nrange=nrange)
    # basemass = 126792+17

    # cav, ints = lr(data, basemass=basemass, lmass=lmass)
    # print(basemass)
    if True:
        data[:, 1] *= 0
        dataset = []
        cv = []
        for s in eng.data.spectra:

            if s.var1 >= 0 and s.var1 <= 200:
                isodat = ud.datachop(s.massdat, basemass-1000, 350000)
                dataset.append(isodat)
                cv.append(s.var1)

                d2 = ud.mergedata(data, isodat)
                data[:, 1] += d2[:, 1]

        coms = []
        ints = []
        for d in dataset:
            com, int, indexes = lr(d, guess=guess, lmass=lmass, nrange=nrange)
            coms.append(com)
            ints.append(int)

        cav = np.average(coms, axis=0, weights=np.array(ints) ** 2)

    return cav, basemass, indexes, data


if True:
    data2 = []
    cavs2 = []
    bmasses2 = []
    for i, f in enumerate(files2):
        cav, bm, indexes, d = calc_e(os.path.join(dir2, f), guess=bm2[i], lmass=lmass, nrange=nrange)
        bmasses2.append(bm)
        cavs2.append(cav)
        data2.append(d)

    cav2 = np.average(cavs2, axis=0)
    stds2 = np.std(cavs2, axis=0)
    print(np.mean(bmasses2))

    plt.scatter(indexes[1:], cav2, c='yellow')
    plt.errorbar(indexes[1:], cav2, yerr=stds2, c='black', capsize=3, capthick=1, elinewidth=1)

    outdat = np.array([indexes[1:], cav2, stds2])
    np.savetxt(outheader + "_mixed_lipav.txt", outdat)

# plt.subplot(211)
# plt.ylim((749, 760))
# f2=interp1d(indexes[1:],cav2,kind='cubic')
# xnew = np.linspace(3, 39, num=1000,endpoint=True)
# plt.plot(indexes[1:],cav2,'o', xnew,f2(xnew),'-')

# plt.hlines((0.5 * popc + 0.5 * SOPE), 0, len(cav))
plt.hlines(POPE, 0, len(cav))
plt.hlines(popg, 0, len(cav))
plt.hlines((0.75 * popg + 0.25 * POPE), 0, len(cav))

# plt.subplot(212)

datavg = []
for i, d in enumerate(data2):
    b2 = bmasses2[i]
    lm = np.mean(cavs2[i])
    boo1 = d[:, 0] > b2 - 500
    dat = d[boo1]
    dat[:, 0] = (dat[:, 0] - b2) / lm
    if i == 0:
        datavg = dat
    else:
        d2 = ud.mergedata(datavg, dat)
        datavg[:, 1] += d2[:, 1]
    plt.plot(dat[:, 0], 760 + 10 * d[boo1, 1] / np.amax(d[boo1, 1]), color="b")

np.savetxt(outheader + "_mixed_datavg.txt", datavg)

# plt.show()

# plt.figure()
# b1=np.mean(bmasses)
# for i,d in enumerate(data):
#    b1 = bmasses[i]
#    boo1 = d[:,0] > b1
#    plt.plot(d[boo1, 0],  d[boo1, 1]/np.amax(d[boo1,1]), label="1"+str(i))

# b2=np.mean(bmasses2)
# for i, d in enumerate(data2):
#    b2 = bmasses2[i]
#    boo1 = d[:,0] > b2
#    plt.plot(d[boo1, 0], d[boo1, 1]/np.amax(d[boo1,1]), label="2"+str(i))
# plt.legend()
plt.show()
exit()

calc_e(os.path.join(dir2, files2[0]))
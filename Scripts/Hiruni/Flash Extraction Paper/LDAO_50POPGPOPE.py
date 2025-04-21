import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec import metaunidec as me

dir2 = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_JEK AqpZ ND Project\\Orbitrap Data\\POPG_POPE\\LDAO"
os.chdir(dir2)
files2 = ud.match_files(dir2, "*.hdf5")

POPG = 748.999
SOPE = 746.050
POPC= 760.076
POPE= 717.996
msp = 29982
#amtb = 126796
aqpz = 98861
basemass = aqpz

lipids = [POPG, POPE]
voltage_range = [120, 180]

lmass = np.average(lipids)
print(lmass)

window = 250
method = 5
nrange = [1, 7]
basemasscorrection = True
outheader = "LDAO"
print(files2)
bm2 = [basemass, basemass, basemass]



def get_basemass(d, guess=126796):
    # Use a simple thing to calculate the base mass from the apo mass peak alone
    isodat = ud.datachop(d, guess - window, guess + window)
    basemass = ud.data_extract(isodat, guess, method, window=window)
    return basemass


def get_basemass2(d, guess=126796, lmass=754.538, nrange=nrange):
    # Calculate the base mass from a linear regression
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


def lr(d, guess=126802, lmass=754.538, nrange=nrange):
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


def calc_e(file, lmass=754.538, guess=126805, nrange=nrange):
    # Launch an engine and import the data
    eng = me.MetaUniDec()
    eng.open(file)
    eng.data.import_grids()

    data = ud.datachop(eng.data.massdat, 0, 350000)

    if ud.isempty(data):
        eng.make_grids()
        eng.data.import_grids()
        data = ud.datachop(eng.data.massdat, guess-1000, 350000)

    # Use Simple Apo Mass
    # basemass = get_basemass(eng.data.massdat, guess=guess)

    # Use Linear Regression
    basemass, lmass = get_basemass2(eng.data.massdat, guess=guess, lmass=lmass, nrange=nrange)
    # basemass = 126792+17

    # cav, ints = lr(data, basemass=basemass, lmass=lmass)
    print(basemass)
    if True:
        data[:, 1] *= 0
        dataset = []
        cv = []
        for s in eng.data.spectra:

            if s.var1 >= voltage_range[0] and s.var1 <= voltage_range[1]:
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

# Calculate Average Lipid Masses and Extracted Data
data2 = []
cavs2 = []
bmasses2 = []
for i, f in enumerate(files2):
    cav, bm, indexes, d = calc_e(os.path.join(dir2, f), guess=bm2[i], lmass=lmass, nrange=nrange)
    bmasses2.append(bm)
    cavs2.append(cav)
    data2.append(d)


# Make Average Mass Plot
cav2 = np.average(cavs2, axis=0)
stds2 = np.std(cavs2, axis=0)
print('Average Base Mass:', np.mean(bmasses2))
plt.scatter(indexes[1:], cav2, c='yellow')
plt.errorbar(indexes[1:], cav2, yerr=stds2, c='black', capsize=3, capthick=1, elinewidth=1)

# Save Output Data
outdat = np.array([indexes[1:], cav2, stds2])
np.savetxt(outheader + "_mixed_lipav.txt", outdat)


# Create Horizontal lines
plt.hlines(lipids[0], 0, len(cav))
plt.hlines(lipids[1], 0, len(cav))
plt.hlines(lmass, 0, len(cav))

# Merging and plotting extracted data
datavg = []
for i, d in enumerate(data2):
    b2 = bmasses2[i]
    lm = np.mean(cavs2[i])
    boo1 = d[:, 0] > b2 - 500
    boo2 = d[:,0] < b2 + lmass * (nrange[1] +2)
    boo1 = np.logical_and(boo1, boo2)
    dat = d[boo1]
    dat[:, 0] = (dat[:, 0] - b2) / lm
    dat[:, 1] /= np.amax(dat[:,1])
    if i == 0:
        datavg = dat
    else:
        d2 = ud.mergedata(datavg, dat)
        datavg[:, 1] += d2[:, 1]
    #plt.plot(dat[:, 0], 760 + 10 * d[boo1, 1] / np.amax(d[boo1, 1]))#, color="b")

# Plot Average
plt.plot(datavg[:,0], 760 + 10 * datavg[:, 1] / np.amax(datavg[:, 1]), color="k")

#Save Data Ouput
np.savetxt(outheader + "_mixed_datavg.txt", datavg)
plt.savefig("Figure.pdf")
plt.show()
exit()
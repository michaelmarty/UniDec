import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec import metaunidec as me

dir = "Z:\\Group Share\\JamesKeener Keener\\AmtB 50PG 50PC\\GC P1_P2_P3"
#dir = "C:\\Users\\Michael Marty\\Desktop"
os.chdir(dir)
files = ud.match_files(dir, "*.hdf5")

def lr(d):
    print("yup")


def calc_e(file):
    eng = me.MetaUniDec()
    eng.open(file)
    eng.data.import_grids()

    data = ud.datachop(eng.data.massdat, 125000, 180000)
    plt.plot(data[:,0],data[:,1])

    dataset = []
    cv = []
    for s in eng.data.spectra:
        if s.var1 > 50:
            isodat = ud.datachop(s.massdat, 125000, 180000)
            dataset.append(isodat)
            plt.plot(isodat[:, 0], isodat[:, 1])
            cv.append(s.var1)

    for d in dataset:
        lr(d)

    plt.show()
    exit()
    indexes = np.arange(0, 50)
    popc = 760.076
    popg = 749.02

    lmass = (popc + popg) / 2.
    basemass = 126796
    window = 200
    method = 4

    isodat = ud.datachop(data, basemass - window, basemass + window)
    basemass = ud.data_extract(isodat, basemass, method, window=window)
    basemass = 126797.63
    print(basemass)

    lmasses = basemass + lmass * indexes


calc_e(os.path.join(dir, files[0]))
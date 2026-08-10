from IsoDec.match import match_charge, cplot
from IsoDec.datatools import get_centroids, get_noise, calc_match_simp, create_isodist, simp_charge
from IsoDec.engine import get_charge_nn, IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
# Set backend to Agg
mpl.use('WxAgg')


def get_charge(data, peakmz, mzwindow=[-1.5, 3.5], eng=None):
    centroids, data = get_centroids(data, peakmz, mzwindow=mzwindow)

    # print("Centroids:", centroids)
    if eng is not None:
        pks = eng.fancy_prediction(centroids)
        for p in pks.peaks:
            isodist = p.isodist
            cplot(isodist, color='g', factor=-1)

    else:
        charge1 = get_charge_nn(centroids)
        isodist1 = create_isodist(peakmz, charge1, data)
        cplot(isodist1, color='g', factor=-1)

    cplot(centroids)
    plt.plot(data[:, 0], data[:, 1])
    plt.xlim(peakmz + mzwindow[0], peakmz + mzwindow[1])
    plt.show()


def isotope_detection(data, mzwindow=1.5, eng=None):
    print(len(data))
    nl = get_noise(data)
    # Chop data below noise level
    data = data[data[:, 1] > nl]
    print(nl, len(data))
    peaks = ud.peakdetect(data, window=50, threshold=0.01)
    # sort peaks
    peaks = np.array(sorted(peaks, key=lambda x: x[1], reverse=True))
    for i, p in enumerate(peaks[6:]):
        print(i, p)
        mz = p[0]
        get_charge(data, mz, eng=eng)


eng = IsoDecEngine()

file = "C:\\Data\\TabbData\\test3.txt"

data = np.loadtxt(file, skiprows=16)

isotope_detection(data, eng=eng)

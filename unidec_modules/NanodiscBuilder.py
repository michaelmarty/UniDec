import numpy as np
import matplotlib.pyplot as plt
import unidec_modules.MassSpecBuilder as msb
from unidec_modules import unidectools as ud
import scipy.stats as stats

belts = {"T0": 22044., "G1": 22101., "T1": 22145., "R1": 22200., "T2": 22246., "MSP1D1(-)": 22044., }
lipids = {"POPC": 760.076, "POPG": 749.02, "PGPC75": 757.312, "PGPC50": 754.548, "Chol": 386.65, "POPE": 717.996,
          "DMPC": 677.9, "TMCL": 1241.6, "CHS": 486.73, "Erg": 396.66, "GM1": 1573, "DMPG": 667, "TOCL": 1458, "SM": 731.081,
          "DPPC":734.039}


def make_nanodiscs(array):
    """
    :param array: Array of different nanodiscs of form:
     N x 6 [belt mass, lipid mass, lipid num, lipid std dev, mass sig, z, zstd dev]
    :return:
    """
    arrays = []
    for n in array:
        bmass = n[0]
        lmass = n[1]
        lnum = n[2]
        lsig = n[3]
        z = n[5]
        zsig = n[6]
        msig = n[4]
        thresh = 4
        dist = np.arange(-int(lsig * thresh), int(lsig * thresh))
        masses = np.array([bmass * 2. + lmass * (lnum + i) for i in dist])
        ints = ud.ndis_std(dist, 0, lsig)

        newarray = np.array([[m, msig, z, zsig, ints[i]] for i, m in enumerate(masses)])
        for a in newarray:
            arrays.append(a)
    arrays = np.array(arrays)
    out, ztab = msb.make_mass_spectrum(arrays, mzrange=(6000, 9000))
    plt.plot(out[:, 0], out[:, 1])
    plt.show()


def make_mixed_nanodiscs(array, plot=True):
    outarray = []
    for n in array:
        arrays = []
        bmass = n[0]
        ratio = n[1]
        lnum = n[2]
        lsig = n[3]
        z = n[5]
        zsig = n[6]
        msig = n[4]
        lmass1 = n[7]
        lmass2 = n[8]
        thresh = 4
        dist = np.arange(-int(lsig * thresh), int(lsig * thresh)) + lnum
        ints = ud.ndis_std(dist, lnum, lsig)

        masses = []
        ivals = []
        for n, d in enumerate(dist):
            nvals = np.arange(0, d + 1.)
            nvals2 = d - nvals

            p1 = stats.binom.pmf(nvals, d, ratio)

            m = np.array([bmass * 2. + lmass1 * nvals[i] + lmass2 * nvals2[i] for i in range(0, len(nvals))])
            ival = np.array([ints[n] * p1[i] for i in range(0, len(nvals))])

            for mass in m:
                masses.append(mass)
            for i in ival:
                ivals.append(i)

        masses = np.array(masses)
        ivals = np.array(ivals)

        print(len(masses))

        b1 = ivals / np.amax(ivals) > 0.01

        masses = masses[b1]
        ivals = ivals[b1]

        newarray = np.array([[m, msig, z, zsig, ivals[i]] for i, m in enumerate(masses)])
        for a in newarray:
            arrays.append(a)
        arrays = np.array(arrays)
        avgmass=np.average(masses)
        #start = np.amin(masses)
        #end = np.amax(masses)
        start=avgmass-10000
        end=avgmass+10000
        out, ztab = msb.make_mass_spectrum(arrays, mzrange=(start, end), zrange=(1, 2), mz_bin_size=1)
        outarray.append(out)

    if plot:
        plt.figure()
        for out in outarray:
            plt.plot(out[:, 0], out[:, 1])
        # print(ud.fft_diff(out, [360, 400])[0])
        # print(ud.fft_diff(out, [750, 770])[0])
        plt.show()

    return outarray


if __name__ == '__main__':
    # array = [[belts["T0"], lipids["PGPC75"], 140, 10, 100, 20, 1],
    #         [belts["T1"], lipids["PGPC75"], 140, 10, 100, 20, 1]]
    # array = [[belts["T0"], lipids["POPC"], 140, 10, 100, 20, 1],
    #         [belts["T0"], lipids["DMPC"], 140, 10, 100, 20, 1]]
    # make_nanodiscs(array)
    array = [[belts["T0"], 0.70, 140, 5, 200, 1, 0.0001, lipids["POPE"], lipids["POPG"]]]  # ,
    # [belts["T0"], 0.05, 140, 10, 100, 20, 1, lipids["POPC"], lipids["Chol"]]]
    make_mixed_nanodiscs(array)

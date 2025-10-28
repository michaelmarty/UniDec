import numpy as np
import matplotlib.pyplot as plt
import unidec.modules.MassSpecBuilder as msb
import unidec.tools as ud
import scipy.stats as stats
import os

lipids = {"POPC": 760.076, "POPG": 749.02, "PGPC75": 757.312, "PGPC25": 751.784, "PGPC50": 754.548, "Chol": 386.65, "POPE": 717.996, "PEPG50": 733.508, "PEPG25": 725.752, "PEPG75": 741.264, "PEPC50": 739.036, "PEPC25": 728.516, "PEPC75": 749.556,
          "DMPC": 677.9, "TMCL": 1241.6, "CHS": 486.73, "Erg": 396.66, "GM1": 1573, "DMPG": 667, "TOCL": 1458, "SM": 731.081,
          "DPPC":734.039}

def make_spectrum(array, plot=True):
    outarray = []
    for n in array:
        arrays = []
        # Array values
        bmass = n[0]
        ratio = n[1]
        lnum_min = n[2]
        lnum_max = n[3]
        msig = n[4]
        z = n[5]
        zsig = n[6]
        lmass1 = n[7]
        lmass2 = n[8]

        # Set up distribution in number of lipids and intensity of each peak
        dist = np.arange(lnum_min, lnum_max)
        lnum_avg = (lnum_min + lnum_max)/2.
        lnum_diff = (lnum_max-lnum_min)
        ints = ud.ndis_std(dist, lnum_avg-5, lnum_diff/3.)

        masses = []
        ivals = []
        # Loop over each bound lipid
        for n, d in enumerate(dist):
            # Number of lipid 1
            nvals = np.arange(0, d + 1.)
            # number of lipid 2
            nvals2 = d - nvals
            # Use binomial distribution for intensities of each substate
            p1 = stats.binom.pmf(nvals, d, ratio)
            # Calcualte masses and intensities of each substate
            m = np.array([bmass + lmass1 * nvals[i] + lmass2 * nvals2[i] for i in range(0, len(nvals))])
            ival = np.array([ints[n] * p1[i] for i in range(0, len(nvals))])
            for mass in m:
                masses.append(mass)
            for i in ival:
                ivals.append(i)

        masses = np.array(masses)
        ivals = np.array(ivals)

        print(len(masses))
        # Remove low abundance peaks
        b1 = ivals / np.amax(ivals) > 0.001

        masses = masses[b1]
        ivals = ivals[b1]
        # Create specra using the msb module
        newarray = np.array([[m, msig, z, zsig, ivals[i]] for i, m in enumerate(masses)])
        for a in newarray:
            arrays.append(a)
        arrays = np.array(arrays)
        avgmass=np.average(masses)
        #start = np.amin(masses)
        #end = np.amax(masses)
        start=avgmass-10000
        end=avgmass+10000
        start = 1000
        end = 20000
        out, ztab = msb.make_mass_spectrum(arrays, mzrange=(start, end), zrange=(1, 25), mz_bin_size=1)
        outarray.append(out)

    # Make the plot
    if plot:
        plt.figure()
        for out in outarray:
            plt.plot(out[:, 0], out[:, 1])
        plt.show()

    return outarray

# Protein Mass
# Ratio of Lipids
# Min # of Lipids
# Max # of Lipids
# Mass FWHM
# Avg. Charge
# Std. Dev Charge
# Lipid 1 mass (pull from dict)
# Lipid 2 mass
array1 = [[126796, 0.25, 0, 10, 100, 18, 1, lipids["POPE"], lipids["POPC"]]]
array2 = [[126796, 0.5, 0, 20, 100, 19, 1, lipids["POPE"], lipids["POPC"]]]
data1 = make_spectrum(array1)
data2 = make_spectrum(array2)
dir = "C:\JEK data copy"
os.chdir(dir)
np.savetxt('enrich.txt', data1[0])
np.savetxt('enrich2.txt',data2[0])

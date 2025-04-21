from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pyperclip
import pandas as pd
from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

# Set backend to Agg
mpl.use('WxAgg')
if False:
    file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
    os.chdir(os.path.dirname(file2))

    reader = ImporterFactory.create_importer(file2)
    spectrum = reader.grab_centroid_data(1)

    spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

    mzrange= [1196, 1201]
    mzrange = None

    if mzrange is not None:
        spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

    os.chdir("D:\\Johnny\\IsoDec\\TestData")
    np.savetxt("21tetd.txt", spectrum)
else:
    os.chdir("D:\\Johnny\\IsoDec\\TestData")
    spectrum = np.loadtxt("21tetd.txt")

verbose = False
knockdown = 5

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.config.knockdown_rounds = knockdown
eng.predmode = 0

timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True)

mc1 = eng.pks

print("Time Python:", time.time() - timestart)

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.config.knockdown_rounds = knockdown
eng.predmode = 3

timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True)

mc2: List[MatchedPeak] = eng.pks

print("Time C:", time.time() - timestart)

shared, u1, u2 = compare_matchedcollections(mc1, mc2, ppmtol=50, objecttocompare="monoiso", maxshift=3)

print("N Peaks:", len(mc1), len(mc2))
print("Shared Peaks: " + str(len(shared)), "Unique Peaks Python: " + str(len(u1)),
      "Unique Peaks C: " + str(len(u2)))

try:
    percent_matched = (len(mc2) - len(u2)) / len(mc2) * 100
    print("Percent of C peaks with a match: " + str(percent_matched))

    percent_clusters_matched = len(shared) / len(mc1) * 100
    print("Percent of Python peaks with a match: " + str(percent_clusters_matched))
except:
    pass


plt.figure()
plot_pks(u1,forcecolor="r", centroids=spectrum, ccolor="k", show=False, plotmass=False, nocentroids=True)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True)
plot_pks(shared, forcecolor="g", show=False, plotmass=False, nocentroids=True)
plt.show()
exit()
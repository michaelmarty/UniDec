from unidec.IsoDec.runtime import IsoDecRuntime
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import *
import time
from unidec.IsoDec.plots import plot_pks

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoDecRuntime(phaseres=4)
eng.config.knockdown_rounds = 5

#file2 = "D:\\Johnny\\GitClones\\mzLib\\mzLib\\Test\\DataFiles\\14kDaProteoformMzIntensityMs1.txt"
# file2 = "C:\\Data\\TabbData\\23_04_21_PEPPI_1B_A_centroid.mzML.gz"
file2 = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control\\23_04_21_PEPPI_1B_B.raw"
file2 = "C:\\Data\\IsoNN\\23_04_21_PEPPI_1B_B.raw"
# file2 = "Z:\\Group Share\\JGP\\PXD045560\\20220709_nLC1000_E_CEW_Isoform-3_MS2_ET20hcD20_Targetted.raw"

os.chdir(os.path.dirname(file2))

scans =[3398]
#scans = [4002]
#scans = [3642]
scans = range(4000, 4500)
#scans = None

#scans = None
timestart = time.time()
reader = eng.process_file(file2, scans=scans)
print("Time:", time.time() - timestart)
print("Time Per Peak: ", (time.time() - timestart) / len(eng.pks.peaks))
#eng.pks.load_pks()
exit()

for s in reader.scans:
    if scans is not None:
        if s not in scans:
            continue
    # Open the scan and get the spectrum
    try:
        spectrum = reader.grab_centroid_data(s)
        #data = ud.datachop(spectrum, 1014, 1035)
        #np.savetxt("C:\\Data\\IsoNN\\test9.txt", data)
        #exit()
    except:
        print("Error Reading Scan", s)
        continue
    # If the spectrum is too short, skip it
    if len(spectrum) < 3:
        continue

    plot_pks(eng.pks, centroids=spectrum, scan=s)
    plt.show()




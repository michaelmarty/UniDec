import pandas as pd
from pyopenms import *
import os
import fnmatch
import Scripts.JGP.xic_xtractor as xic_xtractor
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.IsoDec.datatools as dt

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

# Set backend to Agg
mpl.use('WxAgg')

existing_library = r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\dda_lipid_library.csv"
library_df= pd.read_csv(existing_library, delimiter=",")
masses = library_df["m/z"]
retention_times = library_df["Retention Time"]
names = library_df["Lipid Name"]
new_rts = [[] for i in range(len(retention_times))]

dda_file_directory = "Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\calibrated_mzML"
os.chdir(dda_file_directory)

dda_files = match_files(dda_file_directory, "*.mzML")

min_peak_intensity = 1000
ppm_tol = 10
rt_width = 2

for file in dda_files:
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    print("Creating indexed file")
    ind_file = xic_xtractor.IndexedFile(exp)

    for i in range(len(masses)):
        mass_tol = (masses[i]/1e6) * ppm_tol
        min_rt = retention_times[i] - (rt_width/2)
        max_rt = retention_times[i] + (rt_width/2)
        xic = ind_file.extract_xic(masses[i], mass_tol, min_rt, max_rt)
        peak_index = np.argmax(xic[:,2])
        peak_intensity = xic[peak_index,2]

        if peak_intensity > min_peak_intensity:
            rt_int = np.transpose(np.array([xic[:, 1]/60, xic[:, 2]]))
            peaks = dt.fastpeakdetect(rt_int, window=10, threshold=0.05)
            peaks = [p for p in peaks if p[1] > min_peak_intensity]
            new_rt = -1
            if len(peaks) > 1:
                delta_rts = [abs(p[0] - retention_times[i]) for p in peaks]
                sorted_indices = np.argsort(delta_rts)
                new_rt = peaks[sorted_indices[0]][0]
            elif len(peaks) == 1:
                new_rt = peaks[0][0]

            if new_rt != -1:
                new_rts[i].append(new_rt)
                # plt.plot(xic[:,1]/60, xic[:,2])
                # plt.title(names[i])
                # plt.vlines(retention_times[i], 0, np.max(xic[:,2]), color='r', linestyles='dashed')
                # plt.vlines(new_rt, 0, np.max(xic[:,2]), color='b', linestyles='dashed')
                # plt.show()

rt_averages = [0.0 for i in range(len(new_rts))]
for i in range(len(new_rts)):
    if(len(new_rts[i])>0):
        rt_averages[i] = float(np.average(new_rts[i]))
    else:
        print(names[i], "not found in any file.")
        rt_averages[i] = retention_times[i]

library_df["Retention Time"] = rt_averages
library_df.to_csv("Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\rtcalibrated_dda_library.csv")



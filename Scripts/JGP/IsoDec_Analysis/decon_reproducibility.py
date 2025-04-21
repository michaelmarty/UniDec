from unidec.IsoDec.match import *
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
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

minscans = 10

# toppath = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\"
# toppath = "D:\\Johnny\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\"
# toppath = "D:\\Johnny\\IsoDec\\PXD019247\\isodec_msalign\\8bit\\"
toppath = "D:\\Johnny\\IsoDec\\PXD019247\\isodec_msalign\\4bit\\"

ms1_msalign_files = match_files(toppath, "*ms1.msalign")
ms1_mcs = []
pks_in_masses = []
totalpks = []

for ms1_msalign_file in ms1_msalign_files:
    mc = read_msalign_to_matchedcollection(toppath + ms1_msalign_file, mz_type=None)
    totalpks.append(len(mc.peaks))
    for peak in mc.peaks:
        mc.add_pk_to_masses2(peak, ppmtol=20, rt_tol=3)

    peaks_in_masses = 0
    mc1_indices_to_remove = []
    for mi, mass in enumerate(mc.masses):
        if len(mass.scans) < minscans:
            mc1_indices_to_remove.append(mi)

        #Make a 2D data array of scans and intensities
        data = list(zip(mass.scan_intensities.keys(), mass.scan_intensities.values()))
        data = np.array(data)
        #Sort the data by scan number
        data = data[data[:, 0].argsort()]
        peaks = dt.fastpeakdetect(data, window=minscans, threshold=0.1)
        if len(peaks) == 0:
            mc1_indices_to_remove.append(mi)


    mc.masses = [mass for i, mass in enumerate(mc.masses) if i not in mc1_indices_to_remove]
    #Get the total number of peaks in the masses
    for mass in mc.masses:
        peaks_in_masses += mass.totalpeaks

    pks_in_masses.append(peaks_in_masses)
    ms1_mcs.append(mc)

    if False:
        #Order the masses by decreasing total intensity
        mc.masses = sorted(mc.masses, key=lambda x: x.totalintensity, reverse=True)
        for mass in mc.masses:
            #Order the scans_intensities dictionary by scan number
            scans_intensities = {k: v for k, v in sorted(mass.scan_intensities.items(), key=lambda item: item[0])}
            plt.plot(list(scans_intensities.keys()), list(scans_intensities.values()))
            plt.xlabel("Scan Number", fontsize=14)
            plt.ylabel("Intensity", fontsize=14)
            #Set the minimum y to 0
            plt.ylim(bottom=0)
            plt.title("LC Trace of Mass " + str(mass.monoiso) + " Da")
            plt.show()

# #Plot the LC trace of each matched mass in the top 10 most intense masses
# for ms1_mc in ms1_mcs:
#         #Sort the masses by intensity
#         ms1_mc.masses = sorted(ms1_mc.masses, key=lambda x: x.totalintensity, reverse=True)
#         for i, mass in enumerate(ms1_mc.masses):
#             if i >= 10:
#                 break
#             #Order the scans_intensities dictionary by scan number
#             scans_intensities = {k: v for k, v in sorted(mass.scan_intensities.items(), key=lambda item: item[0])}
#             plt.plot(list(scans_intensities.keys()), list(scans_intensities.values()))
#             plt.show()

matched_mass_percentages = []
matched_masses = []
matched_intensity_percentages = []
total_masses = []

for i in range(len(ms1_mcs)):
    for j in range(len(ms1_mcs)):
        if j <= i:
            continue
        print("Comparing", ms1_msalign_files[i], "to", ms1_msalign_files[j])
        matched1, matched2, matchedint1, matchedint2, nummatched1, nummatched2, n1, n2 = compare_matchedmasses(ms1_mcs[i], ms1_mcs[j], ppmtol=20, rt_tol=5, maxshift=3, f_shared_zs=0.5)
        matched_mass_percentages.append(matched1)
        matched_intensity_percentages.append(matchedint1)
        matched_mass_percentages.append(matched2)
        matched_intensity_percentages.append(matchedint2)
        matched_masses.append(nummatched1)
        matched_masses.append(nummatched2)
        print("Matched Masses: ", nummatched1, nummatched2)
        total_masses.append(n1)
        total_masses.append(n2)



print("Average % Matched Masses: ", np.mean(matched_mass_percentages))
print("Standard Deviation: ", np.std(matched_mass_percentages))
print("Average Matched Intensities: ", np.mean(matched_intensity_percentages))
print("Standard Deviation: ", np.std(matched_intensity_percentages))
print("Average Total masses: ", np.mean(total_masses))
print("Standard Deviation: ", np.std(total_masses))
print("Average Matched Masses: ", np.mean(matched_masses))
print("Standard Deviation: ", np.std(matched_masses))

percent_pks_in_masses = []
for i in range(len(ms1_mcs)):
    percent_pks_in_masses.append((pks_in_masses[i] / totalpks[i]) * 100)
print("Average % Peaks in Masses: ", np.mean(percent_pks_in_masses))
print("Standard Deviation: ", np.std(percent_pks_in_masses))
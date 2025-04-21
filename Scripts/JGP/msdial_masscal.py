import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("D:\\Johnny\\Lipidomics\\e_coli\\2024-07-18_04_lockspray_Ecoli_DDA_TandE-test.txt", delimiter="\t")

msms_matched = df["MS/MS matched"]
experimental_mz = df["Precursor m/z"]
reference_mz = df["Reference m/z"]
raw_areaas = df["Area"]

mz_errors_Th = []
mz_errors_ppm = []
areas = []

#collect all mz errors and intensities for MS/MS-identified lipids

for i in range(len(experimental_mz)):
    if not math.isnan(reference_mz[i]):
        if msms_matched[i] == True:
            mz_error = experimental_mz[i] - reference_mz[i]
            mz_errors_Th.append(mz_error)
            mz_errors_ppm.append((mz_error / reference_mz[i]) * 1e6)
            areas.append(raw_areaas[i])

mz_errors_Th_array = np.array(mz_errors_Th)
mz_errors_ppm_array = np.array(mz_errors_ppm)
areas_array = np.array(areas)
contributions = []


plt.hist(mz_errors_Th_array, bins=50)
plt.show()

sum_area = areas_array.sum()

for i in range(len(mz_errors_Th_array)):
    contributions.append((areas_array[i] / sum_area) * mz_errors_Th_array[i])

print("Weighted average m/z error:", np.array(contributions).sum())

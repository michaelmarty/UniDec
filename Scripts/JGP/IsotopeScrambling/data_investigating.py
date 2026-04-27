import numpy as np
from Scripts.JGP.IsotopeScrambling.MetaMorpheus_tools import *
import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.modules.isotopetools import *
from unidec.IsoDec.datatools import fastnearest
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

h_mass = 1.007276467

def get_psm_isodist(psm, threshold=0.01):
    isodist_intensities = fft_pep_seq_to_dist(psm.base_sequence)
    mono_mz = (float(psm.monoisotopic_mass) + (psm.precursor_charge * h_mass)) / psm.precursor_charge
    mz_diff = mass_diff_c / psm.precursor_charge

    mzs = [0 for i in range(len(isodist_intensities))]
    for i in range(len(mzs)):
        mzs[i] = mono_mz + (i * mz_diff)

    isodist = np.column_stack((mzs, isodist_intensities))

    mask = isodist[:,1] > threshold
    return isodist[mask]

def scale_isodist(spec, isodist, mz_match_tol=0.01, min_matched_peaks=5):
    matched_experimental_intensities = [0.0 for i in range(len(isodist))]
    matched_theoretical_intensities = [0.0 for i in range(len(isodist))]
    matched_mzs = []
    for i in range(len(isodist)):
        nearest_idx = fastnearest(spec[:,0], isodist[i,0])
        if np.abs(spec[nearest_idx,0] - isodist[i,0]) < mz_match_tol:
            matched_mzs.append(spec[nearest_idx,0])
            matched_experimental_intensities[i] = spec[nearest_idx,1]
            matched_theoretical_intensities[i] = isodist[i,1]
    if len(matched_mzs) < min_matched_peaks:
        return 0
    css = calculate_cosinesimilarity(matched_experimental_intensities, matched_theoretical_intensities)

    isodist[:,1] *= max(matched_experimental_intensities)
    return css

def plot_psm_precursor(psm, reader, low_mz_range=0.5, high_mz_range=4):
    if "|" in psm.monoisotopic_mass or psm.qvalue > 0.01 or "K" not in psm.base_sequence:
        return 0,0

    spec = reader.get_single_scan(psm.ms1_scan_num)

    light_isodist = get_psm_isodist(psm)

    k_count = psm.base_sequence.count('K')
    heavy_mz = psm.precursor_mz + ((k_count * 8) / psm.precursor_charge)

    shift = heavy_mz - light_isodist[0, 0]
    heavy_isodist = light_isodist.copy()
    heavy_isodist[:, 0] += shift

    light_css = scale_isodist(spec, light_isodist)
    heavy_css = scale_isodist(spec, heavy_isodist)

    if not light_css or not heavy_css or heavy_css < 0.7 or light_css < 0.7:
        return 0, 0

    # min_mz = psm.precursor_mz - low_mz_range
    # max_mz = heavy_mz + high_mz_range
    #
    # max_intensity = 0
    #
    #
    #
    # plt.hlines(0, min_mz-0.5, max_mz+0.5, colors='k')
    # for i in range(len(spec)):
    #     if min_mz <= spec[i,0] <= max_mz:
    #         if spec[i,1] > max_intensity:
    #             max_intensity = spec[i,1]
    #         plt.vlines(spec[i,0],0,spec[i,1])
    # plt.title("Light CSS:" + str(round(light_css, 3)) + " Heavy CSS:" + str(round(heavy_css,3)))
    # plt.plot(light_isodist[:,0], light_isodist[:,1])
    # plt.plot(heavy_isodist[:,0], heavy_isodist[:,1])
    # plt.xlim(min_mz, max_mz)
    # plt.ylim(0,max_intensity*1.1)
    # plt.show()

    return light_css, heavy_css


os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover")

psm_filepath = r"2026-02-13-16-36-52\Task4-SearchTask\AllPSMs.psmtsv"
raw_filepath = r"2026-02-13-16-36-52\Task2-AveragingTask\190404_Graft_24d_1-calib-averaged.mzML"

psms = load_psm_file(psm_filepath)
psms.sort(key=lambda p: p.total_ion_current, reverse=True)

reader = ImporterFactory.create_importer(raw_filepath)

pairs = 0
light_css_vals = []
heavy_css_vals = []

for psm in psms:
    try:
        light_css, heavy_css = plot_psm_precursor(psm, reader)
        if light_css != 0 and heavy_css != 0:
            pairs += 1
            light_css_vals.append(light_css)
            heavy_css_vals.append(heavy_css)
    except Exception as e:
        print(e)

print("Total identified pairs:", str(pairs))

plt.hist(light_css_vals, bins=20, label="Light", alpha=0.75)
plt.hist(heavy_css_vals, bins=20, label="Heavy", alpha=0.75)
plt.legend()
plt.show()
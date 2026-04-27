import numpy as np
from Scripts.JGP.IsotopeScrambling.MetaMorpheus_tools import *
import os
from unidec.IsoDec.IsoGen.isogenpep import IsoGenPepEngine
from unidec.UniDecImporter.ImportTools import merge_spectra
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.IsoDec.IsoGen.isogenc import *
import matplotlib.pyplot as plt
from Scripts.JGP.LipidQuant.LipidQuant_engine import *
from unidec.IsoDec.datatools import get_all_centroids, remove_noise_cdata
from unidec.tools import datachop
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *
from unidec.modules.unidecstructure import IsoDecConfig
import unidec.IsoDec.match as match
from typing import List, Tuple
from numba import njit
import re
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
from matplotlib.pyplot import rcParams
import matplotlib.patches as patches

iso_mass_diff = 1.0033
proton_mass = 1.007276467

def validate_filepaths(paths):
    for path in paths:
        if os.path.exists(path):
            continue
        else:
            print("Path does not exist: ", path)


def get_psm_isodist(psm, isogenpep_eng=None, heavy=False, threshold=0.001):
    if not heavy:
        isodist_intensities = fft_pep_seq_to_dist(psm.base_sequence)
    else:
        if isogenpep_eng is None:
            isogenpep_eng = IsoGenPepEngine(isolen=16)
        isodist_intensities = isogenpep_eng.predict(psm.base_sequence)

    if "|" in psm.monoisotopic_mass:
        split = psm.monoisotopic_mass.split('|')
        psm.monoisotopic_mass = float(split[0])

    mono_mz = (float(psm.monoisotopic_mass) + (psm.precursor_charge * proton_mass)) / psm.precursor_charge
    mz_diff = iso_mass_diff / psm.precursor_charge

    mzs = [0 for i in range(len(isodist_intensities))]
    for i in range(len(mzs)):
        mzs[i] = mono_mz + (i * mz_diff)

    isodist = np.column_stack((mzs, isodist_intensities))

    above_threshold = isodist[:, 1] > threshold

    # Find last index where intensity > threshold
    if np.any(above_threshold):
        last_valid = np.where(above_threshold)[0][-1]
        mask = np.zeros(len(isodist), dtype=bool)
        mask[:last_valid + 1] = True
    else:
        mask = np.zeros(len(isodist), dtype=bool)

    return isodist[mask]


def get_averaged_spectrum(psm, reader, min_rt, max_rt, scan_orders):
    scan_nums_to_average = []
    scans_to_average = []

    for i in range(len(reader.scans)):
        if min_rt < reader.times[i] < max_rt and scan_orders[i] == 1:
            scan_nums_to_average.append(reader.scans[i])
            scans_to_average.append(reader.grab_centroid_data(reader.scans[i]))

    min_mz = psm.precursor_mz - 5
    max_mz = psm.precursor_mz + 20

    #averaged_scan = reader.get_avg_scan(scan_range=[scan_nums_to_average[0], scan_nums_to_average[len(scan_nums_to_average) - 1]])
    averaged_scan = merge_spectra(scans_to_average, mzbins=0.005)


    centroided = get_all_centroids(averaged_scan, window=10)

    return centroided, averaged_scan


@njit(fastmath=True)
def find_matched_intensities_indices(spec1_mz: np.ndarray, spec1_intensity: np.ndarray, spec2_mz: np.ndarray,
                             max_shift: int, tolerance: float, z: int, peakmz: float):
    """
    Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from low to high m/z).

    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as numpy array. Peak mz values must be ordered.
    spec2_mz:
        Theoretical isotope peak m/z values as numpy array. Peak mz values must be ordered.
    tolerance:
        Peaks will be considered a match when within [tolerance] ppm of the theoretical value.

    Returns
    -------
    matches
        List containing entries of type (centroid intensity)

    """
    mass_diff_c = 1.0033
    if z == 0:
        z = 1

    query_mzs = [float(0) for i in range(len(spec2_mz) + 2 * max_shift)]

    cent_intensities = np.zeros(len(query_mzs))
    cent_indices = [-1 for i in range(len(query_mzs))]

    mono_idx = max_shift

    for i in range(max_shift + 1):
        if i == 0:
            continue
        else:
            query_mzs[max_shift + len(spec2_mz) - 1 + i] = spec2_mz[len(spec2_mz) - 1] + (i * mass_diff_c) / z
            query_mzs[mono_idx - i] = spec2_mz[0] - (i * mass_diff_c) / z

    for i in range(len(spec2_mz)):
        query_mzs[i + max_shift] = spec2_mz[i]

    diff = (peakmz * tolerance * 1e-6)
    for i in range(len(query_mzs)):
        mz = query_mzs[i]
        low_bound = mz - diff
        high_bound = mz + diff
        for j in range(len(spec1_mz)):
            mz2 = spec1_mz[j]
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                continue
            else:
                value = spec1_intensity[j]
                if value > cent_intensities[i]:
                    cent_intensities[i] = spec1_intensity[j]
                    cent_indices[i] = j

    return cent_intensities, cent_indices


def subtract_matched_isodist(spec, maxshift, bestshift, matched_indices, isodist):
    for i in range(len(isodist)):
        spec_index = matched_indices[i + maxshift + bestshift]
        if spec_index >= 0:
            spec_intensity = spec[spec_index, 1]
            if isodist[i,1] > spec_intensity:
                spec[spec_index, 1] = 0
            else:
                spec[spec_index,1] -= isodist[i,1]


def match_isodist(spec, isodist, psm, min_css=0.9, ppm_match_tol=20, no_shift=False):
    #Limit shift for low charges if necessary
    if not no_shift:
        if psm.precursor_charge < 3:
            maxshift = 1
        elif psm.precursor_charge < 6:
            maxshift = 2
        else:
            maxshift = 3
    else:
        maxshift = 0

    cent_intensities, cent_indices = find_matched_intensities_indices(spec[:, 0], spec[:, 1], isodist[:, 0], maxshift,
                                                tolerance=ppm_match_tol, z=psm.precursor_charge, peakmz=psm.precursor_mz)

    if len(np.nonzero(cent_intensities)[0]) == 0:
        return None, 0, None

    norm_factor = max(cent_intensities) / max(isodist[:, 1])
    isodist[:, 1] *= norm_factor

    shiftrange = np.arange(-maxshift, maxshift + 1)

    shifts_scores = np.zeros((len(shiftrange), 2))

    bestshift = -1
    sum = -1
    for i, shift in enumerate(shiftrange):
        s = match.calculate_cosinesimilarity(cent_intensities, isodist[:, 1], shift, maxshift, minusoneaszero=False)
        shifts_scores[i, 0] = shift
        shifts_scores[i, 1] = s
        if s > sum:
            sum = s
            bestshift = shift

    isodist[:,0] += bestshift * (iso_mass_diff / psm.precursor_charge)

    if np.amax(shifts_scores[:,1]) > min_css:
        subtract_matched_isodist(spec, maxshift, bestshift, cent_indices, isodist)

    if maxshift == bestshift:
        matched_intensities = cent_intensities[maxshift+bestshift:]
    else:
        matched_intensities = cent_intensities[maxshift+bestshift:len(isodist)+bestshift+1]

    norm_factor = max(matched_intensities) / max(isodist[:, 1])
    isodist[:, 1] *= norm_factor
    return bestshift, np.amax(shifts_scores[:,1]), matched_intensities


def has_isotope_shape(matched_intensities,
                      isodist_intensities,
                      min_peak_height=0.0):
    diff = np.asarray((matched_intensities/max(matched_intensities)) - (isodist_intensities/max(isodist_intensities)))

    # 1. Must start negative
    if diff[0] >= 0:
        #print("shape does not start negative")
        return False

    # 2. Must have a positive maximum
    peak_idx = np.argmax(diff)
    peak_val = diff[peak_idx]

    if peak_val <= min_peak_height:
        #print("peak height not positive")
        return False

    return True


def extract_dist_for_psm(psm, reader, scan_orders, lipid_quant_config, isogenpep_eng=None):
    #This is the isodist adjusted to match the position in m/z space. No intensity scaling or shifting performed
    light_isodist = get_psm_isodist(psm)

    min_mz = psm.precursor_mz - 10
    max_mz = psm.precursor_mz + 20

    precursor_spectrum = datachop(reader.grab_centroid_data(psm.ms1_scan_num), min_mz, max_mz)
    precursor_spectrum_copy = precursor_spectrum.copy()

    bestshift, light_css, light_matched_intensities = match_isodist(precursor_spectrum, light_isodist, psm)

    top_n = 2
    indices = np.argsort(light_isodist[:, 1])[-top_n:][::-1]
    most_abundant_mzs = light_isodist[indices, 0]

    rt_tol = 2
    min_rt = psm.retention_time - rt_tol
    min_scan_idx = fastnearest(reader.times, min_rt)
    max_rt = psm.retention_time + rt_tol
    max_scan_idx = fastnearest(reader.times, max_rt)

    mask = (np.array(scan_orders[min_scan_idx:max_scan_idx+1]) == 1)

    ppm_tol = 20
    mz_tol = (most_abundant_mzs[0] / 1000000) * ppm_tol

    eics = []
    for i in range(len(most_abundant_mzs)):
        eics.append(reader.get_eic(most_abundant_mzs[i], mz_tol=mz_tol, rt_range=[min_rt, max_rt])[mask])

    int_indices = get_integration_indices(eics[0], eics[1], psm.retention_time, lipid_quant_config)

    averaged_centroided =None
    profile = None
    if int_indices is None:
        print("No averaging window found.")
        averaged_centroided = precursor_spectrum_copy
    else:
        min_avgd_rt = eics[0][int_indices[0], 0]
        max_avgd_rt = eics[0][int_indices[len(int_indices) - 1], 0]
        averaged_centroided, profile = get_averaged_spectrum(psm, reader, min_avgd_rt, max_avgd_rt, scan_orders)
        averaged_centroided = datachop(averaged_centroided, min_mz, max_mz)
        profile = datachop(profile, min_mz, max_mz)

    averaged_centroided_copy = averaged_centroided.copy()

    bestshift, light_css, light_matched_intensities = match_isodist(averaged_centroided, light_isodist, psm)


    if isogenpep_eng is None:
        heavy_isodist = light_isodist.copy()
    else:
        heavy_isodist = get_psm_isodist(psm, heavy=True, isogenpep_eng=isogenpep_eng)

    k_count = psm.base_sequence.count("K")
    k_shift = (k_count * 8.0142) / psm.precursor_charge
    heavy_isodist[:,0] += k_shift

    heavy_bestshift, heavy_css, heavy_matched_intensities = match_isodist(averaged_centroided, heavy_isodist, psm, min_css=0.5, no_shift=True)

    heavy_mono_mz = ((float(psm.monoisotopic_mass) + (psm.precursor_charge * proton_mass)) / psm.precursor_charge) + k_shift

    if heavy_bestshift is None:
        return None



    # for p in averaged_centroided_copy:
    #     plt.vlines(p[0], 0, p[1], color='k')
    # # if profile is not None:
    # #     plt.plot(profile[:,0], profile[:,1], color='b')
    # plt.vlines(heavy_mono_mz, 0, max(averaged_centroided_copy[:,1]), color='r')
    # plt.hlines(0, min(averaged_centroided_copy[:,0]), max(averaged_centroided_copy[:,0]), color='k')
    # plt.plot(light_isodist[:,0], light_isodist[:,1], color='b')
    # plt.ylim(bottom=0)
    # plt.show()

    heavy_match_count = np.count_nonzero(np.array(heavy_matched_intensities))
    matched_frac = heavy_match_count / len(heavy_matched_intensities)


    if light_css < 0.8:
        #print("Failed light css")
        return None
    if heavy_css < 0.6:
        #print("Failed Heavy css")
        return None
    if matched_frac < 0.6:
        #print("Failed matched frac")
        return None
    if not has_isotope_shape(heavy_matched_intensities, heavy_isodist[:,1]):
        return None

    print("file: ", psm.filename)
    print("scan: ", str(psm.ms1_scan_num))

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(5, 5), constrained_layout=True)
    axs[0].plot(heavy_matched_intensities - heavy_isodist[:,1])
    axs[0].hlines(0, 0, 8, colors='k')

    for i in range(len(averaged_centroided_copy)):
        axs[1].vlines(averaged_centroided_copy[i,0], 0, averaged_centroided_copy[i,1])
    axs[1].plot(light_isodist[:,0], light_isodist[:,1], color='r')
    axs[1].plot(heavy_isodist[:,0], heavy_isodist[:,1], color='orange')
    axs[1].vlines(most_abundant_mzs[0], 0, max(light_matched_intensities), colors='r', linestyle='--')
    plt.show()

    return heavy_matched_intensities / max(heavy_matched_intensities)


def process_file_set(psm_filepaths, raw_filepaths, lipid_quant_config, isogenpep_eng=None):
    all_dists = []
    all_seqs = []
    all_full_seqs = []

    for i in range(len(psm_filepaths)):
        filename = os.path.splitext(raw_filepaths[i])[0]
        print("Processing ", filename)

        psms = load_psm_file(psm_filepaths[i])

        reader = ImporterFactory.create_importer(raw_filepaths[i])
        scan_orders = [1 for i in range(len(reader.scans))]
        for i in range(len(reader.scans)):
            scan_orders[i] = reader.get_ms_order(reader.scans[i])

        heavy_dists = []
        seqs = []
        full_seqs = []

        psms_to_process = []
        for psm in psms:
            if 'K' in psm.base_sequence and psm.qvalue < 0.01 and psm.target_decoy_contaminant == "T" and "Homo" in psm.organism:
                psms_to_process.append(psm)

        for i in range(len(psms_to_process)):
            if i % 10 == 0:
                print("Processing peptide " + str(i) + " of " + str(len(psms_to_process)))
            psm = psms_to_process[i]
            heavy_dist = extract_dist_for_psm(psm, reader, scan_orders, config, isogenpep_eng=isogenpep_eng)
            if heavy_dist is not None:
                heavy_dists.append(heavy_dist)
                seqs.append(psm.base_sequence)
                full_seqs.append(psm.full_sequence)

        print("Found " + str(len(heavy_dists)) + " training dists.")
        # Find max length
        max_len = 16

        # Create padded array
        padded_dists = np.zeros((len(heavy_dists), max_len))
        all_dists.extend(heavy_dists)
        all_seqs.extend(seqs)
        all_full_seqs.extend(full_seqs)

        for i, d in enumerate(heavy_dists):
            padded_dists[i, :len(d)] = d
        np.savez_compressed(filename + "_transfer_dists.npz", dists=padded_dists, seqs=seqs, full_seqs=full_seqs)

    return all_dists, all_seqs, all_full_seqs

def add_seqs_to_unique(existing_unique, to_add):
    for seq in to_add:
        if seq not in existing_unique:
            existing_unique.append(seq)
    return

def load_results_from_npzs(raw_filepaths):
    dists = []
    seqs = []
    full_seqs = []
    for raw_file in raw_filepaths:
        filename = os.path.splitext(raw_file)[0]
        data = np.load(filename + "_transfer_dists.npz")
        dists.extend(data['dists'])
        seqs.extend(data['seqs'])
        full_seqs.extend(data['full_seqs'])

    return dists, seqs, full_seqs

def convert_fullseq_to_baseseq(sequence):
    # Remove brackets and everything inside
    clean_text = re.sub(r"\[.*?\]", "", sequence)
    return clean_text


def pad_distributions(dists):
    # Find the maximum length
    max_len = max(len(d) for d in dists)

    # Pad each distribution with zeros at the end
    padded_dists = []
    for d in dists:
        arr = np.array(d)
        if len(arr) < max_len:
            arr = np.pad(arr, (0, max_len - len(arr)), mode='constant')
        padded_dists.append(arr)

    return np.array(padded_dists)  # shape: (n_distributions, max_len)


def average_dists(filtered_dists):
    padded_dists = pad_distributions(filtered_dists)

    mean_dist = padded_dists.mean(axis=0)
    return mean_dist


def find_dist_pairs(d12_results, d24_results, d36_results, d48_results, d60_results):
    unique_full_seqs = []

    # d12_fullseqs = d12_results[2]
    # d24_fullseqs = d24_results[2]
    # d36_fullseqs = d36_results[2]
    # d48_fullseqs = d48_results[2]
    # d60_fullseqs = d60_results[2]

    d12_fullseqs = d12_results[1]
    d24_fullseqs = d24_results[1]
    d36_fullseqs = d36_results[1]
    d48_fullseqs = d48_results[1]
    d60_fullseqs = d60_results[1]

    add_seqs_to_unique(unique_full_seqs, d12_fullseqs)
    add_seqs_to_unique(unique_full_seqs, d24_fullseqs)
    add_seqs_to_unique(unique_full_seqs, d36_fullseqs)
    add_seqs_to_unique(unique_full_seqs, d48_fullseqs)
    add_seqs_to_unique(unique_full_seqs, d60_fullseqs)

    print("Found " + str(len(unique_full_seqs)) + " unique base sequences")

    base_seqs = []
    mean_dists = []

    bigseq = ''

    for fullseq in unique_full_seqs:
        d12_indices = [i for i, val in enumerate(d12_fullseqs) if val == fullseq]
        d24_indices = [i for i, val in enumerate(d24_fullseqs) if val == fullseq]
        d36_indices = [i for i, val in enumerate(d36_fullseqs) if val == fullseq]
        d48_indices = [i for i, val in enumerate(d48_fullseqs) if val == fullseq]
        d60_indices = [i for i, val in enumerate(d60_fullseqs) if val == fullseq]

        base_seq = convert_fullseq_to_baseseq(fullseq)
        print("Base seq: ", base_seq)
        theoretical_dist = fft_pep_seq_to_dist(base_seq, isolen=8)

        dists = []

        for i, d12_ind in enumerate(d12_indices):
            dists.append(d12_results[0][d12_ind])
            #plt.plot(d12_results[0][d12_ind], label="d12_" + str(i))

        for i, d24_ind in enumerate(d24_indices):
            dists.append(d24_results[0][d24_ind])
            #plt.plot(d24_results[0][d24_ind], label="d24_" + str(i))

        for i, d36_ind in enumerate(d36_indices):
            dists.append(d36_results[0][d36_ind])
            #plt.plot(d36_results[0][d36_ind], label="d36_" + str(i))

        for i, d48_ind in enumerate(d48_indices):
            dists.append(d48_results[0][d48_ind])
            #plt.plot(d48_results[0][d48_ind], label="d48_" + str(i))

        for i, d60_ind in enumerate(d60_indices):
            dists.append(d60_results[0][d60_ind])
            #plt.plot(d60_results[0][d60_ind], label="d60_" + str(i))

        # plt.plot(theoretical_dist, color='k', linewidth=2, label="THEORETICAL")
        # plt.legend()
        # plt.show()

        if len(dists) < 1:
            continue
        dists = np.array(dists)

        css_matrix = cosine_similarity(dists)
        avg_similarity = css_matrix.mean(axis=1)
        median_dist = np.median(avg_similarity)
        mad = np.median(np.abs(avg_similarity - median_dist))
        threshold = median_dist - (3 * mad)  # tweak 3 as needed

        mask = avg_similarity >= threshold
        oppo_mask = avg_similarity <= threshold

        filtered_distributions = dists[mask]
        mean_dist = average_dists(filtered_distributions)
        bad_dists = dists[oppo_mask]

        print("Removed " + str(len(dists) - len(filtered_distributions)) + " distributions.")

        # for filtered_dist in filtered_distributions:
        #     plt.plot(filtered_dist, color='b')
        # for bad_dist in bad_dists:
        #     plt.plot(bad_dist, color='r', linestyle='--')
        # plt.plot(theoretical_dist, color='k', linewidth=2, label="THEORETICAL")
        # plt.plot(mean_dist, color='b', linewidth=2, linestyle='--')
        # plt.legend()
        # plt.show()

        base_seqs.append(fullseq)
        mean_dists.append(mean_dist)

        bigseq = bigseq + fullseq


    np.savez_compressed("transfer_learning_dists_" + str(len(base_seqs)) + ".npz",
                        seqs=base_seqs, dists=mean_dists)

    aa_counts = Counter(bigseq)
    aas = []
    counts = []
    for k,v in aa_counts.items():
        aas.append(k)
        counts.append(v)

    plt.bar(aas, counts)
    plt.xlabel('AA')
    plt.ylabel('Count')
    plt.show()


def make_figure_for_psm(psm, reader, scan_orders, lipid_quant_config, engine):
    #This is the isodist adjusted to match the position in m/z space. No intensity scaling or shifting performed
    light_isodist = get_psm_isodist(psm)

    min_mz = psm.precursor_mz - 2
    max_mz = psm.precursor_mz + 10

    precursor_spectrum = datachop(reader.grab_centroid_data(psm.ms1_scan_num), min_mz, max_mz)
    precursor_spectrum_copy = precursor_spectrum.copy()

    bestshift, light_css, light_matched_intensities = match_isodist(precursor_spectrum, light_isodist, psm)

    top_n = 2
    indices = np.argsort(light_isodist[:, 1])[-top_n:][::-1]
    most_abundant_mzs = light_isodist[indices, 0]

    rt_tol = 2
    min_rt = psm.retention_time - rt_tol
    min_scan_idx = fastnearest(reader.times, min_rt)
    max_rt = psm.retention_time + rt_tol
    max_scan_idx = fastnearest(reader.times, max_rt)

    mask = (np.array(scan_orders[min_scan_idx:max_scan_idx+1]) == 1)

    ppm_tol = 20
    mz_tol = (most_abundant_mzs[0] / 1000000) * ppm_tol

    eics = []
    for i in range(len(most_abundant_mzs)):
        eics.append(reader.get_eic(most_abundant_mzs[i], mz_tol=mz_tol, rt_range=[min_rt, max_rt])[mask])

    int_indices = get_integration_indices(eics[0], eics[1], psm.retention_time, lipid_quant_config)

    averaged_centroided =None
    profile = None
    if int_indices is None:
        print("No averaging window found.")
        averaged_centroided = precursor_spectrum_copy
    else:
        min_avgd_rt = eics[0][int_indices[0], 0]
        max_avgd_rt = eics[0][int_indices[len(int_indices) - 1], 0]
        averaged_centroided, profile = get_averaged_spectrum(psm, reader, min_avgd_rt, max_avgd_rt, scan_orders)
        averaged_centroided = datachop(averaged_centroided, min_mz, max_mz)
        profile = datachop(profile, min_mz, max_mz)

    averaged_centroided_copy = averaged_centroided.copy()
    averaged_centroided_copy1 = averaged_centroided.copy()

    bestshift, light_css, light_matched_intensities = match_isodist(averaged_centroided, light_isodist, psm)


    if engine is None:
        heavy_isodist = light_isodist.copy()
    else:
        heavy_isodist = get_psm_isodist(psm, heavy=True, isogenpep_eng=engine)

    original_heavy = light_isodist.copy()

    k_count = psm.base_sequence.count("K")
    k_shift = (k_count * 8.0142) / psm.precursor_charge
    heavy_isodist[:,0] += k_shift
    original_heavy[:,0] += k_shift



    heavy_bestshift, heavy_css, heavy_matched_intensities = match_isodist(averaged_centroided, heavy_isodist, psm, min_css=0.5, no_shift=True)
    og_heavy_bestshift, og_heavy_css, og_heavy_matched_intensities = match_isodist(averaged_centroided_copy1, original_heavy, psm, min_css=0.5, no_shift=True)

    heavy_mono_mz = ((float(psm.monoisotopic_mass) + (psm.precursor_charge * proton_mass)) / psm.precursor_charge) + k_shift

    # Set backend to Agg
    mpl.use('WxAgg')
    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'

    tick_size = 16
    label_size = 20
    linewidth = 3

    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                      '#f781bf', '#a65628', '#984ea3',
                      '#999999', '#e41a1c', '#dede00']

    save_figs = True

    if True:
        light_rect = patches.Rectangle((1001.5, 0), 3, 0.6e7, linewidth=0, facecolor=CB_color_cycle[0], alpha=0.3)
        heavy_rect = patches.Rectangle((1006, 0), 6, 8e7, linewidth=0, facecolor=CB_color_cycle[1], alpha=0.3)

        fig, ax = plt.subplots(figsize=(6,3.5))
        ax.add_patch(light_rect)
        ax.add_patch(heavy_rect)
        for i in range(len(averaged_centroided_copy)):
            plt.vlines(averaged_centroided_copy[i,0], 0, averaged_centroided_copy[i,1], color='k', linewidth=2)
        plt.xlim(1000, 1012)
        plt.ylim(0,8e7)
        plt.xticks([1000, 1002, 1004, 1006, 1008, 1010, 1012], fontsize=tick_size)
        plt.yticks(ticks=[0,4e7,8e7], labels=["0","4e7","8e7"], fontsize=tick_size)
        plt.xlabel('m/z', fontstyle='italic', fontsize=label_size)
        plt.ylabel("Intensity", fontsize=label_size)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.xticks(fontsize=tick_size)
        plt.yticks(fontsize=tick_size)
        plt.tick_params(length=5, width=2)
        plt.tight_layout()
        if save_figs:
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_spectrum.pdf", dpi=300, transparent=True)
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_spectrum.png", dpi=300, transparent=True)
        plt.show()

    if False:
        light_min_mz = 1001.5
        light_max_mz = 1004.5

        fig,ax = plt.subplots(figsize=(3,3))
        for i in range(len(averaged_centroided_copy)):
            plt.vlines(averaged_centroided_copy[i,0], 0, averaged_centroided_copy[i,1], linewidth=2, color='k')
        plt.plot(light_isodist[:, 0], light_isodist[:, 1], color=CB_color_cycle[0], linewidth=2)
        plt.ylim(0, 0.6e7)
        plt.xticks([1001, 1002, 1003, 1004], fontsize=tick_size)
        plt.yticks(ticks=[0, 3e6, 6e6], labels=["0","3e6","6e6"], fontsize=tick_size)
        plt.xlim(light_min_mz, light_max_mz)
        plt.xlabel('m/z', fontstyle='italic', fontsize=label_size)
        plt.ylabel("Intensity", fontsize=label_size)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.tick_params(length=5, width=2)
        plt.tight_layout()
        if save_figs:
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_light.pdf", dpi=300, transparent=True)
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_light.png", dpi=300, transparent=True)
        plt.show()

    if False:
        heavy_min_mz = 1006
        heavy_max_mz = 1012

        fig,ax = plt.subplots(figsize=(3,3))
        for i in range(len(averaged_centroided_copy)):
            plt.vlines(averaged_centroided_copy[i,0], 0, averaged_centroided_copy[i,1], linewidth=2, color='k')
        plt.plot(heavy_isodist[:, 0], heavy_isodist[:, 1], color=CB_color_cycle[1], linewidth=2)
        plt.plot(original_heavy[:,0], original_heavy[:,1], color=CB_color_cycle[2], linestyle='--', linewidth=2)
        plt.ylim(0, 8e7)
        plt.xticks([1006, 1008, 1010, 1012], fontsize=tick_size)
        plt.yticks(ticks=[0, 4e7, 8e7], labels=["0","4e7","8e7"], fontsize=tick_size)
        plt.xlim(heavy_min_mz, heavy_max_mz)
        plt.xlabel('m/z', fontstyle='italic', fontsize=label_size)
        plt.ylabel("Intensity", fontsize=label_size)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.tick_params(length=5, width=2)
        plt.tight_layout()
        if save_figs:
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_heavy.pdf", dpi=300, transparent=True)
            plt.savefig(r"C:\Users\Admin\OneDrive\Desktop\martylab\IsoGen\Figures\scrambled_dist_heavy.png", dpi=300, transparent=True)
        plt.show()

    return


if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover")

    d12_psm_filepaths = [r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_12d_1-calib_Peptides.psmtsv",
                     r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_12d_2-calib_Peptides.psmtsv",
                    r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_12d_3-calib_Peptides.psmtsv",
                     r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_12d_4-calib_Peptides.psmtsv"]

    d12_raw_filepaths = [r"190404_Graft_12d_1.raw",
                         r"190404_Graft_12d_2.raw",
                         r"190404_Graft_12d_3.raw",
                         r"190404_Graft_12d_4.raw"]


    d24_psm_filepaths = [r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_24d_1-calib_Peptides.psmtsv",
                     r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_24d_2-calib_Peptides.psmtsv",
                    r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_24d_3-calib_Peptides.psmtsv",
                     r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_24d_4-calib_Peptides.psmtsv",
                     r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_24d_5-calib_Peptides.psmtsv"]

    d24_raw_filepaths = [r"190404_Graft_24d_1.raw",
                         r"190404_Graft_24d_2.raw",
                         r"190405_Graft_24d_3.raw",
                         r"190405_Graft_24d_4.raw",
                         r"190405_Graft_24d_5.raw"]


    d36_psm_filepaths = [
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_36d_1-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_36d_2-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_36d_3-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_36d_4-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190405_Graft_36d_5-calib_Peptides.psmtsv"]

    d36_raw_filepaths = [r"190405_Graft_36d_1.raw",
                         r"190405_Graft_36d_2.raw",
                         r"190405_Graft_36d_3.raw",
                         r"190405_Graft_36d_4.raw",
                         r"190405_Graft_36d_5.raw"]


    d48_psm_filepaths = [
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d48_1-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d48_2-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d48_3_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d48_4-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d48_5_Peptides.psmtsv"]

    d48_raw_filepaths = [r"190406_Graft_d48_1.raw",
                         r"190406_Graft_d48_2.raw",
                         r"190406_Graft_d48_3.raw",
                         r"190406_Graft_d48_4.raw",
                         r"190406_Graft_d48_5.raw"]


    d60_psm_filepaths = [
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d60_1_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d60_2-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d60_3-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d60_4-calib_Peptides.psmtsv",
        r"FullLightSearch\Task3-SearchTask\Individual File Results\190406_Graft_d60_5-calib_Peptides.psmtsv"]

    d60_raw_filepaths = [r"190406_Graft_d60_1.raw",
                         r"190406_Graft_d60_2.raw",
                         r"190406_Graft_d60_3.raw",
                         r"190406_Graft_d60_4.raw",
                         r"190406_Graft_d60_5.raw"]


    config = LipidQuantConfig()
    config.isorange = [0.1,4]
    config.max_increasing = 5
    config.min_intensity_ratio = 0.4
    config.min_consec_peaks = 1

    #isogenpep_eng = IsoGenPepEngine(isolen=16)
    isogenpep_eng = None

    if True:
        eng = IsoGenPepEngine(isolen=16)

        psm_file = r"FullLightSearch\Task3-SearchTask\Individual File Results\190404_Graft_24d_1-calib_Peptides.psmtsv"
        raw_file = r"190404_Graft_24d_1.raw"

        psms = load_psm_file(psm_file)
        reader = ImporterFactory.create_importer(raw_file)
        scan_orders = [1 for i in range(len(reader.scans))]
        for i in range(len(reader.scans)):
            scan_orders[i] = reader.get_ms_order(reader.scans[i])

        psm = None
        for p in psms:
            if p.ms1_scan_num == 21506:
                psm = p
                break

        print(psm.full_sequence)
        exit()
        make_figure_for_psm(psm, reader, scan_orders, config, eng)
        exit()



    #d12_results = process_file_set(d12_psm_filepaths, d12_raw_filepaths, config, isogenpep_eng)
    d24_results = process_file_set(d24_psm_filepaths, d24_raw_filepaths, config, isogenpep_eng)
    d36_results = process_file_set(d36_psm_filepaths, d36_raw_filepaths, config, isogenpep_eng)
    d48_results = process_file_set(d48_psm_filepaths, d48_raw_filepaths, config, isogenpep_eng)
    d60_results = process_file_set(d60_psm_filepaths, d60_raw_filepaths, config, isogenpep_eng)
    exit()

    d12_results = load_results_from_npzs(d12_raw_filepaths)
    d24_results = load_results_from_npzs(d24_raw_filepaths)
    d36_results = load_results_from_npzs(d36_raw_filepaths)
    d48_results = load_results_from_npzs(d48_raw_filepaths)
    d60_results = load_results_from_npzs(d60_raw_filepaths)

    find_dist_pairs(d12_results, d24_results, d36_results, d48_results, d60_results)


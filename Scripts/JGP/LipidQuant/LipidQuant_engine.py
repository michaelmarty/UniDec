import itertools
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.IsoDec.datatools as dt
from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol, fastpeakdetect
import numpy as np
import time
import os
import fnmatch
from scipy.stats import linregress
from UniDecImporter import ImporterFactory
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
import builtins
from Scripts.JGP.LipidQuant.quantification_results_plotting import *
from numba import njit

#This stops it from printing that it can't find a
#particular term when reading the mzML files.
real_print = builtins.print

def silent_print(*args, **kwargs):
    if not args or not str(args[0]).startswith("term not found"):
        real_print(*args, **kwargs)

builtins.print = silent_print


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

class LipidQuantConfig:
    def __init__(self):
        #XIC creation parameters
        self.ppm_tol = 50
        self.rt_tol = 4
        self.iso_count = 2

        #XIC Filtering parameters
        self.min_peak_intensity = 500
        self.noise_threshold = 100
        self.min_consec_peaks = 3
        self.isorange = [0.1,0.9]
        self.r_thresh = 0.5
        self.max_increasing = 3
        self.min_intensity_ratio = 0.01

        #Normalization parameters
        self.normalize_to_int_stds = True
        self.normalize_to_total = True

class LibraryLipid:
    def __init__(self, name, mass, rt, lipid_class, internal_std, adduct=None, formula=None):
        self.name = name
        self.mass = mass
        self.rt = rt
        self.calibrated_rt = rt
        self.lipid_class = lipid_class
        self.formula = formula
        self.adduct = adduct
        self.internal_std = internal_std == "Y"

class LipidLibrary:
    def __init__(self, filepath, config):
        self.filepath = filepath
        self.lipids, self.lipid_indices = self.parse_library(filepath)
        self.lipid_groups = self.group_lipids_by_rt_mz(config.ppm_tol, config.rt_tol)

    def parse_library(self, filepath):
        lipids = np.array([])
        lipid_indices = {}

        library_df = pd.read_csv(filepath, delimiter=",")
        rts = library_df["Retention Time"]
        names = library_df["Lipid name"]
        adducts = library_df["Adduct"]
        mzs = library_df["m/z"]
        formulas = library_df["Chemical Formula"]
        classes = library_df["Lipid class"]
        standards = library_df["Standard"]

        lipid_index = 0
        for i in range(len(rts)):
            if names[i] in lipid_indices:
                continue
            lipids = np.append(lipids, LibraryLipid(names[i], mzs[i], rts[i], classes[i], standards[i], adducts[i], formulas[i]))
            lipid_indices[names[i]] = lipid_index
            lipid_index += 1

        return lipids, lipid_indices

    def group_lipids_by_rt_mz(self, ppm_tol, rt_tol):
        groups = []

        lipid_masses = np.array([l.mass for l in self.lipids])

        sorted_indices = np.argsort(lipid_masses)
        lipid_masses = lipid_masses[sorted_indices]
        ordered_lipids = self.lipids[sorted_indices]

        grouped_lipids = []
        #Find cases where
        for i in range(len(ordered_lipids)):
            if i in grouped_lipids:
                continue
            mass_tol = (ordered_lipids[i].mass/1e6)*ppm_tol
            within_tol = fastwithin_abstol(lipid_masses, ordered_lipids[i].mass, mass_tol)
            current_group = [ordered_lipids[i]]
            grouped_lipids.append(i)
            for j in range(len(within_tol)):
                if i == within_tol[j] or within_tol[j] in grouped_lipids:
                    continue
                else:
                    delta_rt = abs(ordered_lipids[within_tol[j]].rt - ordered_lipids[i].rt)
                    if abs(ordered_lipids[within_tol[j]].rt - ordered_lipids[i].rt) < rt_tol:
                        current_group.append(ordered_lipids[within_tol[j]])
                        grouped_lipids.append(within_tol[j])

            if len(current_group) == 1:
                groups.append(current_group)
            else:
                filtered_group = []
                names = set()
                for l in current_group:
                    name = l.name
                    if name not in names:
                        names.add(name)
                        filtered_group.append(l)
                groups.append(filtered_group)

        return groups

class SpectraFile:
    def __init__(self, filepath, sample_num, biorep, techrep, samplevol, stdvol, load_importer=True):
        self.filepath = filepath

        if load_importer:
            self.importer = ImporterFactory.create_importer(filepath)
        else:
            self.importer = None
        self.sample_num = sample_num
        self.biorep = biorep
        self.techrep = techrep
        self.fileid = str("S" + str(sample_num) + "_BR" + str(biorep) + "_TR" + str(techrep))

        self.samplevol = float(samplevol)
        self.stdvol = float(stdvol)

def create_spectrafile(filepath):
    return SpectraFile(filepath, 0, 0)

def calibrate_library_rts(file, library, config):
    is_rts = np.array([])
    delta_rts = np.array([])

    for i in range(len(library.lipids)):
        if library.lipids[i].internal_std:
            mz = library.lipids[i].mass
            abs_width = (mz/1e6)*config.ppm_tol
            min_rt = library.lipids[i].rt - (config.rt_tol/2)
            max_rt = library.lipids[i].rt + (config.rt_tol/2)

            is_xic = file.importer.get_eic(mz, abs_width, [min_rt, max_rt])
            max_index = np.argmax(is_xic[:,1])
            if is_xic[max_index,1] > config.min_peak_intensity:
                is_rts = np.append(is_rts,library.lipids[i].rt)
                rt_error = (is_xic[max_index,0]) - library.lipids[i].rt
                delta_rts = np.append(delta_rts, rt_error)
                library.lipids[i].calibrated_rt =is_xic[max_index,0]

    is_rts = np.append(is_rts,0)
    delta_rts = np.append(delta_rts, 0)

    is_rts = np.append(is_rts,25)
    delta_rts = np.append(delta_rts, 0)

    sorted_indices = np.argsort(is_rts)
    is_rts = is_rts[sorted_indices]
    delta_rts = delta_rts[sorted_indices]
    print("Delta RTs:", delta_rts)

    for i in range(len(library.lipids)):
        if library.lipids[i].internal_std:
            continue
        else:
            rt = library.lipids[i].rt
            rt_correction = np.interp(rt, is_rts, delta_rts)
            library.lipids[i].calibrated_rt = rt + rt_correction

def check_xics(xics, config):
    #First check if the mono xic meets intensity threshold
    if np.max(xics[0][:,1]) < config.min_peak_intensity:
        return False

    #Next check if the mono xic meets the consecutive scans threshold
    passes_consec = False
    curr_consec = 0
    for i in range(len(xics[0][:,1])):
        if xics[0][i,1] > config.noise_threshold:
            curr_consec += 1
            if curr_consec >= config.min_consec_peaks:
                passes_consec = True
                break
        else:
            curr_consec = 0

    if not passes_consec:
        return False

    return True

def score_xics(intensities_mono, intensities_iso1):
    m,b,r,p,err = linregress(intensities_mono, intensities_iso1)

    # equation = f"y = {m:.2f}x + {b:.2f}"
    # r_squared = f"R² = {r ** 2:.3f}"
    #
    # plt.scatter(intensities_mono, intensities_iso1, label="Data")
    # plt.plot(intensities_mono, m*intensities_mono + b, label="Linear Fit")
    # plt.text(0.05, 0.95, f"{equation}\n{r_squared}",
    #          transform=plt.gca().transAxes,
    #          fontsize=12, verticalalignment='top',
    #          bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))
    # plt.legend()
    # plt.show()
    return m, r

def get_integration_indices(xic_mono, xic_iso1, rt, config):
    rt_int_mono = np.transpose(np.array([xic_mono[:,0], xic_mono[:,1]]))

    peaks_mono = dt.fastpeakdetect(rt_int_mono, window=3, threshold=0.05)
    possible_integration_indexes = np.empty(len(peaks_mono), dtype=object)

    for i in range(len(peaks_mono)):
        possible_integration_indexes[i] = check_valid_peak(xic_mono, xic_iso1, peaks_mono[i], config)

    delta_rts = [abs(p[0] - rt) for p in peaks_mono]
    sorted_indices = np.argsort(delta_rts).astype(int)
    peaks_mono = peaks_mono[sorted_indices]
    possible_integration_indexes = possible_integration_indexes[sorted_indices]

    for i in range(len(peaks_mono)):
        if len(possible_integration_indexes[i]) > 0:
            return possible_integration_indexes[i]
        else:
            continue
    return None

def check_valid_peak(xic_mono, xic_iso1, peak, config):
    nearest_idx = fastnearest(np.array(xic_mono[:, 0]), peak[0])
    peak_intensity = peak[1]
    adding_peaks_left = True
    adding_peaks_right = True

    left_integration_indices = []
    right_integration_indices = []

    # we expect intensity to go down, so if it increases for more than N datapoints, cutoff before the first increase
    previous_left_intensity = peak_intensity
    previous_right_intensity = peak_intensity

    current_left_increasing = 0
    current_right_increasing = 0

    mono_int = xic_mono[:,1]

    intensity_threshold = peak_intensity * config.min_intensity_ratio

    for j in range(1, len(xic_mono)):
        left_idx = nearest_idx - j
        right_idx = nearest_idx + j

        if adding_peaks_left:
            if left_idx >= 0 and mono_int[left_idx] > intensity_threshold:
                left_integration_indices.append(left_idx)

                # Check for poorly resolved peaks
                if mono_int[left_idx] > previous_left_intensity:
                    current_left_increasing += 1
                    if current_left_increasing >= config.max_increasing:
                        # we remove the recent increasing values
                        for _ in range(current_left_increasing):
                            left_integration_indices.pop()
                        adding_peaks_left = False
                else:
                    current_left_increasing = 0
            else:
                if left_idx >= 0:
                    left_integration_indices.append(left_idx)
                adding_peaks_left = False

        if adding_peaks_right:
            if right_idx < len(xic_mono) and mono_int[right_idx] > intensity_threshold:
                right_integration_indices.append(right_idx)

                if mono_int[right_idx] > previous_right_intensity:
                    current_right_increasing += 1

                    if current_right_increasing >= config.max_increasing:
                        for _ in range(current_right_increasing):
                            right_integration_indices.pop()
                        adding_peaks_right = False
                else:
                    current_right_increasing = 0

            else:
                if right_idx < len(xic_mono):
                    right_integration_indices.append(right_idx)
                adding_peaks_right = False

        if right_idx < len(xic_mono):
            previous_right_intensity = mono_int[right_idx]
        if left_idx >= 0:
            previous_left_intensity = mono_int[left_idx]

        if not adding_peaks_left and not adding_peaks_right:
            break

    # order the integration indices in ascending order
    integration_indices = [nearest_idx]
    integration_indices.extend(left_integration_indices)
    integration_indices.extend(right_integration_indices)
    integration_indices = sorted(integration_indices)

    mono_intensities = np.array([v[1] for v in xic_mono[integration_indices]])
    iso1_intensities = np.array([v[1] for v in xic_iso1[integration_indices]])

    m, r = score_xics(mono_intensities, iso1_intensities)

    sufficient_peaks = len(integration_indices) >= config.min_consec_peaks + 2
    passes_iso_ratio = config.isorange[0] < m < config.isorange[1]
    meets_rthresh = r >= config.r_thresh

    if sufficient_peaks and passes_iso_ratio and meets_rthresh:
        return np.array(integration_indices)
    else:
        return np.array([])

def integrate_xics(xics, integration_indices, lipid):
    areas = []
    summed_area = 0
    max_intensity = np.max(xics[0][:,1])
    for i in range(len(xics)):
        sub_xic = np.vstack(xics[i][integration_indices])
        area = np.trapezoid(sub_xic[:,1], sub_xic[:,0])
        areas.append(area)
        summed_area += area

    # if lipid.internal_std:
    # plt.vlines(xics[0][integration_indices[0], 0], 0, max_intensity, linestyles='dashed', color='r')
    # plt.vlines(xics[0][integration_indices[-1], 0], 0, max_intensity, linestyles='dashed', color='r')
    # plt.vlines(lipid.calibrated_rt, 0, max_intensity, linestyles='dashed', color='b')
    # for i in range(len(xics)):
    #     plt.plot(xics[i][:,0], xics[i][:,1], label="iso"+str(i))
    # plt.legend()
    # plt.title("Lipid: " + lipid.name + " Area: " + str(summed_area))
    # plt.show()
    return summed_area

def quantify_lipid(spectra_file, lipid, config):
    abs_width = (lipid.mass / 1e6) * config.ppm_tol
    min_rt = lipid.calibrated_rt - (config.rt_tol/2)
    max_rt = lipid.calibrated_rt + (config.rt_tol/2)

    masses = [lipid.mass]
    for i in range(1,config.iso_count):
        masses.append(lipid.mass + (i * 1.0033))

    xics = []
    for i in range(len(masses)):
        xic = spectra_file.importer.get_eic(masses[i], abs_width, [min_rt, max_rt])
        xics.append(xic)

    if check_xics(xics, config):
        integration_indices = get_integration_indices(xics[0], xics[1], lipid.calibrated_rt, config)
        if integration_indices is not None:
            area = integrate_xics(xics, integration_indices, lipid)
            return area
        else:
            # print("Failed to quantify", lipid.name, "in", spectra_file.filepath)
            # plt.vlines(lipid.calibrated_rt, 0, np.max(xics[0][:, 1]), linestyles='dashed', color='r')
            # for i in range(len(xics)):
            #     plt.plot(xics[i][:, 0], xics[i][:, 1], label="iso" + str(i))
            # plt.legend()
            # plt.title("Lipid: " + lipid.name)
            # plt.show()
            return 0.0
    else:
        # print("Failed to quantify", lipid.name, "in", spectra_file.filepath)
        # plt.vlines(lipid.calibrated_rt, 0,np.max(xics[0][:,1]), linestyles='dashed', color='r')
        # for i in range(len(xics)):
        #     plt.plot(xics[i][:, 0], xics[i][:, 1], label="iso" + str(i))
        # plt.legend()
        # plt.title("Lipid: " + lipid.name)
        # plt.show()
        return 0.0

def quantify_isomeric_lipids(spectra_file, lipid_group, config):
    areas = [0.0] * len(lipid_group)

    lipid = lipid_group[0]
    min_lipid_rt = min(lipid_group, key=lambda x:x.calibrated_rt).calibrated_rt
    max_lipid_rt = max(lipid_group, key=lambda x:x.calibrated_rt).calibrated_rt

    #They all have the same mass, so we only need to make one set of xics.
    abs_width = (lipid.mass / 1e6) * config.ppm_tol
    min_rt = min_lipid_rt - (config.rt_tol / 2)
    max_rt = max_lipid_rt + (config.rt_tol / 2)

    masses = [lipid.mass]
    for i in range(1, config.iso_count):
        masses.append(lipid.mass + (i * 1.0033))

    xics = []
    for i in range(len(masses)):
        xic = spectra_file.importer.get_eic(masses[i], abs_width, [min_rt, max_rt])
        xics.append(xic)


    if check_xics(xics, config):
        rt_int_mono = np.transpose(np.array([xics[0][:, 0], xics[0][:, 1]]))
        peaks_mono = fastpeakdetect(rt_int_mono)
        possible_integration_indices = np.empty(len(peaks_mono), dtype=object)
        for i in range(len(peaks_mono)):
            possible_integration_indices[i] = check_valid_peak(xics[0], xics[1], peaks_mono[i], config)

        #Remove peaks and integration indices where the peak was not found to be valid
        invalid_integration_indices = []
        for i in range(len(possible_integration_indices)):
            if len(possible_integration_indices[i]) == 0:
                invalid_integration_indices.append(i)

        #Now remove these indices from both peaks_mono and possible_integration_indices
        peaks_mono = [peak for i,peak in enumerate(peaks_mono) if i not in invalid_integration_indices]
        possible_integration_indices = [integ for i, integ in enumerate(possible_integration_indices) if i not in invalid_integration_indices]

        # max_intensity = np.max(xics[0][:,1])
        # for i in range(len(xics)):
        #     plt.plot(xics[i][:, 0], xics[i][:, 1], label="iso"+str(i))


        if len(peaks_mono) == 0:
            print("No peaks detected")
            # plt.show()
            return [0.0 for i in range(len(lipid_group))]

        #In the case where there is only a single valid peak
        #Assign the closest lipid in RT to that peak
        elif len(peaks_mono) == 1:
            delta_rts = [0.0] * len(lipid_group)
            for i in range(len(lipid_group)):
                delta_rts[i] = abs(peaks_mono[0][0] - lipid_group[i].calibrated_rt)

            min_index = np.argmin(delta_rts)
            area = integrate_xics(xics, possible_integration_indices[0], lipid_group[min_index])
            areas[min_index] = area

            # #print the accepted lipid name and RT in green at the max_intensity
            # plt.vlines(lipid_group[min_index].calibrated_rt, 0, max_intensity, linestyles='dashed', color='g')
            # plt.text(lipid_group[min_index].calibrated_rt, max_intensity,
            #          lipid_group[min_index].name + "_"+ str(lipid_group[min_index].calibrated_rt), color='g')
            #
            # for i in range(len(lipid_group)):
            #     if i == min_index:
            #         continue
            #     else:
            #         plt.vlines(lipid_group[i].calibrated_rt, 0, max_intensity, linestyles='dashed', color='r')
            #         plt.text(lipid_group[i].calibrated_rt, max_intensity,
            #                  lipid_group[i].name + "_" + str(lipid_group[i].calibrated_rt), color='r')
            # plt.legend()
            # plt.show()

            return areas

        elif len(peaks_mono) < len(lipid_group):
            print("Fewer peaks than lipids.")
            rt_sorted_lipid_indices = sorted(range(len(lipid_group)), key=lambda i: lipid_group[i].calibrated_rt)
            index_combinations = []
            total_rt_errors = []

            for start in range(len(lipid_group)-len(peaks_mono)):
                combination = rt_sorted_lipid_indices[start:start+len(peaks_mono)]
                index_combinations.append(combination)

            for combo in index_combinations:
                total_rt_error = 0
                for i,lipid_index in enumerate(combo):
                    rt_error = abs(peaks_mono[i][0] - lipid_group[lipid_index].calibrated_rt)
                    total_rt_error += rt_error
                total_rt_errors.append(total_rt_error)

            min_index = np.argmin(total_rt_errors)
            for i in range(len(index_combinations[min_index])):
                lipid_index = index_combinations[min_index][i]
                area = integrate_xics(xics, possible_integration_indices[i], lipid_group[lipid_index])
                areas[lipid_index] = area

            # for i in range(len(index_combinations[min_index])):
            #     lipid_index = index_combinations[min_index][i]
            #     plt.vlines(lipid_group[lipid_index].calibrated_rt, 0, max_intensity, linestyles='dashed', color='g')
            #     plt.text(lipid_group[lipid_index].calibrated_rt, max_intensity,
            #              lipid_group[lipid_index].name + "_" + str(lipid_group[lipid_index].calibrated_rt), color='g')
            #
            # for i in range(len(lipid_group)):
            #     if i in index_combinations[min_index]:
            #         continue
            #     else:
            #         plt.vlines(lipid_group[i].calibrated_rt, 0, max_intensity, linestyles='dashed', color='r')
            #         plt.text(lipid_group[i].calibrated_rt, max_intensity,
            #                  lipid_group[i].name + "_" + str(lipid_group[i].calibrated_rt), color='r')
            # plt.legend()
            # plt.show()

            return areas

        elif len(peaks_mono) == len(lipid_group):
            print("Same number of peaks and lipids.")
            rt_sorted_lipid_indices = sorted(range(len(lipid_group)), key=lambda i: lipid_group[i].calibrated_rt)

            for i in range(len(peaks_mono)):
                lipid_index = rt_sorted_lipid_indices[i]
                area = integrate_xics(xics, possible_integration_indices[i], lipid_group[lipid_index])
                areas[lipid_index] = area

            # for i in range(len(rt_sorted_lipid_indices)):
            #     lipid_index = rt_sorted_lipid_indices[i]
            #     plt.vlines(lipid_group[lipid_index].calibrated_rt, 0, max_intensity, linestyles='dashed', color='g')
            #     plt.text(lipid_group[lipid_index].calibrated_rt, max_intensity,
            #              lipid_group[lipid_index].name + "_" + str(lipid_group[lipid_index].calibrated_rt), color='g')
            #
            # plt.legend()
            # plt.show()
            return areas

        #We only get here if more peaks are valid than we have lipids
        else:
            print("More peaks than lipids.")
            rt_sorted_lipid_indices = sorted(range(len(lipid_group)), key=lambda i: lipid_group[i].calibrated_rt)
            peak_indices = [i for i in range(len(peaks_mono))]
            combinations = list(itertools.combinations(peak_indices, len(lipid_group)))
            rt_errors = [0.0] * len(combinations)
            for i, combo in enumerate(combinations):
                total_rt_error = 0
                for j in range(len(combo)):
                    peak_index = combo[j]
                    rt_error = abs(peaks_mono[peak_index][0] - lipid_group[j].calibrated_rt)
                    total_rt_error += rt_error
                rt_errors[i] = total_rt_error

            min_index = np.argmin(rt_errors)
            for i in range(len(combinations[min_index])):
                peak_index = combinations[min_index][i]
                lipid_index = rt_sorted_lipid_indices[i]
                area = integrate_xics(xics, possible_integration_indices[peak_index], lipid_group[lipid_index])
                areas[lipid_index] = area

            # for i in range(len(rt_sorted_lipid_indices)):
            #     lipid_index = rt_sorted_lipid_indices[i]
            #     plt.vlines(lipid_group[lipid_index].calibrated_rt, 0, max_intensity, linestyles='dashed', color='g')
            #     plt.text(lipid_group[lipid_index].calibrated_rt, max_intensity,
            #              lipid_group[lipid_index].name + "_" + str(lipid_group[lipid_index].calibrated_rt), color='g')
            #
            # for i in range(len(combinations[min_index])):
            #     peak_index = combinations[min_index][i]
            #     peak = peaks_mono[peak_index]
            #     plt.scatter([peak[0]], [peak[1]], s=5, c='g')
            #
            # for i in range(len(peaks_mono)):
            #     if i in combinations[min_index]:
            #         continue
            #     else:
            #         plt.scatter([peaks_mono[i][0]], [peaks_mono[i][1]], s=5, c='r')
            #
            # plt.legend()
            # plt.show()

            return areas

    return np.zeros(len(lipid_group))

def create_spectrafiles_fromcsv(filepath, load_importers=True):
    files_df = pd.read_csv(filepath)

    def make_spectra_file(row):
        return SpectraFile(
            row["Filepath"],
            row["Samplenum"],
            row["Biorep"],
            row["Techrep"],
            row["Samplevol"],
            row["Stdvol"],
            load_importers
        )

    with ThreadPoolExecutor() as executor:
        spectra_files = executor.map(make_spectra_file, files_df.to_dict("records"))

    return list(spectra_files)

def normalize_areas_to_stds(spectra_files, quant_df, config):
    #These will eventually be removed from the quantification dataframe, but not yet.
    standard_class_indices = {}
    for i in range(len(quant_df["Standard"])):
        if quant_df["Standard"][i]:
            standard_class_indices[quant_df["Lipid Class"][i]] = i

    #First level of this is class
    #Second level is BR which contains a list of non-zero areas
    br_standard_areas = {}
    for lipid_class in standard_class_indices:
        br_standard_areas[lipid_class] = {}
        for spectra_file in spectra_files:
            std_area_in_file = quant_df[spectra_file.fileid][standard_class_indices[lipid_class]]
            if std_area_in_file > 0:
                if spectra_file.biorep not in br_standard_areas[lipid_class]:
                    br_standard_areas[lipid_class][spectra_file.biorep] = [std_area_in_file]
                else:
                    br_standard_areas[lipid_class][spectra_file.biorep].append(std_area_in_file)

    br_std_averages = {}
    for lipid_class in br_standard_areas:
        br_std_averages[lipid_class] = {}
        for br in br_standard_areas[lipid_class]:
            br_std_averages[lipid_class][br] = np.mean(br_standard_areas[lipid_class][br])

    std_indices = np.array(standard_class_indices.values())

    for i in range(len(spectra_files)):
        spectra_file = spectra_files[i]
        standardized_areas = np.array([0.0 for i in range(len(quant_df[spectra_file.fileid]))])
        for j in range(len(quant_df[spectra_file.fileid])):
            if quant_df["Lipid Class"][j] in standard_class_indices and j not in std_indices:
                lipid_area = quant_df[spectra_file.fileid][j]
                lipid_class = quant_df["Lipid Class"][j]
                standard_area = quant_df[spectra_file.fileid][standard_class_indices[lipid_class]]
                if standard_area == 0:
                    standard_area = br_std_averages[lipid_class][spectra_file.biorep]
                if config.normalize_to_int_stds:
                    ratio = lipid_area / standard_area
                    v_normalized = ratio * (spectra_file.stdvol / spectra_file.samplevol)
                    standardized_areas[j] = v_normalized
                else:
                    standardized_areas[j] = lipid_area

            elif j in std_indices:
                standardized_areas[j] = 0

        if config.normalize_to_total:
            standardized_areas = standardized_areas/np.sum(standardized_areas)
        quant_df[spectra_file.fileid] = standardized_areas

    std_indices = np.array(standard_class_indices.values())
    quant_df = quant_df.drop(index=std_indices)
    quant_df.reset_index(drop=True, inplace=True)

    return quant_df

def process_file(file_lib_config):
    file, lib, config = file_lib_config
    print("Processing", file.filepath)

    calibrate_library_rts(file, lib, config)

    file_areas = [0.0] * len(lib.lipid_indices)
    for j in range(len(lib.lipid_groups)):
        if len(lib.lipid_groups[j]) == 1:
            area = quantify_lipid(file, lib.lipid_groups[j][0], config)
            index = lib.lipid_indices[lib.lipid_groups[j][0].name]
            file_areas[index] = area
            # print(lib.lipid_groups[j][0].name, "quantity:", area)
        else:
            areas = quantify_isomeric_lipids(file, lib.lipid_groups[j], config)
            for i in range(len(lib.lipid_groups[j])):
                index = lib.lipid_indices[lib.lipid_groups[j][i].name]
                # print(lib.lipid_groups[j][i].name, "quantity:", areas[i])
                file_areas[index] = areas[i]

    # Clean up
    file.importer.close()
    del file.importer.indexed_file

    return file.fileid, file_areas

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    timestart = time.perf_counter()
    config = LipidQuantConfig()
    config.normalize_to_total = True
    config.normalize_to_int_stds = True

    lib = LipidLibrary("Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\Raw\\DDA\\calibrated_mzML\\rt_calibrated_filtered_library.csv", config)
    #spectra_files = create_spectrafiles(spectra_files_directory)

    filenames_id_br_trs_svol_eqvol = {"2024_09_04_AS_46_CS_08_ESIneg_fixed.mzML": [4,1,1,85,45],
                                       "2024_09_04_AS_48_CS_08_ESIneg_fixed.mzML": [4,1,2,85,45],
                                       "2024_09_04_AS_50_CS_08_ESIneg_fixed.mzML": [4,1,3,85,45],

                                       "2024_09_04_AS_52_CS_09_ESIneg_fixed.mzML": [5,3,1,30,15],
                                       "2024_09_04_AS_54_CS_09_ESIneg_fixed.mzML": [5,3,2,30,15],
                                       "2024_09_04_AS_56_CS_09_ESIneg_fixed.mzML": [5,3,3,30,15],

                                       "2024_09_04_AS_142_CS_25_ESIneg_fixed.mzML": [4,2,1,50,25],
                                       "2024_09_04_AS_144_CS_25_ESIneg_fixed.mzML": [4,2,2,50,25],
                                       "2024_09_04_AS_146_CS_25_ESIneg_fixed.mzML": [4,2,3,50,25],

                                       "2024_09_04_AS_148_CS_26_ESIneg_fixed.mzML": [6,1,1,100,50],
                                       "2024_09_04_AS_150_CS_26_ESIneg_fixed.mzML": [6,1,2,100,50],
                                       "2024_09_04_AS_152_CS_26_ESIneg_fixed.mzML": [6,1,3,100,50],

                                       "2024_09_04_AS_118_CS_21_ESIneg_fixed.mzML": [6,2,1,50,25],
                                       "2024_09_04_AS_120_CS_21_ESIneg_fixed.mzML": [6,2,2,50,25],
                                       "2024_09_04_AS_122_CS_21_ESIneg_fixed.mzML": [6,2,3,50,25],

                                       "2024_09_04_AS_136_CS_24_ESIneg_fixed.mzML": [6, 3, 1,30,15],
                                       "2024_09_04_AS_138_CS_24_ESIneg_fixed.mzML": [6, 3, 2,30,15],
                                       "2024_09_04_AS_140_CS_24_ESIneg_fixed.mzML": [6, 3, 3,30,15],

                                       "2024_09_04_AS_100_CS_17_ESIneg_fixed.mzML": [5, 1, 1,40,40],
                                       "2024_09_04_AS_102_CS_17_ESIneg_fixed.mzML": [5, 1, 2,40,40],
                                       "2024_09_04_AS_104_CS_17_ESIneg_fixed.mzML": [5, 1, 3,40,40],

                                       "2024_09_04_AS_88_CS_15_ESIneg_fixed.mzML": [5, 2, 1,50,25],
                                       "2024_09_04_AS_90_CS_15_ESIneg_fixed.mzML": [5, 2, 2,50,25],
                                       "2024_09_04_AS_92_CS_15_ESIneg_fixed.mzML": [5, 2, 3,50,25],

                                       "2024_09_04_AS_76_CS_13_ESIneg_fixed.mzML": [4, 3, 1,30,15],
                                       "2024_09_04_AS_78_CS_13_ESIneg_fixed.mzML": [4, 3, 2,30,15],
                                       "2024_09_04_AS_80_CS_13_ESIneg_fixed.mzML": [4, 3, 3,30,15]
                                       }

    spectra_files_csv = r"C:\Users\MartyLabPC\Desktop\MartyLab\calibrated_mzML\lipidquant_CS_LetAB_local.csv"
    os.chdir(os.path.dirname(spectra_files_csv))

    spectra_files = create_spectrafiles_fromcsv(spectra_files_csv, load_importers=False)

    if True:
        quantification_results = pd.DataFrame()
        quantification_results["Lipid name"] = [lipid.name for lipid in lib.lipids]
        quantification_results["Lipid Class"] = [lipid.lipid_class for lipid in lib.lipids]
        quantification_results["Standard"] = [lipid.internal_std for lipid in lib.lipids]

        args = [(file, lib, config) for file in spectra_files]

        if True:
            with ThreadPoolExecutor(max_workers=2) as executor:
                results = list(executor.map(process_file, args))
        else:
            results = [process_file(f_l_c) for f_l_c in args]

        # Gather results
        for res in results:
            fileid, file_areas = res
            quantification_results[fileid] = file_areas

        quantification_results.to_csv("quantified_lipids.csv", index=False)

    else:
        quantification_results = pd.read_csv("quantified_lipids.csv", delimiter=",")

    normalized_quantification_results = normalize_areas_to_stds(spectra_files, quantification_results, config)
    plot_class_piecharts(spectra_files, normalized_quantification_results)
    compare_skyline_lipidquant(spectra_files, quantification_results)
    exit()
    # skyline_quant_results = load_skyline_results()
    # normalized_skyquant_results = normalize_areas_to_stds(spectra_files, skyline_quant_results, config)
    #plot_indlipid_barcharts(spectra_files, normalized_quantification_results)
    plot_class_piecharts(spectra_files, normalized_skyquant_results)
    timeend = time.perf_counter()
    print("Time to process files:", timeend - timestart)
    exit()




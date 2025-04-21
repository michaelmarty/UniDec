import pandas as pd
from pyopenms import *
import Scripts.JGP.xic_xtractor as xic
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.IsoDec.datatools as dt
from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol
import numpy as np
import time
import os
import fnmatch
from scipy.stats import linregress

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

class PyLineConfig:
    def __init__(self):
        #XIC creation parameters
        self.ppm_tol = 20
        self.rt_tol = 5
        self.iso_count = 3

        #XIC Filtering parameters
        self.min_peak_intensity = 500
        self.noise_threshold = 100
        self.min_consec_peaks = 4
        self.isorange = [0.1,0.9]
        self.r_thresh = 0.5
        self.max_increasing = 3
        self.min_intensity_ratio = 0.01

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
        self.lipids = self.parse_library(filepath)
        self.lipid_groups = self.group_lipids_by_rt_mz(config.ppm_tol, config.rt_tol)

    def parse_library(self, filepath):
        lipids = np.array([])

        library_df = pd.read_csv(filepath, delimiter=",")
        rts = library_df["Retention Time"]
        names = library_df["Lipid Name"]
        adducts = library_df["Adduct"]
        mzs = library_df["m/z"]
        formulas = library_df["Chemical Formula"]
        classes = library_df["Lipid class"]
        standards = library_df["Standard"]

        for i in range(len(rts)):
            lipids = np.append(lipids, LibraryLipid(names[i], mzs[i], rts[i], classes[i], standards[i], adducts[i], formulas[i]))

        return lipids

    def group_lipids_by_rt_mz(self, ppm_tol, rt_tol):
        groups = []

        lipid_masses = np.array([l.mass for l in self.lipids])

        sorted_indices = np.argsort(lipid_masses)
        lipid_masses = lipid_masses[sorted_indices]
        self.lipids = self.lipids[sorted_indices]

        grouped_lipids = []
        #Find cases where
        for i in range(len(self.lipids)):
            if i in grouped_lipids:
                continue
            mass_tol = (self.lipids[i].mass/1e6)*ppm_tol
            within_tol = fastwithin_abstol(lipid_masses, self.lipids[i].mass, mass_tol)
            current_group = [self.lipids[i]]
            grouped_lipids.append(i)
            for j in range(len(within_tol)):
                if i == within_tol[j] or within_tol[j] in grouped_lipids:
                    continue
                else:
                    delta_rt = abs(self.lipids[within_tol[j]].rt - self.lipids[i].rt)
                    if abs(self.lipids[within_tol[j]].rt - self.lipids[i].rt) < rt_tol:
                        current_group.append(self.lipids[within_tol[j]])
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
    def __init__(self, filepath, sample_num, biorep, techrep, samplevol, stdvol):
        self.filepath = filepath
        self.ind_file = None
        self.sample_id = sample_num
        self.biorep = biorep
        self.techrep = techrep
        self.fileid = str("S" + str(sample_num) + "_BR" + str(biorep) + "_TR" + str(techrep))

        self.samplevol = float(samplevol)
        self.stdvol = float(stdvol)

    def load_indexed_file(self):
        exp = MSExperiment()
        MzMLFile().load(self.filepath, exp)
        self.ind_file = xic.IndexedFile(exp)

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



            is_xic = file.ind_file.extract_xic(mz, abs_width, min_rt, max_rt)
            max_index = np.argmax(is_xic[:,2])
            if is_xic[max_index,2] > config.min_peak_intensity:
                is_rts = np.append(is_rts,library.lipids[i].rt)
                rt_error = (is_xic[max_index,0]/60) - library.lipids[i].rt
                delta_rts = np.append(delta_rts, rt_error)
                library.lipids[i].calibrated_rt =is_xic[max_index,0]/60

    is_rts = np.append(is_rts,0)
    delta_rts = np.append(delta_rts, 0)

    is_rts = np.append(is_rts,20)
    delta_rts = np.append(delta_rts, 0)

    sorted_indices = np.argsort(is_rts)
    is_rts = is_rts[sorted_indices]
    delta_rts = delta_rts[sorted_indices]

    for i in range(len(library.lipids)):
        if library.lipids[i].internal_std:
            continue
        else:
            rt = library.lipids[i].rt
            rt_correction = np.interp(rt, is_rts, delta_rts)
            library.lipids[i].calibrated_rt = rt + rt_correction

def check_xics(xics, config):
    #First check if the mono xic meets intensity threshold
    if np.max(xics[0][:,2]) < config.min_peak_intensity:
        return False

    #Next check if the mono xic meets the consecutive scans threshold
    passes_consec = False
    curr_consec = 0
    for i in range(len(xics[0][:,2])):
        if xics[0][i,2] > config.noise_threshold:
            curr_consec += 1
            if curr_consec >= config.min_consec_peaks:
                passes_consec = True
                break
        else:
            curr_consec = 0

    if not passes_consec:
        return False

    return True

def calculate_cosinesimilarity(predicted_dist, expected_dist):
    """
    Calculate the cosine similarity between two distributions.

    :param predicted_dist: The predicted distribution (numpy array).
    :param expected_dist: The expected distribution (numpy array).
    :return: The cosine similarity.
    """

    # Ensure both distributions are the same length
    if len(predicted_dist) != len(expected_dist):
        raise ValueError("Both distributions must have the same length.")

    # Ensure both distributions are numpy arrays
    predicted_dist = np.array(predicted_dist)
    expected_dist = np.array(expected_dist)

    # Normalize the distributions to max 1
    predicted_dist = predicted_dist / np.max(predicted_dist)
    expected_dist = expected_dist / np.max(expected_dist)

    # Calculate the dot product of the two distributions
    dot_product = np.dot(predicted_dist, expected_dist)

    # Calculate the magnitude of each distribution
    magnitude_predicted = np.linalg.norm(predicted_dist)
    magnitude_expected = np.linalg.norm(expected_dist)

    # Calculate the cosine similarity
    cosine_similarity = dot_product / (magnitude_predicted * magnitude_expected)

    return cosine_similarity

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
    rt_int_mono = np.transpose(np.array([xic_mono[:,1], xic_mono[:,2]]))

    peaks_mono = dt.fastpeakdetect(rt_int_mono, window=5, threshold=0.05)
    possible_integration_indexes = np.empty(len(peaks_mono), dtype=object)

    for i in range(len(peaks_mono)):
        nearest_idx = fastnearest(np.array(xic_mono[:,1]), peaks_mono[i][0])
        peak_intensity = peaks_mono[i][1]
        adding_peaks_left = True
        adding_peaks_right = True

        left_integration_indices = []
        right_integration_indices = []

        #we expect intensity to go down, so if it increases for more than N datapoints, cutoff before the first increase
        previous_left_intensity = peak_intensity
        previous_right_intensity = peak_intensity

        current_left_increasing = 0
        current_right_increasing = 0

        for j in range(1, len(xic_mono)):
            left_idx =nearest_idx - j
            right_idx = nearest_idx + j

            if adding_peaks_left:
                if left_idx >= 0 and xic_mono[left_idx,2]/peak_intensity > config.min_intensity_ratio:
                    left_integration_indices.append(left_idx)

                    #Check for poorly resolved peaks
                    if xic_mono[left_idx,2] > previous_left_intensity:
                        current_left_increasing += 1
                        if current_left_increasing >= config.max_increasing:
                            #we remove the recent increasing values
                            for _ in range(current_left_increasing):
                                left_integration_indices.pop()
                            adding_peaks_left = False
                    else:
                        current_left_increasing = 0
                else:
                    if left_idx >=0:
                        left_integration_indices.append(left_idx)
                    adding_peaks_left = False


            if adding_peaks_right:
                if right_idx < len(xic_mono) and xic_mono[right_idx,2]/peak_intensity > config.min_intensity_ratio:
                    right_integration_indices.append(right_idx)

                    if xic_mono[right_idx,2] > previous_right_intensity:
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
                previous_right_intensity = xic_mono[right_idx, 2]
            if left_idx >= 0:
                previous_left_intensity = xic_mono[left_idx, 2]

            if not adding_peaks_left and not adding_peaks_right:
                break

        #order the integration indices in ascending order
        integration_indices = [nearest_idx]
        integration_indices.extend(left_integration_indices)
        integration_indices.extend(right_integration_indices)
        integration_indices = sorted(integration_indices)

        mono_intensities = np.array([v[2] for v in xic_mono[integration_indices]])
        iso1_intensities = np.array([v[2] for v in xic_iso1[integration_indices]])

        m,r = score_xics(mono_intensities, iso1_intensities)

        sufficient_peaks = len(integration_indices) >= config.min_consec_peaks+2
        passes_iso_ratio = config.isorange[0] < m < config.isorange[1]
        meets_rthresh = r >= config.r_thresh

        if sufficient_peaks and passes_iso_ratio and meets_rthresh:
            possible_integration_indexes[i] = np.array(integration_indices)
        else:
            possible_integration_indexes[i] = np.array([])

    delta_rts = [abs(p[0] - (rt*60)) for p in peaks_mono]
    sorted_indices = np.argsort(delta_rts).astype(int)
    peaks_mono = peaks_mono[sorted_indices]
    possible_integration_indexes = possible_integration_indexes[sorted_indices]

    for i in range(len(peaks_mono)):
        if len(possible_integration_indexes[i]) > 0:
            return possible_integration_indexes[i]
        else:
            continue
    return None

def integrate_xics(xics, integration_indices, lipid):
    areas = []
    summed_area = 0
    max_intensity = np.max(xics[0][:,2])
    for i in range(len(xics)):
        sub_xic = np.vstack(xics[i][integration_indices])
        area = np.trapezoid(sub_xic[:,2], sub_xic[:,1])
        areas.append(area)
        summed_area += area

    # plt.vlines(xics[0][integration_indices[0], 1]/60, 0, max_intensity, linestyles='dashed', color='r')
    # plt.vlines(xics[0][integration_indices[-1], 1]/60, 0, max_intensity, linestyles='dashed', color='r')
    # plt.vlines(lipid.calibrated_rt, 0, max_intensity, linestyles='dashed', color='b')
    # for i in range(len(xics)):
    #     plt.plot(xics[i][:,1]/60, xics[i][:,2], label="iso"+str(i))
    # plt.legend()
    # plt.title("Lipid: " + lipid.name + " Area: " + str(summed_area))
    # plt.show()
    return summed_area

def quantify_lipid(spectra_file, lipid, config):
    abs_width = (lipid.mass / 1e6) * config.ppm_tol
    min_rt = lipid.rt - (config.rt_tol/2)
    max_rt = lipid.rt + (config.rt_tol/2)

    masses = [lipid.mass]
    for i in range(1,config.iso_count):
        masses.append(lipid.mass + (i * 1.0033))

    xics = []
    for i in range(len(masses)):
        xic = spectra_file.ind_file.extract_xic(masses[i], abs_width, min_rt, max_rt)
        xics.append(xic)

    if check_xics(xics, config):
        integration_indices = get_integration_indices(xics[0], xics[1], lipid.calibrated_rt, config)
        if integration_indices is not None:
            area = integrate_xics(xics, integration_indices, lipid)
            return area
        else:
            if lipid.internal_std:
                print("Failed to quantify", lipid.name, "in", spectra_file.filepath)
                plt.vlines(lipid.calibrated_rt * 60, 0, np.max(xics[0][:, 2]), linestyles='dashed', color='r')
                for i in range(len(xics)):
                    plt.plot(xics[i][:, 1], xics[i][:, 2], label="iso" + str(i))
                plt.legend()
                plt.title("Lipid: " + lipid.name)
                plt.show()

            return None
    else:
        if lipid.internal_std:
            print("Failed to quantify", lipid.name, "in", spectra_file.filepath)
            plt.vlines(lipid.calibrated_rt*60, 0,np.max(xics[0][:,2]), linestyles='dashed', color='r')
            for i in range(len(xics)):
                plt.plot(xics[i][:, 1], xics[i][:, 2], label="iso" + str(i))
            plt.legend()
            plt.title("Lipid: " + lipid.name)
            plt.show()
        return None

def create_spectrafiles(directory):
    files = match_files(directory, "*.mzML")
    spectra_files = []
    for file in files:
        print("Processing", file)
        entering_file_info = True
        while entering_file_info:
            sample_number = int(input("Enter the sample number: "))
            biorep = int(input("Enter the biorep number: "))
            techrep = int(input("Enter the techrep number: "))
            print("Sample number:", sample_number, "Biorep number:", biorep, "Techrep number:", techrep)
            accepted = input("Is this correct? (Y/N): ")
            if accepted.lower() == "y":
                entering_file_info = False

        spectra_files.append(SpectraFile(file, sample_number, biorep, techrep))

    return spectra_files

def create_spectrafiles_fromdict(dictionary, sample_of_interest=None):
    spectra_files = []
    for file in dictionary:
        if sample_of_interest is not None:
            if dictionary[file][0] != sample_of_interest:
                continue
        spectra_files.append(SpectraFile(file, dictionary[file][0], dictionary[file][1],
                                         dictionary[file][2],dictionary[file][3], dictionary[file][4]))
    return spectra_files

def normalize_areas_to_stds(spectra_files, quant_df):
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
                    br_standard_areas[lipid_class][spectra_file.biorep] = [std_area_in_file/spectra_file.stdvol]
                else:
                    br_standard_areas[lipid_class][spectra_file.biorep].append(std_area_in_file/spectra_file.stdvol)

    br_std_averages = {}
    for lipid_class in br_standard_areas:
        br_std_averages[lipid_class] = {}
        for br in br_standard_areas[lipid_class]:
            br_std_averages[lipid_class][br] = np.mean(br_standard_areas[lipid_class][br])


    for i in range(len(spectra_files)):
        spectra_file = spectra_files[i]
        standardized_areas = np.array([0.0 for i in range(len(quant_df[spectra_file.fileid]))])
        for j in range(len(quant_df[spectra_file.fileid])):
            if quant_df["Lipid Class"][j] in standard_class_indices:
                lipid_area = quant_df[spectra_file.fileid][j] / spectra_file.samplevol
                lipid_class = quant_df["Lipid Class"][j]
                standard_area = quant_df[spectra_file.fileid][standard_class_indices[lipid_class]]/spectra_file.stdvol
                if standard_area == 0:
                    standard_area = br_std_averages[lipid_class][spectra_file.biorep]
                standardized_area = lipid_area / standard_area
                standardized_areas[j] = standardized_area

        standardized_areas = standardized_areas/np.sum(standardized_areas)
        quant_df[spectra_file.fileid] = standardized_areas

    std_indices = np.array(standard_class_indices.values())
    quant_df = quant_df.drop(index=std_indices)

    sample_of_interest = 5
    #Plot the average standardized area of each lipid for each biological replicate
    for i in range(len(quant_df["Lipid Name"])):
        br_areas = {}
        for j in range(len(spectra_files)):
            if spectra_files[j].sample_id == sample_of_interest:
                if spectra_files[j].biorep in br_areas:
                    br_areas[spectra_files[j].biorep].append(quant_df[spectra_files[j].fileid][i])
                else:
                    br_areas[spectra_files[j].biorep] = [quant_df[spectra_files[j].fileid][i]]

        x = list(br_areas.keys())
        y = [np.mean(v) for v in br_areas.values()]

        plt.bar(x, y)
        plt.title(quant_df["Lipid Name"][i])
        plt.show()

    return quant_df

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    timestart = time.perf_counter()
    config = PyLineConfig()

    lib = LipidLibrary("Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\rtcalibrated_dda_library.csv", config)

    spectra_files_directory = "Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\Raw\\MS1\\calibrated_mzML\\"
    os.chdir(spectra_files_directory)

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

    spectra_files = create_spectrafiles_fromdict(filenames_id_br_trs_svol_eqvol)

    if True:
        quantification_results = pd.DataFrame()
        quantification_results["Lipid Name"] = [lipid.name for lipid in lib.lipids]
        quantification_results["Lipid Class"] = [lipid.lipid_class for lipid in lib.lipids]
        quantification_results["Standard"] = [lipid.internal_std for lipid in lib.lipids]

        for i in range(len(spectra_files)):
            print("Processing file:", spectra_files[i].filepath)
            spectra_files[i].load_indexed_file()
            file_areas = [0.0 for _ in range(len(lib.lipids))]
            for j in range(len(lib.lipids)):
                area = quantify_lipid(spectra_files[i], lib.lipids[j], config)
                if area is not None:
                    file_areas[j] = area
            quantification_results[spectra_files[i].fileid] = file_areas
            del spectra_files[i].ind_file

        quantification_results.to_csv("quantified_lipids.csv", index=False)

    else:
        quantification_results = pd.read_csv("quantified_lipids.csv", delimiter=",")

    normalize_areas_to_stds(spectra_files, quantification_results)
    timeend = time.perf_counter()
    print("Time to process files:", timeend - timestart)
    exit()




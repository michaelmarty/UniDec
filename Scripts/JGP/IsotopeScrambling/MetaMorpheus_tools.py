import matplotlib as mpl
import os


class MM_PSMTsv_FileHeader:
    def __init__(self, header_line):
        self.header_dict = {"filename": None,
                            "ms2_scan_num": None,
                            "retention_time": None,
                            "num_experimental_peaks": None,
                            "total_ion_current": None,
                            "ms1_scan_num": None,
                            "precursor_charge": None,
                            "precursor_intensity": None,
                            "precursor_mz": None,
                            "precursor_mass": None,
                            "score": None,
                            "delta_score": None,
                            "notch": None,
                            "base_sequence": None,
                            "full_sequence": None,
                            "ambiguity_level": None,
                            "spectrum_match_count": None,
                            "mods": None,
                            "mods_chemical_formula": None,
                            "num_variable_mods": None,
                            "missed_cleavages": None,
                            "monoisotopic_mass": None,
                            "mass_diff_da": None,
                            "mass_diff_ppm": None,
                            "accession": None,
                            "protein_name": None,
                            "organism": None,
                            "aa_range": None,
                            "target_decoy_contaminant": None,
                            "matched_ion_series": None,
                            "matched_ion_mzs": None,
                            "matched_ion_mass_diff_da": None,
                            "matched_ion_mass_diff_ppm": None,
                            "matched_ion_intensities": None,
                            "matched_ion_count": None,
                            "cumulative_target": None,
                            "cumulative_decoy": None,
                            "qvalue": None
                            }
        header_split = list(header_line.split("\t"))
        try:
            self.header_dict["filename"] = header_split.index("File Name")
            self.header_dict["ms2_scan_num"] = header_split.index("Scan Number")
            self.header_dict["retention_time"] = header_split.index("Scan Retention Time")
            self.header_dict["num_experimental_peaks"] = header_split.index("Num Experimental Peaks")
            self.header_dict["total_ion_current"] = header_split.index("Total Ion Current")
            self.header_dict["ms1_scan_num"] = header_split.index("Precursor Scan Number")
            self.header_dict["precursor_charge"] = header_split.index("Precursor Charge")
            self.header_dict["precursor_intensity"] = header_split.index("Precursor Intensity")
            self.header_dict["precursor_mz"] = header_split.index("Precursor MZ")
            self.header_dict["precursor_mass"] = header_split.index("Precursor Mass")
            self.header_dict["score"] = header_split.index("Score")
            self.header_dict["delta_score"] = header_split.index("Delta Score")
            self.header_dict["notch"] = header_split.index("Notch")
            self.header_dict["base_sequence"] = header_split.index("Base Sequence")
            self.header_dict["full_sequence"] = header_split.index("Full Sequence")
            self.header_dict["ambiguity_level"] = header_split.index("Ambiguity Level")
            self.header_dict["spectrum_match_count"] = header_split.index("Spectrum Match Count")
            self.header_dict["mods"] = header_split.index("Mods")
            self.header_dict["mods_chemical_formula"] = header_split.index("Mods Chemical Formulas")
            self.header_dict["num_variable_mods"] = header_split.index("Num Variable Mods")
            self.header_dict["missed_cleavages"] = header_split.index("Missed Cleavages")
            self.header_dict["monoisotopic_mass"] = header_split.index("Monoisotopic Mass")
            self.header_dict["mass_diff_da"] = header_split.index("Mass Diff (Da)")
            self.header_dict["mass_diff_ppm"] = header_split.index("Mass Diff (ppm)")
            self.header_dict["accession"] = header_split.index("Accession")
            self.header_dict["protein_name"] = header_split.index("Name")
            self.header_dict["organism"] = header_split.index("Organism Name")
            self.header_dict["aa_range"] = header_split.index("Start and End Residues In Full Sequence")
            self.header_dict["target_decoy_contaminant"] = header_split.index("Decoy/Contaminant/Target")
            self.header_dict["matched_ion_series"] = header_split.index("Matched Ion Series")
            self.header_dict["matched_ion_mzs"] = header_split.index("Matched Ion Mass-To-Charge Ratios")
            self.header_dict["matched_ion_mass_diff_da"] = header_split.index("Matched Ion Mass Diff (Da)")
            self.header_dict["matched_ion_mass_diff_ppm"] = header_split.index("Matched Ion Mass Diff (Ppm)")
            self.header_dict["matched_ion_intensities"] = header_split.index("Matched Ion Intensities")
            self.header_dict["matched_ion_count"] = header_split.index("Matched Ion Counts")
            self.header_dict["cumulative_target"] = header_split.index("Cumulative Target")
            self.header_dict["cumulative_decoy"] = header_split.index("Cumulative Decoy")
            self.header_dict["qvalue"] = header_split.index("QValue")
        except Exception as e:
            print(e)
            return


class MM_PSM:
    def __init__(self, line, current_id, parsed_header):
        try:
            split_line = line.split("\t")
            self.global_id = current_id
            self.filename = split_line[parsed_header["filename"]]
            self.ms2_scan_num = int(split_line[parsed_header["ms2_scan_num"]])
            self.retention_time = float(split_line[parsed_header["retention_time"]])
            self.num_experimental_peaks = int(float(split_line[parsed_header["num_experimental_peaks"]]))
            self.total_ion_current = float(split_line[parsed_header["total_ion_current"]])
            self.ms1_scan_num = int(split_line[parsed_header["ms1_scan_num"]])
            self.precursor_charge = int(float(split_line[parsed_header["precursor_charge"]]))
            self.precursor_intensity = float(split_line[parsed_header["precursor_intensity"]])
            self.precursor_mz = float(split_line[parsed_header["precursor_mz"]])
            self.precursor_mass = float(split_line[parsed_header["precursor_mass"]])
            self.score = float(split_line[parsed_header["score"]])
            self.delta_score = float(split_line[parsed_header["delta_score"]])
            self.notch = split_line[parsed_header["notch"]]
            self.base_sequence = split_line[parsed_header["base_sequence"]]
            self.full_sequence = split_line[parsed_header["full_sequence"]]
            self.ambiguity_level = split_line[parsed_header["ambiguity_level"]]
            self.spectrum_match_count = int(split_line[parsed_header["spectrum_match_count"]])
            self.mods = split_line[parsed_header["mods"]]
            self.mods_chemical_formula = split_line[parsed_header["mods_chemical_formula"]]
            self.num_variable_mods = split_line[parsed_header["num_variable_mods"]]
            self.missed_cleavages = split_line[parsed_header["missed_cleavages"]]
            self.monoisotopic_mass = split_line[parsed_header["monoisotopic_mass"]]
            self.mass_diff_da = split_line[parsed_header["mass_diff_da"]]
            self.mass_diff_ppm = split_line[parsed_header["mass_diff_ppm"]]
            self.accession = split_line[parsed_header["accession"]]
            self.protein_name = split_line[parsed_header["protein_name"]]
            self.organism = split_line[parsed_header["organism"]]
            self.aa_range = split_line[parsed_header["aa_range"]]
            self.target_decoy_contaminant = split_line[parsed_header["target_decoy_contaminant"]]
            self.matched_ion_series = split_line[parsed_header["matched_ion_series"]]
            self.matched_ion_mzs = split_line[parsed_header["matched_ion_mzs"]]
            self.matched_ion_mass_diff_da = split_line[parsed_header["matched_ion_mass_diff_da"]]
            self.matched_ion_mass_diff_ppm = split_line[parsed_header["matched_ion_mass_diff_ppm"]]
            self.matched_ion_intensities = split_line[parsed_header["matched_ion_intensities"]]
            self.matched_ion_count = split_line[parsed_header["matched_ion_count"]]
            self.cumulative_target = float(split_line[parsed_header["cumulative_target"]])
            self.cumulative_decoy = float(split_line[parsed_header["cumulative_decoy"]])
            self.qvalue = float(split_line[parsed_header["qvalue"]])
        except Exception as e:
            print(e)
            self.global_id = None
            return

def load_psm_file(filepath):
    psms = []

    current_id = 0
    parsed_header = None
    found_header = False
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            if found_header:
                psm = MM_PSM(line, current_id, parsed_header.header_dict)
                psms.append(psm)
                current_id += 1
            else:
                split = line.split("\t")
                if split[0] == "File Name":
                    parsed_header = MM_PSMTsv_FileHeader(line)
                    found_header = True
    return psms



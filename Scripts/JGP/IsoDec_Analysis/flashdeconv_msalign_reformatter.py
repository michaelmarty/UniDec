import os
import fnmatch
import numpy as np
from Scripts.JGP.IsoDec_Analysis.compare_msaligns import load_isolation_info
import pandas as pd


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def correct_precursor_mz(file):
    with open(file, "r") as f:
        lines = f.readlines()

    current_precursor_mass = 0
    current_precursor_charge = 0
    current_precursor_mz_line = 0
    for i, line in enumerate(lines):
        if "BEGIN IONS" in line:
            current_precursor_mass = 0
            current_precursor_charge = 0
            current_precursor_mz_line = 0
        if "END IONS" in line:
            if current_precursor_mass != 0:
                new_precursor_mz = (current_precursor_mass / current_precursor_charge) + 1.007276467
                lines[current_precursor_mz_line] = "PRECURSOR_MZ=" + str(new_precursor_mz) + "\n"
        if "PRECURSOR_MASS" in line:
            current_precursor_mass = float(line.split("=")[1])
        if "PRECURSOR_CHARGE" in line:
            current_precursor_charge = int(line.split("=")[1])
        if "PRECURSOR_MZ" in line:
            current_precursor_mz_line = i

    with open(file, "w") as f:
        f.writelines(lines)

def correct_ID(file):
    with open(file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "ID" in line and "FRACTION" not in line and "MS_ONE" not in line:
            id = line.split("=")[1]
            lines[i] = "SPECTRUM_ID=" + id
    with open(file, "w") as f:
        f.writelines(lines)

def add_precursor_window(file, file_isolation_info, range=1.2):
    with open(file, "r") as f:
        lines = f.readlines()

    new_lines = np.array(lines)
    current_ms1_scan_num = 0
    current_ms2_scan_num = 0
    added_lines = 0
    for i, line in enumerate(lines):
        if "BEGIN IONS" in line:
            current_ms1_scan_num = 0
        if "SCANS" in line:
            current_ms2_scan_num = int(line.split("=")[1])
        if "MS_ONE_SCAN=" in line:
            current_ms1_scan_num = int(line.split("=")[1])
            if (current_ms1_scan_num not in file_isolation_info.keys() or
                    current_ms2_scan_num not in file_isolation_info[current_ms1_scan_num].keys()):
                #decrement from the current MS2 scan until we find the correct MS1 scan
                new_ms1_scan = current_ms2_scan_num - 1
                if new_ms1_scan in file_isolation_info.keys() and current_ms2_scan_num in file_isolation_info[new_ms1_scan].keys():
                    current_ms1_scan_num = new_ms1_scan
                else:
                    current_ms1_scan_num = 0
                    continue
            iso_mz = file_isolation_info[current_ms1_scan_num][current_ms2_scan_num]
            min = iso_mz - (range/2)
            max = iso_mz + (range/2)
            #now we insert the precursor window information
            begin_line = "PRECURSOR_WINDOW_BEGIN=" + str(min)
            end_line = "PRECURSOR_WINDOW_END=" + str(max)
            new_lines = np.insert(new_lines, i+added_lines, begin_line)
            added_lines += 1
            new_lines = np.insert(new_lines, i +added_lines, "\n")
            added_lines += 1
            new_lines = np.insert(new_lines, i+added_lines, end_line)
            added_lines += 1
            new_lines = np.insert(new_lines, i + added_lines, "\n")
            added_lines += 1

    with open(file, "w") as f:
        f.writelines(new_lines)







toppath = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\FLASHDeconv_nopreceding_ms1s\\"


files = match_files(toppath, "*_ms2.msalign")

for file in files:
    print("Processing file:", file)
    correct_precursor_mz(toppath + file)
exit()
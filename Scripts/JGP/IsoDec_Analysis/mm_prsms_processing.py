import os
import fnmatch
import numpy as np
import pandas as pd

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)


toppath = "D:\\Johnny\\IsoDec\\PXD019247\\MM_ID_inline\\Task1-CalibrateTask\\ReSearch_rep_all_MMs\\Task1-SearchTask\\Individual File Results\\"

psm_files = match_files(toppath, "*PSMs.psmtsv")
proteoforms_files = match_files(toppath, "*Proteoforms.psmtsv")



q_value_threshold = 0.01


for psm_file in psm_files:
    psms_df = pd.read_csv(toppath + psm_file, sep="\t")
    print("Processing file:", psm_file)
    print("Raw PSMs: ", len(psms_df))
    psms_df = psms_df[psms_df["QValue"] < q_value_threshold]
    psms_df  = psms_df[psms_df["Decoy/Contaminant/Target"] == "T"]
    ambiguous_psms = psms_df["Full Sequence"].str.contains(r'\|').sum()
    print("Ambiguous PSMs:", ambiguous_psms)
    print("Filtered PSMs:", len(psms_df))
    avg_matchedfrags = np.mean(psms_df["Matched Ion Counts"])
    print("Average Matched Fragments:", avg_matchedfrags)

for proteoforms_file in proteoforms_files:
    pfms_df = pd.read_csv(toppath + proteoforms_file, sep="\t")
    print("Processing file:", proteoforms_file)
    print("Raw Proteoforms: ", len(pfms_df))
    pfms_df = pfms_df[pfms_df["QValue"] < q_value_threshold]
    pfms_df  = pfms_df[pfms_df["Decoy/Contaminant/Target"] == "T"]
    print("Filtered Proteoforms:", len(pfms_df))
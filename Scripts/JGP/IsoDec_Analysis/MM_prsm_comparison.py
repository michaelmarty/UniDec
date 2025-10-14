import os

from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection_fast
from unidec.IsoDec.datatools import *
import seaborn as sns
import fnmatch
from matplotlib import rcParams
from unidec.IsoDec.match import *

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def parse_full_sequence(full_sequence):
    sequence = ""
    mods = []

    aa_number = 0
    reading_mod = False
    current_mod = ""
    for i in range(len(full_sequence)):
        if full_sequence[i] == "[":
            reading_mod = True
        elif full_sequence[i] == "]":
            reading_mod = False
            mods.append((aa_number, current_mod))
            current_mod = ""

        elif reading_mod:
            current_mod += full_sequence[i]
        else:
            sequence += full_sequence[i]
            aa_number += 1
    return sequence, mods

def compare_mods(mods1, mods2, match_modloc=True):
    if len(mods1) != len(mods2):
        return False
    for i in range(len(mods1)):
        if match_modloc:
            if mods1[i][0] != mods2[i][0]:
                return False
        if mods1[i][1] != mods2[i][1]:
            return False
    return True

def compare_rows(row1, row2, match_scan=True, match_charge=True, match_mods=True, match_modloc=True, allow_ambiguous=True):
    ##Check ambiguity by seeing if the Full Sequence contains a "|"
    psm1_isambiguous = row1["Full Sequence"].find("|") != -1
    psm2_isambiguous = row2["Full Sequence"].find("|") != -1

    ##First check if either psm is ambiguous and if we are allowing ambiguous psms
    if not allow_ambiguous and (psm1_isambiguous or psm2_isambiguous):
        return False

    psm1_scan = int(row1["Scan Number"])
    psm2_scan = int(row2["Scan Number"])
    scan_matched = psm1_scan == psm2_scan
    charge_matched = row1["Precursor Charge"] == row2["Precursor Charge"]

    ##If neither psm is ambiguous, we can compare them directly
    if not psm1_isambiguous and not psm2_isambiguous:
        sequence1, mods1 = parse_full_sequence(row1["Full Sequence"])
        sequence2, mods2 = parse_full_sequence(row2["Full Sequence"])
        seq_matched = sequence1 == sequence2
        mods_matched = compare_mods(mods1, mods2, match_modloc=match_modloc)

        if not match_mods:
            mods_matched = True
        if not match_scan:
            scan_matched = True
        if not match_charge:
            charge_matched = True

        return scan_matched and charge_matched and seq_matched and mods_matched

    elif psm1_isambiguous and not psm2_isambiguous:
        psm1_sequences = row1["Full Sequence"].split("|")
        psm1_seqs, psm1_mods = [], []
        for seq in psm1_sequences:
            seq, mods = parse_full_sequence(seq)
            psm1_seqs.append(seq)
            psm1_mods.append(mods)
        psm2_seq, psm2_mods = parse_full_sequence(row2["Full Sequence"])

        seq_matched = False
        mods_matched = False

        for i in range(len(psm1_sequences)):
            seq_matched = psm1_seqs[i] == psm2_seq
            mods_matched = compare_mods(psm1_mods[i], psm2_mods, match_modloc=match_modloc)

            if not seq_matched:
                continue

            if seq_matched and mods_matched:
                break

            if seq_matched and not match_mods:
                break

        if not match_scan:
            scan_matched = True
        if not match_charge:
            charge_matched = True
        if not match_mods:
            mods_matched = True

        return scan_matched and charge_matched and seq_matched and mods_matched

    elif not psm1_isambiguous and psm2_isambiguous:
        psm2_sequences = row2["Full Sequence"].split("|")
        psm2_seqs, psm2_mods = [], []
        for seq in psm2_sequences:
            seq, mods = parse_full_sequence(seq)
            psm2_seqs.append(seq)
            psm2_mods.append(mods)
        psm1_seq, psm1_mods = parse_full_sequence(row2["Full Sequence"])


        seq_matched = False
        mods_matched = False

        for i in range(len(psm2_sequences)):
            seq_matched = psm2_seqs[i] == psm1_seq
            mods_matched = compare_mods(psm2_mods[i], psm1_mods, match_modloc=match_modloc)

            if not seq_matched:
                continue

            if seq_matched and mods_matched:
                break

            if seq_matched and not match_mods:
                break

        if not match_scan:
            scan_matched = True
        if not match_charge:
            charge_matched = True
        if not match_mods:
            mods_matched = True

        return scan_matched and charge_matched and seq_matched and mods_matched

    elif psm1_isambiguous and psm2_isambiguous:
        psm1_sequences = row1["Full Sequence"].split("|")
        psm1_seqs, psm1_mods = [], []
        for seq in psm1_sequences:
            seq, mods = parse_full_sequence(seq)
            psm1_seqs.append(seq)
            psm1_mods.append(mods)

        psm2_sequences = row2["Full Sequence"].split("|")
        psm2_seqs, psm2_mods = [], []
        for seq in psm2_sequences:
            seq, mods = parse_full_sequence(seq)
            psm2_seqs.append(seq)
            psm2_mods.append(mods)

        seq_matched = False
        mods_matched = False

        for i in range(len(psm1_sequences)):
            for j in range(len(psm2_sequences)):
                seq_matched = psm1_seqs[i] == psm2_seqs[j]
                mods_matched = compare_mods(psm1_mods[i], psm2_mods[j], match_modloc=match_modloc)

                if not seq_matched:
                    continue

                if seq_matched and mods_matched:
                    break

                if seq_matched and not match_mods:
                    break

            if seq_matched and mods_matched:
                break

        if not match_scan:
            scan_matched = True
        if not match_charge:
            charge_matched = True
        if not match_mods:
            mods_matched = True

        return scan_matched and charge_matched and seq_matched and mods_matched

    return False

def find_amb_unamb_pairs(psms1, psms2, psm1_masses, psm2_masses):
    pairs = []
    print("Finding Ambiguous/Unambiguous Pairs")
    print("Total PSMs in File 1: ", len(psms1))
    print("Total PSMs in File 2: ", len(psms2))
    for i,row1 in psms1.iterrows():
        within_tol = fastwithin_abstol(psm2_masses, float(row1["Peptide Monoisotopic Mass"]), 0.1)
        for j in within_tol:
            row2 = psms2.iloc[j]
            if (int(row1["Scan Number"]) == int(row2["Scan Number"]) and
                    int(row1["Precursor Charge"]) == int(row2["Precursor Charge"]) and
                    row1["Protein Accession"] != row2["Protein Accession"]):
                psm1_ambiguous = row1["Full Sequence"].find("|") != -1
                psm2_ambiguous = row2["Full Sequence"].find("|") != -1
                if not psm1_ambiguous and psm2_ambiguous:
                    pairs.append([row1, row2])

    for pair in pairs:
        if "Hydroxylation" in pair[0]["Full Sequence"] or "Oxidation" in pair[0]["Full Sequence"]:
            continue
        else:
            print("Accession:", pair[0]["Protein Accession"])
            print("Scan Number File 1:", pair[0]["Scan Number"])
            print("Scan Number File 2:", pair[1]["Scan Number"])
            print("Protein Name:", pair[0]["Protein Name"])
            print("Unambiguous: ", parse_full_sequence(pair[0]["Full Sequence"])[1])
            amb_seqs = pair[1]["Full Sequence"].split("|")
            for seq in amb_seqs:
                print("Ambiguous: ", parse_full_sequence(seq)[1])
    print("Total Pairs: ", len(pairs))
    return pairs

def load_precursors():
    topdir = "Z:\\Group Share\\JGP\\PXD019247\\precursors\\"
    classic_files = match_files(topdir, "*classic.txt")
    classic_precursors = []
    for file in classic_files:
        classic_df = pd.read_csv(topdir + file, sep='\t')
        masses = classic_df["Precursor Mass (Da)"]
        for mass in masses:
            classic_precursors.append(mass)


    isodec_files = match_files(topdir, "*isodec.txt")
    isodec_precursors = []
    for file in isodec_files:
        isodec_df = pd.read_csv(topdir + file, sep='\t')
        masses = isodec_df["Precursor Mass (Da)"]
        for mass in masses:
            isodec_precursors.append(mass)

    return classic_precursors, isodec_precursors

def compare_prsms(psms1, psm1_peptidemasses, psms2, psm2_peptidemasses,
                  match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename):
    # Get the sorted indices for the peptide masses
    psm1_sorted_indices = np.argsort(psm1_peptidemasses)
    psm2_sorted_indices = np.argsort(psm2_peptidemasses)

    # Reorder the peptide masses and Psms
    psms1 = psms1.iloc[psm1_sorted_indices]
    psm1_peptidemasses = psm1_peptidemasses[psm1_sorted_indices]
    psms2 = psms2.iloc[psm2_sorted_indices]
    psm2_peptidemasses = psm2_peptidemasses[psm2_sorted_indices]

    psms1.reset_index(drop=True, inplace=True)
    psms2.reset_index(drop=True, inplace=True)

    # Now order the psms by peptide mass

    psms1_count = len(psms1)
    psms2_count = len(psms2)


    matched_psms1 = []
    unique_psms1 = []

    matched_psms2 = []
    unique_psms2 = []

    ambiguous_psms1 = []
    ambiguous_psms2 = []

    for i, row1 in psms1.iterrows():
        ambiguous = row1["Full Sequence"].find("|") != -1
        if ambiguous:
            ambiguous_psms1.append(i)
        if not allow_ambiguous and ambiguous:
            continue
        psm1_matched = False
        indices = fastwithin_abstol(psm2_peptidemasses, float(row1["Peptide Monoisotopic Mass"]), 0.1)
        for j in indices:
            row2 = psms2.iloc[j]
            if match_filename:
                if row1["File Name"] != row2["File Name"]:
                    continue
            if compare_rows(row1, row2, match_scan=match_scan, match_charge=match_charge, match_modloc=match_modloc,
                            match_mods=match_mods, allow_ambiguous=allow_ambiguous):
                if i not in matched_psms1:
                    matched_psms1.append(i)
                if j not in matched_psms2:
                    matched_psms2.append(j)
                psm1_matched = True
        if not psm1_matched:
            unique_psms1.append(i)

    # Now populate unique_psms2 with everything not in matched_psms2
    for i, row2 in psms2.iterrows():
        ambiguous = row2["Full Sequence"].find("|") != -1
        if ambiguous:
            ambiguous_psms2.append(i)
        if i in matched_psms2:
            continue
        elif not allow_ambiguous and ambiguous:
            continue
        unique_psms2.append(i)

    print("Total PSMs in File 1: ", psms1_count)
    print("Total PSMs in File 2: ", psms2_count)
    print("Total Unambiguous PSMs in File 1: ", psms1_count - len(ambiguous_psms1))
    print("Total Unambiguous PSMs in File 2: ", psms2_count - len(ambiguous_psms2))
    print("Total Matched PSMs in File 1: ", len(matched_psms1))
    print("Total Matched PSMs in File 2: ", len(matched_psms2))
    print("Total Unique PSMs in File 1: ", len(unique_psms1))
    print("Total Unique PSMs in File 2: ", len(unique_psms2))
    print("Total Ambiguous PSMs in File 1: ", len(ambiguous_psms1))
    print("Total Ambiguous PSMs in File 2: ", len(ambiguous_psms2))
    print("Fraction Ambiguous PSMs in File 1: ", len(ambiguous_psms1) / psms1_count)
    print("Fraction Ambiguous PSMs in File 2: ", len(ambiguous_psms2) / psms2_count)

    if False:
        #manually make a histogram of the proteoform masses
        min_mass = 0
        max_mass = 30000
        bin_width = 1000
        bins = np.arange(min_mass, max_mass, bin_width)
        mmpsms_counts = np.zeros(len(bins) - 1)
        mmprecursors_counts = np.zeros(len(bins) - 1)
        isodecpsms_counts = np.zeros(len(bins) - 1)
        isodecprecursors_counts = np.zeros(len(bins) - 1)

        classic_masses, isodec_masses = load_precursors()

        mm_prsm_masses = []
        isodec_prsm_masses = []
        for i, row in psms1.iterrows():
            mass = float(row["Peptide Monoisotopic Mass"])
            isodec_prsm_masses.append(mass)
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    isodecpsms_counts[j] += 1
                    break

        for i, row in psms2.iterrows():
            mass = float(row["Peptide Monoisotopic Mass"])
            mm_prsm_masses.append(mass)
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    mmpsms_counts[j] += 1
                    break

        for mass in classic_masses:
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    mmprecursors_counts[j] += 1
                    break

        for mass in isodec_masses:
            for j in range(len(bins - 1)):
                if mass >= bins[j] and mass < bins[j + 1]:
                    isodecprecursors_counts[j] += 1
                    break

        # #Perform a mann-whitney U test to compare the precursor distributions and the PrSM distributions.
        # stat_id, p_value_id = ks_2samp(isodec_masses, isodec_prsm_masses)
        # print("Kolmogorov-Smirnov Test Results--IsoDec--")
        # print("Statistic: ", stat_id, " p-value: ", p_value_id)
        #
        # stat_mm, p_value_mm  = ks_2samp(classic_masses, mm_prsm_masses)
        # print("Kolmogorov-Smirnov Test Results--MM-Classic--")
        # print("Statistic: ", stat_mm, " p-value: ", p_value_mm)
        # exit()
        #Normalize counts to the total number of PSMs
        mmpsms_counts /= psms1_count
        isodecpsms_counts /= psms2_count
        mmprecursors_counts /= len(classic_masses)
        isodecprecursors_counts /= len(isodec_masses)

        bar_width = 0.20

        index = np.arange(len(bins[:-1]))

        fig, ax = plt.subplots()

        bar1 = ax.bar(index - 1.5 * bar_width, mmpsms_counts, bar_width, label="MM-PrSMs")
        bar2 = ax.bar(index - 0.5 * bar_width, mmprecursors_counts, bar_width, label="MM-Precursors")
        bar3 = ax.bar(index + 0.5 * bar_width, isodecpsms_counts, bar_width, label="IsoDec-PrSMs")
        bar4 = ax.bar(index + 1.5 * bar_width, isodecprecursors_counts, bar_width, label="Isodec-Precursors")

        ax.set_xlabel("Mass (Da)")
        ax.set_ylabel("Frequency")

        ax.set_xticks(index)
        ax.set_xticklabels(str(bins[i]) for i in range(len(bins[:-1])))

        ax.legend()
        plt.show()

        data = {'Category':[],
                'Value':[]}


        for i in range(len(psm2_peptidemasses)):
            data['Category'].append('MM-PrSMs')
            data['Value'].append(psm2_peptidemasses[i])

        for i in range(len(classic_masses)):
            data['Category'].append('MM-Precursors')
            data['Value'].append(classic_masses[i])

        for i in range(len(psm1_peptidemasses)):
            data['Category'].append('IsoDec-PrSMs')
            data['Value'].append(psm1_peptidemasses[i])

        for i in range(len(isodec_masses)):
            data['Category'].append('IsoDec-Precursors')
            data['Value'].append(isodec_masses[i])


        data = pd.DataFrame(data)

        print("Percent MM Precursors Matched:", (psms2_count / len(classic_masses))*100)
        print("Percent IsoDec Precursors Matched:", (psms1_count / len(isodec_masses))*100)

        plt.figure(figsize=(6,3))
        sns.violinplot(x='Value', y='Category', data=data)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlabel("Mass (Da)", fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylabel("")
        plt.xlim(0,30000)
        plt.tight_layout()
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\precursor_violinplot.png", dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\precursor_violinplot.pdf", dpi=300)
        plt.show()

    return matched_psms1, matched_psms2, unique_psms1, unique_psms2

def compare_prsms_return_uniqueIDs(psms1, psm1_peptidemasses, psms2, psm2_peptidemasses,
                  match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename):
    # Get the sorted indices for the peptide masses
    psm1_sorted_indices = np.argsort(psm1_peptidemasses)
    psm2_sorted_indices = np.argsort(psm2_peptidemasses)

    # Reorder the peptide masses and Psms
    psms1 = psms1.iloc[psm1_sorted_indices]
    psm1_peptidemasses = psm1_peptidemasses[psm1_sorted_indices]
    psms2 = psms2.iloc[psm2_sorted_indices]
    psm2_peptidemasses = psm2_peptidemasses[psm2_sorted_indices]

    psms1.reset_index(drop=True, inplace=True)
    psms2.reset_index(drop=True, inplace=True)

    # Now order the psms by peptide mass

    psms1_count = len(psms1)
    psms2_count = len(psms2)


    matched_psms1 = []
    unique_psms1 = []

    matched_psms2 = []
    unique_psms2 = []

    ambiguous_psms1 = []
    ambiguous_psms2 = []

    for i, row1 in psms1.iterrows():
        ambiguous = row1["Full Sequence"].find("|") != -1
        scan_number1 = row1["Scan Number"]
        if ambiguous:
            ambiguous_psms1.append(i)
        if not allow_ambiguous and ambiguous:
            continue
        psm1_matched = False
        indices = fastwithin_abstol(psm2_peptidemasses, float(row1["Peptide Monoisotopic Mass"]), 0.1)
        for j in indices:
            row2 = psms2.iloc[j]
            scan_number2 = row2["Scan Number"]
            if match_filename:
                if row1["File Name"] != row2["File Name"]:
                    continue
            if compare_rows(row1, row2, match_scan=match_scan, match_charge=match_charge, match_modloc=match_modloc,
                            match_mods=match_mods, allow_ambiguous=allow_ambiguous):
                if (str(scan_number1) + ";" + row1["File Name"]) not in matched_psms1:
                    matched_psms1.append(str(scan_number1) + ";" + str(row1["File Name"]))
                if (str(scan_number2) + ";" + row2["File Name"]) not in matched_psms2:
                    matched_psms2.append(str(scan_number2) + ";" + str(row2["File Name"]))
                psm1_matched = True
        if not psm1_matched:
            unique_psms1.append(str(scan_number1) + ";" + str(row1["File Name"]))

    # Now populate unique_psms2 with everything not in matched_psms2
    for i, row2 in psms2.iterrows():
        ambiguous = row2["Full Sequence"].find("|") != -1
        scan_number2 = row2["Scan Number"]
        if ambiguous:
            ambiguous_psms2.append(i)
        if (str(scan_number2) + ";" + str(row2["File Name"])) in matched_psms2:
            continue
        elif not allow_ambiguous and ambiguous:
            continue
        unique_psms2.append(str(scan_number2) + ";" + str(row2["File Name"]))

    print("Total PSMs in File 1: ", psms1_count)
    print("Total PSMs in File 2: ", psms2_count)
    print("Total Unambiguous PSMs in File 1: ", psms1_count - len(ambiguous_psms1))
    print("Total Unambiguous PSMs in File 2: ", psms2_count - len(ambiguous_psms2))
    print("Total Matched PSMs in File 1: ", len(matched_psms1))
    print("Total Matched PSMs in File 2: ", len(matched_psms2))
    print("Total Unique PSMs in File 1: ", len(unique_psms1))
    print("Total Unique PSMs in File 2: ", len(unique_psms2))
    print("Total Ambiguous PSMs in File 1: ", len(ambiguous_psms1))
    print("Total Ambiguous PSMs in File 2: ", len(ambiguous_psms2))
    print("Fraction Ambiguous PSMs in File 1: ", len(ambiguous_psms1) / psms1_count)
    print("Fraction Ambiguous PSMs in File 2: ", len(ambiguous_psms2) / psms2_count)

    if False:
        #manually make a histogram of the proteoform masses
        min_mass = 0
        max_mass = 30000
        bin_width = 1000
        bins = np.arange(min_mass, max_mass, bin_width)
        mmpsms_counts = np.zeros(len(bins) - 1)
        mmprecursors_counts = np.zeros(len(bins) - 1)
        isodecpsms_counts = np.zeros(len(bins) - 1)
        isodecprecursors_counts = np.zeros(len(bins) - 1)

        classic_masses, isodec_masses = load_precursors()

        mm_prsm_masses = []
        isodec_prsm_masses = []
        for i, row in psms1.iterrows():
            mass = float(row["Peptide Monoisotopic Mass"])
            isodec_prsm_masses.append(mass)
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    isodecpsms_counts[j] += 1
                    break

        for i, row in psms2.iterrows():
            mass = float(row["Peptide Monoisotopic Mass"])
            mm_prsm_masses.append(mass)
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    mmpsms_counts[j] += 1
                    break

        for mass in classic_masses:
            for j in range(len(bins) - 1):
                if mass >= bins[j] and mass < bins[j + 1]:
                    mmprecursors_counts[j] += 1
                    break

        for mass in isodec_masses:
            for j in range(len(bins - 1)):
                if mass >= bins[j] and mass < bins[j + 1]:
                    isodecprecursors_counts[j] += 1
                    break

        # #Perform a mann-whitney U test to compare the precursor distributions and the PrSM distributions.
        # stat_id, p_value_id = ks_2samp(isodec_masses, isodec_prsm_masses)
        # print("Kolmogorov-Smirnov Test Results--IsoDec--")
        # print("Statistic: ", stat_id, " p-value: ", p_value_id)
        #
        # stat_mm, p_value_mm  = ks_2samp(classic_masses, mm_prsm_masses)
        # print("Kolmogorov-Smirnov Test Results--MM-Classic--")
        # print("Statistic: ", stat_mm, " p-value: ", p_value_mm)
        # exit()
        #Normalize counts to the total number of PSMs
        mmpsms_counts /= psms1_count
        isodecpsms_counts /= psms2_count
        mmprecursors_counts /= len(classic_masses)
        isodecprecursors_counts /= len(isodec_masses)

        bar_width = 0.20

        index = np.arange(len(bins[:-1]))

        fig, ax = plt.subplots()

        bar1 = ax.bar(index - 1.5 * bar_width, mmpsms_counts, bar_width, label="MM-PrSMs")
        bar2 = ax.bar(index - 0.5 * bar_width, mmprecursors_counts, bar_width, label="MM-Precursors")
        bar3 = ax.bar(index + 0.5 * bar_width, isodecpsms_counts, bar_width, label="IsoDec-PrSMs")
        bar4 = ax.bar(index + 1.5 * bar_width, isodecprecursors_counts, bar_width, label="Isodec-Precursors")

        ax.set_xlabel("Mass (Da)")
        ax.set_ylabel("Frequency")

        ax.set_xticks(index)
        ax.set_xticklabels(str(bins[i]) for i in range(len(bins[:-1])))

        ax.legend()
        plt.show()

        data = {'Category':[],
                'Value':[]}


        for i in range(len(psm2_peptidemasses)):
            data['Category'].append('MM-PrSMs')
            data['Value'].append(psm2_peptidemasses[i])

        for i in range(len(classic_masses)):
            data['Category'].append('MM-Precursors')
            data['Value'].append(classic_masses[i])

        for i in range(len(psm1_peptidemasses)):
            data['Category'].append('IsoDec-PrSMs')
            data['Value'].append(psm1_peptidemasses[i])

        for i in range(len(isodec_masses)):
            data['Category'].append('IsoDec-Precursors')
            data['Value'].append(isodec_masses[i])


        data = pd.DataFrame(data)

        print("Percent MM Precursors Matched:", (psms2_count / len(classic_masses))*100)
        print("Percent IsoDec Precursors Matched:", (psms1_count / len(isodec_masses))*100)

        plt.figure(figsize=(6,3))
        sns.violinplot(x='Value', y='Category', data=data)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlabel("Mass (Da)", fontsize=16)
        plt.yticks(fontsize=16)
        plt.ylabel("")
        plt.xlim(0,30000)
        plt.tight_layout()
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\precursor_violinplot.png", dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\precursor_violinplot.pdf", dpi=300)
        plt.show()

    return matched_psms1, matched_psms2, unique_psms1, unique_psms2

def get_start_end(range_cell):
    range_cell = range_cell.strip("[]")
    split = range_cell.split("to")
    start = int(split[0])
    end = int(split[1])
    return start, end

def compare_matchedions(psms1, psm1_peptidemasses, psms2, psm2_peptidemasses):
    #I want to be able to compare the matched ions between two sets of psms to determine which ions are unique to each set
    #Additionally, I want to determine which are shared between the two sets of psms
    #I want to see which scan numbers are unique to each set of psms as well

    #First, sort the psms by peptide mass
    psm1_sorted_indices = np.argsort(psm1_peptidemasses)
    psm2_sorted_indices = np.argsort(psm2_peptidemasses)

    psms1 = psms1.iloc[psm1_sorted_indices]
    psm1_peptidemasses = psm1_peptidemasses[psm1_sorted_indices]
    psms2 = psms2.iloc[psm2_sorted_indices]
    psm2_peptidemasses = psm2_peptidemasses[psm2_sorted_indices]

    psms1.reset_index(drop=True, inplace=True)
    psms2.reset_index(drop=True, inplace=True)

    shared_scans = []
    unique_scans1 = []
    unique_scans2 = []

    shared_ions = []
    unique_ions1 = []
    unique_ions2 = []
    coverages1 = []
    coverages2 = []

    for i, row1 in psms1.iterrows():
        scan1 = int(row1["Scan Number"])
        charge1 = int(row1["Precursor Charge"])
        acc1 = row1["Protein Accession"]
        scan1_matched = False
        start, end = get_start_end(row1["Start and End Residues In Protein"])
        indices = fastwithin_abstol(psm2_peptidemasses, float(row1["Peptide Monoisotopic Mass"]), 0.1)
        for i in indices:
            row2 = psms2.iloc[i]
            if compare_rows(row1, row2, match_scan=True, match_charge=True, match_mods=True, match_modloc=True, allow_ambiguous=False):
                shared_scans.append(scan1)
                scan1_matched = True
                #Now compare the matched ions
                #They first need to be split on ';', then each of those needs to be stripped of the brackets and split on ','
                ions1_firstsplit = row1["Matched Ion Series"].split(";")
                ions2_firstsplit = row2["Matched Ion Series"].split(";")
                ions1 = []
                ions2 = []


                for ionsplit in ions1_firstsplit:
                    ions1.extend(ionsplit.strip("[]").split(","))
                for ionsplit in ions2_firstsplit:
                    ions2.extend(ionsplit.strip("[]").split(","))


                length = len(row1["Base Sequence"])
                shared, unique1, unique2, cov1, cov2 = compare_matchedioncells(ions1, ions2, start, end, length)
                shared_ions.append(shared)
                unique_ions1.append(unique1)
                unique_ions2.append(unique2)
                coverages1.append(cov1)
                coverages2.append(cov2)
                break
        if not scan1_matched:
            unique_scans1.append(scan1)

    print("Average Shared Ions:", np.mean(shared_ions))
    print("Average Unique Ions PSM 1:", np.mean(unique_ions1))
    print("Average Unique Ions PSM 2:", np.mean(unique_ions2))
    print("Average Coverage PSM 1:", np.mean(coverages1))
    print("Standard Deviation Coverage PSM 1:", np.std(coverages1))
    print("Average Coverage PSM 2:", np.mean(coverages2))
    print("Standard Deviation Coverage PSM 2:", np.std(coverages2))
    print("Total Shared Scans:", len(shared_scans))

    return shared_scans, shared_ions, unique_ions1, unique_ions2, coverages1, coverages2

def compare_matchedioncells(ions1, ions2, start, end, length):
    shared_ions = []
    unique_ions1 = []
    unique_ions2 = []

    for ion1 in ions1:
        ion1_matched = False
        for ion2 in ions2:
            if ion1 == ion2:
                shared_ions.append(ion1)
                ion1_matched = True
                break
        if not ion1_matched:
            unique_ions1.append(ion1)

    for ion2 in ions2:
        ion2_matched = False
        for ion1 in ions1:
            if ion2 == ion1:
                ion2_matched = True
                break
        if not ion2_matched:
            unique_ions2.append(ion2)

    #From the shared/unique ions, we can determine the shared/unique fragmentation locations
    shared_locs = []
    unique_locs1 = []
    unique_locs2 = []
    ion_types = ["b", "c", "y", "zDot"]

    fragmentation_locs1 = []

    for ion in ions1:
        ionmatched = False
        split = ion.split("+")
        #Now check if the ion contains the ion type
        for ion_type in ion_types:
            if ion_type in split[0]:
                ionmatched = True
                #Remove the ion type from the string
                loc = int(split[0].replace(ion_type, ""))
                if ion_type == "y" or ion_type == "zDot":
                    fragloc = int(length - loc + (start - 1))
                    if fragloc <= 0:
                        break
                    else:
                        fragmentation_locs1.append(fragloc)
                        break
                else:
                    fragmentation_locs1.append(int(loc))
                    break
        if not ionmatched:
            print("Ion not matched:", ion)
            print(split[0])


    fragmentation_locs2 = []
    for ion in ions2:
        ionmatched = False
        split = ion.split("+")
        for ion_type in ion_types:
            if ion_type in split[0]:
                ionmatched = True
                loc = int(split[0].replace(ion_type, ""))
                if ion_type == "y" or ion_type == "zDot":
                    fragloc = int(length - loc + (start - 1))
                    if fragloc <= 0:
                        break
                    else:
                        fragmentation_locs2.append(fragloc)
                        break
                else:
                    fragmentation_locs2.append(int(loc))
                    break
        if not ionmatched:
            print("Ion not matched:", ion)
            print(split[0])

    sharedlocs = []
    uniquelocs1 = []
    uniquelocs2 = []

    for loc1 in fragmentation_locs1:
        loc1_matched = False
        for loc2 in fragmentation_locs2:
            if loc1 == loc2:
                sharedlocs.append(loc1)
                loc1_matched = True
                break
        if not loc1_matched:
            uniquelocs1.append(loc1)

    for loc2 in fragmentation_locs2:
        loc2_matched = False
        for loc1 in fragmentation_locs1:
            if loc2 == loc1:
                loc2_matched = True
                break
        if not loc2_matched:
            uniquelocs2.append(loc2)


    locs1 = set(fragmentation_locs1)
    locs2 = set(fragmentation_locs2)

    coverage1 = (len(locs1) / (length - 1)) * 100
    coverage2 = (len(locs2) / (length - 1)) * 100

    return len(sharedlocs), len(uniquelocs1), len(uniquelocs2), coverage1, coverage2

def get_coverage(ions, start, length):
    fragmentation_locs = []
    ion_types = ["b", "c", "y", "zDot"]
    for ion in ions:
        ionmatched = False
        split = ion.split("+")
        for ion_type in ion_types:
            if ion_type in split[0]:
                ionmatched = True
                loc = int(split[0].replace(ion_type, ""))
                if ion_type == "y" or ion_type == "zDot":
                    fragloc = int(length - loc + (start - 1))
                    if fragloc <= 0:
                        break
                    else:
                        fragmentation_locs.append(fragloc)
                        break
                else:
                    fragmentation_locs.append(int(loc))
                    break
        if not ionmatched:
            print("Ion not matched:", ion)
            print(split[0])

    locs = set(fragmentation_locs)
    coverage = (len(locs) / (length - 1)) * 100

    return coverage

def assess_decon_accuracy(psms, msalignfiles_directory, matchedscans=None):
    percents_matched = {}
    coverages = {}

    matched_scan_numbers = []
    matched_filenames = []

    if matchedscans is not None:
        for i in range(len(matchedscans)):
            split = matchedscans[i].split(";")
            matched_scan_numbers.append(int(split[0]))
            matched_filenames.append(split[1])

    ms2_msalign_files = match_files(msalignfiles_directory, "*ms2.msalign")

    msalign_dict = {}
    for file in ms2_msalign_files:
        #Get the filename without extension
        filename = file.split('.')[0].replace("_ms2", "")
        print("Reading MSAlign File:", filename)
        fullpath = msalignfiles_directory + file
        msalign_dict[filename] = read_msalign_to_matchedcollection_fast(fullpath, mz_type="peak")

    for i, row in psms.iterrows():
        file = row["File Name"]
        filename = file.split('.')[0].replace("_ms2", "")
        scan = int(row["Scan Number"])
        #Need to check of this combination of scan and file exists in the matched scans
        #but only if matchedscans is not None
        prsm_id = str(scan) + ";" + str(file)
        if matchedscans is not None and prsm_id not in matchedscans:
            continue
        if file in msalign_dict:
            ms2_mc = msalign_dict[filename]
            # Get every peak in the matched collection that has the same scan number
            peaks = [peak for peak in ms2_mc.peaks if peak.scan == scan]

            mi_int_split = row["Matched Ion Intensities"].split(";")
            mi_intensities = []
            for type in mi_int_split:
                mi_intensities.extend(type.strip("[]").split(","))

            mi_intensities = [float(intensity.split(":")[1]) for intensity in mi_intensities]
            mi_intensities = set(mi_intensities)

            perc_matched = (len(mi_intensities) / len(peaks)) * 100
            percents_matched[str(scan) + ";" + str(file)] = perc_matched

            #Get the coverage of the matched ions
            start, end = get_start_end(row["Start and End Residues In Protein"])
            length = len(row["Base Sequence"])
            ions_firstsplit = row["Matched Ion Series"].split(";")
            ions = []


            for ionsplit in ions_firstsplit:
                ions.extend(ionsplit.strip("[]").split(","))

            coverage = get_coverage(ions, start, length)
            coverages[str(scan) + ";" + str(file)] = coverage
        else:
            print("File: ", filename,  " not found.")

    return percents_matched, coverages

def prepare_prsms(prsm_filepath, qval_cutoff, allow_ambiguous):
    prsms = pd.read_csv(prsm_filepath, sep="\t")
    prsms = prsms[prsms["QValue"] <= qval_cutoff]
    prsms = prsms[prsms["Peptide Monoisotopic Mass"].astype(str).str.find("|") == -1]
    if not allow_ambiguous:
        prsms = prsms[prsms["Full Sequence"].astype(str).str.find("|") == -1]
    prsms.reset_index(drop=True, inplace=True)
    peptidemasses = prsms["Peptide Monoisotopic Mass"]
    peptidemasses = np.array(peptidemasses.astype(float))
    sorted_indices = np.argsort(peptidemasses)
    prsms = prsms.iloc[sorted_indices]
    peptidemasses = peptidemasses[sorted_indices]
    new_names = np.array([])
    for filename in prsms["File Name"]:
        new_name = filename.split('.')[0].replace("_ms2", "")
        new_names = np.append(new_names, new_name)
    prsms["File Name"] = new_names
    return prsms, peptidemasses




# Set backend to Agg
mpl.use('WxAgg')

#Cannot match_modloc without match_mods
qval_cutoff = 0.01
match_modloc = True
match_mods = True
match_scan = True
match_charge = True
allow_ambiguous = True

#Used when comparing AllPsms files.
match_filename = True

if False:
    topfd_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\topFD_msalign\TFD\Task2-SearchTask\AllPSMs.psmtsv"
    FLASHDeconv_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\FLASHDeconv_msalign\FLD\Task2-SearchTask\AllPSMs.psmtsv"
    isodec_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\msalign_for_toppic\8bin\ID_8\Task2-SearchTask\AllPSMs.psmtsv"

    topfd_prsms, topfd_peptidemasses = prepare_prsms(topfd_psms_file, qval_cutoff, allow_ambiguous)
    FLASHDeconv_prsms, FLASHDeconv_peptidemasses = prepare_prsms(FLASHDeconv_psms_file, qval_cutoff, allow_ambiguous)
    isodec_prsms, isodec_peptidemasses = prepare_prsms(isodec_psms_file, qval_cutoff, allow_ambiguous)


    topfd_scans = topfd_prsms["Scan Number"].astype(int)
    topfd_files = topfd_prsms["File Name"]
    topfd_scans_files = []
    for i in range(len(topfd_scans)):
        topfd_scans_files.append(str(topfd_scans[i]) + "_" + topfd_files[i])

    FLASHDeconv_scans = FLASHDeconv_prsms["Scan Number"].astype(int)
    FLASHDeconv_files = FLASHDeconv_prsms["File Name"]
    FLASHDeconv_scans_files = []
    for i in range(len(FLASHDeconv_scans)):
        FLASHDeconv_scans_files.append(str(FLASHDeconv_scans[i]) + "_" + FLASHDeconv_files[i])


    isodec_scans = isodec_prsms["Scan Number"].astype(int)
    isodec_files = isodec_prsms["File Name"]
    isodec_scans_files = []
    for i in range(len(isodec_scans)):
        isodec_scans_files.append(str(isodec_scans[i]) + "_" + isodec_files[i])

    #Get the overlapped scans and files between the three files


    #Get the overlapped scans between the three files
    uniquescans_topfd = [i for i in topfd_scans_files if i not in FLASHDeconv_scans_files and i not in isodec_scans_files]
    uniquescans_flashdeconv = [i for i in FLASHDeconv_scans_files if i not in topfd_scans_files and i not in isodec_scans_files]
    uniquescans_isodec = [i for i in isodec_scans_files if i not in topfd_scans_files and i not in FLASHDeconv_scans_files]

    topfd_flashdeconv_scans = [i for i in topfd_scans_files if i in FLASHDeconv_scans_files and i not in isodec_scans_files]
    topfd_isodec_scans = [i for i in topfd_scans_files if i in isodec_scans_files and i not in FLASHDeconv_scans_files]
    flashdeconv_isodec_scans = [i for i in FLASHDeconv_scans_files if i in isodec_scans_files and i not in topfd_scans_files]

    scan_overlap_all = [i for i in topfd_scans_files if i in FLASHDeconv_scans_files and i in isodec_scans_files]

    print("Total Scans in TopFD: ", len(topfd_scans))
    print("Total Scans in FLASHDeconv: ", len(FLASHDeconv_scans))
    print("Total Scans in IsoDec: ", len(isodec_scans))

    print("Scans Unique to TopFD: ", len(uniquescans_topfd))
    print("Scans Unique to FLASHDeconv: ", len(uniquescans_flashdeconv))
    print("Scans Unique to IsoDec: ", len(uniquescans_isodec))
    print("Scan Overlap between TopFD and FLASHDeconv, but not IsoDec: ", len(topfd_flashdeconv_scans))
    print("Scan Overlap between TopFD and IsoDec, but not FLASHDeconv: ", len(topfd_isodec_scans))
    print("Scan Overlap between FLASHDeconv and IsoDec, but not TopFD: ", len(flashdeconv_isodec_scans))
    print("Scan Overlap between all three: ", len(scan_overlap_all))

    #First get overlap between each pair of files
    print("Comparing TopFD and FLASHDeconv")
    topfd_vs_FLASHDeconv = compare_prsms(topfd_prsms, topfd_peptidemasses, FLASHDeconv_prsms, FLASHDeconv_peptidemasses,
                                         match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing TopFD and IsoDec")
    topfd_vs_isodec = compare_prsms(topfd_prsms, topfd_peptidemasses, isodec_prsms, isodec_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing FLASHDeconv and IsoDec")
    FLASHDeconv_vs_isodec = compare_prsms(FLASHDeconv_prsms, FLASHDeconv_peptidemasses, isodec_prsms, isodec_peptidemasses,
                                          match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)



    #Needed regions
    #Unique to TopFD
    #Unique to FLASHDeconv
    #Unique to IsoDec
    #Overlap between TopFD and FLASHDeconv, but not IsoDec
    #Overlap between TopFD and IsoDec, but not FLASHDeconv
    #Overlap between FLASHDeconv and IsoDec, but not TopFD
    #Overlap between all three

    unique_tfd = [i for i in topfd_vs_FLASHDeconv[2] if i in topfd_vs_isodec[2]]
    unique_fd = [i for i in topfd_vs_FLASHDeconv[3] if i in FLASHDeconv_vs_isodec[2]]
    unique_id = [i for i in topfd_vs_isodec[3] if i in FLASHDeconv_vs_isodec[3]]

    overlap_tfd_fd = [i for i in topfd_vs_FLASHDeconv[0] if i in topfd_vs_isodec[2]]
    overlap_tfd_id = [i for i in topfd_vs_isodec[0] if i in topfd_vs_FLASHDeconv[2]]
    overlap_fd_id = [i for i in FLASHDeconv_vs_isodec[0] if i in topfd_vs_FLASHDeconv[2]]

    overlap_all = [i for i in topfd_vs_FLASHDeconv[0] if i in topfd_vs_isodec[0]]

    print("Unique to TopFD: ", len(unique_tfd))
    print("Unique to FLASHDeconv: ", len(unique_fd))
    print("Unique to IsoDec: ", len(unique_id))
    print("Overlap between TopFD and FLASHDeconv, but not IsoDec: ", len(overlap_tfd_fd))
    print("Overlap between TopFD and IsoDec, but not FLASHDeconv: ", len(overlap_tfd_id))
    print("Overlap between FLASHDeconv and IsoDec, but not TopFD: ", len(overlap_fd_id))
    print("Overlap between all three: ", len(overlap_all))
    exit()

if False:
    isodec_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\msalign_for_toppic\8bin\ID_8\Task2-SearchTask\AllPSMs.psmtsv"
    topfd_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\topFD_msalign\TFD\Task2-SearchTask\AllPSMs.psmtsv"
    fld_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\FLASHDeconv_msalign\FLD\Task2-SearchTask\AllPSMs.psmtsv"

    topfd_prsms, topfd_peptidemasses = prepare_prsms(topfd_psms_file, qval_cutoff, allow_ambiguous)
    isodec_prsms, isodec_peptidemasses = prepare_prsms(isodec_psms_file, qval_cutoff, allow_ambiguous)
    fld_prsms, fld_peptidemasses = prepare_prsms(fld_psms_file, qval_cutoff, allow_ambiguous)

    topfd_scans = topfd_prsms["Scan Number"].astype(int)
    topfd_files = topfd_prsms["File Name"]
    topfd_scans_files = []
    for i in range(len(topfd_scans)):
        filename = topfd_files[i].split(".")[0]
        filename = filename.replace("_ms2", '')
        topfd_scans_files.append(str(topfd_scans[i]) + "_" + filename)

    isodec_scans = isodec_prsms["Scan Number"].astype(int)
    isodec_files = isodec_prsms["File Name"]
    isodec_scans_files = []
    for i in range(len(isodec_scans)):
        filename = isodec_files[i].split(".")[0]
        filename = filename.replace("_ms2", '')
        isodec_scans_files.append(str(isodec_scans[i]) + "_" + filename)

    fld_scans = fld_prsms["Scan Number"].astype(int)
    fld_files = fld_prsms["File Name"]
    fld_scans_files = []
    for i in range(len(fld_scans)):
        filename = fld_files[i].split(".")[0]
        filename = filename.replace("_ms2", '')
        fld_scans_files.append(str(fld_scans[i]) + "_" + filename)

    #Get the overlapped scans and files between each combination of the three files
    unique_fld = [i for i in fld_scans_files if i not in topfd_scans_files and i not in isodec_scans_files]
    unique_tfd = [i for i in topfd_scans_files if i not in fld_scans_files and i not in isodec_scans_files]
    unique_id = [i for i in isodec_scans_files if i not in topfd_scans_files and i not in fld_scans_files]

    fld_tfd_not_id = [i for i in fld_scans_files if i in topfd_scans_files and i not in isodec_scans_files]
    fld_id_not_tfd = [i for i in fld_scans_files if i in isodec_scans_files and i not in topfd_scans_files]
    tfd_id_not_fld = [i for i in topfd_scans_files if i in isodec_scans_files and i not in fld_scans_files]

    all_three = [i for i in topfd_scans_files if i in isodec_scans_files and i in fld_scans_files]

    print("Total scans in TopFD: ", len(topfd_scans_files))
    print("Total scans in IsoDec: ", len(isodec_scans_files))
    print("Total scans in FLASHDeconv: ", len(fld_scans_files))
    print("Unique to TopFD: ", len(unique_tfd))
    print("Unique to FLD: ", len(unique_fld))
    print("Unique to IsoDec: ", len(unique_id))
    print("TopFD and FLASHDeconv but not IsoDec: ", len(fld_tfd_not_id))
    print("TopFD and IsoDec but not FLASHDeconv: ", len(tfd_id_not_fld))
    print("IsoDec and FLASHDeconv but not TopFD: ", len(fld_id_not_tfd))
    print("All Three", len(all_three))
    exit()

if False:
    mm_id_psm_file = "D:\\Johnny\\IsoDec\\PXD019247\\MM_isodec_inline_noAvg\\Task3-SearchTask\\AllPSMs.psmtsv"
    mm_classic_psm_file = "D:\\Johnny\\IsoDec\\PXD019247\\MM_Classic_noAvg\\Task3-SearchTask\\AllPSMs.psmtsv"

    # mm_id_psm_file = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\isodec_8bit_c\\Task2-SearchTask\\AllPSMs.psmtsv"
    # mm_classic_psm_file = "Z:\\Group Share\\JGP\\PXD019247\\Raw_Test\\Task2-SearchTask\\AllPSMs.psmtsv"

    mm_id_psms = pd.read_csv(mm_id_psm_file, sep="\t")
    mm_classic_psms = pd.read_csv(mm_classic_psm_file, sep="\t")

    mm_id_psms = mm_id_psms[mm_id_psms["Decoy/Contaminant/Target"] == "T"]
    mm_classic_psms = mm_classic_psms[mm_classic_psms["Decoy/Contaminant/Target"] == "T"]

    mm_id_psms = mm_id_psms[mm_id_psms["QValue"] <= qval_cutoff]
    mm_classic_psms = mm_classic_psms[mm_classic_psms["QValue"] <= qval_cutoff]

    mm_id_psms = mm_id_psms[mm_id_psms["Peptide Monoisotopic Mass"].astype(str).str.find("|") == -1]
    mm_classic_psms = mm_classic_psms[mm_classic_psms["Peptide Monoisotopic Mass"].astype(str).str.find("|") == -1]

    mm_id_psms.reset_index(drop=True, inplace=True)
    mm_classic_psms.reset_index(drop=True, inplace=True)

    mm_id_peptidemasses = mm_id_psms["Peptide Monoisotopic Mass"]
    mm_classic_peptidemasses = mm_classic_psms["Peptide Monoisotopic Mass"]

    mm_id_peptidemasses = np.array(mm_id_peptidemasses.astype(float))
    mm_classic_peptidemasses = np.array(mm_classic_peptidemasses.astype(float))

    mm_id_sorted_indices = np.argsort(mm_id_peptidemasses)
    mm_classic_sorted_indices = np.argsort(mm_classic_peptidemasses)

    mm_id_psms = mm_id_psms.iloc[mm_id_sorted_indices]
    mm_id_peptidemasses = mm_id_peptidemasses[mm_id_sorted_indices]
    mm_classic_psms = mm_classic_psms.iloc[mm_classic_sorted_indices]
    mm_classic_peptidemasses = mm_classic_peptidemasses[mm_classic_sorted_indices]

    comp = compare_prsms(mm_id_psms, mm_id_peptidemasses, mm_classic_psms, mm_classic_peptidemasses,
                            match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)
    exit()

if False:
    qval_binwidth = 0.01
    qval_min = 0
    qval_max = qval_cutoff

    qval_bins = []
    qval_binmins = np.arange(qval_min, qval_max, qval_binwidth)
    for i in range(len(qval_binmins)):
        qval_bins.append((qval_binmins[i], qval_binmins[i] + qval_binwidth))


    fracs_matched = [0 for i in range(int((qval_max - qval_min) / qval_binwidth))]

    for bin in qval_bins:
        bin_min, bin_max = bin
        matched_psms_inbin = 0
        total_psms_inbin = 0
        for i,row1 in psms1.iterrows():
            if row1["QValue"] >= bin_min and row1["QValue"] < bin_max:
                total_psms_inbin += 1
                if i in matched_psms1:
                    matched_psms_inbin += 1
        if total_psms_inbin == 0:
            continue
        frac_matched = matched_psms_inbin / total_psms_inbin
        bin_index = int((bin_min - qval_min) / qval_binwidth)
        fracs_matched[bin_index] = frac_matched

    for i in range(len(qval_bins)):
        print("QValue Range: ", qval_bins[i])
        print("Fraction of PSMs Matched: ", fracs_matched[i])

    #Now convert from indices to rows
    matched_psms1 = [psms1.iloc[i] for i in matched_psms1]
    matched_psms2 = [psms2.iloc[i] for i in matched_psms2]
    unique_psms1 = [psms1.iloc[i] for i in unique_psms1]
    unique_psms2 = [psms2.iloc[i] for i in unique_psms2]

    plt.plot(qval_binmins, fracs_matched)
    plt.xlabel("QValue Min")
    plt.ylabel("Fraction of PSMs Matched")
    plt.title("QValue Min vs Fraction of PSMs Matched")
    plt.show()

if False:
    id8_prsmfile = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\isodec_8bit_c\\Task2-SearchTask\\AllPSMs.psmtsv"
    id8_msalign_directory = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\"
    id8_prsms, id8_peptidemasses = prepare_prsms(id8_prsmfile, qval_cutoff, allow_ambiguous)

    topfd_prsmfile = "Z:\\Group Share\\JGP\\PXD019247\\topFD_msalign\\TopFD\\Task2-SearchTask\\AllPSMs.psmtsv"
    topfd_msalign_directory = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\"
    topfd_prsms, topfd_peptidemasses = prepare_prsms(topfd_prsmfile, qval_cutoff, allow_ambiguous)

    print("Comparing IsoDec 8-bit and TopFD")
    id8_vs_topfd = compare_prsms_return_uniqueIDs(id8_prsms, id8_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                                  match_mods, match_modloc, match_scan, match_charge, allow_ambiguous,
                                                  match_filename)

    print("Shared:", len(id8_vs_topfd[0]))
    print("Unique to IsoDec:", len(id8_vs_topfd[2]))
    print("Unique to TopFD:", len(id8_vs_topfd[3]))
    exit()

if False:
    id4_prsmfile = "Z:\\Group Share\\JGP\\PXD019247\\isodec_4bit_c_msalign\\isodec_4bit_c\\Task2-SearchTask\\AllPSMs.psmtsv"
    id4_msalign_directory = "Z:\\Group Share\\JGP\\PXD019247\\isodec_4bit_c_msalign\\"
    id4_prsms, id4_peptidemasses = prepare_prsms(id4_prsmfile, qval_cutoff)

    id8_prsmfile = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\isodec_8bit_c\\Task2-SearchTask\\AllPSMs.psmtsv"
    id8_msalign_directory = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\"
    id8_prsms, id8_peptidemasses = prepare_prsms(id8_prsmfile, qval_cutoff)

    flashdeconv_prsmfile = "D:\\Johnny\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\FLASHDeconv\\Task2-SearchTask\\AllPSMs.psmtsv"
    flashdeconv_msalign_directory = "D:\\Johnny\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\"
    flashdeconv_prsms, flashdeconv_peptidemasses = prepare_prsms(flashdeconv_prsmfile, qval_cutoff)

    topfd_prsmfile = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\TopFD\\Task2-SearchTask\\AllPSMs.psmtsv"
    topfd_msalign_directory = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\"
    topfd_prsms, topfd_peptidemasses = prepare_prsms(topfd_prsmfile, qval_cutoff)

    #First get overlap between each pair of files
    print("Comparing IsoDec 4-bit and IsoDec 8-bit")
    id4_vs_id8 = compare_prsms(id4_prsms, id4_peptidemasses, id8_prsms, id8_peptidemasses,
                                         match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing IsoDec 4-bit and FLASHDeconv")
    id4_vs_flashdeconv = compare_prsms(id4_prsms, id4_peptidemasses, flashdeconv_prsms, flashdeconv_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing IsoDec 4-bit and TopFD")
    id4_vs_topfd = compare_prsms(id4_prsms, id4_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing IsoDec 8-bit and FLASHDeconv")
    id8_vs_flashdeconv = compare_prsms(id8_prsms, id8_peptidemasses, flashdeconv_prsms, flashdeconv_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing IsoDec 8-bit and TopFD")
    id8_vs_topfd = compare_prsms(id8_prsms, id8_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing FLASHDeconv and TopFD")
    flashdeconv_vs_topfd = compare_prsms(flashdeconv_prsms, flashdeconv_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    #Now get the prsms that are matched between all 4 files
    unique_id4 = [i for i in id4_vs_id8[2] if i in id4_vs_flashdeconv[2] and i in id4_vs_topfd[2]]
    unique_id8 = [i for i in id4_vs_id8[3] if i in id8_vs_flashdeconv[2] and i in id8_vs_topfd[2]]
    unique_flashdeconv = [i for i in id4_vs_flashdeconv[3] if i in id8_vs_flashdeconv[3] and i in flashdeconv_vs_topfd[2]]
    unique_topfd = [i for i in id4_vs_topfd[3] if i in id8_vs_topfd[3] and i in flashdeconv_vs_topfd[3]]

    overlap_id4_id8 = [i for i in id4_vs_id8[0] if i in id4_vs_flashdeconv[2] and i in id4_vs_topfd[2]]
    overlap_id4_flashdeconv = [i for i in id4_vs_flashdeconv[0] if i in id4_vs_id8[2] and i in id4_vs_topfd[2]]
    overlap_id4_topfd = [i for i in id4_vs_topfd[0] if i in id4_vs_id8[2] and i in id4_vs_flashdeconv[2]]
    overlap_id8_flashdeconv = [i for i in id8_vs_flashdeconv[0] if i in id4_vs_id8[3] and i in id8_vs_topfd[2]]
    overlap_id8_topfd = [i for i in id8_vs_topfd[0] if i in id4_vs_id8[3] and i in id8_vs_flashdeconv[2]]
    overlap_flashdeconv_topfd = [i for i in flashdeconv_vs_topfd[0] if i in id4_vs_flashdeconv[2] and i in id8_vs_flashdeconv[2]]

    overlap_all = [i for i in id4_vs_id8[0] if i in id4_vs_flashdeconv[0] and i in id4_vs_topfd[0]]

    print("Unique to IsoDec 4-bit: ", len(unique_id4))
    print("Unique to IsoDec 8-bit: ", len(unique_id8))
    print("Unique to FLASHDeconv: ", len(unique_flashdeconv))
    print("Unique to TopFD: ", len(unique_topfd))
    print("Overlap between IsoDec 4-bit and IsoDec 8-bit, but not FLASHDeconv or TopFD: ", len(overlap_id4_id8))
    print("Overlap between IsoDec 4-bit and FLASHDeconv, but not IsoDec 8-bit or TopFD: ", len(overlap_id4_flashdeconv))
    print("Overlap between IsoDec 4-bit and TopFD, but not IsoDec 8-bit or FLASHDeconv: ", len(overlap_id4_topfd))
    print("Overlap between IsoDec 8-bit and FLASHDeconv, but not IsoDec 4-bit or TopFD: ", len(overlap_id8_flashdeconv))
    print("Overlap between IsoDec 8-bit and TopFD, but not IsoDec 4-bit or FLASHDeconv: ", len(overlap_id8_topfd))
    print("Overlap between FLASHDeconv and TopFD, but not IsoDec 4-bit or IsoDec 8-bit: ", len(overlap_flashdeconv_topfd))
    print("Overlap between all four: ", len(overlap_all))

    pmatched_id4, coverages_id4 = assess_decon_accuracy(id4_prsms, id4_msalign_directory, matchedscans=overlap_all)
    pmatched_id8, coverages_id8 = assess_decon_accuracy(id8_prsms, id8_msalign_directory, matchedscans=overlap_all)
    pmatched_flashdeconv, coverages_flashdeconv = assess_decon_accuracy(flashdeconv_prsms, flashdeconv_msalign_directory, matchedscans=overlap_all)
    pmatched_topfd, coverages_topfd = assess_decon_accuracy(topfd_prsms, topfd_msalign_directory, matchedscans=overlap_all)

    #print all of the values
    print("Average Percent Matched IsoDec 4-bit:", np.mean(pmatched_id4))
    print("Standard Deviation Percent Matched IsoDec 4-bit:", np.std(pmatched_id4))
    print("Average Percent Matched IsoDec 8-bit:", np.mean(pmatched_id8))
    print("Standard Deviation Percent Matched IsoDec 8-bit:", np.std(pmatched_id8))
    print("Average Percent Matched FLASHDeconv:", np.mean(pmatched_flashdeconv))
    print("Standard Deviation Percent Matched FLASHDeconv:", np.std(pmatched_flashdeconv))
    print("Average Percent Matched TopFD:", np.mean(pmatched_topfd))
    print("Standard Deviation Percent Matched TopFD:", np.std(pmatched_topfd))

    print("Average Coverage IsoDec 4-bit:", np.mean(coverages_id4))
    print("Standard Deviation Coverage IsoDec 4-bit:", np.std(coverages_id4))
    print("Average Coverage IsoDec 8-bit:", np.mean(coverages_id8))
    print("Standard Deviation Coverage IsoDec 8-bit:", np.std(coverages_id8))
    print("Average Coverage FLASHDeconv:", np.mean(coverages_flashdeconv))
    print("Standard Deviation Coverage FLASHDeconv:", np.std(coverages_flashdeconv))
    print("Average Coverage TopFD:", np.mean(coverages_topfd))
    print("Standard Deviation Coverage TopFD:", np.std(coverages_topfd))

if True:
    id8_prsmfile = isodec_psms_file = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\msalign_for_toppic\8bin\ID_8\Task2-SearchTask\AllPSMs.psmtsv"
    id8_msalign_directory = "Z:\\Group Share\\JGP\\PXD019247\\isodec_8bit_c_msalign\\"
    id8_prsms, id8_peptidemasses = prepare_prsms(id8_prsmfile, qval_cutoff, allow_ambiguous)

    flashdeconv_prsmfile = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\FLASHDeconv_msalign\FLD\Task2-SearchTask\AllPSMs.psmtsv"
    flashdeconv_msalign_directory = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\"
    flashdeconv_prsms, flashdeconv_peptidemasses = prepare_prsms(flashdeconv_prsmfile, qval_cutoff, allow_ambiguous)

    topfd_prsmfile = r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\topFD_msalign\TFD\Task2-SearchTask\AllPSMs.psmtsv"
    topfd_msalign_directory = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\topFD_msalign\\"
    topfd_prsms, topfd_peptidemasses = prepare_prsms(topfd_prsmfile, qval_cutoff, allow_ambiguous)

    #First get overlap between each pair of files
    print("Comparing IsoDec 8-bit and FLASHDeconv")
    id8_vs_flashdeconv = compare_prsms_return_uniqueIDs(id8_prsms, id8_peptidemasses, flashdeconv_prsms, flashdeconv_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing IsoDec 8-bit and TopFD")
    id8_vs_topfd = compare_prsms_return_uniqueIDs(id8_prsms, id8_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    print("Comparing FLASHDeconv and TopFD")
    flashdeconv_vs_topfd = compare_prsms_return_uniqueIDs(flashdeconv_prsms, flashdeconv_peptidemasses, topfd_prsms, topfd_peptidemasses,
                                    match_mods, match_modloc, match_scan, match_charge, allow_ambiguous, match_filename)

    unique_id8 = [i for i in id8_vs_flashdeconv[2] if i in id8_vs_topfd[2]]
    unique_flashdeconv = [i for i in id8_vs_flashdeconv[3] if i in flashdeconv_vs_topfd[2]]
    unique_topfd = [i for i in id8_vs_topfd[3] if i in flashdeconv_vs_topfd[3]]

    overlap_id8_flashdeconv = [i for i in id8_vs_flashdeconv[0] if i in id8_vs_topfd[2]]
    overlap_id8_topfd = [i for i in id8_vs_topfd[0] if i in id8_vs_flashdeconv[2]]
    overlap_flashdeconv_topfd = [i for i in flashdeconv_vs_topfd[0] if i in id8_vs_flashdeconv[3]]

    overlap_all = [i for i in id8_vs_flashdeconv[0] if i in id8_vs_topfd[0]]


    print("Unique to IsoDec 8-bit: ", len(unique_id8))
    print("Unique to FLASHDeconv: ", len(unique_flashdeconv))
    print("Unique to TopFD: ", len(unique_topfd))
    print("Overlap between IsoDec 8-bit and FLASHDeconv, but not TopFD: ", len(overlap_id8_flashdeconv))
    print("Overlap between IsoDec 8-bit and TopFD, but not FLASHDeconv: ", len(overlap_id8_topfd))
    print("Overlap between FLASHDeconv and TopFD, but not IsoDec 8-bit: ", len(overlap_flashdeconv_topfd))
    print("Overlap between all three: ", len(overlap_all))

    pmatched_id8, coverages_id8 = assess_decon_accuracy(id8_prsms, id8_msalign_directory, matchedscans=overlap_all)
    pmatched_flashdeconv, coverages_flashdeconv = assess_decon_accuracy(flashdeconv_prsms, flashdeconv_msalign_directory, matchedscans=overlap_all)
    pmatched_topfd, coverages_topfd = assess_decon_accuracy(topfd_prsms, topfd_msalign_directory, matchedscans=overlap_all)

    id8_v_flashdeconv_delta_pmatched = []
    id8_v_topfd_delta_pmatched = []

    id8_v_flashdeconv_delta_coverage = []
    id8_v_topfd_delta_coverage = []

    #Iterate through the keys in pmatched_id8 dictionary
    for key in pmatched_id8:
        #Get the value from pmatched_flashdeconv dictionary with the same key
        if key in pmatched_flashdeconv:
            id8_v_flashdeconv_delta_pmatched.append(pmatched_id8[key] - pmatched_flashdeconv[key])
        else:
            print("Key not found in pmatched_flashdeconv:", key)
        if key in pmatched_topfd:
            id8_v_topfd_delta_pmatched.append(pmatched_id8[key] - pmatched_topfd[key])
        else:
            print("Key not found in pmatched_topfd:", key)

        if key in coverages_flashdeconv:
            id8_v_flashdeconv_delta_coverage.append(coverages_id8[key] - coverages_flashdeconv[key])
        if key in coverages_topfd:
            id8_v_topfd_delta_coverage.append(coverages_id8[key] - coverages_topfd[key])

    print("Average Delta Percent Matched IsoDec 8-bit - FLASHDeconv:", np.mean(id8_v_flashdeconv_delta_pmatched))
    print("Standard Deviation Delta Percent Matched IsoDec 8-bit - FLASHDeconv:", np.std(id8_v_flashdeconv_delta_pmatched))
    print("Average Delta Percent Matched IsoDec 8-bit - TopFD:", np.mean(id8_v_topfd_delta_pmatched))
    print("Standard Deviation Delta Percent Matched IsoDec 8-bit - TopFD:", np.std(id8_v_topfd_delta_pmatched))
    print("Average Delta Coverage IsoDec 8-bit - FLASHDeconv:", np.mean(id8_v_flashdeconv_delta_coverage))
    print("Standard Deviation Delta Coverage IsoDec 8-bit - FLASHDeconv:", np.std(id8_v_flashdeconv_delta_coverage))
    print("Average Delta Coverage IsoDec 8-bit - TopFD:", np.mean(id8_v_topfd_delta_coverage))
    print("Standard Deviation Delta Coverage IsoDec 8-bit - TopFD:", np.std(id8_v_topfd_delta_coverage))

    print("Average Percent Matched IsoDec 8-bit:", np.mean(list(pmatched_id8.values())))
    print("Standard Deviation Percent Matched IsoDec 8-bit:", np.std(list(pmatched_id8.values())))
    print("Average Percent Matched FLASHDeconv:", np.mean(list(pmatched_flashdeconv.values())))
    print("Standard Deviation Percent Matched FLASHDeconv:", np.std(list(pmatched_flashdeconv.values())))
    print("Average Percent Matched TopFD:", np.mean(list(pmatched_topfd.values())))
    print("Standard Deviation Percent Matched TopFD:", np.std(list(pmatched_topfd.values())))

    print("Average Coverage IsoDec 8-bit:", np.mean(list(coverages_id8.values())))
    print("Standard Deviation Coverage IsoDec 8-bit:", np.std(list(coverages_id8.values())))
    print("Average Coverage FLASHDeconv:", np.mean(list(coverages_flashdeconv.values())))
    print("Standard Deviation Coverage FLASHDeconv:", np.std(list(coverages_flashdeconv.values())))
    print("Average Coverage TopFD:", np.mean(list(coverages_topfd.values())))
    print("Standard Deviation Coverage TopFD:", np.std(list(coverages_topfd.values())))









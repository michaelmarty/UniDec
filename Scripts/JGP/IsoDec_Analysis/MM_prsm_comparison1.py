from Scripts.JGP.IsoDec_Workflow.toppic_prsm_comparison import PTM, PrSM_Collection, sort_collections
import numpy as np
import re
import unidec.IsoDec.datatools as dt
import itertools
import bisect
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams

mod_pattern = r'\[([^\]]*)\]'

class MMTsvFileHeader:
    def __init__(self, header_line):
        self.header_dict = {"filename": None,
                            "scan_num": None,
                            "retention_time": None,
                            "precursor_charge": None,
                            "precursor_mass": None,
                            "accession": None,
                            "aa_range": None,
                            "db_seq": None,
                            "full_seq": None,
                            "pfm_mass": None,
                            "qvalue": None,
                            "target_decoy": None,
                            "score": None
                            }
        header_split = list(header_line.split("\t"))
        try:
            self.header_dict["filename"] = header_split.index("File Name")
            self.header_dict["scan_num"] = header_split.index("Scan Number")
            self.header_dict["retention_time"] = header_split.index("Scan Retention Time")
            self.header_dict["precursor_charge"] = header_split.index("Precursor Charge")
            self.header_dict["precursor_mass"] = header_split.index("Precursor Mass")
            self.header_dict["accession"] = header_split.index("Accession")
            self.header_dict["aa_range"] = header_split.index("Start and End Residues In Full Sequence")
            self.header_dict["db_seq"] = header_split.index("Base Sequence")
            self.header_dict["full_seq"] = header_split.index("Full Sequence")
            self.header_dict["pfm_mass"] = header_split.index("Monoisotopic Mass")
            self.header_dict["qvalue"] = header_split.index("QValue")
            self.header_dict["target_decoy"] = header_split.index("Decoy")
            self.header_dict["score"] = header_split.index("Score")
        except Exception as e:
            print(e)
            return None

class MM_PrSM:
    def __init__(self, line, current_id, parsed_header):
        try:
            if "|" in line:
                return
            split_line = line.split("\t")
            self.global_id = current_id
            filename = split_line[parsed_header["filename"]].split('.')[0]
            self.filename = filename.replace("_ms2", "")
            self.ms2_scan_num = int(split_line[parsed_header["scan_num"]])
            self.unique_id = self.filename + "_" + str(self.ms2_scan_num)
            self.rt = float(split_line[parsed_header["retention_time"]])
            self.precursor_charge = int(float(split_line[parsed_header["precursor_charge"]]))
            self.precursor_mass = float(split_line[parsed_header["precursor_mass"]])
            self.accession = split_line[parsed_header["accession"]].split('|')[0]
            self.parse_aa_range(split_line[parsed_header["aa_range"]])
            self.db_seq = split_line[parsed_header["db_seq"]]
            self.pfm_mass = float(split_line[parsed_header["pfm_mass"]])
            self.is_ambiguous = False #This immediately gets fixed in full sequence parsing
            self.ptm_set = None
            self.possible_ptm_sets = np.array([])
            self.parse_fullseqs(split_line[parsed_header["full_seq"]])
            self.qvalue = float(split_line[parsed_header["qvalue"]])
            self.target_decoy = split_line[parsed_header["target_decoy"]]
            self.score = float(split_line[parsed_header["score"]])
        except Exception as e:
            print(e)
            self.global_id = None
            return



    def parse_aa_range(self, aa_range):
        if "|" in aa_range:
            aa_range_split = aa_range.split('|')
            aa_range = aa_range_split[0]

        #First remove the brackets
        aa_range = aa_range.replace('[','').replace(']','')
        aa_split = aa_range.split('to')
        self.aa_start = int(aa_split[0])
        self.aa_end = int(aa_split[1])

    def parse_fullseqs(self, fullseq_cell):
        fullseqs = []
        if "|" in fullseq_cell:
            fullseqs = fullseq_cell.split("|")
            self.is_ambiguous = True


        if self.is_ambiguous:
            self.possible_ptm_sets = [np.array([]) for i in range(len(fullseqs))]
            for f, fullseq in enumerate(fullseqs):
                matches = list(re.finditer(mod_pattern, fullseq))
                positions = np.array([])
                added_chars = 0
                for i, match in enumerate(matches):
                    positions = np.append(positions, match.start() - added_chars)
                    added_chars += len(str(match.group(1))) + 2


                ptms = np.array([])
                for i,match in enumerate(matches):
                    ptm = PTM(str(match.group(1)), "metamorpheus", position=positions[i])
                    ptms = np.append(ptms, ptm)

                self.possible_ptm_sets[f] = ptms
            return
        else:
            self.possible_ptm_sets = None
            matches = list(re.finditer(mod_pattern, fullseq_cell))
            positions = np.array([])
            added_chars = 0
            if len(list(matches)) == 0:
                self.ptm_set = np.array([])
            else:
                for i, match in enumerate(matches):
                    positions = np.append(positions, match.start() - added_chars)
                    added_chars += len(str(match.group(1))) + 2


                ptms = np.array([])
                for i, match in enumerate(matches):
                    ptm = PTM(str(match.group(1)), "metamorpheus", position=positions[i])
                    ptms = np.append(ptms, ptm)

                self.ptm_set = ptms
                return

def compare_ptm_sets(ptm_set1, ptm_set2, match_mod_localization):
    if ptm_set1 is None or ptm_set2 is None:
        return False

    if len(ptm_set1) != len(ptm_set2):
        return False

    for i in range(len(ptm_set1)):
        if ptm_set1[i].name != ptm_set2[i].name:
            return False
        if ptm_set1[i].type != ptm_set2[i].type:
            return False
        if match_mod_localization and ptm_set1[i].position != ptm_set2[i].position:
            return False

    return True

def compare_prsms(prsm1, prsm2, match_mod_identity, match_mod_localization):
    if prsm1.filename != prsm2.filename:
        return False

    matched_scan_num = False
    matched_precursor_charge = False
    matched_accessions = False
    matched_aa_range = False
    matched_ptm_sets = False

    if prsm1.ms2_scan_num == prsm2.ms2_scan_num:
        matched_scan_num = True
    if prsm1.precursor_charge == prsm2.precursor_charge:
        matched_precursor_charge = True
    if prsm1.accession == prsm2.accession:
        matched_accessions = True
    if prsm1.aa_start == prsm2.aa_start and prsm1.aa_end == prsm2.aa_end:
        matched_aa_range = True

    if match_mod_identity:
        if prsm1.is_ambiguous and prsm2.is_ambiguous:
            for i in range(len(prsm1.possible_ptm_sets)):
                for j in range(len(prsm2.possible_ptm_sets)):
                    if compare_ptm_sets(prsm1.possible_ptm_sets[i], prsm2.possible_ptm_sets[j], match_mod_localization):
                        matched_ptm_sets = True
                        break
        elif prsm1.is_ambiguous and not prsm2.is_ambiguous:
            for i in range(len(prsm1.possible_ptm_sets)):
                if compare_ptm_sets(prsm1.possible_ptm_sets[i], prsm2.ptm_set, match_mod_localization):
                    matched_ptm_sets = True
                    break
        elif not prsm1.is_ambiguous and prsm2.is_ambiguous:
            for j in range(len(prsm2.possible_ptm_sets)):
                if compare_ptm_sets(prsm1.ptm_set, prsm2.possible_ptm_sets[j], match_mod_localization):
                    matched_ptm_sets = True
                    break
        elif not prsm1.is_ambiguous and not prsm2.is_ambiguous:
            if compare_ptm_sets(prsm1.ptm_set, prsm2.ptm_set, match_mod_localization):
                matched_ptm_sets = True
    else:
        matched_ptm_sets = True

    return matched_scan_num and matched_precursor_charge and matched_accessions and matched_aa_range and matched_ptm_sets

def create_prsm_collections(engines, prsms, match_mod_identity=True, match_mod_localization=True,
                            mod_matchtol=0.1, max_missed_monos=3):
    collections = np.array([])
    collection_masses = np.array([])

    for i, engine in enumerate(engines):
        for prsm in prsms[i]:
            if len(collections) == 0:
                new_collection = PrSM_Collection(prsm, engine)
                collections = np.append(collections, new_collection)
                collection_masses = np.append(collection_masses, prsm.pfm_mass)
            else:
                possible_pfm_masses = [prsm.pfm_mass]
                for k in range(max_missed_monos):
                    if k == 0:
                        continue
                    else:
                        possible_pfm_masses.append(prsm.pfm_mass + (k * 1.0033))
                        possible_pfm_masses.append(prsm.pfm_mass - (k * 1.0033))
                indices_within_tol = []
                for mass in possible_pfm_masses:
                    indices_within_tol.extend(dt.fastwithin_abstol(collection_masses, mass, 0.6))

                nearest = bisect.bisect_left(collection_masses, prsm.pfm_mass)
                if len(indices_within_tol) == 0:
                    # Create a new collection
                    new_collection = PrSM_Collection(prsm, engine)
                    # Insert the mass into the collection_masses array such that the array remains ordered
                    collection_masses = np.insert(collection_masses, nearest, prsm.pfm_mass)
                    collections = np.insert(collections, nearest, new_collection)
                else:
                    added_to_existing = False
                    for index in indices_within_tol:
                        if compare_prsms(prsm, collections[index].prsms[0], match_mod_identity, match_mod_localization):
                            collections[index].prsms.append(prsm)
                            collections[index].engines.append(engine)
                            added_to_existing = True
                            break
                    if not added_to_existing:
                        new_collection = PrSM_Collection(prsm, engine)
                        collection_masses = np.insert(collection_masses, nearest, prsm.pfm_mass)
                        collections = np.insert(collections, nearest, new_collection)
    return collections

def load_prsms(engines, filepaths, allow_ambiguous=False, qvalue_cutoff=0.01):
    prsms = [[] for i in range(len(engines))]

    for i, file in enumerate(filepaths):
        current_id = 0
        parsed_header = {}
        found_header = False
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if found_header:
                    prsm = MM_PrSM(line, current_id, parsed_header.header_dict)
                    if prsm.global_id is None or prsm.qvalue > qvalue_cutoff:
                        continue
                    if allow_ambiguous:
                        prsms[i].append(prsm)
                        current_id += 1
                    if not allow_ambiguous and prsm.is_ambiguous:
                        continue
                    if not allow_ambiguous and not prsm.is_ambiguous:
                        prsms[i].append(prsm)
                        current_id += 1

                else:
                    split = line.split("\t")
                    if split[0] == "File Name":
                        parsed_header = MMTsvFileHeader(line)
                        found_header = True
    return prsms

def load_psm_file(filepath, allow_ambiguous=False, qvalue_cutoff=0.01):
    psms = []

    current_id = 0
    parsed_header = {}
    found_header = False
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            if found_header:
                psm = MM_PrSM(line, current_id, parsed_header.header_dict)
                if psm.global_id is None or psm.qvalue > qvalue_cutoff:
                    continue
                if allow_ambiguous:
                    psms.append(psm)
                    current_id += 1
                if not allow_ambiguous and psm.is_ambiguous:
                    continue
                if not allow_ambiguous and not psm.is_ambiguous:
                    psms.append(psm)
                    current_id += 1

            else:
                split = line.split("\t")
                if split[0] == "File Name":
                    parsed_header = MMTsvFileHeader(line)
                    found_header = True
    return psms


def get_prsm_scan_overlap(engines, prsms):
    print("Getting Scan Overlap...")
    scan_files = [[] for i in range(len(engines))]
    for i, engine in enumerate(engines):
        for prsm in prsms[i]:
            prsm_id = prsm.filename + "_" + str(prsm.ms2_scan_num)
            scan_files[i].append(prsm_id)

    possible_index_combos = []
    for r in range(1, len(engines) + 1):
        possible_index_combos.extend(list(itertools.combinations(range(len(engines)), r)))


    combo_vals = [[] for i in range(len(possible_index_combos))]

    for i, combo in enumerate(possible_index_combos):
        code = np.zeros(len(engines))
        for v in combo:
            code[v] = 1

        included = [[] for v in code if v == 1]
        includes_added = 0
        for j in range(len(code)):
            if code[j] == 1:
                includes_added += 1
                included[includes_added - 1] = scan_files[j]


        excluded_set = set()
        for j in range(len(code)):
            if code[j] == 0:
                excluded_set.update(scan_files[j])

        included_sets = [set(sublist) for sublist in included]
        included_common = reduce(lambda x,y:x.intersection(y), included_sets)
        combo_vals[i] = included_common - excluded_set
        current_engines = [engines[ind] for ind,v in enumerate(code) if int(v) == 1]
        print("Engine(s):", current_engines, "Overlap:", len(combo_vals[i]))

    print("\n")

def group_prsms_by_file(engine, prsms):
    file_prsms = {}
    for prsm in prsms:
        if prsm.filename not in file_prsms:
            file_prsms[prsm.filename] = [prsm]
        else:
            file_prsms[prsm.filename].append(prsm)
    ambiguous_counts = []
    unambiguous_counts = []

    for filename in file_prsms:
        unambiguous = len([p for p in file_prsms[filename] if p.is_ambiguous == False])
        unambiguous_counts.append(unambiguous)
        ambiguous = len([p for p in file_prsms[filename] if p.is_ambiguous == True])
        ambiguous_counts.append(ambiguous)

    print("Engine:", engine)
    print("Average Unambiguous PrSMs per File:", np.mean(unambiguous_counts))
    print("Standard Deviation of Unambiguous PrSMs per File:", np.std(unambiguous_counts))
    print("Average Ambiguous PrSMs per File:", np.mean(ambiguous_counts))
    print("Standard Deviation of Ambiguous PrSMs per File:", np.std(ambiguous_counts))
    print("\n")
    return

def get_decoy_score_dists(engines, prsms):
    FLD_score_cutoff = 8.131
    TFD_score_cutoff = 11.036
    ID8_score_cutoff = 16.013

    targets_above_FLD_cutoff = [0 for i in range(len(engines))]
    targets_above_TFD_cutoff = [0 for i in range(len(engines))]
    targets_above_ID8_cutoff = [0 for i in range(len(engines))]

    score_cutoffs = [ID8_score_cutoff, TFD_score_cutoff, FLD_score_cutoff]
    colors = ['darkblue', 'red', 'orange']

    decoy_scores = [[] for i in range(len(engines))]
    target_scores = [[] for i in range(len(engines))]
    for i, engine in enumerate(engines):
        for prsm in prsms[i]:
            if prsm.target_decoy == "Y":
                decoy_scores[i].append(prsm.score)
            if prsm.target_decoy == "N":
                target_scores[i].append(prsm.score)
                if prsm.score > FLD_score_cutoff:
                    targets_above_FLD_cutoff[i] += 1
                if prsm.score > TFD_score_cutoff:
                    targets_above_TFD_cutoff[i] += 1
                if prsm.score > ID8_score_cutoff:
                    targets_above_ID8_cutoff[i] += 1

    for i in range(len(engines)):
        print("Engine:", engines[i])
        print("Targets Above FLD Cutoff:", targets_above_FLD_cutoff[i])
        print("Targets Above TFD Cutoff:", targets_above_TFD_cutoff[i])
        print("Targets Above ID8 Cutoff:", targets_above_ID8_cutoff[i])
        print("\n")

    min_score = 0
    max_score = 150
    bin_width = 1

    bin_mins = np.arange(min_score, max_score + bin_width, bin_width)
    bin_maxs = bin_mins + bin_width
    bin_centres = (bin_mins + bin_maxs) / 2
    decoy_bin_counts = [np.zeros(len(bin_mins)) for i in range(len(engines))]
    target_bin_counts = [np.zeros(len(bin_mins)) for i in range(len(engines))]

    for i in range(len(engines)):
        for score in decoy_scores[i]:
            bin_index = int(np.floor((score - min_score) / bin_width))
            decoy_bin_counts[i][bin_index] += 1
        for score in target_scores[i]:
            bin_index = int(np.floor((score - min_score) / bin_width))
            target_bin_counts[i][bin_index] += 1




    for i in range(len(engines)):
        plt.plot(bin_centres, decoy_bin_counts[i], label=engines[i], color=colors[i])
        plt.vlines(score_cutoffs[i], 0, 150, linestyles="dashed", color=colors[i])
    plt.xlabel("Score", fontsize=18)
    plt.ylabel("Count", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,150)
    plt.xlim(0,50)
    plt.legend()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    # plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\Decoy_Score_Distributions.png",
    #             dpi=300, transparent=True)
    # plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\Decoy_Score_Distributions.pdf",
    #     dpi=300, transparent=True)
    plt.show()

    for i in range(len(engines)):
        plt.plot(bin_centres, target_bin_counts[i], label=engines[i], color=colors[i])
        plt.vlines(score_cutoffs[i], 0, 150, linestyles="dashed", color=colors[i])
    plt.xlabel("Score", fontsize=18)
    plt.ylabel("Count", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,150)
    plt.legend()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    # plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\Target_Score_Distributions.png",
    #             dpi=300, transparent=True)
    # plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\Target_Score_Distributions.pdf",
    #     dpi=300, transparent=True)
    plt.show()




if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42
    rcParams['font.family'] = 'Arial'

    decon_engines = ["IsoDec 8-bin", "TopFD", "FLASHDeconv"]
    psm_filepaths = [ r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\msalign_for_toppic\8bin\ID_8\Task2-SearchTask\AllPSMs.psmtsv",
                     r"Z:\Group Share\JGP\PXD019247\TD_Search\topFD_msalign\TFD\Task2-SearchTask\AllPSMs.psmtsv",
                     r"C:\Users\MartyLabPC\Desktop\MartyLab\IsoDec\PXD019247\FLASHDeconv_msalign\FLD\Task2-SearchTask\AllPSMs.psmtsv"]


    prsms = load_prsms(decon_engines, psm_filepaths, allow_ambiguous=False, qvalue_cutoff=0.01)
    get_decoy_score_dists(decon_engines, prsms)
    for i in range(len(decon_engines)):
        group_prsms_by_file(decon_engines[i], prsms[i])
    get_prsm_scan_overlap(decon_engines, prsms)

    print("Getting PrSM Overlap...")
    for i in range(len(prsms)):
        print("Engine:", decon_engines[i], "PrSMs:", len(prsms[i]))
    collections = create_prsm_collections(decon_engines, prsms, match_mod_identity=True, match_mod_localization=True)
    sorted_collections = sort_collections(decon_engines, collections)
    for key in sorted_collections.keys():
        print("Engine(s):", key, "PrSMs:", len(sorted_collections[key]))




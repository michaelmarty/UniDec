import fnmatch
import numpy as np
import os
import unidec.IsoDec.datatools as dt
import itertools
import bisect
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib as mpl
from unidec.tools import within_abs


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

class ToppicTsvFileHeader:
    def __init__(self, header_line):
        self.header_dict = {"filename": None,
                       "scan_num": None,
                       "retention_time": None,
                       "precursor_charge": None,
                       "precursor_mass": None,
                       "accession": None,
                       "aa_start": None,
                       "aa_end": None,
                       "db_seq": None,
                       "pfm_seq": None,
                       "pfm_mass": None,
                       "nterm": None,
                       "fixed_mods": None,
                       "unexpected_mods": None,
                       "matchedpeaks_count": None,
                       "matchedfrags_count": None,
                       "evalue": None}

        header_split = list(header_line.split("\t"))
        try:
            self.header_dict["filename"] = header_split.index("Data file name")
            self.header_dict["scan_num"] = header_split.index("Scan(s)")
            self.header_dict["retention_time"] = header_split.index("Retention time")
            self.header_dict["precursor_charge"] = header_split.index("Charge")
            self.header_dict["precursor_mass"] = header_split.index("Precursor mass")
            self.header_dict["accession"] = header_split.index("Protein accession")
            self.header_dict["aa_start"] = header_split.index("First residue")
            self.header_dict["aa_end"] = header_split.index("Last residue")
            self.header_dict["db_seq"] = header_split.index("Database protein sequence")
            self.header_dict["pfm_seq"] = header_split.index("Proteoform")
            self.header_dict["pfm_mass"] = header_split.index("Proteoform mass")
            self.header_dict["nterm"] = header_split.index("Protein N-terminal form")
            self.header_dict["fixed_mods"] = header_split.index("Fixed PTMs")
            self.header_dict["unexpected_mods"] = header_split.index("unexpected modifications")
            self.header_dict["matchedpeaks_count"] = header_split.index("#matched peaks")
            self.header_dict["matchedfrags_count"] = header_split.index("#matched fragment ions")
            self.header_dict["evalue"] = header_split.index("E-value")
        except Exception as e:
            print(e)
            return None

class PTM:
    def __init__(self, split, type, position=None):
        try:
            if type == "fixed":
                newsplit = split.split(':')
                self.name = newsplit[0]
                self.position = int(newsplit[1].replace('[','').replace(']',''))
                self.position_range = None
                self.ambiguous_localization = False
                self.type = type
            elif type == "unexpected":
                newsplit = split.split(':')
                self.name = newsplit[0]
                self.mass = float(newsplit[0])
                if '-' in newsplit[1]:
                    newsplit[1] = newsplit[1].replace('[','').replace(']','')
                    pos_split = newsplit[1].split('-')
                    self.ambiguous_localization = True
                    self.position = None
                    self.position_range = [int(pos_split[0]), int(pos_split[1])]
                else:
                    self.position = int(newsplit[1].replace('[', '').replace(']', ''))
                    self.ambiguous_localization = False
                    self.position_range = None
            elif type == "metamorpheus":
                newsplit = split.split(':')
                self.type = type
                self.mod_type = newsplit[0]
                newsplit2 = newsplit[1].split(" on ")
                self.name = newsplit2[0]
                if position is not None:
                    self.position = position



            else:
                raise ValueError("Invalid modification type", type)
        except ValueError as e:
            print(e)

    def __str__(self):
        if self.ambiguous_localization:
            return self.name + ":" + str(self.position_range[0]) + "-" + str(self.position_range[1])
        else:
            return self.name + ":" + str(self.position)

class PTM_Set:
    def __init__(self, fixed_mods, unexpected_mods):
        fixed_mods_split = fixed_mods.split(';')
        self.fixed_mods = []
        for mod in fixed_mods_split:
            if mod == '':
                continue
            fixed_mod = PTM(mod, "fixed")
            self.fixed_mods.append(fixed_mod)

        unexpected_mods_split = unexpected_mods.split(';')
        self.unexpected_mods = []
        for mod in unexpected_mods_split:
            if mod == '':
                continue
            unexpected_mod = PTM(mod, "unexpected")
            self.unexpected_mods.append(unexpected_mod)


    def __str__(self):
        fixed_mods_str = "Fixed Mods:"
        for mod in self.fixed_mods:
            fixed_mods_str += str(mod) + ";"
        fixed_mods_str = fixed_mods_str[:-1]

        unexpected_mods_str = "Unexpected Mods:"
        for mod in self.unexpected_mods:
            unexpected_mods_str += str(mod) + ";"
        unexpected_mods_str = unexpected_mods_str[:-1]

        return fixed_mods_str + " " + unexpected_mods_str

class Toppic_PrSM:
    def __init__(self, line, current_id, parsed_header):
        try:
            split_line = line.split("\t")
            self.global_id = current_id
            filename = os.path.basename(split_line[parsed_header["filename"]])
            self.filename = filename.split('.')[0].replace('_ms2','')
            self.ms2_scan_num = int(split_line[parsed_header["scan_num"]])
            self.unique_id = self.filename + "_" + str(self.ms2_scan_num)
            self.rt = float(split_line[parsed_header["retention_time"]])/60 #Puts this in minutes
            self.precursor_charge = int(split_line[parsed_header["precursor_charge"]])
            self.precursor_mass = float(split_line[parsed_header["precursor_mass"]])
            self.accession = split_line[parsed_header["accession"]].split('|')[1]
            self.aa_start = int(split_line[parsed_header["aa_start"]])
            self.aa_end = int(split_line[parsed_header["aa_end"]])
            self.db_seq = split_line[parsed_header["db_seq"]]
            self.pfm_seq = split_line[parsed_header["pfm_seq"]]
            self.pfm_mass = float(split_line[parsed_header["pfm_mass"]])
            self.nterm = split_line[parsed_header["nterm"]]
            self.ptm_set = PTM_Set(split_line[parsed_header["fixed_mods"]], split_line[parsed_header["unexpected_mods"]])
            self.matchedpeaks_count = int(split_line[parsed_header["matchedpeaks_count"]])
            self.matchedfrags_count = int(split_line[parsed_header["matchedfrags_count"]])
            self.evalue = float(split_line[parsed_header["evalue"]])
        except Exception as e:
            self.global_id = None
            return
        
class PrSM_Collection:
    def __init__(self, initial_prsm, engine):
        self.prsms = [initial_prsm]
        self.engines = [engine]

def compare_ptm_locs(ptm1, ptm2):
    #Handle the case where both are localized
    if ptm1.ambiguous_localization is False and ptm2.ambiguous_localization is False:
        if ptm1.position == ptm2.position:
            return True
        else:
            return False

    #Handle the case where ptm1 is ambiguously localized
    if ptm1.ambiguous_localization is True and ptm2.ambiguous_localization is False:
        if ptm2.position > ptm1.position_range[0] and ptm2.position < ptm1.position_range[1]:
            return True
        else:
            return False

    #Handle the case where just ptm2 is ambiguously localized
    if ptm1.ambiguous_localization is False and ptm2.ambiguous_localization is True:
        if ptm1.position > ptm2.position_range[0] and ptm1.position < ptm2.position_range[1]:
            return True
        else:
            return False

    #Handle the case where both are ambiuguous
    if ptm1.ambiguous_localization is True and ptm2.ambiguous_localization is True:
        #Check if the ranges overlap at all.
        if ptm1.position_range[0] <= ptm2.position_range[1] and ptm2.position_range[0] <= ptm1.position_range[1]:
            return True
        else:
            return False

def compare_ptm_sets(set1, set2, match_mod_localization, mod_matchtol, max_missed_monos=2):
    if len(set1.fixed_mods) != len(set2.fixed_mods) or len(set1.unexpected_mods) != len(set2.unexpected_mods):
        return False

    #We know that the lengths match, if they're both 0, return true
    if len(set1.fixed_mods) == 0 and len(set1.unexpected_mods) == 0:
        return True

    fixed_mods_matched = False

    #Now we know that both have the same length and are not 0.
    set1_fixedmods_matched = []
    set2_fixedmods_matched = []
    for i, mod1 in enumerate(set1.fixed_mods):
        for j, mod2 in enumerate(set2.fixed_mods):
            if j in set2_fixedmods_matched:
                continue
            matched_mod_identity = False
            matched_mod_localization = False
            if mod1.name == mod2.name:
                matched_mod_identity = True
            if match_mod_localization:
                matched_mod_localization = compare_ptm_locs(mod1, mod2)
            else:
                matched_mod_localization = True

            if matched_mod_identity and matched_mod_localization:
                set1_fixedmods_matched.append(i)
                set2_fixedmods_matched.append(j)

    if (len(set1_fixedmods_matched) == len(set1.fixed_mods) and
            len(set2_fixedmods_matched) == len(set2.fixed_mods)):
        fixed_mods_matched = True

    unexpected_mods_matched = False

    #Now compare the unexpected mods
    set1_unexpectedmods_matched = []
    set2_unexpectedmods_matched = []
    for i, mod1 in enumerate(set1.unexpected_mods):
        for j, mod2 in enumerate(set2.unexpected_mods):
            if j in set2_unexpectedmods_matched:
                continue
            matched_mod_identity = False
            matched_mod_localization = False

            possible_mod1_masses = []
            possible_mod2_masses = []
            for k in range(max_missed_monos):
                if k == 0:
                    possible_mod1_masses.append(mod1.mass)
                    possible_mod2_masses.append(mod2.mass)
                else:
                    possible_mod1_masses.append(mod1.mass + (k * 1.0033))
                    possible_mod2_masses.append(mod2.mass + (k * 1.0033))
                    possible_mod1_masses.append(mod1.mass - (k * 1.0033))
                    possible_mod2_masses.append(mod2.mass - (k * 1.0033))


            if any(within_abs(m1, m2, mod_matchtol) for m1 in possible_mod1_masses for m2 in possible_mod2_masses):
                matched_mod_identity = True
            if match_mod_localization:
                matched_mod_localization = compare_ptm_locs(mod1, mod2)
            else:
                matched_mod_localization = True

            if matched_mod_identity and matched_mod_localization:
                set1_unexpectedmods_matched.append(i)
                set2_unexpectedmods_matched.append(j)


    if (len(set1_unexpectedmods_matched) == len(set1.unexpected_mods) and
        len(set2_unexpectedmods_matched) == len(set2.unexpected_mods)):
        unexpected_mods_matched = True

    return fixed_mods_matched and unexpected_mods_matched

def compare_prsms(prsm1, prsm2, match_mod_identity, match_mod_localization, mod_matchtol):
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
        if (compare_ptm_sets(prsm1.ptm_set, prsm2.ptm_set, match_mod_localization, mod_matchtol) and
            prsm1.nterm == prsm2.nterm):
            matched_ptm_sets = True
        else:
            matched_ptm_sets = False
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
                    #Create a new collection
                    new_collection = PrSM_Collection(prsm, engine)
                    #Insert the mass into the collection_masses array such that the array remains ordered
                    collection_masses = np.insert(collection_masses, nearest, prsm.pfm_mass)
                    collections = np.insert(collections, nearest, new_collection)
                else:
                    added_to_existing = False
                    for index in indices_within_tol:
                        if compare_prsms(prsm, collections[index].prsms[0], match_mod_identity, match_mod_localization, mod_matchtol):
                            collections[index].prsms.append(prsm)
                            collections[index].engines.append(engine)
                            added_to_existing = True
                            break
                    if not added_to_existing:
                        new_collection = PrSM_Collection(prsm, engine)
                        collection_masses = np.insert(collection_masses, nearest, prsm.pfm_mass)
                        collections = np.insert(collections, nearest, new_collection)
    return collections

def sort_collections(engines, collections):
    engine_combos = []
    for i in range(len(engines)):
        combinations = itertools.combinations(engines, i+1)
        engine_combos.extend(combinations)

    collections_dict = {}
    for combo in engine_combos:
        collections_dict[tuple(sorted(combo))] = []

    for collection in collections:
        for key in collections_dict.keys():
            if key == tuple(sorted(collection.engines)):
                collections_dict[key].append(collection)
                break

    return collections_dict

def load_prsms(directories, allow_unexpected_mods=True, allow_ambiguous_localization=True):
    prsm_files = [[] for i in range(len(decon_engines))]
    prsms = [[] for i in range(len(decon_engines))]

    for i, directory in enumerate(directories):
        prsm_files[i] = match_files(directory, "*_prsm.tsv")
        current_id = 0
        parsed_header = {}
        for file in prsm_files[i]:
            found_header = False
            with open(directory + "\\" + file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if found_header:
                        prsm = Toppic_PrSM(line, current_id, parsed_header.header_dict)
                        if prsm.global_id is None:
                            continue
                        else:
                            if allow_unexpected_mods:
                                if any(p.ambiguous_localization for p in prsm.ptm_set.unexpected_mods) and not allow_ambiguous_localization:
                                    continue
                                else:
                                    prsms[i].append(prsm)
                                    current_id += 1
                            else:
                                if len(prsm.ptm_set.unexpected_mods) == 0:
                                    prsms[i].append(prsm)
                                    current_id += 1
                    else:
                        split = line.split("\t")
                        if split[0] == "Data file name":
                            parsed_header = ToppicTsvFileHeader(line)
                            found_header = True
    return prsms

if __name__ == "__main__":
    decon_engines = ["IsoDec 8-bin", "TopFD", "FLASHDeconv"]

    directories = [r"Z:\Group Share\JGP\PXD019247\Laptop_searches\IsoDec\isodec_8bit_c_msalign",
                   r"Z:\Group Share\JGP\PXD019247\Laptop_searches\IsoDec\topFD_msalign",
                   r"Z:\Group Share\JGP\PXD019247\Laptop_searches\IsoDec\FLASHDeconv_msalign"]

    prsms = load_prsms(directories, allow_unexpected_mods=True, allow_ambiguous_localization=False)
    isodec_prsms = prsms[0]
    isodec_scan_files = []
    for prsm in isodec_prsms:
        isodec_scan_files.append(str(prsm.ms2_scan_num) + "_" + prsm.filename)

    topfd_prsms = prsms[1]
    topfd_scan_files = []
    for prsm in topfd_prsms:
        topfd_scan_files.append(str(prsm.ms2_scan_num) + "_" + prsm.filename)

    fld_prsms = prsms[2]
    fld_scan_files = []
    for prsm in fld_prsms:
        fld_scan_files.append(str(prsm.ms2_scan_num) + "_" + prsm.filename)


    id_only = [i for i in isodec_scan_files if i not in topfd_scan_files and i not in fld_scan_files]
    tfd_only = [i for i in topfd_scan_files if i not in fld_scan_files and i not in isodec_scan_files]
    fld_only = [i for i in fld_scan_files if i not in topfd_scan_files and i not in isodec_scan_files]
    id_and_tfd = [i for i in isodec_scan_files if i in topfd_scan_files and i not in fld_scan_files]
    id_and_fld = [i for i in isodec_scan_files if i in fld_scan_files and i not in topfd_scan_files]
    tfd_and_fld = [i for i in fld_scan_files if i in topfd_scan_files and i not in isodec_scan_files]
    all_three = [i for i in fld_scan_files if i in topfd_scan_files and i in isodec_scan_files]

    print("Scan Overlap...")
    print("ID Only:", len(id_only))
    print("TopFD Only:", len(tfd_only))
    print("FLASHDeconv Only:", len(fld_only))
    print("ID and TopFD:", len(id_and_tfd))
    print("ID and FLASHDeconv:", len(id_and_fld))
    print("TopFD and FLASHDeconv:", len(tfd_and_fld))
    print("All Three:", len(all_three))

    collections = create_prsm_collections(decon_engines, prsms, match_mod_identity=False, match_mod_localization=False)
    print("Total Collections:", len(collections))
    sorted_collections = sort_collections(decon_engines, collections)
    for key in sorted_collections.keys():
        print("Engine(s):", key, "PrSMs:", len(sorted_collections[key]))

    # for coll in sorted_collections[tuple(["IsoDec 4-bin"])]:
    #     prsm = coll.prsms[0]
    #     print("File name:", prsm.filename, "Scan number:", prsm.ms2_scan_num, "Accession:", prsm.accession,
    #           "AA Range:", prsm.aa_start, "-", prsm.aa_end, "PTM Set:", str(prsm.ptm_set))








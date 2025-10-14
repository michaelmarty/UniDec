import os
import fnmatch
import numpy as np
from unidec.IsoDec.match import *
import unidec.IsoDec.datatools as dt
from unidec.UniDecImporter.ImporterFactory import *
from unidec.IsoDec.msalign_export import *
from Scripts.JGP.IsoDec_Workflow.toppic_prsm_comparison import *
from unidec.tools import within_ppm


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def load_ms1_msaligns(directories, engines, make_masses=False, mass_tol=20, rt_tol=3, min_scans=10, max_scans=400, min_pfm_mass=4000):
    msalign_files = {}
    matched_collections = {}
    for engine in engines:
        matched_collections[engine] = {}
    for i, directory in enumerate(directories):
        msalign_files[engines[i]] = match_files(directory, "*ms1.msalign")
        for j, msalign_file in enumerate(msalign_files[engines[i]]):
            msalign_filename = msalign_file.split(".")[0]
            msalign_filename = msalign_filename.replace("_ms1", "")
            print(msalign_filename)
            mc = read_msalign_to_matchedcollection_fast(directory + msalign_file, engines[i], load_precursors=False, mz_type=None, msorder=1)
            if make_masses:
                for peak in mc.peaks:
                    mc.add_pk_to_masses(peak, rt_tol=rt_tol)
                mass_indices_to_remove = []
                for mi, mass in enumerate(mc.masses):
                    if len(mass.scans) < min_scans or len(mass.scans) > max_scans or mass.monoiso < min_pfm_mass:
                        mass_indices_to_remove.append(mi)
                        continue

                    data = list(zip(mass.scan_intensities.keys(), mass.scan_intensities.values()))
                    data = np.array(data)
                    # Sort the data by scan number
                    data = data[data[:, 0].argsort()]
                    peaks = dt.fastpeakdetect(data, window=int(min_scans/2), threshold=0.25)
                    intensities = np.array([float(val) for val in mass.scan_intensities.values()])

                    count = int(0.3 * (len(mass.scans)))

                    max_sn = np.max(intensities) / np.mean(np.sort(intensities)[:count])

                    # plt.plot(mass.scan_intensities.keys(), intensities)
                    # plt.title(max_sn)
                    # plt.show()

                    if len(peaks) == 0 and max_sn < 2:
                        mass_indices_to_remove.append(mi)

                mc.masses = [mass for i, mass in enumerate(mc.masses) if i not in mass_indices_to_remove]
            matched_collections[engines[i]][msalign_filename] = mc

    return matched_collections

def load_ms2_msaligns(directories, engines):
    msalign_files = {}
    matched_collections = {}
    for engine in engines:
        matched_collections[engine] = {}
    for i, directory in enumerate(directories):
        msalign_files[engines[i]] = match_files(directory, "*ms2.msalign")
        for j, msalign_file in enumerate(msalign_files[engines[i]]):
            msalign_filename = msalign_file.split(".")[0]
            msalign_filename = msalign_filename.replace("_ms2", "")
            print("Loading", msalign_filename, "from engine", engines[i])
            mc = read_msalign_to_matchedcollection_fast(directory + msalign_file, engines[i],
                                                        load_precursors=True, mz_type=None, msorder=2)
            matched_collections[engines[i]][msalign_filename] = mc

    return matched_collections

def load_ms1_ms2_scandict(rawfile_directory):
    rawfiles = match_files(rawfile_directory, "*.raw")
    ms1_ms2_scandict = {}


    for rawfile in rawfiles:
        rawfile_name = rawfile.split(".")[0]
        ms1_ms2_scandict[rawfile_name] = {}
        reader = ImporterFactory.create_importer(rawfile_directory + rawfile)
        current_ms1_scan = -1
        for i in range(1, len(reader.scans)+1):
            if reader.get_ms_order(i) == 1:
                current_ms1_scan = i
                ms1_ms2_scandict[rawfile_name][current_ms1_scan] = []
            else:
                ms1_ms2_scandict[rawfile_name][current_ms1_scan].append(i)
        keys_to_remove = []
        for key in ms1_ms2_scandict[rawfile_name]:
            if len(ms1_ms2_scandict[rawfile_name][key]) == 0:
                keys_to_remove.append(key)
        for key in keys_to_remove:
            del ms1_ms2_scandict[rawfile_name][key]
    return ms1_ms2_scandict

def extract_ms1_scan_isolation_info(directory, rawfile):
    reader = ImporterFactory.create_importer(directory + rawfile)
    ms1_scan_nums = np.array([])
    ms2_scan_nums = np.array([])
    mzs = np.array([])
    previous_ms1 = 0
    for i in range(1, len(reader.scans)+1):
        if reader.get_ms_order(i) != 2:
            previous_ms1 = int(i)
            continue
        else:
            mz, width = reader.get_isolation_mz_width(i)
            ms1_scan_nums = np.append(ms1_scan_nums, previous_ms1)
            ms2_scan_nums = np.append(ms2_scan_nums, int(i))
            mzs = np.append(mzs, mz)
    np.savez_compressed("Z:\\Group Share\\JGP\\PXD019247\\" + rawfile.split('.')[0] + "_isolation_info.npz",
                        ms1_scans=ms1_scan_nums, ms2_scans=ms2_scan_nums, mzs=mzs)
    return

def load_isolation_info():
    directory = "Z:\\Group Share\\JGP\\PXD019247\\"

    isolation_info = {}

    npz_files = match_files(directory, "*_isolation_info.npz")
    for npz_file in npz_files:
        npz = np.load(directory + npz_file)
        ms1_scans = npz["ms1_scans"]
        ms2_scans = npz["ms2_scans"]
        mzs = npz["mzs"]

        filename = os.path.splitext(os.path.basename(npz_file))[0]
        filename = filename.replace("_isolation_info", "")
        isolation_info[filename] = {}
        for i in range(len(ms1_scans)):
            if ms1_scans[i] not in isolation_info[filename]:
                isolation_info[filename][ms1_scans[i]] = {ms2_scans[i]: mzs[i]}
            else:
                isolation_info[filename][ms1_scans[i]][ms2_scans[i]] = mzs[i]
    return isolation_info

def get_precursors(engines, ms1_mcs, isolation_info):
    precursors = {}
    for engine in engines:
        precursors[engine] = {}
    for filename in isolation_info:
        for engine in engines:
            precursors[engine][filename] = {}

    for engine in engines:
        for filename in isolation_info:
            if filename not in ms1_mcs[engine]:
                continue
            ms1_scan_dict = sort_by_scan_order(ms1_mcs[engine][filename],1)
            for key1 in isolation_info[filename]:
                for ms2_scan in isolation_info[filename][key1]:
                    precursor_min = isolation_info[filename][key1][ms2_scan] - 1
                    precursor_max = isolation_info[filename][key1][ms2_scan] + 1
                    curr_precursors = findprecursors(precursor_min, precursor_max, key1, ms1_scan_dict, max_precursors=1)
                    if len(curr_precursors) == 0:
                        continue
                    else:
                        newprec  = Precursor(filename, engine, key1, ms2_scan, curr_precursors[0].monoiso, curr_precursors[0].z)
                        precursors[engine][filename][ms2_scan] = newprec

    return precursors

def get_precursors_from_mcs(engines, ms2_mcs, isolation_info):
    precursors = {}
    for engine in engines:
        precursors[engine] = {}
    for filename in isolation_info:
        for engine in engines:
            precursors[engine][filename] = {}

    for filename in isolation_info:
        for engine in engines:
            if filename not in ms2_mcs[engine]:
                continue
            for precursor in ms2_mcs[engine][filename].precursors:
                precursors[engine][filename][precursor.ms2_scan] = precursor
    return precursors

def group_precursors(precursors, engines, filenames):
    precursor_collections = {}
    for file in filenames:
        precursor_collections[file] = {}
    for engine in engines:
        for file in filenames:
            if file not in precursors[engine]:
                continue
            for ms2_scan in precursors[engine][file]:
                if ms2_scan not in precursor_collections[file]:
                    collection = PrecursorCollection(file, engine, precursors[engine][file][ms2_scan])
                    precursor_collections[file][ms2_scan] = collection
                else:
                    precursor_collections[file][ms2_scan].engines.append(engine)
                    if engine in precursor_collections[file][ms2_scan].eng_prec_dict:
                        print("Error! Attempting to add two precursors from the same engine!")
                        continue
                    else:
                        precursor_collections[file][ms2_scan].eng_prec_dict[engine] = precursors[engine][file][ms2_scan]


    for file in filenames:
        if file not in precursor_collections:
            continue
        for scan in precursor_collections[file]:
            precursor_collections[file][scan].process_collection()


    return precursor_collections

def sort_collections(engines, collections):
    engine_combos = []
    for i in range(len(engines)):
        combinations = itertools.combinations(engines, i+1)
        engine_combos.extend(combinations)

    collections_dict = {}
    for combo in engine_combos:
        collections_dict[tuple(sorted(combo))] = []

    for file in collections:
        for ms2_scan in collections[file]:
            prec_collection = collections[file][ms2_scan]
            for key in collections_dict.keys():
                if tuple(sorted(key)) == tuple(sorted(prec_collection.engines)):
                    collections_dict[key].append(prec_collection)
                    break

    return collections_dict

def process_all_collections(collections):
    all_three_matched_z = 0
    two_of_three_matched_z = 0
    one_of_three_matched_z = 0
    two_of_two_matched_z = 0
    one_of_two_matched_z = 0
    one_of_one_matched_z = 0

    all_three_matched_masses = 0
    two_of_three_matched_masses = 0
    one_of_three_matched_masses = 0
    two_of_two_matched_masses = 0
    one_of_two_matched_masses = 0
    one_of_one_matched_masses = 0
    for file in collections:
        for ms2_scan in collections[file]:
            prec_collection = collections[file][ms2_scan]

            if prec_collection.matched_zs == 3:
                all_three_matched_z += 1
            elif prec_collection.matched_zs == 2:
                engines = len(prec_collection.engines)
                if engines == 2:
                    two_of_two_matched_z += 1
                if engines == 3:
                    two_of_three_matched_z += 1
            else:
                engines = len(prec_collection.engines)
                if engines == 1:
                    one_of_one_matched_z += 1
                if engines == 2:
                    one_of_two_matched_z += 1
                if engines == 3:
                    one_of_three_matched_z += 1

            if prec_collection.matched_masses == 3:
                all_three_matched_masses += 1
            elif prec_collection.matched_masses == 2:
                engines = len(prec_collection.engines)
                if engines == 2:
                    two_of_two_matched_masses += 1
                if engines == 3:
                    two_of_three_matched_masses += 1
            else:
                engines = len(prec_collection.engines)
                if engines == 1:
                    one_of_one_matched_masses += 1
                if engines == 2:
                    one_of_two_matched_masses += 1
                if engines == 3:
                    one_of_three_matched_masses += 1

    print("All three matched Z's:", all_three_matched_z)
    print("Two of three matched Z's:", two_of_three_matched_z)
    print("One of three matched Z's:", one_of_three_matched_z)
    print("Two of two matched Z's:", two_of_two_matched_z)
    print("One of two matched Z's:", one_of_two_matched_z)
    print("One of one matched Z's:", one_of_one_matched_z)

    print("All three matched masses:", all_three_matched_masses)
    print("Two of three matched masses:", two_of_three_matched_masses)
    print("One of three matched masses:", one_of_three_matched_masses)
    print("Two of two matched masses:", two_of_two_matched_masses)
    print("One of two matched masses:", one_of_two_matched_masses)
    print("One of one matched masses:", one_of_one_matched_masses)

    return

def compare_matchedmasses(mass1, mass2, mass_tol_ppm=20, rt_tol=3, min_z_overlap=0.2):
    matched_rt = False
    pass_z_overlap  = False
    matched_filename = False

    delta_rt = abs(mass1.apexrt - mass2.apexrt)
    if delta_rt <= rt_tol:
        matched_rt = True

    #Get the overlap of mass1.zs and mass2.zs
    matched_zs = set(mass1.zs).intersection(set(mass2.zs))
    if len(matched_zs) >= min_z_overlap * len(mass1.zs):
        pass_z_overlap = True



    return matched_rt and pass_z_overlap

def group_masses(engines, filenames, ms1_mcs, max_missed_monos=3, mass_tol_ppm=20, rt_tol=3, min_z_overlap=0.5):
    mass_collections = np.array([])
    mass_collection_masses = np.array([])

    for engine in engines:
        for filename in filenames:
            if filename not in ms1_mcs[engine]:
                continue
            for mass in ms1_mcs[engine][filename].masses:
                if len(mass_collection_masses) == 0:
                    new_mass_collection = MassCollection(filename, engine, mass)
                    mass_collections = np.append(mass_collections, new_mass_collection)
                    mass_collection_masses = np.append(mass_collection_masses, mass.monoiso)
                else:
                    possible_feature_masses = [mass.monoiso]
                    for k in range(max_missed_monos):
                        if k == 0:
                            continue
                        else:
                            possible_feature_masses.append(mass.monoiso + (k * 1.0033))
                            possible_feature_masses.append(mass.monoiso - (k * 1.0033))

                    indices_within_masstol = []
                    mass_tol = (mass_tol_ppm * mass.monoiso) / 1000000
                    for monoiso in possible_feature_masses:
                        indices_within_masstol.extend(dt.fastwithin_abstol(mass_collection_masses, monoiso, mass_tol))

                    nearest = bisect.bisect_left(mass_collection_masses, mass.monoiso)
                    added_to_existing = False

                    for i in indices_within_masstol:
                        if (filename == mass_collections[i].filename and
                                compare_matchedmasses(mass, mass_collections[i].masses[0], mass_tol_ppm, rt_tol, min_z_overlap)):
                            mass_collections[i].masses = np.append(mass_collections[i].masses,  mass)
                            mass_collections[i].engines = np.append(mass_collections[i].engines, engine)
                            mass_collections[i].eng_mass_dict[engine] = mass
                            added_to_existing = True
                            break

                    if not added_to_existing:
                        new_mass_collection = MassCollection(filename, engine, mass)
                        mass_collections = np.insert(mass_collections, nearest, new_mass_collection)
                        mass_collection_masses = np.insert(mass_collection_masses, nearest, mass.monoiso)
    return mass_collections

def get_scan_collections(engines, filenames, reported_ms2s):
    scan_collections = {}
    for filename in filenames:
        scan_collections[filename] = {}

    for engine in engines:
        for filename in filenames:
            for ms2_scan in reported_ms2s[engine][filename]:
                if scan_collections[filename].get(ms2_scan) is None:
                    scan_collections[filename][ms2_scan] = ScanCollection(filename, engine, ms2_scan)
                else:
                    scan_collections[filename][ms2_scan].engines.append(engine)
    return scan_collections

def group_scan_collections(engines, collections):
    engine_combos = []
    for i in range(len(engines)):
        combinations = itertools.combinations(engines, i + 1)
        engine_combos.extend(combinations)

    collections_dict = {}
    for combo in engine_combos:
        collections_dict[tuple(sorted(combo))] = []

    for file in collections:
        for ms2_scan in collections[file]:
            engines_id = tuple(sorted(collections[file][ms2_scan].engines))
            collections_dict[engines_id].append(collections[file][ms2_scan])
    return collections_dict


class Precursor:
    def __init__(self, filename, engine, ms1_scan, ms2_scan, precursor_mass, precursor_charge, precursor_intensity=None, precursor_mz=None):
        self.filename = filename
        self.engine = engine
        self.ms1_scan = ms1_scan
        self.ms2_scan = ms2_scan
        self.precursor_mass = precursor_mass
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        self.precursor_intensity = precursor_intensity

class PrecursorCollection:
    '''
    Stores a set of precursors (separated by engine).
    Importantly, only take s
    '''
    def __init__(self, filename, engine, precursor):
        self.filename = filename
        self.engines = [engine]
        self.eng_prec_dict = {engine: precursor}
        self.matched_zs = 0
        self.matched_masses = 0


    def process_collection(self, mass_tol_ppm=20, max_missed_monos=3):
        #First compare the z's:
        z_matches_dict = {}
        for engine in self.eng_prec_dict:
            if self.eng_prec_dict[engine].precursor_charge not in z_matches_dict:
                z_matches_dict[self.eng_prec_dict[engine].precursor_charge] = 1
            else:
                z_matches_dict[self.eng_prec_dict[engine].precursor_charge] += 1

        mass_matches_dict = {}
        for engine in self.eng_prec_dict:
            if len(mass_matches_dict) == 0:
                mass_matches_dict[self.eng_prec_dict[engine].precursor_mass] = 1
            else:
                possible_masses = [self.eng_prec_dict[engine].precursor_mass]
                for i in range(1,max_missed_monos+1):
                    possible_masses.append(self.eng_prec_dict[engine].precursor_mass + (i*1.0033))
                    possible_masses.append(self.eng_prec_dict[engine].precursor_mass - (i * 1.0033))
                for key in mass_matches_dict.keys():
                    if any(within_ppm(m1, key, mass_tol_ppm) for m1 in possible_masses):
                        mass_matches_dict[key] += 1
                        break
        self.matched_zs = max(z_matches_dict.values())
        self.matched_masses = max(mass_matches_dict.values())

class ScanCollection:
    def __init__(self, filename, engine, scan):
        self.filename = filename
        self.scan = scan
        self.engines = [engine]


class MassCollection:
    def __init__(self, filename, engine, mass):
        self.filename = filename
        self.engines = np.array([engine])
        self.eng_mass_dict = {engine: mass}
        self.masses = np.array([mass])



if __name__ == "__main__":
    engines = ["IsoDec 8-bin", "TopFD", "FLASHDeconv"]

    directories = ["C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\msalign_for_toppic\\8bin\\",
                   "Z:\\Group Share\\JGP\PXD019247\\topFD_msalign\\",
                   "Z:\\Group Share\\JGP\\PXD019247\\FLASHDeconv_msalign\\"]


    if False:
        rawfiledirectory = "Z:\\Group Share\\JGP\\PXD019247\\"
        rawfiles = match_files(rawfiledirectory, "*.raw")
        for rawfile in rawfiles:
            extract_ms1_scan_isolation_info(rawfiledirectory, rawfile)

    isolation_info = load_isolation_info()
    print("Loading MS1s...")
    ms1_matched_collections = load_ms1_msaligns(directories, engines, make_masses=False)
    print("Loading MS2s...")
    ms2_matched_collections = load_ms2_msaligns(directories, engines)
    reported_ms2s = {}
    for engine in engines:
        reported_ms2s[engine] = {}
        for filename in isolation_info:
            reported_ms2s[engine][filename] = set()
            for pk in ms2_matched_collections[engine][filename].peaks:
                reported_ms2s[engine][filename].add(pk.scan)

    #now print the total number of unique scans present in each engine
    for engine in engines:
        total_scans = 0
        for filename in isolation_info:
            total_scans += len(reported_ms2s[engine][filename])
        print(engine, "total scans:", total_scans)
    scan_collections = get_scan_collections(engines, isolation_info.keys(), reported_ms2s)
    scan_collections_dict = sort_collections(engines, scan_collections)
    for key in scan_collections_dict:
        print("Engine(s):", key, "total scans:", len(scan_collections_dict[key]))




    if False:
        mass_collections = group_masses(engines, isolation_info.keys(), ms1_matched_collections)
        print("Total number of mass collections:", len(mass_collections))
        for coll in  mass_collections:
            for e in coll.eng_mass_dict:
                m = coll.eng_mass_dict[e]
                m.scans_intensities = {k: v for k, v in sorted(m.scan_intensities.items(), key=lambda item: item[0])}
                plt.plot(m.scans_intensities.keys(), m.scans_intensities.values(), label=e)
            plt.legend()
            plt.xlabel("Scan number")
            plt.ylabel("Intensity")
            plt.title(coll.masses[0].monoiso)
            plt.show()

        exit()

    #Because of differences in the output, this has become necessary, and I'm livid...
    #I really wanted to pull precursors from each output in the same way
    id8_precursors = get_precursors(["IsoDec 8-bin"], ms1_matched_collections, isolation_info)
    tfd_precursors = get_precursors(["TopFD"], ms1_matched_collections, isolation_info)
    fld_precursors = get_precursors_from_mcs(["FLASHDeconv"], ms2_matched_collections, isolation_info)
    precursors = {}

    precursors["IsoDec 8-bin"] = id8_precursors["IsoDec 8-bin"]
    precursors["TopFD"] = tfd_precursors["TopFD"]
    precursors["FLASHDeconv"] = fld_precursors["FLASHDeconv"]

    groups = group_precursors(precursors, engines, isolation_info.keys())
    process_all_collections(groups)
    total_sum = 0
    for group in groups:
        total_sum += len(groups[group])
    print("Total precursors:", total_sum)
    sorted_collections = sort_collections(engines, groups)
    for key in sorted_collections.keys():
        print("Engine(s):", key, "Total Scans:", len(sorted_collections[key]))

    #What actually needs to happen?
    '''
    The goal is to determine:
    1) In how many cases do all 3 tools report a precursor across all 6 files.
    2) In how many cases do all 3 tools report a precursor, and agree on precursor charge across all 6 files.
    3) In how many cases do all 3 tools report a precursor, agree on the charge, and agree on the mass across all 6 files.
    '''
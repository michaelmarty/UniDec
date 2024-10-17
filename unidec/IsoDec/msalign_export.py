import time
import os


def sort_by_scan_order(matched_collection, order):
    scan_dict = {}
    timestart = time.perf_counter()
    for p in matched_collection.peaks:
        if p.ms_order == order:
            if p.scan not in scan_dict.keys():
                scan_dict[p.scan] = [p]
            else:
                scan_dict[p.scan].append(p)
    print("Time to sort by scan", time.perf_counter() - timestart)
    return scan_dict


def get_ms1_scan_num_id(ms1_scan_dict, k, v):
    sorted_ms1_scan_nums = list(dict(sorted(ms1_scan_dict.items())).keys())
    for i in range(len(sorted_ms1_scan_nums)):
        if i == len(sorted_ms1_scan_nums) - 1:
            return i, sorted_ms1_scan_nums[i]
        elif sorted_ms1_scan_nums[i] > k:
            return i - 1, sorted_ms1_scan_nums[i - 1]


def findprecursors(precursor_min, precursor_max, precursor_scanNum, ms1_scan_dict, max_precursors=None):
    result = []
    for k, v in ms1_scan_dict.items():
        if k == precursor_scanNum:
            for p in v:
                if p.mz >= precursor_min and p.mz <= precursor_max:
                    result.append(p)
    if max_precursors is not None:
        result = sorted(result, key=lambda x: x.matchedintensity, reverse=True)[:max_precursors]
    return result

def findprecursors_noms1(precursor_min, precursor_max, pks, max_precursors=None):
    result = []
    for p in pks:
        if p.mz >= precursor_min and p.mz <= precursor_max:
            result.append(p)
    if max_precursors is not None:
        result = sorted(result, key=lambda x: x.matchedintensity, reverse=True)[:max_precursors]
    return result


def write_ms1_msalign(ms1_scan_dict, ms2_scan_dict, file):
    name, _ = os.path.splitext(os.path.basename(file))
    file = file.replace('\\', '/')
    file = file.replace(".raw", ".mzML")
    with open(name + "_ms1.msalign", "w") as f:
        id = 0

        f.write("#IsoDec Pirates\n")
        f.write("#Timestamp: " + str(time.time()) + "\n")
        f.write("####################### Parameters ######################\n")
        f.write("#File name:\t" + file + "\n")
        f.write("#Faims data:\tNo\n")
        f.write("#Faims voltage:\tN/A\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Number of MS1 scans:\t" + str(-1) + "\n")
        else:
            f.write("#Number of MS1 scans:\t" + str(len(ms1_scan_dict)) + "\n")
        f.write("#Number of MS2 scans:\t" + str(len(ms2_scan_dict)) + "\n")
        f.write("#Spectral data type:\tcentroid\n")
        f.write("#Peak error tolerance:\t0.01\n")
        f.write("#MS1 signal/noise ratio:\t3\n")
        f.write("#MS/MS signal/noise ratio:\t1\n")
        f.write("#Thread number:\t1\n")
        f.write("#Default precursor window:\t3 m/z\n")
        f.write("#Activation type:\tFILE\n")
        f.write("#Use MS-Deconv score:\tNo\n")
        f.write("#Miss MS1 spectra:\tYes\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Max scan number:\t" + str(-1) + "\n")
        else:
            f.write("#Min scan number:\t" + str(list(ms1_scan_dict.keys())[0]) + "\n")
        f.write("#Use single scan noise level:\tNo\n")
        f.write("#ECScore cutoff:\t0.5\n")
        f.write("#Additional feature search:\tNo\n")
        f.write("#Generate Html files:\tNo\n")
        f.write("#Do final filtering:\tYes\n")
        f.write("#Version:\t1.0.0\n")
        f.write("####################### Parameters ######################\n\n")

        for k, v in ms1_scan_dict.items():
            f.write("BEGIN IONS\n")
            f.write("FILE NAME=" + file + "\n")
            f.write("SPECTRUM ID=" + str(id) + "\n")
            f.write("TITLE=Scan_" + str(k) + "\n")
            f.write("RETENTION_TIME=" + str(v[0].rt * 60) + "\n")
            f.write("SCANS=" + str(k) + "\n")
            f.write("LEVEL=" + str(1) + "\n")
            for p in v:
                for monoiso in p.monoisos:
                    f.write(str(monoiso) + "\t" + str(p.mz) + "\t" + str(p.z) + "\t" + str(1) + "\n")
            f.write("END IONS\n\n")
            id += 1


def write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, file, act_type, max_precursors=None):
    starttime = time.perf_counter()
    file = file.replace('\\', '/')
    name, _ = os.path.splitext(os.path.basename(file))
    with open(name + "_ms2.msalign", "w") as f:
        id = 0
        feature_id = 0

        f.write("#IsoDec Pirates\n")
        f.write("#Timestamp: " + str(time.time()) + "\n")
        f.write("####################### Parameters ######################\n")
        f.write("#File name:\t" + file + "\n")
        f.write("#Faims data:\tNo\n")
        f.write("#Faims voltage:\tN/A\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Number of MS1 scans:\t" + str(-1) + "\n")
        else:
            f.write("#Number of MS1 scans:\t" + str(len(ms1_scan_dict)) + "\n")
        f.write("#Number of MS2 scans:\t" + str(len(ms2_scan_dict)) + "\n")
        f.write("#Spectral data type:\tcentroid\n")
        f.write("#Peak error tolerance:\t0.01\n")
        f.write("#MS1 signal/noise ratio:\t3\n")
        f.write("#MS/MS signal/noise ratio:\t1\n")
        f.write("#Thread number:\t1\n")
        f.write("#Default precursor window:\t3 m/z\n")
        f.write("#Activation type:\tFILE\n")
        f.write("#Use MS-Deconv score:\tNo\n")
        f.write("#Miss MS1 spectra:\tYes\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Max scan number:\t" + str(-1) + "\n")
        else:
            f.write("#Min scan number:\t" + str(list(ms1_scan_dict.keys())[0]) + "\n")
        f.write("#Use single scan noise level:\tNo\n")
        f.write("#ECScore cutoff:\t0.5\n")
        f.write("#Additional feature search:\tNo\n")
        f.write("#Generate Html files:\tNo\n")
        f.write("#Do final filtering:\tYes\n")
        f.write("#Version:\t1.0.0\n")
        f.write("####################### Parameters ######################\n\n")


        for k, v in ms2_scan_dict.items():
            mz, width = reader.get_isolation_mz_width(k)
            precursor_min = mz - width / 2
            precursor_max = mz + width / 2

            if (len(ms1_scan_dict) == 0):
                ms1_id, ms1_scan = -1, -1
                precursors = findprecursors_noms1(precursor_min, precursor_max, v, max_precursors)
            else:
                ms1_id, ms1_scan = get_ms1_scan_num_id(ms1_scan_dict, k, v)
                precursors = findprecursors(precursor_min, precursor_max, ms1_scan, ms1_scan_dict, max_precursors)

            if precursors == []:
                continue
            else:
                f.write("BEGIN IONS\n")
                f.write("FILE_NAME=" + file.replace('\\', '/') + "\n")
                f.write("SPECTRUM_ID=" + str(id) + "\n")
                f.write("TITLE=Scan_" + str(k) + "\n")
                f.write("SCANS=" + str(k) + "\n")
                f.write("RETENTION_TIME=" + str(v[0].rt * 60) + "\n")
                f.write("LEVEL=" + str(2) + "\n")
                f.write("MS_ONE_ID=" + str(ms1_id) + "\n")
                f.write("MS_ONE_SCAN=" + str(ms1_scan) + "\n")
                f.write("PRECURSOR_WINDOW_BEGIN=" + str(precursor_min) + "\n")
                f.write("PRECURSOR_WINDOW_END=" + str(precursor_max) + "\n")
                f.write("ACTIVATION=" + act_type + "\n")
                f.write("PRECURSOR_MZ=")
                for i in range(len(precursors)):
                    if i == len(precursors) - 1:
                        f.write(str(precursors[i].mz) + "\n")
                    else:
                        f.write(str(precursors[i].mz) + ":")
                f.write("PRECURSOR_CHARGE=")
                for i in range(len(precursors)):
                    if i == len(precursors) - 1:
                        f.write(str(precursors[i].z) + "\n")
                    else:
                        f.write(str(precursors[i].z) + ":")
                f.write("PRECURSOR_MASS=")
                for i in range(len(precursors)):
                    if i == len(precursors) - 1:
                        f.write(str(precursors[i].monoiso) + "\n")
                    else:
                        f.write(str(precursors[i].monoiso) + ":")
                f.write("PRECURSOR_INTENSITY=")
                for i in range(len(precursors)):
                    if i == len(precursors) - 1:
                        f.write(str(precursors[i].matchedintensity) + "\n")
                    else:
                        f.write(str(precursors[i].matchedintensity) + ":")
                f.write("PRECURSOR_FEATURE_ID=")
                for i in range(len(precursors)):
                    if i == len(precursors) - 1:
                        f.write(str(feature_id) + "\n")
                        feature_id += 1
                    else:
                        f.write(str(feature_id) + ":")
                        feature_id += 1
                for p in v:
                    f.write(str(p.monoiso) + "\t" + str(p.matchedintensity) + "\t" + str(p.z) + "\t" + str(1) + "\n")
                f.write("END IONS\n\n")
                id += 1


def write_msalign_files(eng, reader, file):
    ms1_scan_dict = sort_by_scan_order(eng.pks, 1)
    ms2_scan_dict = sort_by_scan_order(eng.pks, 2)

    write_ms1_msalign(ms1_scan_dict, file)
    write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, file, act_type="HCD")


def export_msalign(pks, reader, filename, act_type="HCD"):
    ms1_scan_dict = msalign.sort_by_scan_order(pks, 1)

    ms2_scan_dict = msalign.sort_by_scan_order(pks, 2)

    msalign.write_ms1_msalign(ms1_scan_dict, filename)
    msalign.write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, filename, act_type=act_type)

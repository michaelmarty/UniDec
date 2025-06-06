import time
import os
import unidec.IsoDec.msalign_export as msalign


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


def write_ms1_msalign(ms1_scan_dict, ms2_scan_dict, file, config):
    with open(file + "_ms1.msalign", "w") as f:
        id = 0

        f.write("#IsoDec\n")
        f.write("#Timestamp: " + str(time.time()) + "\n")
        f.write("####################### Parameters ######################\n")
        f.write("#File name:\t" + file + "\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Number of MS1 scans:\t" + str(-1) + "\n")
        else:
            f.write("#Number of MS1 scans:\t" + str(len(ms1_scan_dict)) + "\n")
        f.write("#Number of MS2 scans:\t" + str(len(ms2_scan_dict)) + "\n")
        f.write("#Adduct mass:\t"+str(config.adductmass)+"\n")
        f.write("#Peak Detection Window (peak count):\t"+str(config.peakwindow)+"\n")
        f.write("#Peak Detection Threshold:\t" + str(config.peakthresh) + "\n")
        f.write("#Phase Resolution:\t"+str(config.phaseres)+"\n")
        f.write("#Peak Match Tolerance (ppm):\t"+str(config.matchtol)+"\n")
        f.write("#Minimum peaks:\t"+str(config.minpeaks)+"\n")
        f.write("#Cosine Similarity Score Threhsold:\t"+str(config.css_thresh)+"\n")
        f.write("#Minimum Experimental Peak Area Explained:\t"+str(config.minareacovered)+"\n")
        f.write("#Maximum Monoisotopic Shift:\t"+str(config.maxshift)+"\n")
        f.write("#Peak Cluster Window (m/z):\t"+str(config.mzwindow[0]) + "-" + str(config.mzwindow[1])+"\n")
        f.write("#Activation type:\tFILE\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Max scan number:\t" + str(-1) + "\n")
        else:
            f.write("#Min scan number:\t" + str(list(ms1_scan_dict.keys())[0]) + "\n")
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
                    f.write(str(monoiso) + "\t" + str(p.matchedintensity) + "\t" + str(p.z) + "\t" + str(1) + "\n")
            f.write("END IONS\n\n")
            id += 1


def write_ms2_msalign(ms2_scan_dict, ms1_scan_dict, reader, file, config, act_type="HCD",
                      max_precursors=None, report_multiple_monoisos=True):
    with open(file + "_ms2.msalign", "w") as f:
        id = 0
        feature_id = 0

        f.write("#IsoDec\n")
        f.write("#Timestamp: " + str(time.time()) + "\n")
        f.write("####################### Parameters ######################\n")
        f.write("#File name:\t" + file + "\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Number of MS1 scans:\t" + str(-1) + "\n")
        else:
            f.write("#Number of MS1 scans:\t" + str(len(ms1_scan_dict)) + "\n")
        f.write("#Number of MS2 scans:\t" + str(len(ms2_scan_dict)) + "\n")
        f.write("#Adduct mass:\t" + str(config.adductmass) + "\n")
        f.write("#Peak Detection Window (peak count):\t" + str(config.peakwindow) + "\n")
        f.write("#Peak Detection Threshold:\t" + str(config.peakthresh) + "\n")
        f.write("#Phase Resolution:\t" + str(config.phaseres) + "\n")
        f.write("#Peak Match Tolerance (ppm):\t" + str(config.matchtol) + "\n")
        f.write("#Minimum peaks:\t" + str(config.minpeaks) + "\n")
        f.write("#Cosine Similarity Score Threhsold:\t" + str(config.css_thresh) + "\n")
        f.write("#Minimum Experimental Peak Area Explained:\t" + str(config.minareacovered) + "\n")
        f.write("#Maximum Monoisotopic Shift:\t" + str(config.maxshift) + "\n")
        f.write("#Peak Cluster Window (m/z):\t" + str(config.mzwindow[0]) + "-" + str(config.mzwindow[1]) + "\n")
        f.write("#Activation type:\tFILE\n")
        if len(ms1_scan_dict) == 0:
            f.write("#Max scan number:\t" + str(-1) + "\n")
        else:
            f.write("#Min scan number:\t" + str(list(ms1_scan_dict.keys())[0]) + "\n")
        f.write("#Version:\t1.0.0\n")
        f.write("####################### Parameters ######################\n\n")

        for k, v in ms2_scan_dict.items():
            mz, width = 0,0

            if (len(ms1_scan_dict) == 0):
                ms1_id, ms1_scan = -1, -1
            else:
                ms1_id, ms1_scan = get_ms1_scan_num_id(ms1_scan_dict, k, v)

            try:
                mz, width = reader.get_isolation_mz_width(k)
                precursor_min = mz - width / 2
                precursor_max = mz + width / 2
            except Exception as e:
                print("Error getting precursor info: " + str(e))
                precursor_min = 0
                precursor_max = 0

            if (len(ms1_scan_dict) == 0):
                precursors = findprecursors_noms1(precursor_min, precursor_max, v, max_precursors)
            else:
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
                    mono_mz = (precursors[i].monoiso / precursors[i].z) + 1.007276466812
                    if i == len(precursors) - 1:
                        f.write(str(mono_mz) + "\n")
                    else:
                        f.write(str(mono_mz) + ":")

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
                    f.write("PRECURSOR_FEATURE_ID=0" + "\n")


            for p in v:
                if report_multiple_monoisos:
                    for monoiso in p.monoisos:
                        f.write(str(round(monoiso, 4)) + "\t" + str(p.matchedintensity) + "\t" + str(p.z) + "\n")
                else:
                    f.write(str(p.monoiso) + "\t" + str(p.matchedintensity) + "\t" + str(p.z) + "\t" + str(1) + "\n")
            f.write("END IONS\n\n")
            id += 1



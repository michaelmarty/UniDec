import math
import os

import numpy as np
import pandas as pd

from unidec import tools as ud
from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol
from unidec.IsoDec.match import MatchedCollection, MatchedPeak, get_estimated_monoiso, remove_multiple_monoisos
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist


# Will generate list of matchedpeak objects from text file
# Just specify the desired field with charge
# Text file should be in the format: field(monoiso or mz) " " charge
def read_manual_annotations(path=None, delimiter=' ', data=None):
    if path is None:
        path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt"
    mc = MatchedCollection()
    peaks = []
    with open(path, "r") as file:
        for line in file:
            if line.strip():
                row = line.strip().split(delimiter, 1)
                peak = round(float(row[0]), 2)
                currCharge = row[1].strip()
                z = MatchedPeak(mz=peak, z=int(currCharge))
                peakmass = (peak - 1.007276467) * int(currCharge)
                z.monoiso = get_estimated_monoiso(peakmass)
                z.monoisos = [z.monoiso]
                isodist = fast_calc_averagine_isotope_dist(z.monoiso, currCharge)
                if data is not None:
                    # Find nearest peak in data
                    mz = isodist[np.argmax(isodist[:, 1]), 0]
                    startindex = fastnearest(data[:, 0], mz)
                    # endindex = fastnearest(data[:, 0], mz + window[1])

                    dmax = data[startindex, 1]

                else:
                    dmax = 10000

                for i in range(len(isodist)):
                    isodist[i, 1] *= dmax

                z.isodist = isodist

                mc.add_peak(z)
                peaks.append(z)
    return mc


def get_unique_matchedions(coll1, coll2):
    shared = []
    unique1 = []
    unique2 = []

    coll2_matchedindices = []
    for p1 in range(len(coll1.peaks)):
        if not pd.isnull(coll1.peaks[p1].matchedion):
            matched = False
            for p2 in range(len(coll2.peaks)):
                if (not pd.isnull(coll2.peaks[p2].matchedion) and
                        coll1.peaks[p1].matchedion == coll2.peaks[p2].matchedion and
                        coll1.peaks[p1].z == coll2.peaks[p2].z):
                    coll2_matchedindices.append(p2)
                    matched = True

            if matched:
                shared.append(coll1.peaks[p1])
            else:
                unique1.append(coll1.peaks[p1])

    for p2 in range(len(coll2.peaks)):
        if p2 not in coll2_matchedindices and not pd.isnull(coll2.peaks[p2].matchedion):
            unique2.append(coll2.peaks[p2])

    return shared, unique1, unique2


def compare_matched_ions(coll1, coll2, other_alg=None):
    shared1 = []
    unique1 = []

    shared2 = []
    unique2 = []

    for p1 in coll1.peaks:
        matched = False
        for p2 in coll2.peaks:
            if set(p1.matchedions) & set(p2.matchedions):
                if p1 not in shared1:
                    shared1.append(p1)
                matched = True
        if not matched and len(p1.matchedions) > 0:
            unique1.append(p1)

    for p2 in coll2.peaks:
        matched = False
        for p1 in coll1.peaks:
            if set(p1.matchedions) & set(p2.matchedions):
                if p2 not in shared2:
                    shared2.append(p2)
                matched = True
        if not matched and len(p2.matchedions) > 0:
            unique2.append(p2)

    other_name = "Other"
    if other_alg is not None:
        other_name = other_alg

    print("% Shared IsoDec:", len(shared1) / (len(shared1) + len(unique1)) * 100)
    print("% Shared", other_name, ":", len(shared2) / (len(shared2) + len(unique2)) * 100)
    print("% Unique IsoDec:", len(unique1) / (len(shared1) + len(unique1)) * 100)
    print("% Unique", other_name, ":", len(unique2) / (len(shared2) + len(unique2)) * 100)

    print("Number shared IsoDec:", len(shared1))
    print("Number shared", other_name, ":", len(shared2))
    print("Number unique IsoDec:", len(unique1))
    print("Number unique", other_name, ":", len(unique2))

    # Create matched collections of shared and unique peaks
    shared1 = MatchedCollection().add_peaks(shared1)
    unique1 = MatchedCollection().add_peaks(unique1)
    shared2 = MatchedCollection().add_peaks(shared2)
    unique2 = MatchedCollection().add_peaks(unique2)

    return shared1, unique1, shared2, unique2


def compare_matchedmasses(coll1, coll2, ppmtol=20, maxshift=3, rt_tol=2, f_shared_zs=0.5):
    # make a list of the masses within the list of matched masses contained by each collection
    monoisos1 = np.array([m.monoiso for m in coll1.masses])
    monoisos2 = np.array([m.monoiso for m in coll2.masses])

    sum_intensity1 = np.sum(np.array([m.totalintensity for m in coll1.masses]))
    sum_intensity2 = np.sum(np.array([m.totalintensity for m in coll2.masses]))

    # now order the monoisos and apply that ordering to the masses
    sorted_indices1 = np.argsort(monoisos1)
    sorted_indices2 = np.argsort(monoisos2)
    monoisos1 = monoisos1[sorted_indices1]
    monoisos2 = monoisos2[sorted_indices2]
    coll1.masses = [coll1.masses[i] for i in sorted_indices1]
    coll2.masses = [coll2.masses[i] for i in sorted_indices2]

    matched_coll1 = []
    matched_coll2 = []

    unique_coll1 = []
    unique_coll2 = []

    # Initially we'll store the indices of the matches/unique masses
    for i in range(len(monoisos1)):
        foundmatch = False
        within_tol = fastwithin_abstol(monoisos2, monoisos1[i], int(math.ceil(maxshift + 0.5)))
        if len(within_tol) > 0:
            # Order the within_tol by descending matched intensity
            within_tol = sorted(within_tol, key=lambda x: coll2.masses[x].totalintensity, reverse=True)
            for j in within_tol:
                # ensure that the retention times are close enough
                if abs(coll1.masses[i].apexrt - coll2.masses[j].apexrt) < rt_tol:
                    # ensure that the monoisos match within the ppmtol allowing the maxshift
                    possible_monoisos = [monoisos1[i] + x * 1.0033 for x in range(-maxshift, maxshift + 1)]
                    if any([ud.within_ppm(m, monoisos2[j], ppmtol) for m in possible_monoisos]):
                        # ensure that the z values overlap enough
                        z_overlap = np.intersect1d(coll1.masses[i].zs, coll2.masses[j].zs)
                        min_frac = len(z_overlap) / max(len(coll1.masses[i].zs), len(coll2.masses[j].zs))
                        if min_frac >= f_shared_zs:
                            matched_coll1.append(i)
                            foundmatch = True
                            break
        if not foundmatch:
            unique_coll1.append(i)

    for i in range(len(monoisos2)):
        if i not in matched_coll2:
            foundmatch = False
            within_tol = fastwithin_abstol(monoisos1, monoisos2[i], int(math.ceil(maxshift + 0.5)))
            if len(within_tol) > 0:
                for j in within_tol:
                    # ensure that the retention times are close enough
                    if abs(coll1.masses[j].apexrt - coll2.masses[i].apexrt) < rt_tol:
                        # ensure that the monoisos match within the ppmtol allowing the maxshift
                        possible_monoisos = [monoisos2[i] + x * 1.0033 for x in range(-maxshift, maxshift + 1)]
                        if any([ud.within_ppm(m, monoisos1[j], ppmtol) for m in possible_monoisos]):
                            # ensure that the z values overlap enough
                            z_overlap = np.intersect1d(coll1.masses[j].zs, coll2.masses[i].zs)
                            min_frac = len(z_overlap) / max(len(coll1.masses[j].zs), len(coll2.masses[i].zs))
                            if min_frac >= f_shared_zs:
                                matched_coll2.append(i)
                                foundmatch = True
                                break
            if not foundmatch:
                unique_coll2.append(i)

    # Now we'll convert the indices to the actual masses
    matched_coll1 = [coll1.masses[i] for i in matched_coll1]
    matched_coll2 = [coll2.masses[i] for i in matched_coll2]
    unique_coll1 = [coll1.masses[i] for i in unique_coll1]
    unique_coll2 = [coll2.masses[i] for i in unique_coll2]

    matched_intensity1 = np.sum(np.array([m.totalintensity for m in matched_coll1]))
    matched_intensity2 = np.sum(np.array([m.totalintensity for m in matched_coll2]))

    percent_matched1 = (len(matched_coll1) / len(coll1.masses)) * 100
    percent_matched2 = (len(matched_coll2) / len(coll2.masses)) * 100
    percent_intensity1 = (matched_intensity1 / sum_intensity1) * 100
    percent_intensity2 = (matched_intensity2 / sum_intensity2) * 100
    print("% Masses Matched Collection 1:", (len(matched_coll1) / len(coll1.masses)) * 100)
    print("% Masses Matched Collection 2:", (len(matched_coll2) / len(coll2.masses)) * 100)
    print("% Matched Intensity Collection 1:", (matched_intensity1 / sum_intensity1) * 100)
    print("% Matched Intensity Collection 2:", (matched_intensity2 / sum_intensity2) * 100)

    return percent_matched1, percent_matched2, percent_intensity1, percent_intensity2, len(matched_coll1), len(
        matched_coll2), len(monoisos1), len(monoisos2)


def compare_annotated(l1, l2, ppmtol, maxshift):
    unique_annotated = []
    unique_experimental = []
    shared = []

    annotated_mzs = []
    annotated_charges = []
    annotated_indices = []
    for i in range(len(l1.peaks)):
        annotated_mzs.append(l1[i].mz)
        annotated_charges.append(l1[i].z)
        annotated_indices.append(i)
        currentshift = 1
        while currentshift <= maxshift:
            negative_mz = l1[i].mz - (1 / l1[i].z) * currentshift
            annotated_mzs.append(negative_mz)
            annotated_charges.append(l1[i].z)
            annotated_indices.append(i)
            positive_mz = l1[i].mz + (1 / l1[i].z) * currentshift
            annotated_mzs.append(positive_mz)
            annotated_charges.append(l1[i].z)
            annotated_indices.append(i)
            currentshift += 1

    # Sort annotated mzs and apply that sort to the other lists
    annotated_mzs, annotated_charges, annotated_indices = zip(
        *sorted(zip(annotated_mzs, annotated_charges, annotated_indices)))

    matched_annotated_indices = []
    for i in range(len(l2.peaks)):
        foundmatch = False
        curr_mz = l2[i].mz
        curr_tol = (curr_mz / 1e6) * ppmtol
        matched_indices = fastwithin_abstol(np.array(annotated_mzs), curr_mz, curr_tol)
        for index in matched_indices:
            matched_charge = annotated_charges[index]
            if matched_charge == l2[i].z:
                shared.append(l2[i])
                matched_annotated_indices.append(annotated_indices[index])
                foundmatch = True
                break
        if not foundmatch:
            unique_experimental.append(l2[i])

    unique_annotated = MatchedCollection().add_peaks(
        [l1[i] for i in range(len(l1.peaks)) if i not in matched_annotated_indices])
    shared = MatchedCollection().add_peaks(shared)
    unique_experimental = MatchedCollection().add_peaks(unique_experimental)

    return shared, unique_annotated, unique_experimental


def read_msalign_to_matchedcollection(file, data=None, mz_type="monoiso", mass_type="monoiso"):
    """
    Read an msalign file to a MatchedCollection object
    :param file: Path to msalign file
    :return: MatchedCollection object
    """

    mc = MatchedCollection()

    current_scan = 1
    current_rt = 0
    scan_peaks = []
    with open(file, "r") as f:
        for line in f:
            # print(line)
            if line[0] == "#":
                continue
            elif "BEGIN IONS" in line:
                continue
            elif "END IONS" in line:
                grouped_peaks = mc.group_peaks(scan_peaks)
                mc.peaks.extend(grouped_peaks)
                scan_peaks = []
                continue
            elif "SCANS" in line:
                current_scan = int(line.split("=")[1])
                continue
            elif "RETENTION_TIME" in line:
                current_rt = float(line.split("=")[1]) / 60
                continue
            elif "=" not in line and len(line) > 1:
                if mass_type == "monoiso":
                    split = line.split("\t")
                    mz = (float(split[0]) + (1.007276467 * int(split[2]))) / int(split[2])
                    pk = MatchedPeak(int(split[2]), mz)
                    pk.monoisos = [float(split[0])]
                    pk.scan = current_scan
                    pk.monoiso = float(split[0])
                    pk.matchedintensity = float(split[1])
                    pk.rt = current_rt
                    scan_peaks.append(pk)
                elif mass_type == "peak":
                    split = line.split("\t")

                    estimated_monoiso = get_estimated_monoiso(float(split[0]))
                    mz = (estimated_monoiso + (1.007276467 * int(split[2]))) / int(split[2])
                    pk = MatchedPeak(int(split[2]), mz)
                    pk.monoisos = [estimated_monoiso]
                    pk.scan = current_scan
                    pk.monoiso = estimated_monoiso
                    pk.matchedintensity = float(split[1])
                    pk.rt = current_rt
                    scan_peaks.append(pk)

                isodist = fast_calc_averagine_isotope_dist(pk.monoiso, pk.z)
                if data is not None:
                    # Find nearest peak in data
                    mz = isodist[np.argmax(isodist[:, 1]), 0]
                    startindex = fastnearest(data[:, 0], mz)
                    # endindex = fastnearest(data[:, 0], mz + window[1])

                    dmax = data[startindex, 1]

                else:
                    dmax = 10000

                for i in range(len(isodist)):
                    isodist[i, 1] *= dmax
                if mz_type == "monoiso":
                    pk.mz = isodist[np.argmax(isodist[:, 1]), 0]
                pk.isodist = isodist
    print("Loaded", len(mc.peaks), "peaks from msalign file")
    return mc


def read_fd_tsv_to_matchedcollection(file, data=None, mz_type="monoiso"):
    mc = MatchedCollection()

    with open(file, "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            else:
                split = line.split('\t')
                monoiso = float(split[6])
                min_z = int(split[8])
                max_z = int(split[9])
                # snr = float(split[14])
                # if snr < 2:
                #     continue
                for z in range(min_z, max_z + 1):
                    mz = (monoiso + (1.007276467 * z)) / z
                    pk = MatchedPeak(z, mz)
                    pk.monoisos = [monoiso]
                    pk.scan = int(split[2])
                    pk.monoiso = monoiso
                    pk.matchedintensity = float(split[7])
                    pk.rt = float(split[3])

                    isodist = fast_calc_averagine_isotope_dist(pk.monoiso, pk.z)
                    if data is not None:
                        # Find nearest peak in data
                        mz = isodist[np.argmax(isodist[:, 1]), 0]
                        startindex = fastnearest(data[:, 0], mz)
                        # endindex = fastnearest(data[:, 0], mz + window[1])

                        dmax = data[startindex, 1]

                    else:
                        dmax = 10000

                    for i in range(len(isodist)):
                        isodist[i, 1] *= dmax
                    if mz_type == "monoiso":
                        pk.mz = isodist[np.argmax(isodist[:, 1]), 0]
                    pk.isodist = isodist
                    mc.add_peak(pk)

    print("Loaded", len(mc.peaks), "peaks from FLASHDeconv .tsv file")
    return mc


def peak_mz_z_df_to_matchedcollection(df, data=None):
    """
    Convert a pandas dataframe of peak mzs (col1) and zs (col2) to a MatchedCollection object
    :param df: Pandas dataframe with columns mz, z
    :param data: Optional data to match to as 2D numpy array [m/z, intensity]
    :return: MatchedCollection object
    """
    mc = MatchedCollection()
    for i, row in df.iterrows():
        pk = MatchedPeak(z=int(row['z']), mz=row['peakmz'])
        pk.monoiso = (pk.mz - 1.007276467) * pk.z
        pk.monoisos = [pk.monoiso]
        mc.add_peak(pk)
        isodist = fast_calc_averagine_isotope_dist(pk.monoiso, pk.z)

        if data is not None:
            # Find nearest peak in data
            mz = isodist[np.argmax(isodist[:, 1]), 0]
            startindex = fastnearest(data[:, 0], mz)
            # endindex = fastnearest(data[:, 0], mz + window[1])

            dmax = data[startindex, 1]


        else:
            dmax = 10000

        for i in range(len(isodist)):
            isodist[i, 1] *= dmax
        pk.mz = isodist[np.argmax(isodist[:, 1]), 0]
        pk.isodist = isodist
    print("Loaded", len(mc.peaks), "peaks from dataframe")
    return mc


def compare_matchedcollections(coll1, coll2, ppmtol=50, objecttocompare="monoisos", maxshift=3, ignorescan=False):
    unique1 = []
    shared = []
    unique2 = []

    diffs = []

    if not ignorescan:
        scans = np.unique([p.scan for p in coll1.peaks])
        scans2 = np.unique([p.scan for p in coll2.peaks])
        scans = np.union1d(scans, scans2)
    else:
        scans = [0]
    for s in scans:
        if not ignorescan:
            # Pull out peaks for each individual scan.
            coll1_sub = MatchedCollection().add_peaks([p for p in coll1.peaks if p.scan == s])
            coll2_sub = MatchedCollection().add_peaks([p for p in coll2.peaks if p.scan == s])
        else:
            coll1_sub = coll1
            coll2_sub = coll2

        # Sort the peak lists and the masses by monoisotopic mass
        coll1_sub.peaks = sorted(coll1_sub.peaks, key=lambda x: x.monoiso)
        coll2_sub.peaks = sorted(coll2_sub.peaks, key=lambda x: x.monoiso)
        masses1 = np.array([p.monoiso for p in coll1_sub.peaks])
        masses2 = np.array([p.monoiso for p in coll2_sub.peaks])
        # print(masses1, masses2)
        coll2_matched_indices = []

        for i in range(len(coll1_sub.peaks)):
            foundmatch = False

            within_tol = fastwithin_abstol(masses2, masses1[i], int(math.ceil(maxshift + 0.5)))
            if len(within_tol) > 0:
                for j in within_tol:
                    if coll1_sub.peaks[i].z == coll2_sub.peaks[j].z:
                        # now check if a monoiso is within the ppmtol
                        for m1 in coll1_sub.peaks[i].monoisos:
                            for m2 in coll2_sub.peaks[j].monoisos:

                                if ud.within_ppm(m1, m2, ppmtol):

                                    coll2_matched_indices.append(j)
                                    foundmatch = True
                                    diffs.append((m1 - m2) / m1 * 1e6)

                                else:
                                    if maxshift > 0:
                                        diff = np.abs(m1 - m2)
                                        # print(m1, m2, diff)
                                        if diff < maxshift * 1.1:
                                            mm_count = int(round(m1 - m2))
                                            newm2 = m2 + mm_count * 1.0033
                                            if ud.within_ppm(m1, newm2, ppmtol):
                                                coll2_matched_indices.append(j)
                                                foundmatch = True
                                                diffs.append((m1 - newm2) / m1 * 1e6)

            if foundmatch:
                shared.append(coll1_sub.peaks[i])
            else:
                unique1.append(coll1_sub.peaks[i])

        coll2_unique_matchedinds = list(set(coll2_matched_indices))
        for i in range(len(coll2_sub.peaks)):
            if i not in coll2_matched_indices:
                unique2.append(coll2_sub.peaks[i])

    shared = MatchedCollection().add_peaks(shared)
    unique1 = MatchedCollection().add_peaks(unique1)
    unique2 = MatchedCollection().add_peaks(unique2)

    print("Avg PPM Diff:", np.mean(diffs))

    print("Shared:", len(shared.peaks), "Unique 1", len(unique1.peaks), "Unique2:", len(unique2.peaks))

    return shared, unique1, unique2


def read_msalign_to_matchedcollection_fast(file, decon_engine=None, data=None, mz_type="monoiso", msorder=1):
    mc = MatchedCollection()

    msalign_filename = os.path.splitext(os.path.basename(file))[0]
    msalign_filename = msalign_filename.replace("_ms2", "")

    scan_lines = []
    with (open(file, "r") as f):
        adding_to_scan = False
        for line in f:
            if "BEGIN IONS" in line:
                adding_to_scan = True
                current_scan = []
            elif "END IONS" in line:
                scan_lines.append(current_scan)
                adding_to_scan = False
            elif adding_to_scan:
                current_scan.append(line)
            else:
                continue
    for lineset in scan_lines:
        parse_line_set(mc, msalign_filename, decon_engine, lineset, data=data, scan_order=msorder)
    for peak in mc.peaks:
        peak.ms_order = msorder
    print("Loaded", len(mc.peaks), "peaks from msalign file")
    return mc


def parse_line_set(mc, filename, decon_engine, lines, data=None, scan_order=2):
    scan = None
    current_rt = None
    current_precursor_mass = None
    current_precursor_charge = None
    current_precursor_intensity = None
    current_precursor_mz = None
    current_ms1_scan = None
    current_ms2_scan = None
    ionstart = None
    rawpeaks = []
    for i, line in enumerate(lines):
        if "SCANS" in line:
            scan = int(line.split("=")[1])
            current_ms2_scan = int(line.split("=")[1])
        if "RETENTION_TIME" in line:
            current_rt = float(line.split("=")[1]) / 60
        if "PRECURSOR_MASS" in line:
            try:
                mass_split = line.split("=")
                if ":" not in mass_split[1]:
                    current_precursor_mass = float(mass_split[1])
                else:
                    newsplit = mass_split[1].split(":")
                    current_precursor_mass = float(newsplit[0])
            except Exception as e:
                current_precursor_mass = -1
        if "PRECURSOR_CHARGE" in line:
            try:
                zsplit = line.split("=")
                if ":" not in zsplit[1]:
                    current_precursor_charge = int(line.split("=")[1])
                else:
                    newsplit = zsplit[1].split(":")
                    current_precursor_charge = int(newsplit[0])
            except Exception as e:
                current_precursor_charge = -1
        if "PRECURSOR_INTENSITY" in line:
            try:
                intsplit = line.split("=")
                if ":" not in intsplit[1]:
                    current_precursor_intensity = float(line.split("=")[1])
                else:
                    newsplit = intsplit[1].split(":")
                    current_precursor_intensity = float(newsplit[0])
            except Exception as e:
                current_precursor_intensity = 0

        if "PRECURSOR_MZ" in line:
            try:
                mzsplit = line.split("=")
                if ":" not in mzsplit[1]:
                    current_precursor_mz = float(line.split("=")[1])
                else:
                    newsplit = mzsplit[1].split(":")
                    current_precursor_mz = float(newsplit[0])
            except Exception as e:
                current_precursor_mz = -1

        if "MS_ONE_SCAN" in line:
            current_ms1_scan = int(line.split("=")[1])
        if "=" not in line:
            ionstart = i
            break
    if ionstart is None:
        return
    for i in range(ionstart, len(lines)):
        line = lines[i]
        split = line.split("\t")
        mz = (float(split[0]) + (1.007276467 * int(split[2]))) / int(split[2])
        pk = MatchedPeak(int(split[2]), mz)
        pk.monoisos = [float(split[0])]
        pk.scan = scan
        pk.monoiso = float(split[0])
        pk.matchedintensity = float(split[1])
        pk.rt = current_rt
        pk.ms_order = scan_order

        isodist = fast_calc_averagine_isotope_dist(pk.monoiso, pk.z)
        if data is not None:
            # Find nearest peak in data
            mz = isodist[np.argmax(isodist[:, 1]), 0]
            startindex = fastnearest(data[:, 0], mz)
            # endindex = fastnearest(data[:, 0], mz + window[1])

            dmax = data[startindex, 1]

        else:
            dmax = 10000

        for i in range(len(isodist)):
            isodist[i, 1] *= dmax
        pk.isodist = isodist

        pk.mz = isodist[np.argmax(isodist[:, 1]), 0]

        rawpeaks.append(pk)

    mc.peaks.extend(mc.group_peaks(rawpeaks))
    return


def df_to_matchedcollection(df, monoiso="Monoisotopic Mass", peakmz="Most Abundant m/z", peakmass="Most Abundant Mass",
                            scan="Scan", z="Charge", intensity="Abundance", ion="Ion", checkmultiplemonoisos=False):
    """
    Convert a pandas dataframe to a MatchedCollection object
    :param df: Pandas dataframe with columns mz, intensity, z, scan, rt, monoiso, peakmass, avgmass
    :return: MatchedCollection object
    """
    mc = MatchedCollection()
    for i, row in df.iterrows():
        pk = MatchedPeak(row[z], row[peakmz])
        pk.monoiso = row[monoiso]
        pk.monoisos = np.array([pk.monoiso])

        if scan is not None:
            pk.scan = row[scan]
        else:
            pk.scan = -1
        if peakmass is not None:
            pk.peakmass = row[peakmass]
        else:
            pk.peakmass = row[monoiso]
        if ion is not None:
            pk.matchedion = row[ion]
        else:
            pk.matchedion = None
        if pk.matchedion == math.nan:
            pk.matchedion = None
        if intensity is not None:
            pk.matchedintensity = row[intensity]
        else:
            pk.matchedintensity = 1

        isodist = create_isodist2(pk.monoiso, pk.z, pk.matchedintensity)
        pk.isodist = isodist
        mc.add_peak(pk)

    if checkmultiplemonoisos:
        mc = remove_multiple_monoisos(mc)

    print("Loaded", len(mc.peaks), "peaks from dataframe")
    return mc


def create_isodist2(monoiso, charge, maxval, adductmass=1.007276467):
    charge = float(charge)
    isodist = fast_calc_averagine_isotope_dist(monoiso, charge=charge)
    isodist[:, 1] *= maxval
    return isodist

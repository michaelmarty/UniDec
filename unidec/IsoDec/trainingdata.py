import numpy as np
import os
import unidec.tools as ud
from unidec.IsoDec.datatools import *
from copy import deepcopy
import matchms
import pickle as pkl
import time
from unidec.IsoDec.match import create_isodist, find_matches
from typing import List, Tuple
from unidec.UniDecImporter.ImporterFactory import ImporterFactory, recognized_types

@njit(fastmath=True)
def match_peaks(centroids: np.array, isodist: np.array, tol: float = 5.0) -> Tuple[List[int], List[int]]:
    """
    matchingpeaks = matchms.similarity.spectrum_similarity_functions.find_matches(centroids[:, 0],
                                                                                  isodist[:, 0], tol)
    matchedindexes = [match[0] for match in matchingpeaks]
    isomatches = [match[1] for match in matchingpeaks]"""
    matchedindexes, isomatches = find_matches(centroids[:, 0], isodist[:, 0], tol)
    return matchedindexes, isomatches

def match_charge(centroids, peakmz, charge_range=[1, 50], tol=0.01, threshold=0.85, nextbest=0.8, baseline=0.1,
                 silent=True):
    startime = time.perf_counter()
    topcent = deepcopy(centroids)

    b1 = centroids[:, 1] > baseline * np.amax(centroids[:, 1])
    centroids = centroids[b1]
    if len(centroids) < 3:
        return 0, [], [], []

    target = matchms.Spectrum(mz=centroids[:, 0], intensities=centroids[:, 1], metadata={"precursor_mz": peakmz})
    cosine_greedy = matchms.similarity.CosineGreedy(tolerance=tol)

    ztab = np.arange(charge_range[0], charge_range[1])
    charge_guess = simp_charge(centroids, silent=True)
    # sort so that things close to charge_guess come first in ztab
    ztab = np.array(sorted(ztab, key=lambda x: np.abs(x - charge_guess)))

    for i, z in enumerate(ztab):
        isodist = create_isodist(peakmz, z, centroids)

        matchindexes, isomatches = match_peaks(topcent, isodist, tol)
        if len(matchindexes) < 3:
            continue

        spectrum2 = matchms.Spectrum(mz=isodist[:, 0], intensities=isodist[:, 1],
                                     metadata={"charge": z, "precursor_mz": peakmz})

        score = cosine_greedy.pair(target, spectrum2)['score']

        if score > threshold:

            if not silent:
                print("Matched:", z, score)
                endtime = time.perf_counter()
                print("Time:", endtime - startime)
            return z, isodist, matchindexes, isomatches
    if not silent:
        print("No Match")
        endtime = time.perf_counter()
        print("Time:", endtime - startime)
    return 0, [], [], []

def is_valid(file):
    """
    Check file to see if extension is recognized by UniDec
    :param file: File name
    :return: True/False whether file is valid UniDec data file
    """
    extension = "." + file.split(".")[-1]
    if extension.lower() in recognized_types:
        return True
    else:
        return False


def process_file(file, overwrite=False, peakdepth=10, maxpeaks=None, onedropper=0.9):
    """
    Process raw m/z data into centroid pkl files.
    :param file: File name
    :param overwrite: Whether to overwrite existing pkl files
    :param peakdepth: Maximum number of peaks to pick from each spectrum
    :param maxpeaks: Maximum total number of peaks to pick per file
    :param onedropper: Percentage of 1+ charge states to drop
    :return: List of charge states picked for histogram plots
    """
    # Initialize the lists and time
    starttime = time.perf_counter()
    zdat = []
    # Check that the file is valid
    if is_valid(file):
        print(file)
    else:
        print("Unknown File Type", file)
        return []

    # Check whether the output file already exists and if so whether to skip this one
    fileheader = os.path.splitext(os.path.basename(file))[0]
    outfile = fileheader + "_" + str(peakdepth) + ".pkl"
    if os.path.isfile(outfile) and not overwrite:
        print("Pkl found:", outfile)
        return []

    # Get importer and check it
    reader = ImporterFactory.create_importer(file)
    try:
        print("N Scans:", np.amax(reader.scans))
    except Exception as e:
        print("Could not open:", file)
        return []

    # Loop over all scans
    good_centroids = []
    for s in reader.scans:
        # Open the scan and get the spectrum
        try:
            spectrum = reader.grab_scan_data(s)
        except:
            print("Error Reading Scan", s)
            continue
        # If the spectrum is too short, skip it
        if len(spectrum) < 3:
            continue
        # Find the peaks and print progress
        peaks = isotope_finder(spectrum)
        if s % 10 == 0:
            print(s, len(spectrum), len(peaks), len(good_centroids))

        # For each peak, try to find the charge state
        for i, p in enumerate(peaks[:peakdepth]):
            # sort peaks by m/z
            peaks = np.array(sorted(peaks, key=lambda x: x[0]))
            # Get the centroids around the peak
            centroids, d = get_centroids(spectrum, p[0], mzwindow=[-1.5, 3.5])
            # sort centroids by m/z
            if len(centroids) < 3:
                continue
            centroids = np.array(sorted(centroids, key=lambda x: x[0]))
            # Try a simple calculation of the charge state
            testz = simp_charge(centroids, silent=True)
            # If simple charge state calculation is not good, ignore it.
            # Should only be possible if peaks are further apart than they should be
            if testz == 0:
                continue
            # If it is a 1+ charge state, randomly drop some
            elif testz == 1:
                r = np.random.rand()
                if r < onedropper:
                    continue

            # Carefully determine the charge state and get out the matched peaks
            charge, isodist, matched, isomatched = match_charge(centroids, p[0])
            # If it's good, add it to the list
            if charge != 0:
                zdat.append(charge)
                good_centroids.append([centroids, charge, matched, isodist, isomatched, p[0], s])
                # print([centroids, charge, matched, isodist, isomatched])

        # If we have reached the maximum number of peaks, stop
        if maxpeaks is not None:
            if len(good_centroids) > maxpeaks > 0:
                break
    # Dump to pkl file
    with open(outfile, "wb") as f:
        pkl.dump(good_centroids, f)

    print("Time:", time.perf_counter() - starttime, len(zdat))
    return zdat


def process_dir(directory):
    """
    Process directory to turn raw files into centroid pkl files
    :param directory: Target directory
    :return: List of all charges found for histogram plots
    """
    os.chdir(directory)
    files = os.listdir(directory)
    print(files)
    zdatall = []
    for file in files:
        zdatfile = process_file(file)
        for z in zdatfile:
            zdatall.append(z)
    return zdatall


if __name__ == "__main__":
    starttime = time.perf_counter()
    directory = "Z:\Group Share\JGP\PXD027650"
    directory = "Z:\\Group Share\\JGP\\PXD041357"
    directory = "Z:\\Group Share\\JGP\\PXD042298"
    directory = "Z:\\Group Share\\JGP\\PXD042921"
    # os.chdir(directory)
    # file = "20230222_Easy1200_20min_120K_10us_20-80ACN_Iso1.raw"
    # file = "20220709_nLC1000_E_CEW_Isoform-3_MS2_ET20hcD20_Targetted.raw"
    # zdat = process_file(file, overwrite=True)
    process_dir(directory)
    print("Time:", time.perf_counter() - starttime)
    pass

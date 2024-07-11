import numpy as np
import os
import unidec.tools as ud
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import match_charge
import pickle as pkl
import time


def is_valid(file):
    extension = "." + file.split(".")[-1]
    if extension.lower() in ud.known_extensions:
        return True
    else:
        return False


def process_file(file, overwrite=False, peakdepth=100):
    starttime = time.perf_counter()
    zdat = []
    if is_valid(file):
        print(file)
    else:
        print("Unknown File Type", file, extension)
        return []

    fileheader = os.path.basename(file).split(".")[0]
    outfile = fileheader + "_" + str(peakdepth) + ".pkl"

    if os.path.isfile(outfile) and not overwrite:
        print("Pkl found:", outfile)
        return []

    reader = ud.get_importer(file)

    try:
        print("N Scans:", np.amax(reader.scans))
    except Exception as e:
        print("Could not open:", file)
        return []

    good_centroids = []

    for s in reader.scans:
        try:
            spectrum = reader.grab_scan_data(s)
        except:
            print("Error Reading Scan", s)
            continue
        if len(spectrum) < 3:
            continue
        peaks = isotope_finder(spectrum)
        print(s, len(spectrum), len(peaks))

        for i, p in enumerate(peaks[:peakdepth]):
            try:
                centroids, d = get_centroids(spectrum, p[0], mzwindow=[-1.5, 3.5])
                charge, isodist, matched, isomatched = match_charge(centroids, p[0])
                if charge != 0:
                    zdat.append(charge)
                    good_centroids.append([centroids, charge, matched, isodist, isomatched, p[0], s])
                    # print([centroids, charge, matched, isodist, isomatched])
            except Exception as e:
                print("Error:", e)
                print("Error in Scan:", s, p, i)
                pass

    with open(outfile, "wb") as f:
        pkl.dump(good_centroids, f)

    print("Time:", time.perf_counter() - starttime, len(zdat))
    return zdat


def process_dir(directory):
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
    pass

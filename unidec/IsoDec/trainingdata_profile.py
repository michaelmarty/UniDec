import os
import pickle
from itertools import groupby
import matplotlib.pyplot as plt
from unidec.IsoDec.runtime import IsoDecRuntime
import time
import numpy as np

def run_profile_processing(file_path, scan_range=None, mz_range=None, save=False):
    eng = IsoDecRuntime()
    eng.config.css_thresh = 0.95
    eng.config.knockdown_rounds = 1
    eng.config.minareacovered = 0.5
    eng.config.maxshift = 1

    start_time = time.perf_counter()
    reader = eng.process_file(file_path, scans=scan_range, mz_range=mz_range, verbose=False)

    eng.pks.remove_duplicate_peaks()
    print("Number of peaks:", len(eng.pks))
    get_profile_from_pks(reader, eng.pks)
    end_time = time.perf_counter()
    print(f"Time to get profile data: {end_time - start_time:.2f} seconds")

    if save:
        outfile = os.path.splitext(file_path)[0] + "_profile_data.pkl"
        with open(outfile, "wb") as outfile_handle:
            pickle.dump(eng.pks, outfile_handle, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Profile data saved to {outfile}")

    return eng

def get_profile_from_pks(reader, pks):
    """Populate peak profile data while loading each scan only once."""
    peaks_by_scan = sorted(pks, key=lambda peak: peak.scan)
    for scan, scan_peaks in groupby(peaks_by_scan, key=lambda peak: peak.scan):
        scan_data = reader.get_single_scan(scan)
        if scan_data is None:
            print("No scan data for scan", scan)
            continue

        for peak in scan_peaks:
            peak.pull_from_data(scan_data)

def process_dir(directory, check_existing=True):
    total_peaks = 0
    total_files = 0
    for file in os.listdir(directory):
        if file.endswith(".raw"):
            full_path = os.path.join(directory, file)
            # Check existing
            if check_existing:
                outfile = os.path.splitext(full_path)[0] + "_profile_data.pkl"
                if os.path.exists(outfile):
                    print(f"Profile data already exists for {file}, skipping.")
                    continue

            eng = run_profile_processing(full_path, scan_range=None, mz_range=None, save=True)
            total_peaks += len(eng.pks)
            total_files += 1

    print(f"Processed {total_files} files in {directory}, total peaks: {total_peaks}")
    return total_files, total_peaks


def process_dir_list(topdir, dirlist, check_existing=True):
    t1 = time.perf_counter()
    total_peaks = 0
    total_files = 0
    for subdir in dirlist:
        directory = os.path.join(topdir, subdir)
        if os.path.isdir(directory):
            files, peaks = process_dir(directory, check_existing=check_existing)
            total_files += files
            total_peaks += peaks
            print(f"{time.perf_counter() - t1:.2f} seconds")
        else:
            print(f"Directory {directory} does not exist.")

    print(f"Processed {total_files} files in total, total peaks: {total_peaks}")
    t2 = time.perf_counter()
    print(f"Total time: {t2 - t1:.2f} seconds")


def load_training_data(topdir, dirlist):
    training_data = []
    for subdir in dirlist:
        directory = os.path.join(topdir, subdir)
        if os.path.isdir(directory):
            for file in os.listdir(directory):
                if file.endswith("_profile_data.pkl"):
                    full_path = os.path.join(directory, file)
                    print(f"Loading training data from {full_path}")
                    with open(full_path, "rb") as infile:
                        pks = pickle.load(infile)
                        for peak in pks:
                            training_data.append(peak.z)
        else:
            print(f"Directory {directory} does not exist.")
    print(f"Loaded training data from {len(training_data)} peaks.")
    print(f"Max charge state: {max(training_data)}")
    return training_data

def plot_z_dist(training_data, topdir=None):
    z_values = [z for z in training_data]
    plt.hist(z_values, bins=np.arange(min(z_values), max(z_values) + 1) - 0.5, density=True)
    plt.xlabel("Charge State (z)")
    plt.ylabel("Frequency")
    plt.title("Distribution of Charge States in Training Data")
    outfile = os.path.join(topdir, "charge_state_distribution") if topdir else "charge_state_distribution"
    plt.savefig(outfile + ".png", dpi=300, transparent=True)
    plt.savefig(outfile + ".pdf", transparent=True)
    plt.show()


data_dirs = ["MSV000090488",
             "MSV000091923",
             "PXD045560",
             "PXD046651",
             "PXD027650",
             "PXD041357",
             "PXD042921",
             "PXD069439"
             ]



if __name__ == "__main__":
    topdir= r"Z:\Group Share\JGP"

    td = load_training_data(topdir, data_dirs)
    plot_z_dist(td, topdir)
    exit()


    process_dir_list(topdir, data_dirs, check_existing=True)
    exit()

    # directory = r"Z:\Group Share\JGP\PXD069439"
    #
    # process_dir(directory)
    # exit()

    file = r"TDP_LH_F3.raw"


    full_path = os.path.join(directory, file)

    scan_range = [2400, 2450]
    scan_range = None
    mz_range = [628, 631]
    mz_range = None

    eng = run_profile_processing(full_path, scan_range=scan_range, mz_range=mz_range, save=True)

    # # Sort peaks by mz
    # eng.pks.peaks = sorted(eng.pks.peaks, key=lambda peak: peak.mz)
    # # # Filter by charge state 6
    # eng.pks.peaks = [p for p in eng.pks.peaks if p.z == 2 or p.z == 6]
    #
    # print("Number of peaks:", len(eng.pks))

    norm=False
    offset = False
    for i, p in enumerate(eng.pks):
        # Random color
        color = (np.random.rand(), np.random.rand(), np.random.rand())
        # print(p.mz, p.z, p.monoiso)
        if norm:
            p.profile_data[:, 1] /= p.peakint
            p.centroids[:, 1] /= p.peakint

        if offset:
            p.profile_data[:, 1] += i
            p.centroids[:, 1] += i

        plt.plot(p.profile_data[:, 0], p.profile_data[:, 1], color = color)
        plt.plot(p.centroids[:,0], p.centroids[:,1], "o", color = color)
        plt.plot()
    plt.title("Profile Data for Peaks")
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.show()
    exit()

import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy.special import wofz
from sympy.integrals.risch import NonElementaryIntegral

from unidec.IsoDec.datatools import fastpeakdetect, get_fwhm_peak, get_all_centroids
from unidec.IsoDec import datatools as dt
from unidec.UniDecImporter import ImporterFactory
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.match import *


def generate_data(num_samples=1000, intensity_scale=100):
    mz_values = np.random.uniform(0, 4000, num_samples)
    intensity_values = np.zeros_like(mz_values)

    mean_values = [200, 600, 1150, 2000, 3000]
    std_dev_values = [50, 80, 110, 160, 210]

    xmin, xmax = 0, np.max(mz_values)
    x = np.linspace(xmin, xmax, 1000)

    gaussians = []
    for mean, std_dev in zip(mean_values, std_dev_values):
        p = gaussian(x, mean, std_dev)
        alpha = 0.000001
        decay_factor = np.exp(-alpha * (x - mean) ** 2)
        gaussians.append(p * decay_factor)

    for p in gaussians:
        rand_range_val = np.random.uniform(50, 100) * intensity_scale
        p = p / np.max(p) * rand_range_val
        intensity_values += np.interp(mz_values, x, p)

    dat = np.column_stack((mz_values, intensity_values))
    dat = dat[dat[:, 0].argsort()]

    return dat


def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def plot_data(dat, label='Datax', color='b', subplot_num=111, under_gen=False, theo_matches=None, centroided=False,
              cap = None):
    # plt.subplot(subplot_num)
    mz_values = dat[:, 0]
    intensity_values = dat[:, 1]

    if centroided:
        for mz, intensity in zip(mz_values, intensity_values):
            plt.vlines(mz, 0, intensity, color=color, linewidth=1.5, alpha = 0.4)
    else:
        plt.plot(mz_values, intensity_values, color=color, linewidth=1.5)

    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.xlim(0, np.max(mz_values))

    if under_gen:
        if cap is not None:
            plt.ylim(-cap, cap)
        theo_mz_values = theo_matches[:, 0]
        theo_intensity_values = theo_matches[:, 1]
        for mz, intensity in zip(theo_mz_values, theo_intensity_values):
            plt.vlines(mz, 0, -intensity, color='r', linewidth=1.5)
        x_zeros = np.linspace(0, 4000)
        y_zeros = np.zeros_like(x_zeros)
        plt.plot(x_zeros, y_zeros, color='black', linewidth=1.5)

    elif cap is not None:
        plt.ylim(0, cap)
    else:
        max_intensity = np.max(intensity_values)
        plt.ylim(np.min(dat[:,1]), max_intensity)

    plt.xlabel("m/z", fontsize=16)
    plt.ylabel("Intensity", fontsize=14)
    plt.title(f"{label} Peaks", fontsize=16)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.tight_layout(pad=5.0)


def plot_dummy_generated():
    profile = generate_data()

    plt.figure(figsize=(18, 12))

    plot_data(profile, label='Generated Data', subplot_num=221)
    centroided = dt.get_all_centroids(profile)
    plot_data(centroided, label='Centroided Data (Full)', color='g', subplot_num=222)

    indices = set([i for i in range(len(centroided))])
    theo_matches = np.zeros((len(centroided), 2))
    for i in range(len(centroided)):
        idx = random.randint(0, len(centroided) - 1)
        if idx in indices:
            theo_matches[i, 0] = centroided[idx, 0]
            theo_matches[i, 1] = centroided[idx, 1]

            indices.remove(idx)

    plot_data(centroided, label='Matched', color='g', subplot_num=223, under_gen=True,
              theo_matches=theo_matches)

    plt.tight_layout(pad=5.0)

    updated_dat = dt.subtract_matched_centroid_range(profile, theo_matches, centroided)
    # plot_data(centroided, color='g', subplot_num=224, annotate=True)
    plot_data(updated_dat, label='Missed Peaks', color='b', subplot_num=224)
    plt.show()


def centroid_delete_visualize(path, scan=None):

    plt.figure(figsize=(18, 12))
    importer = ImporterFactory.create_importer(path)

    if scan != None:
        print("FETCHING DATA FROM SCAN NUMBER: ", scan)
        profile = importer.get_single_scan(scan)
    else:
        print("FETCHING AVERAGE SCAN")
        profile = importer.get_avg_scan()
    centroided = dt.get_all_centroids(profile)
    plot_data(profile, label='Generated Data', subplot_num=221)
    plot_data(centroided, label='Centroided Data (Full)', color='g', subplot_num=222)
    indices = set([i for i in range(len(centroided))])
    theo_matches = np.zeros((len(centroided), 2))
    for i in range(len(centroided)//2):
        idx = random.randint(0, len(centroided) - 1)
        if idx in indices:
            theo_matches[i, 0] = centroided[idx, 0]
            theo_matches[i, 1] = centroided[idx, 1]
            indices.remove(idx)
    plot_data(centroided, label='Matched', color='g', subplot_num=223, under_gen=True,
              theo_matches=theo_matches)
    plt.tight_layout(pad=10.0)

    updated_dat = dt.subtract_matched_centroid_range(profile, theo_matches, centroided)
    # plot_data(centroided, color='g', subplot_num=224, annotate=True)
    plot_data(updated_dat, label='Missed Peaks', color='b', subplot_num=224)
    plt.show()


def grab_eng_centroids(centroided):
    eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
    eng.config.activescan = 1
    eng.config.activescanrt = -1
    eng.config.activescanorder = 2
    eng.batch_process_spectrum(centroided, centroided=True)
    all_centroids = []
    for mc in eng.pks.peaks:
        for centroid in mc.centroids:
            all_centroids.append([centroid[0], centroid[1]])

    stk = np.array(all_centroids)
    return stk


def apply_gaussian_to_ranges(dat, unresolved_masses, unresolved_maxcharges, unresolved_distwidths, unresolved_heights,
                             all_charges, resolution):
    mz_min = 1000
    x_axis = dat[:, 0][60000:]
    intensities = np.zeros_like(x_axis)
    for i in range(len(unresolved_masses)):
        for charge in unresolved_maxcharges:
            zdist = gaussian(all_charges, charge, unresolved_distwidths[i])
            for j in range(len(zdist)):
                if zdist[j] < 0.01:
                    continue
                mz = unresolved_masses[i] / all_charges[j]
                if mz < mz_min:
                    continue
                sigma = unresolved_masses[i] / resolution
                gauss = gaussian(x_axis, mz, sigma, unresolved_heights[i] * zdist[j])
                intensities += gauss
    #Slice and seperate from roughly 744 mz
    dat[:, 1][60000:] = intensities
    return dat


def etd_edits():
    mpl.use('WxAgg')
    plt.figure(figsize=(18, 12))
    path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\profile_mzml\\1500_scans_200K_16 fills-qb1.mzML"
    importer = ImporterFactory.create_importer(path)
    dat = importer.get_single_scan(1)
    # np.savetxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\1500 scans profile mzml data.txt", dat)
    # return dat
    unresolved_masses = [57650, 63420, 96500, 234570]
    unresolved_maxcharges = [104, 72, 83, 146]
    unresolved_distwidths = [2.5, 2.5, 2.5, 2.25]
    unresolved_heights = [100000, 140000, 150000, 200000]
    resolution = 100000
    all_charges = np.arange(1, 150)


    updated_dat = apply_gaussian_to_ranges(dat, unresolved_masses, unresolved_maxcharges, unresolved_distwidths,
                                           unresolved_heights, all_charges=all_charges, resolution=resolution)
    np.savetxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\Updated 1500 scans mzml merged data.txt", updated_dat)
    return updated_dat



def mzml_scans_comp():
    # dat = etd_edits()

    dat = np.loadtxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\Updated 1500 scans mzml merged data.txt")
    dat = dat[:60000]
    highest = np.max(dat[:, 1])
    # dat = np.loadtxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\Updated 1500 scans mzml merged data.txt")
    # plot_data(dat, label='Control(with profile mixed)', color='g', subplot_num=221, centroided=False)
    centroided = dt.get_all_centroids(dat, background=0, moving_average_smoothing=0)
    print("DIMENSION OF CENTROIDED", centroided.shape)
    # plot_data(centroided, label='Centroided', color='g', subplot_num=222, centroided=True)
    theoretical_matches = grab_eng_centroids(centroided)

    plt.figure(figsize=(18, 12))
    plot_data(centroided, label='Matched', color='g', subplot_num=221, under_gen=True,
              theo_matches=theoretical_matches, centroided=True)
    remainder, matched_count = dt.subtract_matched_centroid_range(dat, theoretical_matches, centroided)
    plt.figure(figsize=(18, 12))
    plot_data(remainder, label='Missed', color='g', subplot_num=222, centroided=False, cap=highest)
    total_count = len(dat)
    accuracy = (matched_count / total_count) * 100
    print(f"Matched count (values zeroed out): {matched_count}")
    print(f"Length of profile (total points): {total_count}")
    print(f"Total points missed: {total_count - matched_count}")
    print(f"Accuracy: {accuracy:.2f}%")
    plt.show()


def txt_scans_comp():
    # path = "Z:\\Group Share\\JGP\\DataForJoe\\Volker\\Protein+cov-binder_#212#919.txt"
    path = "Z:\\Group Share\\JGP\\DataForJoe\\20240718_BHT_1mgmLGroEL_Static_2.RAW"
    importer = ImporterFactory.create_importer(path)
    dat = importer.get_avg_scan()
    # centroided = dt.get_all_centroids(dat, background=0, moving_average_smoothing=0)
    # theoretical_matches = grab_eng_centroids(centroided)
    centroided = fastpeakdetect(dat)
    plt.plot(dat[:, 0], dat[:, 1], color='b', linewidth=1.5)
    plt.show()
    # remainder, matched_count = dt.subtract_matched_centroid_range(dat, theoretical_matches, centroided)
    # plt.figure(figsize=(18, 12))
    # plot_data(remainder, label='Missed', color='g', subplot_num=224, centroided=False)
    # plt.show()
    # total_count = len(dat)
    # accuracy = (matched_count / total_count) * 100
    # print(f"Matched count (values zeroed out): {matched_count}")
    # print(f"Length of profile (total points): {total_count}")
    # print(f"Total points missed: {total_count - matched_count}")
    # print(f"Accuracy: {accuracy:.2f}%")


def multi_gaussian(x, *params):
    n_gaussians = len(params) // 3
    result = np.zeros_like(x)
    for i in range(n_gaussians):
        A = params[i * 3]
        mu = params[i * 3 + 1]
        sigma = params[i * 3 + 2]
        result += gaussian(x, A, mu, sigma)
    return result


def dynamic_window_size(mz_values, intensity_values, peak_position, default_window=10, min_window=2, max_window=20):
    """
    Dynamically calculate the window size around the peak based on intensity profile.
    """
    peak_index = np.argmin(np.abs(mz_values - peak_position))
    left_index = max(0, int(peak_index - default_window))
    right_index = min(len(mz_values), int(peak_index + default_window))

    intensity_region = intensity_values[left_index:right_index]
    peak_width = np.max(intensity_region) - np.min(intensity_region)

    window_size = max(min_window, min(max_window, peak_width))
    return window_size


def multi_gaussian_fit(data, peaks):
    mz = data[:, 0]
    intensity = data[:, 1]

    mz_spacing = 0.005
    mz_width = mz[-1] - mz[0]
    num_mz = int(mz_width / mz_spacing)
    fit_mz = np.linspace(mz[0], mz[-1], num_mz)

    initial_guess = []
    print(len(peaks))
    for p in peaks:
        # plt.vlines(p[0], 0, p[1], color='r', linestyle='--')
        initial_guess.append(p[1])
        initial_guess.append(p[0])
        window_size = dynamic_window_size(mz, intensity, p[1])
        print(window_size)
        sigma_guess = window_size

        initial_guess.append(sigma_guess)


        params, _ = curve_fit(multi_gaussian, xdata=mz, ydata=intensity, p0=initial_guess, maxfev = 50000)


    centroid_mz = []
    centroid_intensities = []
    for i in range(0, len(params), 3):
        centroid_mz.append(params[i + 1])
        centroid_intensities.append(params[i])

    return centroid_mz, centroid_intensities, params, fit_mz


def gen_dat(path):
    importer = ImporterFactory.create_importer(path)
    dat = importer.get_avg_scan()
    return dat


def trim_file():
    dat = gen_dat("Z:\\Group Share\\JGP\\DataForJoe\\20240718_BHT_1mgmLGroEL_Static_2.RAW")
    dat = dat[11126:11300]
    plt.plot(dat[:, 0], dat[:, 1], color='b', linewidth=1.5)
    plt.show()
    np.savetxt("Z:\\Group Share\\JGP\\DataForJoe\\joe_ripped groel 2 scan.txt", dat)


# Plot the data and fit
def plot_fit(data, centroid_mz, centroid_intensities, params, fit_mz):
    mz = data[:, 0]
    intensity = data[:, 1]

    plt.plot(mz, intensity, label="Original Data", color='g', linewidth=1.5)

    plt.plot(fit_mz, multi_gaussian(fit_mz, *params), label="Multi-Gaussian Fit", color='r', linewidth=1.0)

    for i in range(len(centroid_mz)):
        plt.vlines(centroid_mz[i], 0, centroid_intensities[i]>0, color='r', linestyle='--')


    plt.ylim(np.min(params), np.max(params)+1000)
    plt.xlim(np.min(mz), np.max(mz))
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.show()


if __name__ == "__main__":
    # dat = gen_dat("Z:\\Group Share\\JGP\\DataForJoe\\20240718_BHT_1mgmLGroEL_Static_2.RAW")
    mzml_scans_comp()
    exit()
    dat = gen_dat("Z:\\Group Share\\JGP\\DataForJoe\\joe_ripped groel 2 scan.txt")

    cen = get_all_centroids(dat)
    exit()
    peaks = dt.fastpeakdetect(dat, window=3)
    plt.vlines(peaks[:, 0], 0, peaks[:, 1], color='g')
    cmz, cint, params, fit_mz = multi_gaussian_fit(dat, peaks)
    plot_fit(dat, cmz, cint, params, fit_mz)





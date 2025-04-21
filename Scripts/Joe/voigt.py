import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.optimize import curve_fit
from scipy.optimize import minimize

from unidec.tools import calc_FWHM
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.IsoDec.datatools import fastpeakdetect
from unidec.IsoDec import datatools as dt

def voigt(x, mu, sigma, gamma, height=1):
    if sigma == 0 or gamma == 0:
        return 0
    return height * np.real(wofz(((x - mu) + 1j * gamma) / (sigma * np.sqrt(2))))


def voigt_fwhm(sigma, gamma):
    try:
        return 2 * np.sqrt(2 * np.log(2)) * sigma * np.sqrt(1 + (gamma / sigma) ** 2)
    except ZeroDivisionError:
        return 0

def estimate_sigma_gamma(fwhm, initial_guess=None):
    if initial_guess is None:
        initial_guess = [fwhm / 2, fwhm / 2]

    def objective(params):
        sigma, gamma = params
        return (voigt_fwhm(sigma, gamma) - fwhm) ** 2

    result = minimize(objective, initial_guess, bounds=[(0, None), (0, None)])
    sigma, gamma = result.x
    return sigma, gamma


def multi_voigt(x, *params):
    n_voigts = len(params) // 4
    result = np.zeros_like(x)
    for i in range(n_voigts):
        mu = params[i * 4]
        sigma = params[i * 4 + 1]
        gamma = params[i * 4 + 2]
        height = params[i * 4 + 3]
        result += voigt(x, mu, sigma, gamma, height)
    return result


def isolate_peak(data, peak_position, window_size=1):
    peak_idx = np.abs(data[:, 0] - peak_position).argmin()



    left_idx = max(peak_idx - window_size, 0)
    right_idx = min(peak_idx + window_size, len(data) - 1)

    isolated_data = data[left_idx:right_idx]
    return isolated_data


def multi_voigt_fit(data, peaks, window_size=1):
    mz = data[:, 0]
    intensity = data[:, 1]
    initial_guess = []

    for peak in peaks:
        peak_position = peak[0]
        isolated_data = isolate_peak(data, peak_position, window_size)
        fwhm, indices = calc_FWHM(peak[1], isolated_data)
        sigma, gamma = estimate_sigma_gamma(fwhm)
        initial_guess.append(peak_position)
        initial_guess.append(max(sigma,1e-6))
        initial_guess.append(max(gamma, 1e-6))

        initial_guess.append(peak[1])
    print(initial_guess)
    try:
        params, _ = curve_fit(multi_voigt, mz, intensity, p0=initial_guess, maxfev=1000)
    except Exception as e:
        print("Curve fitting failed.")
        print(e)
        return 0, []
    return params





def split_into_groups(data, peaks, window=5):
    groups = []
    for peak in peaks:
        peak_position = peak[0]
        isolated_data = isolate_peak(data, peak_position, window)
        if len(isolated_data) < 4:
            print(f"Skipping peak at {peak_position} due to insufficient data.")
            continue
        groups.append(isolated_data)
    return groups



def find_closest_mz(original_mz, mz_value):
    return np.abs(original_mz - mz_value).argmin()



def fix_overlap(original_intensity, group):
    result = np.zeros_like(original_intensity)
    for i, intensity in enumerate(original_intensity):
        closest_idx = find_closest_mz(group[:, 0], original_intensity[i])
        result[i] = group[closest_idx, 1]
    return result



#Find the values in the groups' overlapped ranges that are CLOSEST to the original intensities from the original data
#Then store the mu from those values as the centroid of the range


if __name__ == '__main__':
    mpl.use('WxAgg')
    path = ("Z:\\Group Share\\JGP\\DataForJoe\\joe_ripped groel 2 scan.txt")
    #path = "Z:\\Group Share\\JGP\\DataForJoe\\dat.txt"
    importer = ImporterFactory.create_importer(path)
    dat = importer.get_single_scan(1)
    # plt.plot(dat[:, 0], dat[:, 1], color='g')
    peaks = fastpeakdetect(dat, window=3)
    # groups = get_peak_groups(dat[:, 0], dat[:, 1], peaks)



    # plt.plot(dat[:, 0], dat[:, 1], color='g')
    # curves = multi_voigt_fit(dat, peaks)
    # fitted_curve = multi_voigt(dat[:, 0], *curves)

    # for i in range(0, len(curves), 4):
    #     peak_position = curves[i]
    #     peak_intensity = curves[i + 3]
    #     plt.vlines(peak_position, 0, peak_intensity, color='r', linestyle='-')
    r = ['r','g','b']
    all_groups = []
    groups = split_into_groups(dat, peaks, window=20)
    ptr = 1
    for idx, group in enumerate(groups):
        idx = idx%3
        peaks = fastpeakdetect(group, window=3)
        curves = multi_voigt_fit(group, peaks)
        try:
            fitted_curve = multi_voigt(group[:, 0], *curves)
            plt.plot(group[:, 0],fitted_curve, color='r')
            all_groups.append((ptr, fitted_curve, curves[0]))
            ptr+=1
            plt.plot(group[:, 0], np.zeros_like(group[:,0]),  color=r[idx])

        except Exception as e:
            print("Failed to fit curve,", curves)

    for i in all_groups:
        print(i)

    plt.show()
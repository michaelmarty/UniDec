import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from plots import cplot
import matplotlib
import math

data_dirs = ["MSV000090488",
             "MSV000091923",
             "PXD045560",
             "PXD046651",
             "PXD027650",
             "PXD041357",
             "PXD042921",
             "PXD069439"
             ]


def make_wavelengths(n_wavelengths: int, lambda_min: float = 0.001, lambda_max: float = 10000) -> np.ndarray:
    """
    Generates a list of wavelengths following the Casanovo formula:
    ꟛ_i = (ꟛ_min / 2π) * (ꟛ_max / ꟛ_min)^(i / (d_sin - 1))
    From: https://colab.research.google.com/github/Alexander-Sol/MS_spectra_encoding/blob/Colab/03_Casanovo.ipynb#scrollTo=787ec355
    """
    if n_wavelengths < 1:
        raise ValueError("n_wavelengths must be at least 1")
    if lambda_min <= 0 or lambda_max <= lambda_min:
        raise ValueError("lambda_min and lambda_max must satisfy 0 < lambda_min < lambda_max")

    # A single sine/cosine pair should use the longest wavelength so that it
    # spans the full supported m/z range without wrapping.
    if n_wavelengths == 1:
        return np.array([lambda_max / (2 * np.pi)])

    wavelengths = np.zeros(n_wavelengths)
    for i in range(n_wavelengths):
        wavelengths[i] = (lambda_min / (2 * np.pi)) * (lambda_max / lambda_min) ** (i / (n_wavelengths - 1))
    return wavelengths


# Define positional encoding function used by Casanovo
def positional_encoding(m_z, d_model=16):
    """
    Encode a single m/z value into a d_model-dimensional vector.

    If d_model is even, then the wavelengths are shared between sine and cosine.
    If d_model is odd, then we will have one more sine wavelength than cosine wavelength.

    d_sin: The first ⌈d_model / 2⌉ dimensions represent the sine encoding
    d_cos: The last d_model - d_sin dimensions represent the cosine encoding
    """
    encoding = np.zeros(d_model)
    d_sin = math.ceil(d_model / 2)
    d_cos = d_model - d_sin
    sin_wavelengths = make_wavelengths(d_sin)
    cos_wavelengths = sin_wavelengths if d_cos == d_sin else make_wavelengths(d_cos)
    for d in range(d_sin):
        wavelength = sin_wavelengths[d]
        encoding[d] = np.sin(m_z / wavelength)  # First half: sine
    for d in range(d_sin, d_model):
        wavelength = cos_wavelengths[d - d_sin]
        encoding[d] = np.cos(m_z / wavelength)  # Second half: cosine
    return encoding

def decode_mz(encoding, lambda_min=0.001, lambda_max=10000):
    """
    Decode a positional encoding back into an m/z value in [0, lambda_max].

    The longest-wavelength sine/cosine pair completes one cycle over the
    supported m/z range.  ``atan2`` recovers its phase without the quadrant
    ambiguity introduced by applying ``arcsin`` and ``arccos`` separately.
    """
    d_model = len(encoding)
    d_sin = math.ceil(d_model / 2)
    d_cos = d_model - d_sin
    if d_cos == 0:
        raise ValueError("encoding must contain at least one sine/cosine pair")

    sin_wavelengths = make_wavelengths(d_sin, lambda_min, lambda_max)
    cos_wavelengths = make_wavelengths(d_cos, lambda_min, lambda_max)

    # Both wavelength sequences end at lambda_max / (2*pi), so their final
    # components form a matched pair even when d_model is odd.
    if not np.isclose(sin_wavelengths[-1], cos_wavelengths[-1]):
        raise ValueError("encoding does not contain a matched sine/cosine pair")

    phase = np.arctan2(encoding[d_sin - 1], encoding[-1])
    if phase < 0:
        phase += 2 * np.pi
    return phase * sin_wavelengths[-1]

def dist_encoding(dist : np.ndarray, npoints=4000, d_model=256):
    # Encode the distribution using positional encoding
    mz_encodings = np.zeros((npoints, d_model))
    if len(dist) > npoints:
        dist = dist[:npoints]
    for i in range(len(dist)):
        mz_encodings[i] = positional_encoding(dist[i, 0], d_model=d_model)

    int_encodings = np.zeros(npoints)
    int_encodings[:len(dist)] = dist[:, 1] / np.amax(dist[:, 1])

    return mz_encodings, int_encodings

def decode_dist(mz_encodings, int_encodings):
    # Decode the distribution using positional encoding
    npoints, d_model = mz_encodings.shape
    dist = np.zeros((npoints, 2))
    for i in range(npoints):
        dist[i, 0] = decode_mz(mz_encodings[i])
        dist[i, 1] = int_encodings[i]

    # Drop all where first dimension is 0 (these are the padded values)
    dist = dist[dist[:, 0] != 0]
    return dist

def calc_sigma(data):
    """
    https://arxiv.org/pdf/1907.07241
    """
    # s = np.trapezoid(data[:,1], data[:,0])
    s = 0
    x1, y1 = data[0]
    for d in data[1:]:
        x2, y2 = d
        s += (y1 + y2) * (x2 - x1) / 2.0
        x1, y1 = x2, y2

    # diffs = np.diff(data[:,0])
    # ys = data[:-1,1]
    # ys2 = data[1:,1]
    # ys = (ys + ys2) / 2.0  # Average the y values
    # sum = np.sum(ys * diffs)
    return s/(np.sqrt(2 * np.pi) * np.amax(data[:, 1]))


def ndis_std(x: np.ndarray, mid: float, sig: float, a: float = 1.0) -> np.array:
    """
    Normal Gaussian function normalized to the max of 1.
    :param x: x values
    :param mid: Mean of Gaussian
    :param sig: Standard Deviation
    :param a: Maximum amplitude (default is 1)
    :return: Gaussian distribution at x values
    """
    # x = np.array(x).astype(float)
    return a * np.exp(-(x - mid) * (x - mid) / (2.0 * sig * sig))


def dist_fit_to_profile(isodist, profile, window = 10):
    x_output = profile[:,0]
    y_output = np.zeros_like(x_output)

    sigmas = []
    starts = []
    ends = []
    ratios = []
    intensities = []

    for i, (mz, intensity) in enumerate(isodist):
        # Isolate local data near the current m/z value
        index = np.argmin(np.abs(x_output - mz))

        if index > 0 and index < len(x_output) - 1:
            local_int = profile[index, 1]
            index_start = max(0, index - window)
            index_end = min(len(profile), index + window)

            local_data = profile[index_start:index_end]

            b1 = local_data[:,1] > 0.25 * local_int
            local_data2 = local_data[b1]

            local_sigma = calc_sigma(local_data2)
            sigmas.append(local_sigma)
            starts.append(index_start)
            ends.append(index_end)
            ratios.append(local_int / intensity)
            intensities.append(intensity)
        else:
            starts.append(-1)
            ends.append(-1)
    # Weighted average of sigmas based on intensity
    global_sigma = np.average(sigmas, weights=intensities)
    avg_ratio = np.average(ratios, weights=intensities)

    print(global_sigma)

    for i, (mz, intensity) in enumerate(isodist):
            index_start = starts[i]
            index_end = ends[i]

            if index_start == -1 or index_end == -1:
                continue

            local_data = profile[index_start:index_end]

            new_gauss = ndis_std(local_data[:,0], mz, global_sigma, a=intensity * avg_ratio)

            y_output[index_start:index_end] += new_gauss

    return np.transpose(np.vstack((x_output, y_output)))



if __name__ == "__main__":
    # mz=500
    # enc = positional_encoding(mz, d_model=16)
    # print("Encoding for m/z =", mz, ":", enc)
    # decoded_mz = decode_mz(enc)
    # print("Decoded m/z from encoding:", decoded_mz)
    #
    #
    # exit()
    # matplotlib.use("WxAgg")
    topdir = r"Z:\Group Share\JGP"
    file_path = r"Z:\Group Share\JGP\PXD069439\TDP_HH_F12_profile_data.pkl"

    # Load file_path data
    with open(file_path, "rb") as f:
        pks = pickle.load(f)

    # plens = []
    # clens = []
    # pdensity = []
    # cdensity = []
    # for p in pks:
    #     plens.append(len(p.profile_data))
    #     clens.append(len(p.centroids))
    #     pdensity.append(np.mean(np.diff(p.profile_data[:, 0])))
    #     cdensity.append(np.mean(np.diff(p.centroids[:, 0])))
    #
    # print("Max profile length:", max(plens))
    # print("Max centroid length:", max(clens))
    #
    # print("Average profile length:", np.mean(plens))
    # print("Average centroid length:", np.mean(clens))
    #
    # print("Average profile density:", np.mean(pdensity))
    # print("Average centroid density:", np.mean(cdensity))

    for peak in pks:
        fit_data = dist_fit_to_profile(peak.isodist, peak.profile_data, window=10)
        plt.figure(figsize=(10, 5))

        # cplot(peak.centroids)
        plt.plot(peak.profile_data[:, 0], peak.profile_data[:, 1], label="Profile Data", color="k")
        plt.plot(fit_data[:, 0], fit_data[:, 1], label="Fitted Data", color="r")
        # cplot(peak.isodist, factor=1)

        plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.tools as ud
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import *
import pickle as pkl
import time
import torch
import matplotlib as mpl
from math import floor
from numba import jit, njit, prange

try:
    mpl.use("WxAgg")
except:
    pass

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


@njit(fastmath=True)
def decimal_to_other(decimal, base=12):
    """
    Convert a decimal number to a dodecimal number or other base.
    :param decimal: float, decimal number
    :param base: int, base to convert to
    :return: 4 element list of each digit in the new base
    """
    db = decimal * base
    d1 = floor(db)
    d2 = floor((db - d1) * base)
    d3 = floor(((db - d1) * base - d2) * base)
    # d4 = floor((((db - d1) * base - d2) * base - d3) * base)

    return [d1 / base, d2 / base, d3 / base]  # , d4 / base]


@njit(fastmath=True)
def get_digit(number, n):
    """
    Get the nth digit of a number in a given base
    :param number: Float number
    :param n: Digit to get
    :return: Digit normalized to 0 to 1
    """
    a = number // 10. ** n % 10.
    return a / 10.


@njit(fastmath=True)
def encode_mz_to_single(mz):
    """
    Encode the mz data point as a single digit normalized to 0 to 1.
    :param mz: mz value
    :return: digit
    """
    if mz > 100000:
        return 1
    else:
        return mz / 100000.


'''
def mz_to_bit(mz):
    """
    Simple function to get digits for each float.
    :param mz: Float of m/z value
    :return: list of 16 digits scaled to between 0 and 1 starting with 1e5 digit.
    """
    intpart = int(mz)
    decpart = mz - intpart
    intpartdigits = [get_digit(intpart, i) for i in range(6)]
    decimaldigits = [get_digit(decpart * 1e10, i) for i in range(10)]
    # combine and reverse digitis
    digits = np.concatenate([intpartdigits[::-1], decimaldigits[::-1]])
    return digits
'''


@njit(fastmath=True)
def mz_to_bit(mz, maxlen=16):
    """
    Convert m/z to a bit array of digits in different bases.

    TODO: Include actual digit values when finalized

    :param mz: m/z value float
    :param maxlen: Max length of the bit array
    :return: List of digits scaled to 0 to 1 in different bases
    """
    if mz == 0:
        return np.zeros(maxlen)

    # Separate integer part and decimal part
    intpart = int(mz)
    decpart = mz - intpart
    floatpart = decpart
    decpart = round(decpart, 4)

    # Get digits for integer part and decimal part in base 10
    intpartdigits = [get_digit(intpart, i) for i in prange(4)]
    decimaldigits = [get_digit(decpart * 1e4, i) for i in prange(4)]

    # combine digits together into a list
    if maxlen > 16:
        l = maxlen
    else:
        l = 16

    digits = np.zeros(l, dtype=float)
    digits[2:6] = intpartdigits[::-1]
    digits[6:10] = decimaldigits[::-1]
    # Sacrifice the 100k digit for the decimal part
    digits[0] = encode_mz_to_single(mz)
    digits[1] = floatpart

    # Get digits for decimal part in other bases
    dd12 = decimal_to_other(decpart, base=12)
    digits[10:13] = dd12

    if maxlen >= 16:
        dd7 = decimal_to_other(decpart, base=7)
        digits[13:16] = dd7
    if maxlen >= 19:
        dd11 = decimal_to_other(decpart, base=11)
        digits[16:19] = dd11
    else:
        return digits
    if maxlen >= 22:
        dd13 = decimal_to_other(decpart, base=13)
        digits[19:22] = dd13
    if maxlen >= 25:
        dd15 = decimal_to_other(decpart, base=15)
        digits[22:25] = dd15
    if maxlen >= 28:
        dd17 = decimal_to_other(decpart, base=17)
        digits[25:28] = dd17
    if maxlen >= 31:
        dd19 = decimal_to_other(decpart, base=19)
        digits[28:31] = dd19
    # if maxlen >= 34:
    #    dd23 = decimal_to_other(decpart, base=23)
    #    digits[31:34] = dd23

    return digits


@njit(fastmath=True)
def bit_matrix(x, maxlen=16):
    """
    Create a matrix of bit arrays for a list of m/z values
    :param x: list of m/z values
    :param maxlen: Max lengths of the bit arrays
    :return: Matrix of bit arrays (maxlen x maxlen)
    """
    z = np.empty((maxlen, maxlen))
    for i in prange(maxlen):
        z[i] = mz_to_bit(x[i], maxlen=maxlen)
    return z


@njit(fastmath=True)
def diff_matrix(x):
    """
    Create a matrix of differences between m/z values. Faster than numpy approach of subtract.outer
    Set all above 1 to 0 and all below 0 to 1/abs(diff)/50.
    This is a simple way to encode the differences and use the full range of the matrix.
    :param x: m/z values
    :return: matrix of differences (maxlen x maxlen)
    """
    lx = len(x)
    z = np.empty((lx, lx))
    for i in prange(lx):
        for j in prange(lx):
            val = (x[i] - x[j])
            if val < 0:
                val = (1 / abs(val)) / 50.
            if val > 1:
                val = 0
            z[i, j] = val
    return z


@njit(fastmath=True)
def encode_isodist(isodist, maxlen=16):
    """
    Encode an isotope distribution into a 3 channel matrix.
    :param isodist: Isotope distribution (m/z, intensity)
    :param maxlen: Max length of the matrix
    :return: Encoded matrix (3 x maxlen x maxlen), indexes of the original isodist after sorting
    """
    indexes = np.arange(len(isodist))
    # sort isodist by intensity and remove lower values
    sortindex = np.argsort(isodist[:, 1])[::-1]
    isodist = isodist[sortindex]
    indexes = indexes[sortindex]

    # Fill to correct length, pad with 0s or -1s
    isodist2 = np.zeros((maxlen, 2))
    indexes2 = np.zeros(maxlen) - 1
    i1 = isodist[:maxlen]
    isodist2[:len(i1)] = i1
    indexes2[:len(i1)] = indexes[:maxlen]
    isodist = isodist2

    # Set up distaince matrix
    # distmatrix = np.subtract.outer(isodist[:, 0], isodist[:, 0])
    distmatrix = diff_matrix(isodist[:, 0])

    # Create the weights matrix
    weightsmatrix = np.outer(isodist[:, 1], isodist[:, 1])
    weightsmatrix /= np.amax(weightsmatrix)

    # Create the digit matrix
    digitmatrix = bit_matrix(isodist[:, 0], maxlen=maxlen)

    # Put it all together
    emat = np.empty((3, maxlen, maxlen), dtype=float)
    emat[0] = digitmatrix
    emat[1] = distmatrix
    emat[2] = weightsmatrix

    return emat, indexes2


@njit(fastmath=True)
def encode_phase(centroids, maxz=50, phaseres=8):
    """
    Encode the charge phases for a set of centroids

        Old Code:
        rescale = centroids[:, 0] * 2 * np.pi * (i + 1) / mass_diff_c
        y = np.sin(rescale)
        z = np.cos(rescale)
        phase = (np.arctan2(z, y) / (2 * np.pi)) % 1

    :param centroids: Centroids (m/z, intensity)
    :param maxz: Maximum charge state to calculate
    :param phaseres: Resolution of phases to encode in number of bins
    :return: Charge phase histogram (maxz x phaseres)
    """
    phases = np.zeros((maxz, phaseres))
    rescale = centroids[:, 0] / mass_diff_c
    for i in range(maxz):
        # phase = (((rescale * (i + 1)) % 1) + 0.25) % 1
        phase = (rescale * (i + 1)) % 1  # Note, this is a much simpler implementation, needs different model
        phaseindexes = np.floor(phase * phaseres)
        for j in range(len(centroids)):
            phases[i, int(phaseindexes[j])] += centroids[j, 1]
    phases /= np.amax(phases)
    return phases


@njit(fastmath=True)
def encode_phase_all(centroids, peaks, lowmz=-1.5, highmz=5.5, phaseres=8):
    """
    Work on speeding this up
    :param centroids:
    :param peaks:
    :param lowmz:
    :param highmz:
    :return:
    """
    goodpeaks = []
    outcentroids = []
    indexes = []
    indexvalues = np.arange(len(centroids))

    for i, p in enumerate(peaks):
        peakmz = p[0]
        # Find all centroids in the neighborhood of the peak
        start = fastnearest(centroids[:, 0], peakmz + lowmz)
        end = fastnearest(centroids[:, 0], peakmz + highmz) + 1

        if i < len(peaks)-1:
            nextpeak = peaks[i+1][0]
            nextindex = fastnearest(centroids[:, 0], nextpeak)
            if end >= nextindex:
                end = nextindex - 1

        if end - start < 3:
            continue
        c = centroids[start:end]
        goodpeaks.append(p)
        outcentroids.append(c)
        indexes.append(indexvalues[start:end])

    emats = [encode_phase(c, phaseres=phaseres) for c in outcentroids]
    return emats, goodpeaks, outcentroids, indexes


@njit(fastmath=True)
def charge_phase_calculator(centroids, maxz=50, phaseres=16, remove_harmonics=True):
    """
    Calculate the charge phases for a set of centroids
    :param centroids: Centroids (m/z, intensity)
    :param maxz: Maximum charge state to calculate
    :return: Charge phases (maxz x maxlen)
    """
    phases = np.zeros((maxz, phaseres))
    indexes = np.zeros((maxz, len(centroids)))
    for i in range(1, maxz + 1):
        rescale = centroids[:, 0] * 2 * np.pi * i / mass_diff_c
        y = np.sin(rescale)
        z = np.cos(rescale)
        phase = (np.arctan2(z, y) / (2 * np.pi)) % 1
        phaseindexes = np.floor(phase * phaseres)
        indexes[i - 1] = phaseindexes
        for j in range(len(centroids)):
            phases[i - 1, int(phaseindexes[j])] += centroids[j, 1]

    maxes = np.array([np.max(p) for p in phases])

    if remove_harmonics:
        peakbool = maxes > np.amax(maxes) * .95
        if len(peakbool) > 1:
            # Find first peak that is true
            top_z = np.argmax(peakbool) + 1
        else:
            top_z = np.argmax(maxes) + 1
    else:
        top_z = np.argmax(maxes) + 1

    best_phase = np.argmax(phases[top_z - 1])
    mask = indexes[top_z - 1] == best_phase
    return top_z, mask


def decode_emat(emat):
    """
    Decode an encoded matrix back into an isotope distribution
    :param emat: Encoded matrix (3 x maxlen x maxlen)
    :return: Isotope distribution (m/z, intensity)
    """
    digitmatrix = emat[0]
    distmatrix = emat[1]
    weightsmatrix = emat[2]

    # Decode the digit matrix
    mzs = np.zeros(len(digitmatrix))
    for i in range(len(digitmatrix)):
        mzs[i] = digitmatrix[i, 0] * 100000

    # Decode the weights matrix
    ints = np.zeros(len(weightsmatrix))
    for i in range(len(weightsmatrix)):
        ints[i] = np.sqrt(weightsmatrix[i, i])

    # b1 = ints == 0
    # ints[b1] = np.amax(ints) * 0.05
    b1 = ints > 0
    mzs = mzs[b1]
    ints = ints[b1]

    return np.vstack([mzs, ints]).T


def plot_emat(emat):
    """
    Simple plot to view the encoded matrix
    :param emat: Encoded matrix
    :return: None
    """
    plt.subplot(221)
    plt.imshow(emat[0], cmap='gray', aspect='auto')
    plt.subplot(222)
    plt.imshow(emat[1], cmap='gray', aspect='auto')
    plt.subplot(223)
    plt.imshow(emat[2], cmap='gray', aspect='auto')
    plt.subplot(224)
    # Color plot for rgb channels
    rgb = emat.transpose(1, 2, 0)
    rgb *= 255 / np.amax(rgb)
    rgb = rgb.astype(np.uint8)
    plt.imshow(rgb, aspect='auto')
    plt.show()


def save_encoding(data, outfile):
    emat = [d[0] for d in data]
    centroids = np.array([d[1] for d in data], dtype=object)
    z = [d[2] for d in data]
    print("Saving to:", outfile)
    np.savez_compressed(outfile, emat=emat, centroids=centroids, z=z)


def encode_dir(pkldir, outdir=None, name="medium", maxfiles=None, plot=False, **kwargs):
    startime = time.perf_counter()
    training = []
    test = []
    zdist = []

    files = ud.match_files_recursive(pkldir, ".pkl")

    if maxfiles is not None:
        files = files[:maxfiles]

    print("Files:", files)

    for file in files:
        if "bad_data" in file:
            continue

        tr, te, zd = encode_phase_file(file, **kwargs)

        training.extend(tr)
        test.extend(te)
        zdist.extend(zd)

    # Write out everything
    print(len(training), len(test))
    if outdir is not None:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        os.chdir(outdir)

    save_encoding(training, "training_data_" + name + ".npz")
    save_encoding(test, "test_data_" + name + ".npz")

    # torch.save(training, "training_data_" + name + ".pth")
    # torch.save(test, "test_data_" + name + ".pth")

    endtime = time.perf_counter()
    print("Time:", endtime - starttime, len(zdist))
    if plot:
        plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
        plt.show()


@njit(fastmath=True)
def encode_noise(peakmz: float, intensity: float, maxlen=16):
    mznoise = np.random.uniform(peakmz - 1.5, peakmz + 3.5, maxlen)
    intnoise = np.abs(np.random.normal(0, intensity, maxlen))

    ndata = np.empty((maxlen, 2))
    ndata[:, 0] = mznoise
    ndata[:, 1] = intnoise

    emat = encode_phase(ndata)
    return emat, ndata, 0


@njit(fastmath=True)
def encode_double(centroid, centroid2, maxdist=1.5, minsep=0.1, intmax=0.2, phaseres=8):
    peakmz = centroid[np.argmax(centroid[:, 1]), 0]
    peakmz2 = centroid2[np.argmax(centroid2[:, 1]), 0]

    # Move peak2 to a random spot around peak1
    shift = minsep + np.random.uniform(-1, 1) * maxdist
    centroid2[:, 0] += peakmz + shift - peakmz2

    # Set into to random value between 0.05 and 0.95% of the max
    centroid2[:, 1] *= np.random.uniform(0.05, intmax) * np.amax(centroid[:, 1]) / np.amax(centroid2[:, 1])

    mergedc = np.empty((len(centroid) + len(centroid2), 2))
    mergedc[:len(centroid)] = centroid
    mergedc[len(centroid):] = centroid2

    emat = encode_phase(mergedc, phaseres=phaseres)

    return emat, mergedc


def encode_phase_file(file, maxlen=8, save=True, outdir="C:\\Data\\IsoNN\\multi", name="medium",
                      onedropper=0.95):
    training = []
    test = []
    zdist = []
    # Load pkl file
    if "peaks.pkl" in file:
        return training, test, zdist

    try:
        with open(file, "rb") as f:
            centroids = pkl.load(f)
    except Exception as e:
        print("Error Loading File:", file, e)
        return training, test, zdist
    # Encode each centroid
    print("File:", file, len(centroids))
    for c in centroids:
        centroid = c[0]
        z = c[1]
        if z == 1:
            # toss it out with a 80% chance
            r = np.random.rand()
            if r < onedropper:
                continue
        zdist.append(z)
        emat = encode_phase(centroid, phaseres=maxlen)
        # emat.astype(np.float32)
        emat = torch.as_tensor(emat, dtype=torch.float32)

        # randomly sort into training and test data
        r = np.random.rand()
        if r < 0.9:
            training.append((emat, centroid, z))
        else:
            test.append((emat, centroid, z))
    return training, test, zdist


data_dirs = ["MSV000090488",
             "MSV000091923",
             "RibosomalPfms_Td_Control",
             "PXD019247",
             "PXD045560",
             "PXD046651",
             "PXD027650",
             "PXD041357",
             "PXD042921"
             ]

small_data_dirs = ["MSV000090488",
                   "PXD019247",
                   "PXD045560",
                   ]

if __name__ == "__main__":
    starttime = time.perf_counter()
    # outdir = "C:\\Data\\IsoNN\\multi"
    # file = "C:\\Data\\TabbData\\20170105_L_MaD_ColC4_Ecoli20161108-10_01_10.pkl"
    # file = "Z:\\Group Share\\JGP\\MSV000090488\\20220207_UreaExtracted_SW480_C4_RPLC_01.pkl"
    # directory = "Z:\\Group Share\\JGP\\MSV000090488\\"
    # directory = "Z:\\Group Share\\JGP\\MSV000091923"
    # directory = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control"
    # directory = "Z:\\Group Share\\JGP\\PXD019247"
    # os.chdir(directory)
    # directory = "C:\\Data\\TabbData\\"
    # outdir = "C:\\Data\\IsoNN\\training\\MSV000090488\\"

    for d in data_dirs:
        topdir = os.path.join("Z:\\Group Share\\JGP", d)
        outdir = os.path.join("C:\\Data\\IsoNN\\training", d)
        print("Directory:", topdir, "Outdir:", outdir)
        encode_dir(topdir, maxlen=2, name="phase2", onedropper=0.0, maxfiles=None,
                   outdir=outdir)
    # encode_multi_file(file, maxlen=32, nfeatures=2, save=True, name="small32x2")
    print("Time:", time.perf_counter() - starttime)

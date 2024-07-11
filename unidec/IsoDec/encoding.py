import numpy as np
import matplotlib.pyplot as plt
import os
import unidec.tools as ud
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import match_charge, cplot
import pickle as pkl
import time
import torch

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


def get_digit(number, n):
    return number // 10 ** n % 10


def mz_to_bit(mz):
    intpart = int(mz)
    decpart = mz - intpart
    intpartdigits = np.array([get_digit(intpart, i) for i in range(6)]) / 10.
    decimaldigits = np.array([get_digit(decpart * 1e10, i) for i in range(10)]) / 10.
    # combine and reverse digitis
    digits = np.concatenate([intpartdigits[::-1], decimaldigits[::-1]])
    return digits


def bit_to_mz(bit):
    intpart = 0
    decpart = 0
    for i in range(6):
        intpart += bit[i] * 10 ** (5 - i)
    for i in range(10):
        decpart += bit[i + 6] * 10 ** (9 - i)
    return intpart + decpart / 1e10


def encode_isodist(isodist, matched=None, maxlen=16):
    indexes = np.arange(len(isodist))
    # sort isodist by intensity and remove lower values
    sortindex = np.argsort(isodist[:, 1])[::-1]
    isodist = isodist[sortindex]
    indexes = indexes[sortindex]

    if matched is not None:
        matched = matched[sortindex]

    if len(isodist) > maxlen:
        isodist = isodist[:maxlen]
        indexes = indexes[:maxlen]
        if matched is not None:
            matched = matched[:maxlen]
    else:
        # Append zeros to get to len maxlen
        isodist = np.vstack([isodist, np.zeros((maxlen - len(isodist), 2))])
        indexes = np.hstack([indexes, np.zeros(maxlen - len(indexes))])
        if matched is not None:
            matched = np.hstack([matched, np.zeros(maxlen - len(matched))])

    # Set up distaince and weighting matrices
    distmatrix = np.zeros((maxlen, maxlen))
    weightsmatrix = np.zeros((maxlen, maxlen))
    for j in range(maxlen):
        for k in range(maxlen):
            distmatrix[j, k] = isodist[j, 0] - isodist[k, 0]
            weightsmatrix[j, k] = isodist[j, 1] * isodist[k, 1]

    b1 = distmatrix > 1
    b2 = distmatrix < 0
    distmatrix[b1 | b2] = 0

    digitmatrix = np.zeros((maxlen, maxlen))
    for j in range(maxlen):
        digitmatrix[j] = mz_to_bit(isodist[j, 0])

    weightsmatrix /= np.amax(weightsmatrix)

    # Put it all together
    emat = np.array([digitmatrix, distmatrix, weightsmatrix])
    if np.amax(emat) > 1:
        print(emat)
        exit()
    if np.amin(emat) < 0:
        print(emat)
        exit()
    if matched is not None:
        return emat, indexes, matched
    else:
        return emat, indexes


def plot_emat(emat):
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


def line_plot(emat):
    digits = emat[0]
    for d in digits:
        mz = bit_to_mz(d)
        plt.vlines(mz, 0, 1)
    plt.show()


def encode_dir(pkldir, type=0, outdir="C:\\Data\\IsoNN\\multi", name="medium"):
    startime = time.perf_counter()
    training = []
    test = []
    zdist = []

    files = ud.match_files_recursive(pkldir, ".pkl")
    print("Files:", files)

    for file in files:
        if type == 0:
            tr, te, zd = encode_file(file)
        elif type == 1:
            tr, te, zd = encode_multi_file(file, save=False)
        else:
            raise ValueError("Invalid Type")
        training.extend(tr)
        test.extend(te)
        zdist.extend(zd)

    # Write out everything
    print(len(training), len(test))
    os.chdir(outdir)
    torch.save(training, "training_data_" + name + ".pth")
    torch.save(test, "test_data_" + name + ".pth")

    endtime = time.perf_counter()
    print("Time:", endtime - starttime, len(centroids))

    plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
    plt.show()


def encode_file(file):
    training = []
    test = []
    zdist = []
    # Load pkl file
    try:
        with open(file, "rb") as f:
            centroids = pkl.load(f)
    except Exception as e:
        print("Error Loading File:", file, e)
        return
    # Encode each centroid
    print("File:", file, len(centroids))
    for c in centroids:
        centroid = c[0]
        z = c[1]
        if z == 1:
            # toss it out with a 80% chance
            r = np.random.rand()
            if r < 0.8:
                continue
        zdist.append(z)
        emat, indexes = encode_isodist(centroid)
        t = torch.tensor(emat, dtype=torch.float32)
        # randomly sort into training and test data
        r = np.random.rand()
        if r < 0.9:
            training.append((t, z))
        else:
            test.append((t, z))
    return training, test, zdist


def encode_multi_individual(c, maxlen=16, nfeatures=3):
    centroid = c[0]
    z = c[1]
    matched = np.array(c[2])
    bmatched = np.zeros(len(centroid))
    bmatched[matched] = 1
    isodist = c[3]
    isomatched = c[4]
    try:
        peakmz = c[5]
    except:
        peakmz = centroid[np.argmax(centroid[:, 1]), 0]
    if z == 1:
        # toss it out with a 80% chance
        r = np.random.rand()
        if r < 0.8:
            return 0, 0, 0, 0

    emat, indexes, matched2 = encode_isodist(centroid, matched=bmatched, maxlen=maxlen)
    t = torch.tensor(emat, dtype=torch.float32)

    target = np.zeros((nfeatures, maxlen))
    target[0] = matched2
    return t, target, z, peakmz


def encode_double(c, c2, maxlen=16, maxdist=1.5, minsep=2, nfeatures=3):
    # Check if the peaks are close enough together
    centroid = c[0]
    centroid2 = c2[0]
    try:
        peakmz2 = c2[5]
    except:
        peakmz2 = centroid2[np.argmax(centroid2[:, 1]), 0]

    try:
        peakmz = c[5]
    except:
        peakmz = centroid[np.argmax(centroid[:, 1]), 0]

    if np.abs(peakmz - peakmz2) > maxdist:
        return 0, 0, 0

    # Filter out peaks that are the same charge state and too close together
    z = c[1]
    z2 = c2[1]
    if z == z2:
        maxdist = minsep / z
        if np.abs(peakmz - peakmz2) < maxdist:
            return 0, 0, 0

    matched = np.array(c[2])
    bmatched = np.zeros(len(centroid))
    bmatched[matched] = 1
    isodist = c[3]
    isomatched = c[4]

    matched2 = np.array(c2[2])
    bmatched2 = np.zeros(len(centroid2))
    bmatched2[matched2] = 1
    isodist2 = c2[3]
    isomatched2 = c2[4]

    cm = centroid[bmatched.astype(bool)]
    cm2 = centroid2[bmatched2.astype(bool)]

    a1 = np.zeros(len(cm))
    a2 = np.ones(len(cm2))

    mergedc = np.vstack([cm, cm2])
    mergeda = np.hstack([a1, a2])
    if len(mergedc) > maxlen:
        # remove the lowest intensity peaks
        sortindex = np.argsort(mergedc[:, 1])[::-1]
        mergedc = mergedc[sortindex]
        mergedc = mergedc[:maxlen]
        mergeda = mergeda[sortindex]
        mergeda = mergeda[:maxlen]

    emat, indexes, mergeda = encode_isodist(mergedc, matched=mergeda, maxlen=maxlen)
    t = torch.tensor(emat, dtype=torch.float32)

    target = np.zeros((nfeatures, maxlen))

    target[0] = mergeda == 0
    target[1] = mergeda == 1

    # cplot(mergedc)
    # plt.show()

    return t, target, 1


def encode_multi_file(file, maxlen=16, nfeatures=3, save=True):
    training = []
    test = []
    zdist = []
    # Load pkl file
    try:
        with open(file, "rb") as f:
            centroids = pkl.load(f)
    except Exception as e:
        print("Error Loading File:", file, e)
        return
    # Encode each centroid
    print("File:", file, len(centroids))
    mzvals = []
    for c in centroids:
        t, target, z, mz = encode_multi_individual(c, maxlen=maxlen, nfeatures=nfeatures)
        if z == 0:
            mzvals.append(0)
            continue
        zdist.append(z)
        mzvals.append(mz)
        # randomly sort into training and test data
        r = np.random.rand()
        if r < 0.9:
            training.append((t, target))
        else:
            test.append((t, target))
    print("Phase 1", len(training), len(test))
    for i, c in enumerate(centroids):
        if mzvals[i] == 0:
            continue
        for j, c2 in enumerate(centroids):
            if i <= j:
                continue
            if mzvals[j] == 0:
                continue
            if np.abs(mzvals[i] - mzvals[j]) > 1.5:
                continue
            t, target, z = encode_double(c, c2, maxlen=maxlen, nfeatures=nfeatures)
            if z == 0:
                continue

            r = np.random.rand()
            if r < 0.9:
                training.append((t, target))
            else:
                test.append((t, target))
    print("Phase 2", len(training), len(test))
    if save:
        # Write out everything
        print(len(training), len(test))
        torch.save(training, "exp_training_data_medium.pth")
        torch.save(test, "exp_test_data_medium.pth")

    return training, test, zdist


if __name__ == "__main__":
    os.chdir("C:\\Data\\IsoNN\\multi")
    file = "C:\\Data\\TabbData\\LPT_CYTO_GFPb.pkl"
    file = "Z:\\Group Share\\JGP\\MSV000090488\\20220207_UreaExtracted_SW480_C4_RPLC_01.pkl"
    directory = "Z:\\Group Share\\JGP\\MSV000090488\\"
    encode_dir(directory, type=1, name="large")
    #encode_multi_file(file)
    pass

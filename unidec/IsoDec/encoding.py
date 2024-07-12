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

try:
    mpl.use("WxAgg")
except:
    pass

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
        sequence = mz_to_bit(isodist[j, 0])
        if len(sequence) > maxlen:
            sequence = sequence[:maxlen]
        if len(sequence) < maxlen:
            sequence = np.hstack([sequence, np.zeros(maxlen - len(sequence))])
        digitmatrix[j] = sequence

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


def encode_dir(pkldir, type=0, outdir="C:\\Data\\IsoNN\\multi", name="medium", **kwargs):
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
            tr, te, zd = encode_multi_file(file, save=False, **kwargs)
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
    print("Time:", endtime - starttime, len(zdist))

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


def encode_multi_individual(c, maxlen=16, nfeatures=3, onedropper=0.8):
    centroid = c[0]
    z = c[1]
    if z == 1:
        # toss it out with a 80% chance
        r = np.random.rand()
        if r < onedropper:
            return 0, 0, 0, 0

    matched = np.array(c[2])

    if len(matched) < 3:
        return 0, 0, 0, 0

    bmatched = np.zeros(len(centroid))
    bmatched[matched] = 1
    isodist = c[3]
    isomatched = c[4]
    try:
        peakmz = c[5]
    except:
        peakmz = centroid[np.argmax(centroid[:, 1]), 0]

    if False:
        cplot(centroid)
        cplot(centroid[bmatched.astype(bool)], color='g', factor=-1)
        plt.show()

    emat, indexes, matched2 = encode_isodist(centroid, matched=bmatched, maxlen=maxlen)
    t = torch.tensor(emat, dtype=torch.float32)

    target = np.zeros((nfeatures, maxlen))
    target[0] = matched2
    return t, target, z, peakmz


def encode_double(c, c2, maxlen=16, maxdist=1.5, minsep=0.1, nfeatures=3):
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

    # Filter out peaks that are the same charge state and too close together
    z = c[1]
    z2 = c2[1]

    matched = np.array(c[2])
    bmatched = np.zeros(len(centroid))
    bmatched[matched] = 1

    matched2 = np.array(c2[2])
    bmatched2 = np.zeros(len(centroid2))
    bmatched2[matched2] = 2

    # Move peak2 to a random spot around peak1
    shift = minsep + np.random.uniform(-1, 1) * maxdist
    centroid2[:, 0] += peakmz + shift - peakmz2

    # Set into to random value between 0.05 and 0.95% of the max
    centroid2[:, 1] *= np.random.uniform(0.05, 0.95) * np.amax(centroid[:, 1]) / np.amax(centroid2[:, 1])

    mergedc = np.vstack([centroid, centroid2])
    mergeda = np.hstack([bmatched, bmatched2])

    if len(mergedc) > maxlen:
        # remove the lowest intensity peaks
        sortindex = np.argsort(mergedc[:, 1])[::-1]
        mergedc = mergedc[sortindex]
        mergedc = mergedc[:maxlen]
        mergeda = mergeda[sortindex]
        mergeda = mergeda[:maxlen]
    elif len(mergedc) < maxlen:
        # Append zeros to get to len maxlen
        mergedc = np.vstack([mergedc, np.zeros((maxlen - len(mergedc), 2))])
        mergeda = np.hstack([mergeda, np.zeros(maxlen - len(mergeda))])

    emat, indexes = encode_isodist(mergedc, maxlen=maxlen)
    t = torch.tensor(emat, dtype=torch.float32)

    indexes = indexes.astype(int)

    target = np.zeros((nfeatures, maxlen))

    ma = mergeda == 1
    ma2 = mergeda == 2

    mergedc = mergedc[indexes]
    target[0] = ma[indexes]
    target[1] = ma2[indexes]

    if False:
        cplot(mergedc[mergedc[:, 1] > 0])

        cm1 = centroid[matched]
        cm22 = centroid2[matched2]
        cm = mergedc[target[0]]
        cm2 = mergedc[target[1]]

        cplot(cm, color='g', factor=-1)
        cplot(cm1, color='b', factor=-1)
        cplot(cm2, color='cyan', factor=-1, base=np.amax(centroid[:, 1]) * -1)

        plt.show()

    return t, target, 1


def encode_multi_file(file, maxlen=16, nfeatures=3, save=True, outdir="C:\\Data\\IsoNN\\multi", name="medium",
                      onedropper=0.9):
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
    ints = []
    for c in centroids:
        t, target, z, mz = encode_multi_individual(c, maxlen=maxlen, nfeatures=nfeatures, onedropper=onedropper)
        if z == 0:
            mzvals.append(0)
            ints.append(0)
            continue
        zdist.append(z)
        mzvals.append(mz)
        ints.append(np.amax(c[0][:, 1]))
        # randomly sort into training and test data
        r = np.random.rand()
        if r < 0.9:
            training.append((t, target))
        else:
            test.append((t, target))
    print("Phase 1", len(training), len(test))
    count = 0

    lp1 = len(training)
    maxcount = lp1 * 2
    for i, c in enumerate(centroids):
        if mzvals[i] == 0:
            continue
        if count > maxcount:
            break
        found = False
        count2 = 0
        while not found and count2 < 10:

            j = np.random.randint(0, len(centroids))
            c2 = centroids[j]

            if i == j:
                count2 += 1
                continue
            if mzvals[j] == 0:
                count2 += 1
                continue

            t, target, z = encode_double(c, c2, maxlen=maxlen, nfeatures=nfeatures)

            if z == 0:
                count2 += 1
                continue

            found = True
            r = np.random.rand()
            if r < 0.9:
                training.append((t, target))
                count += 1
            else:
                test.append((t, target))
    print("Phase 2", len(training), len(test))
    if save:
        os.chdir(outdir)
        # Write out everything
        print(len(training), len(test))
        torch.save(training, "training_data_" + name + ".pth")
        torch.save(test, "test_data_" + name + ".pth")

    return training, test, zdist


if __name__ == "__main__":
    starttime = time.perf_counter()
    os.chdir("C:\\Data\\IsoNN\\multi")
    file = "C:\\Data\\TabbData\\20170105_L_MaD_ColC4_Ecoli20161108-10_01_10.pkl"
    # file = "Z:\\Group Share\\JGP\\MSV000090488\\20220207_UreaExtracted_SW480_C4_RPLC_01.pkl"
    directory = "Z:\\Group Share\\JGP\\MSV000090488\\"
    directory = "Z:\\Group Share\\JGP\\MSV000091923"
    #directory = "C:\\Data\\TabbData\\"
    encode_dir(directory, type=1, maxlen=32, nfeatures=2, name="large32x2", onedropper=0.95)
    #encode_multi_file(file, maxlen=32, nfeatures=2, save=True, name="small32x2")
    print("Time:", time.perf_counter() - starttime)

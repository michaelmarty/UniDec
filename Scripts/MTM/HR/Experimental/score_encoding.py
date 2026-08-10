import numpy as np
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.plots import cplot
from unidec.IsoDec.match import create_isodist, find_matches
import pickle as pkl
from copy import deepcopy
import time
import unidec.tools as ud

def create_vectors(centroids, isodist, matched, isomatched, isolen=16):
    if len(matched) != len(isomatched):
        print("ERROR: Length mismatch")
        raise ValueError
    pred_vector = deepcopy(isodist[:, 1])
    exp_vector = np.zeros_like(pred_vector)
    exp_vector[isomatched] = centroids[matched, 1]

    # Pad or cut to isolen
    if len(pred_vector) < isolen:
        pred_vector = np.pad(pred_vector, (0, isolen - len(pred_vector)), mode="constant")
        exp_vector = np.pad(exp_vector, (0, isolen - len(exp_vector)), mode="constant")
    elif len(pred_vector) > isolen:
        pred_vector = pred_vector[:isolen]
        exp_vector = exp_vector[:isolen]

    # Normalize both
    if np.amax(pred_vector) > 0:
        pred_vector /= np.amax(pred_vector)

    if np.amax(exp_vector) > 0:
        exp_vector /= np.amax(exp_vector)
    # Diff vector
    diff_vector = pred_vector - exp_vector

    # concatenate exp and diff vectors
    combined_vector = np.concatenate([exp_vector, diff_vector])
    return combined_vector


def vec_from_assignment(centroids, z, isolen=16):
    isodist = create_isodist(centroids[np.argmax(centroids[:, 1]), 0], z, centroids)
    matchedindexes, isomatches = find_matches(centroids, isodist, 5.0)
    return create_vectors(centroids, isodist, matchedindexes, isomatches, isolen=isolen)


def get_dummy_charge_state(true_charge):
    fake_charge = -1
    while fake_charge == -1 or fake_charge == true_charge:
        # fake_charge = np.random.randint(1, 5000000)
        fake_charge = int(round(1 + np.random.random() ** 4 * 49))
    if fake_charge == true_charge:
        raise ValueError("Fake charge is the same as true charge")
    return fake_charge

def create_dummy_vector(centroids, true_charge, peakmz, isolen=16):
    fake_charge = get_dummy_charge_state(true_charge)

    isodist = create_isodist(float(peakmz), int(fake_charge), centroids)

    matchedindexes, isomatches = find_matches(centroids, isodist, 5.0)

    # cplot(isodist, color="g")

    dv = create_vectors(centroids, isodist, matchedindexes, isomatches, isolen=isolen)
    return dv

def encode_file(path, testfrac=0.1, isolen=16):
    print("File path:", path)
    dataset = np.load(path, allow_pickle=True)

    training = []
    test = []

    for data in dataset:
        if len(data) == 7:
            centroids = data[0]
            charge = data[1]
            matched = data[2]
            isodist = data[3]
            isomatched = data[4]
            peakmz = data[5]
            # scan = data[6]
        elif len(data) == 5:
            centroids = data[0]
            charge = data[1]
            matched = data[2]
            isodist = data[3]
            isomatched = data[4]
            peakmz = centroids[np.argmax(centroids[:, 1]), 0]
            # scan = -1
        else:
            print("ERROR:", len(data))
            raise ValueError

        cv = create_vectors(centroids, isodist, matched, isomatched, isolen=isolen)
        dv = create_dummy_vector(centroids, charge, peakmz, isolen)

        if np.random.rand() < testfrac:
            test.append((cv, 1))
            test.append((dv, 0))
        else:
            training.append((cv, 1))
            training.append((dv, 0))


        # cplot(centroids, isodist=isodist, zfactor=-1)
        # plt.show()

    return training, test

def save_encoding(data, outfile):
    vecs = [d[0] for d in data]
    res = [d[1] for d in data]
    print("Saving to:", outfile)
    np.savez_compressed(outfile, vecs=vecs, res=res)


def encode_dir(pkldir, outdir=None, name="medium", maxfiles=None, **kwargs):
    startime = time.perf_counter()
    training = []
    test = []

    files = ud.match_files_recursive(pkldir, ".pkl")

    if maxfiles is not None:
        files = files[:maxfiles]

    print("Files:", files)

    for file in files:
        if "bad_data" in file:
            continue

        tr, te = encode_file(file, **kwargs)

        training.extend(tr)
        test.extend(te)

    # Write out everything
    print(len(training), len(test))
    if outdir is not None:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        os.chdir(outdir)
    else:
        os.chdir(pkldir)

    save_encoding(training, "training_data_" + name + ".npz")
    save_encoding(test, "test_data_" + name + ".npz")

    print("Time:", time.perf_counter() - startime)



if __name__ == "__main__":
    path = "Z:\\Group Share\\JGP\\MSV000090488\\20220207_UreaExtracted_SW480_C4_RPLC_01.pkl"
    outpath = "Z:\\Group Share\\JGP\\MSV000090488\\short.pkl"

    dir = "Z:\\Group Share\\JGP\\MSV000090488\\"

    encode_dir(dir, maxfiles=1, name="score1")

    # dataset = np.load(path, allow_pickle=True)
    # shotdataset = dataset[:500:100]
    #
    # # np.save(outpath, shotdataset, allow_pickle=True)
    # pkl.dump(shotdataset, open(outpath, "wb"))

    #encode_file(outpath)

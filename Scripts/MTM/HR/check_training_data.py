import time
from unidec.IsoDec.encoding import data_dirs, small_data_dirs
import os
import numpy as np
import matplotlib.pyplot as plt
from unidec.IsoDec.plots import cplot

smalldir2 = ["PXD027650"]

topdirectory = "C:\\Data\\IsoNN\\training"

def pull_centroids():

    dirs = [os.path.join(topdirectory, d) for d in data_dirs]
    # dirs = [os.path.join(topdirectory, d) for d in smalldir2]

    zs = []
    cs = []

    for d in dirs:
        training_path = "phase83"
        ext = ".npz"
        training_path = "training_data_" + training_path + ext
        training_path = os.path.join(d, training_path)
        td = np.load(training_path, allow_pickle=True)
        training_data = [td["emat"], td["z"], list(td["centroids"])]
        print(d, len(training_data[0]))
        zvals = training_data[1]
        cvals = training_data[2]

        for i, z in enumerate(zvals):
            if z < 40:
                continue
            print("Z:", z)
            zs.append(z)
            cs.append(cvals[i])

    print(len(zs))

    os.chdir(topdirectory)
    # Write out zs and cs
    zs = np.array(zs)
    zs.tofile("zs.txt", sep="\n")
    cs = np.array(cs, dtype=object)
    np.save("cs.npy", cs, allow_pickle=True)

def check_highz():
    os.chdir(topdirectory)
    zs = np.loadtxt("zs.txt")
    cs = np.load("cs.npy", allow_pickle=True)
    zs = zs.astype(int)

    ncorrect = 0
    print(len(zs))

    for i, z in enumerate(zs):
        c = cs[i]
        c = np.array(c)
        plt.figure(figsize=(12, 8))
        cplot(c, z=z)
        peakmz = c[np.argmax(c[:, 1]), 0]
        plt.xlim(peakmz-0.4, peakmz+0.4)
        plt.show()
        # Query user
        x = input("Is this correct? (y/n): ")
        if x == "y":
            ncorrect += 1

    print("Correct:", ncorrect)
    print("Total:", len(zs))
    print("Percent Correct:", ncorrect / len(zs))

def pull_random_centroids():
    dirs = [os.path.join(topdirectory, d) for d in data_dirs]

    total = 5868491
    # Pick 100 random values from the total
    n= 100
    zindexes = np.random.choice(total, n, replace=False)
    print(zindexes)

    zs = []
    cs = []
    count = 0

    for d in dirs:
        training_path = "phase83"
        ext = ".npz"
        training_path = "training_data_" + training_path + ext
        training_path = os.path.join(d, training_path)
        td = np.load(training_path, allow_pickle=True)
        # training_data = [td["emat"], td["z"], list(td["centroids"])]
        print(d, len(td["z"]))

        indexes = zindexes - count
        indexes = indexes[indexes > 0]
        indexes = indexes[indexes < len(td["z"])]

        for i in indexes:
            print("Z:", td["z"][i])
            zs.append(td["z"][i])
            cs.append(td["centroids"][i])

        count += len(td["z"])

    print(count)
    os.chdir(topdirectory)

    # Write out zs and cs
    zs = np.array(zs)
    zs.tofile("zs_rand.txt", sep="\n")
    cs = np.array(cs, dtype=object)
    np.save("cs_rand.npy", cs, allow_pickle=True)

def check_random():
    os.chdir(topdirectory)
    zs = np.loadtxt("zs_rand.txt")
    cs = np.load("cs_rand.npy", allow_pickle=True)
    zs = zs.astype(int)

    ncorrect = 0
    print(len(zs))

    for i, z in enumerate(zs):
        if z != 2:
            continue
        c = cs[i]
        c = np.array(c)
        plt.figure(figsize=(12, 8))
        cplot(c, z=z, zfactor=-1)
        plt.title("+" + str(z))
        peakmz = c[np.argmax(c[:, 1]), 0]
        # plt.xlim(peakmz-0.4, peakmz+0.4)
        plt.show()
        # Query user
        x = input("Is this correct? (y/n): ")
        if x == "y":
            ncorrect += 1

    print("Correct:", ncorrect)
    print("Total:", len(zs))
    print("Percent Correct:", ncorrect / len(zs))

if __name__ == "__main__":
    # pull_centroids()
    # check_highz()
    # pull_random_centroids()
    check_random()
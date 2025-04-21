import os
import pickle as pkl
import matchms
import time
import numpy as np
import matplotlib.pyplot as plt

starttime = time.perf_counter()

dir = "Z:\\Group Share\\JGP\\MSV000091923\\raw"
dir = "Z:\\Group Share\\JGP\\MSV000090488"
dir = "Z:\\Group Share\\JGP\\MSV000091693"

os.chdir(dir)
files = os.listdir(dir)
print(files)
zdat = []
for file in files:
    ftime = time.perf_counter()
    extension = "." + file.split(".")[-1]
    if extension == ".pkl":
        print(file)
        with open(file, "rb") as f:
            centroids = pkl.load(f)
        # Encode each centroid

        fixed = []

        print(len(centroids))

        for c in centroids:
            centroid = c[0]
            z = c[1]
            centmatch = c[2]
            isodist = c[3]
            isomatch = c[4]

            zdat.append(z)

            matchingpeaks = matchms.similarity.spectrum_similarity_functions.find_matches(centroid[:, 0],
                                                                                          isodist[:, 0], 0.01)
            matchedindexes = [match[0] for match in matchingpeaks]
            isomatches = [match[1] for match in matchingpeaks]

            if len(matchedindexes) == 0:
                print("No Matches")
                continue

            outlist = [centroid, z, matchedindexes, isodist, isomatches]
            fixed.append(outlist)

        print("Time:", time.perf_counter()-ftime, len(fixed))

        if len(fixed) != len(centroids):
            print("Error")
            continue

        with open(file, "wb") as f:
            pkl.dump(fixed, f)


plt.hist(zdat, bins=np.arange(0.5, 50.5, 1))
plt.show()
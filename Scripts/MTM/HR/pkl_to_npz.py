import pickle as pkl
import numpy as np
import os
import torch

topdir = "C:\\Data\\IsoNN\\training\\"

for root, dirs, files in os.walk(topdir):
    for file in files:
        if file.endswith(".pth"):

            infile = os.path.join(root, file)
            print(infile)
            outfile = os.path.splitext(infile)[0] + ".npz"
            if os.path.isfile(outfile) or "bad_data" in infile:
                continue
            data = torch.load(infile)

            emat = [d[0] for d in data]
            centroids = np.array([d[1] for d in data], dtype=object)
            z = [d[2] for d in data]

            np.savez_compressed(outfile, emat=emat, centroids=centroids, z=z)

exit()

import numpy as np
import pickle as pkl
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.encoding import *
from unidec.IsoDec.match import cplot

topdir = "C:\\Data\\IsoNN\\training"
directory = small_data_dirs[2]
directory = os.path.join(topdir, directory)
path = os.path.join(directory, "bad_data.pkl")

with open(path, "rb") as f:
    bad_data = pkl.load(f)
    print(f"Loaded {len(bad_data)} bad data points from {path}")
dtype = 1
for i, bd in enumerate(bad_data):
    centroids = bd[0]
    z = int(bd[1])
    pred_z = int(bd[2])
    mask = None
    pred_mask = None
    print(z, pred_z)
    plt.subplot(121)
    plt.title("Training")
    cplot(centroids, mask=mask, z=z)
    plt.subplot(122)
    plt.title("Prediction")
    cplot(centroids, mask=pred_mask, z=pred_z)
    plt.yticks([])
    plt.show()

from IsoDec.encoding import *
import torch
import pickle as pkl
import time
import unidec.tools as ud
import os
import numpy as np

# starttime = time.perf_counter()
#
# os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
#
# outdir = "C:\\Data\\IsoNN\\exp_centroids"
# os.chdir(outdir)
#
# pkldir = "Z:\\Group Share\\JGP"
#
# training = []
# test = []
# zdist = []
#
# files = ud.match_files_recursive(pkldir, ".pkl")
# print("Files:", files)
#
#
# for file in files:
#     # Load pkl file
#     try:
#         with open(file, "rb") as f:
#             centroids = pkl.load(f)
#     except Exception as e:
#         print("Error Loading File:", file, e)
#         continue
#     # Encode each centroid
#     print("File:", file, len(centroids))
#     for c in centroids:
#         centroid = c[0]
#         z = c[1]
#         if z == 1:
#             # toss it out with a 80% chance
#             r = np.random.rand()
#             if r < 0.8:
#                 continue
#         zdist.append(z)
#         emat, indexes = encode_isodist(centroid)
#         t = torch.tensor(emat, dtype=torch.float32)
#         # randomly sort into training and test data
#         r = np.random.rand()
#         if r < 0.9:
#             training.append((t, z))
#         else:
#             test.append((t, z))
#
# # Write out everything
# print(len(training), len(test))
# torch.save(training, "exp_training_data_xlarge.pth")
# torch.save(test, "exp_test_data_xlarge.pth")
#
# endtime = time.perf_counter()
# print("Time:", endtime - starttime, len(centroids))
#
# plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
# plt.show()

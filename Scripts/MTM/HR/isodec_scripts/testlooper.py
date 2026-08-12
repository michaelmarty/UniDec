import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
from IsoDec.engine import IsoDecEngine
import time

eng = IsoDecEngine()

file = "Z:\\Group Share\\JGP\\MSV000091923\\20220414_Ecoli_ProtMix_C_R1.pkl"

with open(file, "rb") as f:
    data = pkl.load(f)

os.chdir("C:\\Data\\IsoNN\\training")

correct = 0
incorrect = 0
total = len(data)
starttime = time.perf_counter()
baddata = []
for i, d in enumerate(data):
    centroid = d[0]
    centroid = np.array(centroid)
    z = d[1]

    zpred = eng.phase_predictor(centroid)
    if i % 1000 == 99:
        print(i, total, correct, incorrect, correct / i)
        print(z, zpred)

    if z == zpred:
        correct += 1
    else:
        incorrect += 1
        baddata.append([centroid, z, zpred])

with open("bad_data.pkl", "wb") as f:
    pkl.dump(baddata, f)

print(correct, incorrect, correct / total * 100, incorrect / total * 100)
print("Done. Time:", time.perf_counter() - starttime, "Time per centroid:", (time.perf_counter() - starttime) / total)

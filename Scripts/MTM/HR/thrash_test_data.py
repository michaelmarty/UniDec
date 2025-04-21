from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.encoding import data_dirs, small_data_dirs, extract_centroids
import os
from unidec.IsoDec.datatools import fastnearest
import numpy as np
import matplotlib.pyplot as plt
from unidec.IsoDec.plots import cplot

smalldir2 = ["PXD027650"]

topdirectory = "C:\\Data\\IsoNN\\training"

if __name__ == "__main__":
    dirs = [os.path.join(topdirectory, d) for d in data_dirs]

    eng = IsoDecEngine(phaseres=8, use_wrapper=False)

    correct = 0
    total = 0

    nn_correct = 0

    for d in dirs[0:1]:
        training_path = "phase83"
        ext = ".npz"
        test_path = "test_data_" + training_path + ext
        test_path = os.path.join(d, test_path)
        td = np.load(test_path, allow_pickle=True)
        # training_data = [td["emat"], td["z"], list(td["centroids"])]
        print(d, len(td["z"]))

        for i, z in enumerate(td["z"]):
            # Random test to only do this 0.1% of the time
            if np.random.rand() > 0.01:
                continue
            # if z == 1:
            #     continue

            c = td["centroids"][i]
            peakmz = c[np.argmax(c[:, 1]), 0]
            lowmz = -1.05
            highmz = 2.05
            # Find all centroids in the neighborhood of the peak
            start = fastnearest(c[:, 0], peakmz + lowmz)
            end = fastnearest(c[:, 0], peakmz + highmz) + 1
            c = c[start:end]

            datathresh = 0.05

            if datathresh > 0:
                b1 = c[:, 1] > (datathresh * np.amax(c[:, 1]))
                new_c = c[b1]
                if len(new_c) < 3:
                    continue
                c = new_c

            predicted_z = eng.thrash_predictor(c)[0]

            nn_predicted_z = eng.phase_predictor(c)[0]

            total += 1
            if predicted_z == z:
                correct += 1

            if nn_predicted_z == z:
                nn_correct += 1

            # print("Z:", z)
            # print("Predicted Z:", predicted_z)
            #
            # if predicted_z == z:
            #     print("Correct")
            # else:
            #     print("Incorrect")
            #     plt.subplot(121)
            #     cplot(c, z=predicted_z)
            #     plt.subplot(122)
            #     cplot(c, z=z)
            #     plt.show()

    print("Total:", total)
    print("Thrash Correct:", correct)

    print("Thrash Accuracy:", correct / total)

    print("NN Correct:", nn_correct)
    print("NN Accuracy:", nn_correct / total)
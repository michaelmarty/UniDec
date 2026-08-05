import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
from unidec.IsoDec.runtime import IsoDecRuntime
import time
import numpy as np

def get_profile_from_pks(reader, pks):
    for p in pks:
        p.pull_from_reader(reader)

if __name__ == "__main__":
    directory = r"C:\\Data\\PXD069439"
    file = r"TDP_LH_F3.raw"


    full_path = os.path.join(directory, file)

    eng = IsoDecRuntime()
    eng.config.css_thresh = 0.9
    eng.config.knockdown_rounds = 1
    scan_range = [2400, 2450]
    scan_range = None
    reader = eng.process_file(full_path, scans=scan_range)
    start_time = time.perf_counter()
    get_profile_from_pks(reader, eng.pks)
    end_time = time.perf_counter()
    print(f"Time to get profile data: {end_time - start_time:.2f} seconds")

    for p in eng.pks:
        # Random color
        color = (np.random.rand(), np.random.rand(), np.random.rand())
        plt.plot(p.profile_data[:, 0], p.profile_data[:, 1], color = color)
        plt.plot(p.centroids[:,0], p.centroids[:,1], "o", color = color)
        plt.plot()
    plt.title("Profile Data for Peaks")
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.show()
    exit()




    tic = importer.get_tic()

    plt.plot(tic[:, 0], tic[:, 1])
    plt.title("TIC")
    plt.xlabel("Time (min)")
    plt.ylabel("Intensity")
    plt.show()

    for s in [1]:
        centroid = importer.grab_centroid_data(s)
        spectrum = importer.get_single_scan(s)


        print(len(centroid), len(spectrum))
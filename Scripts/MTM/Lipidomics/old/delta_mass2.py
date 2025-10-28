import pandas as pd
import os
import matplotlib.pyplot as plt
import unidec.tools as ud
import numpy as np

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

#file = "Mz_0_2022729846.txt"
mzML_files = ud.match_files(toppath, "*.mzml", exclude="fixed")

for file in mzML_files:
    txtfile = file[:-5] + "_fixed2.txt"

    df = pd.read_table(txtfile, header=0)

    df = df[df["MS/MS matched"]]

    avgmz = df["Precursor m/z"].to_numpy()
    refmz = df["Reference m/z"].to_numpy()

    deltamz = avgmz-refmz
    avgd = np.mean(deltamz)
    print(avgd)

    fitd = np.polyfit(avgmz, deltamz, 1)
    fitdat = avgmz * fitd[0] + fitd[1]


    plt.subplot(121)
    plt.hist(deltamz)
    plt.vlines(avgd, 0, 100)
    plt.subplot(122)
    plt.plot(avgmz, deltamz, linestyle=" ", marker="o")
    plt.plot(avgmz, fitdat)

    plt.show()
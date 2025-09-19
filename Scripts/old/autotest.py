import sys

sys.path.append("C:\\Python\\UniDec")
import numpy as np
import matplotlib.pyplot as plt
from UniDec import unidec
import os
import fnmatch

toppath = "C:\\Python\\UniDec\\TestSpectra"


def get_files(directory):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*.txt'):
            print(file)
            path = os.path.join(directory, file)
            files.append(path)

    return files


def continuous_peak_width(data, nwindow, psfun):
    starts = np.arange(data[0, 0], data[len(data) - 1, 0], len(data) / nwindow)
    ends = starts + 2 * len(data) / nwindow
    widths = []
    plt.figure()
    for i, start in enumerate(starts):
        cdata = limit_data(data, start, ends[i])
        if np.amax(cdata[:, 1]) > 0.01 * np.amax(data[:, 1]):
            # mask=ndis(cdata[:,0],cdata[len(data)/nwindow,0],len(data)/nwindow)
            # cdata[:,1]*=mask
            try:
                fwhm, junk, maxval = auto_peak_width(cdata, psfun)
            except:
                fwhm = 0
                maxval = 0
            widths.append([maxval, fwhm])
            plt.plot(cdata[:, 0], cdata[:, 1] + i / float(nwindow))
    widths = np.array(widths)
    widths = widths[widths[:, 0] > 0]
    widths = widths[widths[:, 1] > 0]
    fit, fitdat = poly_fit(widths, 2)

    cpw = data
    boo1 = cpw[:, 0] < np.amin(widths[:, 0])
    boo2 = cpw[:, 0] > np.amax(widths[:, 0])
    cpw[:, 1] = np.poly1d(fit)(cpw[:, 0])
    cpw[boo1, 1] = fitdat[0]
    cpw[boo2, 1] = fitdat[len(widths) - 1]

    plt.figure()
    plt.plot(widths[:, 0], widths[:, 1], marker="o")
    plt.plot(widths[:, 0], fitdat)
    plt.plot(cpw[:, 0], cpw[:, 1])


files = get_files(toppath)

eng = unidec.UniDec()

for file in files[::-1]:
    eng.open_file(file)
    eng.process_data()
    continuous_peak_width(eng.data.data2, 10, 0)
    plt.show()

    # plt.figure()
    # plt.plot(eng.data.data2[:, 0], eng.data.data2[:, 1])
    # plt.show()
    eng.get_auto_peak_width()

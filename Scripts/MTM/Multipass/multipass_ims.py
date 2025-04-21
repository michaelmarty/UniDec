import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud


def run_fft(data, mintime=0, maxtime=60):
    ftdata = ud.rfft(data)
    b1 = ftdata[:, 0] > 0

    if mintime == 0:
        mintime = 1e-12

    b3 = ftdata[:, 0] < 1/mintime
    b2 = ftdata[:, 0] > 1/maxtime
    b1 *= b2 * b3
    ftdata = ftdata[b1]
    ftdata[:, 0] = 1. / ftdata[:, 0]
    return ftdata


if __name__ == "__main__":
    path = "C:\Data\Varun\STP_z14_allPasses_combined_raw.csv"

    data = np.loadtxt(path, delimiter=",")

    ftdata = ud.fft(data)
    b1 = ftdata[:, 0] > 0
    #b1 = ftdata[:, 0] > 0.02
    #b2 = ftdata[:, 0] < 0.06
    #b1 *= b2
    ftdata = ftdata[b1]

    maxindex = np.argmax(ftdata[:, 1])
    maxtime = 1 / ftdata[maxindex, 0]
    print(maxtime)

    plt.figure()

    plt.subplot(121)
    plt.title(path)
    plt.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]))
    plt.xlabel("Time (ms)")
    plt.subplot(122)
    plt.plot(data[:, 0] - 12, data[:, 1] / np.amax(data[:, 1]))
    plt.plot(1 / ftdata[:, 0], ftdata[:, 1])
    plt.xlabel("Time per Cycle (ms)")
    plt.text(maxtime, ftdata[maxindex, 1], str(maxtime))
    #plt.xlim(np.amin(1 / ftdata[:, 0]), np.amax(1 / ftdata[:, 0]))
    plt.vlines(maxtime, 0, ftdata[maxindex, 1], linestyles='--', color="k")
    plt.show()

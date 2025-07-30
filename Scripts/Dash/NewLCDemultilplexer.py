import os
import numpy as np
import matplotlib.pyplot as plt

def nonlinearize(datatop, num_compressed, gap=5):
    """
    Compress the data in a simple and nonlinear way.
    :param datatop: Data array (N x 2)
    :param num_compressed:
    :return:
    """
    if num_compressed == 0:
        return datatop
    else:
        num_compressed = int(num_compressed)
        gap = int(gap)
        return np.array([np.mean(datatop[index:index + num_compressed-gap], axis=0) for index in
                         range(0, len(datatop), num_compressed)])

def calc_start_index(data, numcompressed):
    maxy = np.max(data[:, 1])
    ymaxindex = np.argmax(data[:, 1])
    threshold = 0.95 * maxy
    while ymaxindex > 0 and data[ymaxindex, 1] > threshold:
        ymaxindex -= 1
    return (ymaxindex+1) % numcompressed

def parse_rg_data(data, numcompressed, gap=5):
    startindex = calc_start_index(data, numcompressed)
    pdata = nonlinearize(data[startindex:], numcompressed, gap=gap)

    rb = pdata[:, 1] < -100
    gb = np.logical_not(rb)

    gdata = pdata[gb]
    rdata = pdata[rb]

    rdata[:, 1] = rdata[:, 1] * -1  # Flip the red channel to positive values
    rdata[:, 1] = rdata[:, 1] - np.min(rdata[:, 1])  # Subtract the minimum to avoid negative values

    gdata[:, 1] = gdata[:, 1] - np.min(gdata[:, 1])  # Subtract the minimum to avoid negative values
    return rdata, gdata


dir = "Z:\\Group Share\\Dash\\UIB data test"
file = "001-D1F-A4-MSP1D1sfGFP(+) UIB1B_data.txt"

for file in os.listdir(dir):
    if not file.endswith(".txt"):
        continue
    print("Processing file:", file)


    data = np.loadtxt(os.path.join(dir, file))
    numcompressed = 50
    rdata, gdata = parse_rg_data(data, numcompressed, gap=10)

    plt.plot(rdata[:, 0], rdata[:, 1], label='Red Channel', color='red')
    plt.plot(gdata[:, 0], gdata[:, 1], label='Green Channel', color='green')
    plt.show()
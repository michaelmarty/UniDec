import os
import numpy as np
import matplotlib.pyplot as plt
from unidec.tools import integrate

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

def calc_start_index(data, num_compressed):
    adata = np.abs(data[:, 1])
    maxy = np.max(adata)
    ymaxindex = np.argmax(adata)
    threshold = 0.95 * maxy
    while ymaxindex > 0 and np.abs(data[ymaxindex, 1]) > threshold:
        ymaxindex -= 1
    return (ymaxindex+1) % num_compressed

def parse_rg_data(data, num_compressed, gap=5, negcolor="green", ramp=0):
    startindex = calc_start_index(data, num_compressed)
    print("Start Index:", startindex)

    parseddata = nonlinearize(data[startindex:], num_compressed, gap=gap)

    maxarg = np.argmax(parseddata[:, 1])
    minarg = np.argmin(parseddata[:, 1])

    print("Max Arg:", maxarg, "Min Arg:", minarg)

    nindex = minarg % 2
    pindex = maxarg % 2

    if nindex == pindex:
        nindex += 1 ### MANUALLY ADJUST IT HERE
        print("WARNING!!!: Negative and positive indices are the same, adjusting positive index.")
        # raise ValueError("Negative and positive indices are the same, check the data.")

    ndata = parseddata[nindex::2]  # Negative data
    pdata = parseddata[pindex::2]  # Positive data

    #invert ndata
    ndata[:, 1] = ndata[:, 1] * -1  # Flip the negative channel to positive values

    if negcolor == "green":
        gdata = ndata
        rdata = pdata
    else:
        gdata = pdata
        rdata = ndata

    rdata[:, 1] = rdata[:, 1] - np.min(rdata[:, 1])  # Subtract the minimum to avoid negative values
    gdata[:, 1] = gdata[:, 1] - np.min(gdata[:, 1])  # Subtract the minimum to avoid negative values

    if ramp > 0:
        rdata[:, 1] = rdata[:, 1] * ramp

    return rdata, gdata, parseddata

def calc_percents(data, range1 = [2.2, 5], range2=[6.5,9]):
    """
    Calculate the percentage of data in the specified ranges.
    :param data: Data array (N x 2)
    :param range1: First range to calculate percentage
    :param range2: Second range to calculate percentage
    :return: Tuple of percentages in the two ranges
    """

    s1, _ = integrate(data, range1[0], range1[1])
    s2, _ = integrate(data, range2[0], range2[1])

    total = s1 + s2

    if total==0:
        percent1 = 0
        percent2 = 0
    else:
        percent1 = s1 / total * 100
        percent2 = s2 / total * 100

    return percent1, percent2


dir = "Z:\\Group Share\\Dash\\UIB data test"
file = "001-D1F-A4-MSP1D1sfGFP(+) UIB1B_data.txt"
dir = "Z:\\Group Share\\Dash\\Current Data\\Agilent LC Data\\LXEngine\\20250729 TD UIB_FLDTest_MultiplexBoost 2025-07-30 16-36-56"
file = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\20250729 TD UIB_FLDTest_MultiplexBoost 2025-07-30 16-36-56\009-D1F-A6-MSP1D1sfGFP(+)-MSP1D1mCherry(.D\processed\UIB1B_data.txt"
file = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\20251030 LX\20251029 Control 2025-10-29 13-01-55\sfGFP\006-D1F-C1-d1gfp_disc.D\Text files\UIB1B_rawdata.txt"
file = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\Channel 2 Tests\Tracer test\20251111_fld validation 2025-11-12 12-53-29\005-D1F-B1-Trace Mix_2.D\Text files\UIB1B_rawdata.txt"
#file = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\Channel 2 Tests\Tracer test\20251111_fld validation 2025-11-12 12-53-29\003-NV-blank.D\Text files\UIB1B_rawdata.txt"
file = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\20251209_sfGFP_mCherryXL 2025-12-09 12-08-53\Mix\007-D1F-E3-sfGFP(+)_mCherryXL(-).D\Text files\UIB1B_rawdata.txt"


for file in [file]:#os.listdir(dir):
    if not file.endswith(".txt"):
        continue
    print("Processing file:", file)


    data = np.loadtxt(os.path.join(dir, file))
    plt.plot( data[:,0], data[:, 1], label='Original Data', color='black', marker = 'o', markersize=1, alpha=0.1)
    num_compressed = 50
    rdata, gdata, pdata = parse_rg_data(data, num_compressed, gap=30, negcolor="red", ramp=1)

    r1 = [2.5, 5]
    r2 = [6.5, 9]

    rps = np.round(calc_percents(gdata, r1, r2))
    gps = np.round(calc_percents(rdata, r1, r2))

    print("Red Channel Percent Unbound:", rps[0], "Bound:", rps[1])
    print("Green Channel Percent Unbound:", gps[0], "Bound:", gps[1])

    # plt.plot(data[:, 0], data[:, 1], label='Parsed Data', color='black', alpha=0.5)
    # plt.plot(pdata[:, 0], pdata[:, 1], label='Compressed Data', color='blue', alpha=0.5)

    plt.plot(rdata[:, 0], rdata[:, 1], label='Red Channel', color='red')
    plt.plot(gdata[:, 0], gdata[:, 1], label='Green Channel', color='green')

    # Add vertical lines for the ranges
    plt.axvline(r1[0], color='grey', linestyle='--', label=f'Red Range {r1[0]}-{r1[1]}')
    plt.axvline(r1[1], color='grey', linestyle='--')
    plt.axvline(r2[0], color='blue', linestyle='--', label=f'Blue Range {r2[0]}-{r2[1]}')
    plt.axvline(r2[1], color='blue', linestyle='--')

    plt.show()
    exit()
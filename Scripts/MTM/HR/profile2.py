import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
import scipy.optimize as opt

def check_tol(x, y, ppmtol):
    diff = np.abs(x - y)
    diffcut = x * ppmtol * 1e-6
    if diff < diffcut:
        return True
    else:
        return False


def fancydiff(data,isodist, ppmtol = 10, sort=False):
    if sort:
        data = data[np.argsort(data[:, 0])]
        isodist = isodist[np.argsort(isodist[:, 0])]
    # Otherwise assume sorted

    matchedintensity = np.array([data[ud.nearest(data[:, 0], i), 1] for i in isodist[:, 0]])
    b1 = matchedintensity  > 0
    isodist = isodist[b1]

    iind = 0

    diffs = []
    ints = []

    isolow = isodist[iind]
    if iind < len(isodist) - 1:
        isohigh = isodist[iind + 1]
    else:
        isohigh = [isolow[0] * 2, 0]

    i1 = ud.nearest(data[:, 0], isolow[0])
    i2 = ud.nearest(data[:, 0], isohigh[0])
    norm1 = data[i1, 1]
    norm2 = data[i2, 1]

    for i, d in enumerate(data):
        if d[1] == 0:
            continue

        if iind == 0 and d[0] <= isolow[0]:
            if check_tol(d[0], isolow[0], ppmtol):
                diff = isolow[0] - d[0]
                diffs.append(diff)
                ints.append(d[1] / norm1)
            else:
                continue
        elif isolow[0] < d[0] <=isohigh[0]:
            if check_tol(d[0], isolow[0], ppmtol):
                diff = isolow[0] - d[0]
                diffs.append(diff)
                ints.append(d[1] / norm1)
            elif check_tol(d[0], isohigh[0], ppmtol):
                diff = isohigh[0] - d[0]
                diffs.append(diff)
                ints.append(d[1] / norm2)
            else:
                continue
        else:
            while isohigh[0] < d[0]:
                iind += 1
                isolow = isodist[iind]
                if iind < len(isodist) - 1:
                    isohigh = isodist[iind + 1]
                else:
                    isohigh = [isolow[0] * 2, 0]

                i1 = ud.nearest(data[:, 0], isolow[0])
                i2 = ud.nearest(data[:, 0], isohigh[0])
                norm1 = data[i1, 1]
                norm2 = data[i2, 1]

            if check_tol(d[0], isolow[0], ppmtol):
                diff = isolow[0] - d[0]
                diffs.append(diff)
                ints.append(d[1] / norm1)
            elif check_tol(d[0], isohigh[0], ppmtol):
                diff = isohigh[0] - d[0]
                diffs.append(diff)
                ints.append(d[1] / norm2)
            else:
                continue

            continue

    return np.array(diffs), np.sqrt(np.abs(np.log(np.array(ints))))

def vfit(x, m1):
    m2 = np.abs(m1)
    m1 = -np.abs(m1)

    if x < 0:
        return m1 * x
    elif x == 0:
        return 0
    else:
        return m2 * x

def vfit_array(x, m1):
    return np.array([vfit(i, m1) for i in x])

def fit_vfit(x, y):
    """
    Fit the vfit function to the data
    :param x: x data
    :param y: y data
    :return: fit parameters
    """
    # Initial guess for m1 and b
    initial_guess = [-500]

    # Fit the vfit function to the data
    params, _ = opt.curve_fit(vfit_array, x, y, p0=initial_guess)

    return params

def linear_interpolation(d1, d2, y):
    """
    Perform linear interpolation between two points.
    :param d1: First data point (x1, y1)
    :param d2: Second data point (x2, y2)
    :param y: y-coordinate where to interpolate
    :return: interpolated x value at y
    """
    x1, y1 = d1
    x2, y2 = d2
    if y1 == y2:
        return (x1 + x2) / 2
    # Calculate the slope (m) and intercept (b) of the line
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    # Calculate the x value at the given y
    x = (y - b) / m
    return x

def find_fwhm(data, index, wfactor=4):
    """
    Find the full width at half maximum (FWHM) of a peak in the data.
    :param data: Data array with m/z and intensity
    :param index: Index of the peak in the data
    :param ppmtol: PPM tolerance for window to search for FWHM
    :return: FWHM value
    """
    mz = data[index, 0]
    intensity = data[index, 1]

    if intensity == 0:
        return -1, -1, -1, -1, -1

    half_max = intensity / 2.0

    # Find left side
    left_index = index
    while left_index > 0 and data[left_index, 1] >= half_max:
        left_index -= 1

    # Find right side
    right_index = index
    while right_index < len(data) - 1 and data[right_index, 1] >= half_max:
        right_index += 1

    leftmz = linear_interpolation(data[left_index], data[left_index + 1], half_max)
    rightmz = linear_interpolation(data[right_index - 1], data[right_index], half_max)
    # print("Left MZ:", leftmz, "Right MZ:", rightmz, "Peak MZ:", mz)
    leftwidth = mz - leftmz
    rightwidth = rightmz - mz
    fwhm = leftwidth + rightwidth

    left_index = ud.nearest(data[:, 0], mz-wfactor*leftwidth)
    right_index = ud.nearest(data[:, 0], mz+wfactor*rightwidth)

    if leftwidth == 0 or rightwidth == 0:
        # If either width is zero, we cannot calculate FWHM properly
        return -1, -1, -1, -1, -1

    # If leftwidth is far apart from rightwidth, we can assume the peak is not symmetric
    # if 0.5 < leftwidth / rightwidth < 2:
    #     return -1, -1, -1, -1, -1

    return fwhm, leftwidth, rightwidth, left_index, right_index

def fast_fwhm(data, peaks, sort=False, wfactor=4):
    """
    Fast FWHM calculation for a set of peaks.
    :param data: Data array with m/z and intensity
    :param peaks: List of peak m/z values
    :return: FWHM values for each peak
    """
    if sort:
        data = data[np.argsort(data[:, 0])]
        peaks = peaks[np.argsort(peaks[:, 0])]
    # Otherwise assume sorted

    matchedindexes = np.array([ud.nearest(data[:, 0], i) for i in peaks[:, 0]])

    fwhmlist = []
    for i, p in enumerate(peaks):
        if p[1] == 0:
            fwhmlist.append([-1, -1, -1, -1, -1])
            continue

        index = matchedindexes[i]
        fwhm = find_fwhm(data, index, wfactor=wfactor)

        fwhmlist.append(fwhm)

    return np.array(fwhmlist)

def remove_peaks(data, fwhmlist):
    """
    Remove peaks from the data based on FWHM values.
    :param data: Data array with m/z and intensity
    :param fwhmlist: List of FWHM values for each peak
    :return: Data with peaks removed
    """
    b1 = np.ones(len(data), dtype=bool)
    for i, fwhm in enumerate(fwhmlist):
        if fwhm[0] == -1:
            continue
        s = int(fwhm[3])
        e = int(fwhm[4])
        b1[s:e] = False

    return data[b1]

os.chdir("C:\\Data\\Yuri\\")

data = np.loadtxt("processed.txt")

isodist = np.loadtxt("isodist.txt")

print(len(data), len(isodist))

fwhms = fast_fwhm(data, isodist, wfactor=4)
print("FWHM:", len(fwhms), "Peaks:", len(isodist))

leftoverdata = remove_peaks(data, fwhms)
plt.plot(data[:, 0], data[:, 1], label='Processed Data')
plt.plot(leftoverdata[:, 0], leftoverdata[:, 1], label='Leftover Data', linestyle='--')
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.show()

# b1 = fwhms[:, 0] > -1
#
# plt.plot(isodist[b1,0], fwhms[b1, 0], marker="o", linestyle=" ", label="FWHM")
# plt.show()

exit()
d, i = fancydiff(data, isodist)
b1 = i < 1
b2 = d < 0
b1 = b1 & b2
d = d[b1]
i = i[b1]

print("Diffs:", len(d), "Ints:", len(i))
fit = np.polyfit(i, d, 1)
print("Fit:", fit)

plt.plot(i, d, marker="o", linestyle=" ", alpha=0.01)

x = np.sort(i)
plt.plot(x, x*fit[0] + fit[1], label="vfit", color='r')
plt.show()

#
# plt.plot(data[:, 0], data[:, 1], label='Processed Data')
# plt.plot(isodist[:, 0], isodist[:, 1], label='Isodist', linestyle='--')
# plt.xlabel('m/z')
# plt.show()
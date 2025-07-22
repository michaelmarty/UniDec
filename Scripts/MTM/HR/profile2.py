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

        if diff > 3:
            print(diff, i, d, isolow, isohigh)

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

os.chdir("C:\\Data\\Yuri\\")

data = np.loadtxt("processed.txt")

isodist = np.loadtxt("isodist.txt")

print(len(data), len(isodist))

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
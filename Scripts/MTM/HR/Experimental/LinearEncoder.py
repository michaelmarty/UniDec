import numpy as np
import matplotlib.pyplot as plt
from unidec.IsoDec.plots import *
# import torch
from numba import njit
import time

# Example data
example = np.array([[5.66785531e+02, 1.47770838e+06],
                    [5.67057354e+02, 1.54980838e+06],
                    [5.67507468e+02, 5.21600520e+07],
                    [5.67708173e+02, 8.35557760e+07],
                    [5.67908401e+02, 7.28264240e+07],
                    [5.68060254e+02, 1.87337225e+06],
                    [5.68108674e+02, 4.35435520e+07],
                    [5.68239256e+02, 3.88155375e+06],
                    [5.68309390e+02, 2.05468060e+07],
                    [5.68509951e+02, 7.18109250e+06],
                    [5.68707871e+02, 2.30373500e+06],
                    [5.69150563e+02, 1.57598062e+06],
                    [5.69243121e+02, 1.96390440e+07],
                    [5.69334393e+02, 6.82677120e+07],
                    [5.69425337e+02, 1.22867432e+08],
                    [5.69516492e+02, 1.45702336e+08],
                    [5.69607541e+02, 1.20801936e+08],
                    [5.69698595e+02, 1.06786072e+08],
                    [5.69789906e+02, 6.56232960e+07],
                    [5.69881208e+02, 3.41013880e+07],
                    [5.69972168e+02, 1.70930360e+07],
                    [5.70063432e+02, 9.17621100e+06],
                    [5.70699369e+02, 1.96462650e+06]])



@njit(fastmath=True)
def linear_encoder(data, maxz=50, mass_diff_c=1.0033, nq=2):
    deltax = mass_diff_c / (maxz**2 * nq)
    mzlen = data[-1, 0] - data[0, 0]
    xmax = mass_diff_c * mzlen
    xaxis = np.arange(0, xmax, deltax)
    lx = len(xaxis)

    xdata = data[:, 0]
    xdata -= xdata[0]
    xdata /= mass_diff_c

    factor = lx / xmax
    xindexes = np.floor(xdata * factor)

    outdata = np.zeros((lx, 2))
    outdata[:, 0] = xaxis
    for i, d in enumerate(data):
        index = xindexes[i]
        outdata[int(index), 1] += d[1]

    return outdata


def matrix_expander(evec, maxz=50, nq=2):
    lv = len(evec)
    output = np.zeros((len(evec), maxz))
    output[:, 0] = evec[:, 1]
    kernel = np.zeros((lv, maxz))
    for z in range(maxz):

        spacing = round(maxz**2 * nq / (z+1))
        kernel[0, z] = 1/5.
        kernel[spacing, z] = 1/5.
        kernel[lv - spacing, z] = 1/5.
        kernel[spacing * 2, z] = 1/5.
        kernel[lv-2*spacing, z] = 1/5.

    ftk = np.fft.rfft2(kernel)
    fto = np.fft.rfft2(output)
    output = np.fft.irfft2(ftk * fto)

    return output


evec = linear_encoder(example)

emat = matrix_expander(evec)

plt.imshow(emat)
plt.show()

exit()


file = "C:\\Data\\IsoNN\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.txt"
data = np.loadtxt(file)

evec = linear_encoder(data)

emat = matrix_expander(evec)

plt.imshow(emat, aspect=1/500000)
plt.show()
exit()
st = time.perf_counter()
edat = linear_encoder(data)
print("Encoding Time:", time.perf_counter() - st, len(edat))

# cplot(example)
cplot(edat, color="b")
plt.show()

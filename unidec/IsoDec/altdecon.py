import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
import matplotlib as mpl

mpl.use('WxAgg')

example = np.array(
    [
        [5.66785531e02, 1.47770838e06],
        [5.67057354e02, 1.54980838e06],
        [5.67507468e02, 5.21600520e07],
        [5.67708173e02, 8.35557760e07],
        [5.67908401e02, 7.28264240e07],
        [5.68060254e02, 1.87337225e06],
        [5.68108674e02, 4.35435520e07],
        [5.68239256e02, 3.88155375e06],
        [5.68309390e02, 2.05468060e07],
        [5.68509951e02, 7.18109250e06],
        [5.68707871e02, 2.30373500e06],
        [5.69150563e02, 1.57598062e06],
        [5.69243121e02, 1.96390440e07],
        [5.69334393e02, 6.82677120e07],
        [5.69425337e02, 1.22867432e08],
        [5.69516492e02, 1.45702336e08],
        [5.69607541e02, 1.20801936e08],
        [5.69698595e02, 1.06786072e08],
        [5.69789906e02, 6.56232960e07],
        [5.69881208e02, 3.41013880e07],
        [5.69972168e02, 1.70930360e07],
        [5.70063432e02, 9.17621100e06],
        [5.70699369e02, 1.96462650e06],
    ]
)

import numpy as np
from scipy import signal

MAX_CHARGE = 50


def gen_thrash_arrays(centroids, startpad=10):
    max_mz = centroids[-1, 0]
    min_mz = centroids[0, 0]

    additional = 4 - (max_mz - min_mz)
    if additional < 0:
        additional = 0
    else:
        max_mz += additional / 2
        min_mz -= additional / 2

    num_l = int((max_mz - min_mz) * MAX_CHARGE * 8)

    intx = np.linspace(min_mz, max_mz, num_l)
    if len(intx) <= 1:
        return None, None, None, None, None, None

    linear_data = ud.lintegrate(centroids, intx, fastmode=True)
    if len(linear_data) <= 1:
        return None, None, None, None, None, None

    corry = signal.fftconvolve(linear_data[:, 1], linear_data[:, 1][::-1], mode='same')
    maxposition = np.argmax(corry)
    ac = corry[maxposition + startpad:]


    acx = intx - intx[0]
    acx = 1 / acx[startpad:len(ac) + startpad]

    fft = np.fft.fft(linear_data[:, 1])
    fft = np.abs(fft)
    fft = fft[startpad:len(ac) + startpad]

    fftx = np.fft.fftfreq(len(linear_data[:, 1]), d=(linear_data[1, 0] - linear_data[0, 0]))
    fftx = fftx[startpad:len(ac) + startpad]
    b1 = fftx < MAX_CHARGE
    fft = fft[b1]
    fftx = fftx[b1]

    fft /= np.max(fft)
    ac2 = ud.lintegrate(np.transpose([acx, ac]), fftx)


    acx, ac = ac2[:, 0], ac2[:, 1]
    ac /= np.max(ac)
    mul = ac * fft
    return linear_data, ac, fft, mul, fftx, acx

def thrash_predict(centroids):
    linear_data, ac, fft, mul, fftx, acx = gen_thrash_arrays(centroids)
    if linear_data is None or len(mul) == 0:
        return 0
    maxpos = np.argmax(mul)

    charge = acx[maxpos]
    return round(charge)


if __name__ == "__main__":
    example = np.loadtxt("Z:\Group Share\JGP\MockData\mockdata_centroids.csv", delimiter=",")
    example = ud.datachop(example, 326,327)
    for i in range(len(example)):
        plt.plot([example[i, 0], example[i, 0]], [0, example[i, 1]], color="black")
    plt.show()



    raw_data = example
    charge_state = thrash_predict(raw_data)
    print(charge_state)

    linear_data, ac, fft, mul, fftx, acx = gen_thrash_arrays(raw_data)

    plt.subplot(121)
    plt.plot(linear_data[:, 0], linear_data[:, 1])
    plt.subplot(122)
    plt.plot(acx, mul, label="Mul")
    plt.plot(acx, ac-1, label="AC")
    plt.plot(fftx, fft-2, label="FFT")
    plt.legend()
    plt.show()


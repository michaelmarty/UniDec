import os
import matplotlib.pyplot as plt
import numpy as np

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

os.chdir("Z:\\Group Share\\Dash")

file = "UIB data.txt"

data = np.loadtxt(file)
data[:, 1] = data[:, 1] - np.min(data[:, 1])  # Subtract the minimum to avoid negative values

pdata = nonlinearize(data[6:], 50)

rdata = pdata[::2]
gdata = pdata[1::2]
print("Rsum:", np.sum(rdata[:, 1]), "Gsum:", np.sum(gdata[:, 1]))

plt.plot(rdata[:, 0], rdata[:, 1], marker="o", color="r")
plt.plot(gdata[:, 0], gdata[:, 1], marker="o", color="g")
plt.show()
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from unidec.IsoDec.match import create_isodist

try:
    mpl.use("WxAgg")
except:
    pass
def cplot(centroids, color='r', factor=1, base=0, mask=None, mfactor=-1, mcolor="g", z=0, zcolor="b", zfactor=1):
    """
    Simple script to plot centroids
    :param centroids: Centroid array with m/z in first column and intensity in second
    :param color: Color
    :param factor: Mutiplicative factor for intensity. -1 will set below the axis
    :param base: Base of the lines. Default is 0. Can be adjusted to shift up or down.
    :return: None
    """
    #if centroids is not None and len(centroids) > 0:
    #    plt.hlines(0, np.amin(centroids[:, 0]), np.amax(centroids[:, 0]), color="k")

    if mask is not None:
        if len(centroids) > len(mask):
            mask = np.append(mask, np.zeros(len(centroids) - len(mask)))
        else:
            mask = mask[:len(centroids)]
        for c in centroids[mask.astype(bool)]:
            plt.vlines(c[0], base, base + mfactor * c[1], color=mcolor)

    if z != 0:
        isodist = create_isodist(centroids[np.argmax(centroids[:, 1]), 0], z, centroids)
        for c in isodist:
            plt.vlines(c[0], base, base + zfactor * c[1], color=zcolor, linewidth=3)

    if centroids is not None:
        for c in centroids:
            plt.vlines(c[0], base, base + factor * c[1], color=color)


def plot_pks(pks, data=None, centroids=None, scan=-1, show=False, labelz=True, title=None):
    plt.subplot(121)
    if data is not None:
        plt.plot(data[:, 0], data[:, 1])

    if centroids is not None:
        cplot(centroids)

    for p in pks.peaks:
        if scan == -1 or p.scan == scan:
            color = p.color
            isodist = p.isodist
            plt.subplot(121)
            cplot(isodist, color=color, factor=-1)
            centroids = p.centroids
            peakmz = p.mz
            cplot(centroids)
            plt.subplot(122)
            massdist = p.massdist
            cplot(massdist, color=color)
            mass = p.avgmass

            if labelz:
                try:
                    plt.text(mass, np.amax(centroids[:, 1]) * 1.05, str(p.z), color=color)
                except:
                    try:
                        plt.text(mass, np.amax(isodist[:, 1]) * 1.05, str(p.z), color=color)
                    except:
                        pass

    if title is not None:
        plt.suptitle(title)
    if show:
        plt.show()
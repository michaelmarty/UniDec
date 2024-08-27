import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from unidec.IsoDec.match import create_isodist
import time

try:
    mpl.use("WxAgg")
except:
    pass

def fast_vlines(centroids, color, base, factor):
    """
    Fast vlines plotter
    :param centroids: Centroids array
    :param color: Color
    :param base: Base of the lines
    :param factor: Factor to multiply intensity by
    :return: None
    """
    xpairs = np.transpose([centroids[:, 0], centroids[:, 0]])
    ypairs = np.transpose([base + np.zeros(len(centroids)), base + factor * centroids[:, 1]])
    xlist = []
    ylist = []
    for xends, yends in zip(xpairs, ypairs):
        xlist.extend(xends)
        xlist.append(None)
        ylist.extend(yends)
        ylist.append(None)
    plt.plot(xlist, ylist, color=color)

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
        fast_vlines(centroids[mask.astype(bool)], mcolor, base, mfactor)
        #for c in centroids[mask.astype(bool)]:
        #    plt.vlines(c[0], base, base + mfactor * c[1], color=mcolor)

    if z != 0:
        isodist = create_isodist(centroids[np.argmax(centroids[:, 1]), 0], z, centroids)
        fast_vlines(isodist, zcolor, base, zfactor)
        #for c in isodist:
        #    plt.vlines(c[0], base, base + zfactor * c[1], color=zcolor, linewidth=3)

    if centroids is not None:
        fast_vlines(centroids, color, base, factor)
        #for c in centroids:
        #    plt.vlines(c[0], base, base + factor * c[1], color=color)


def plot_pks(pks, data=None, centroids=None, scan=-1, show=False, labelz=True, title=None, ccolor="r", plotmass=True,
             zcolor=False, zcolormap="nipy_spectral", forcecolor=None, nocentroids=False):

    if plotmass:
        plt.subplot(121)
    if data is not None:
        plt.plot(data[:, 0], data[:, 1])

    if centroids is not None:
        cplot(centroids, color=ccolor)

    # Turn off top and right axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")

    if zcolor:
        cmap = plt.get_cmap(zcolormap)

    print(len(pks.peaks))
    for p in pks.peaks:
        if scan == -1 or p.scan == scan:
            if forcecolor is not None:
                color = forcecolor
            elif zcolor:
                rescale = ((p.z-1)**0.5)/(50**0.5)
                color = cmap(rescale)
            else:
                color = p.color
            isodist = p.isodist
            if plotmass:
                plt.subplot(121)
            cplot(isodist, color=color, factor=-1)
            #plt.text(p.mz, np.amax(isodist[:, 1]) * 1.05, (str(p.peakmass) + str(p.mz) + str(p.z)), color=color)
            if not nocentroids:
                cplot(p.centroids, color=ccolor)
            if plotmass:
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


def plot_zdist(eng):
    zdata = [x[1] for x in eng.training_data]
    plt.hist(zdata, bins=range(0, 50))
    plt.show()


if __name__ == "__main__":
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
    starttime = time.perf_counter()
    for i in range(50):
        cplot(example)
    print("Time to plot 50", time.perf_counter()-starttime)
    plt.show()
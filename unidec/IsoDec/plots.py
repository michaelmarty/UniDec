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

# Global variable to control dragging
zoom_levels = []
start_x = None
start_y = None


def on_scroll(event):
    ax = plt.gca()

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if event.inaxes is not None:
        mouse_x = event.xdata
        mouse_y = event.ydata

        zoom_factor = 0.7
        if event.button == 'up':
            zoom_levels.append((xlim, ylim))

            new_x_range = (xlim[1] - xlim[0]) * (1 - zoom_factor)
            new_y_range = (ylim[1] - ylim[0]) * (1 - zoom_factor)
            ax.set_xlim([mouse_x - new_x_range / 2, mouse_x + new_x_range / 2])
            ax.set_ylim([mouse_y - new_y_range / 2, mouse_y + new_y_range / 2])

        elif event.button == 'down':
            if zoom_levels:
                ax.set_xlim(zoom_levels[-1][0])
                ax.set_ylim(zoom_levels[-1][1])
                zoom_levels.pop()

        plt.draw()


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


def cplot(centroids, color='r', factor=1, base=0, mask=None, mfactor=-1, mcolor="g", z=0, zcolor="b", zfactor=1, isodist=None):
    """
    Simple script to plot centroids
    :param centroids: Centroid array with m/z in first column and intensity in second
    :param color: Color
    :param factor: Mutiplicative factor for intensity. -1 will set below the axis
    :param base: Base of the lines. Default is 0. Can be adjusted to shift up or down.
    :return: None
    """
    # if centroids is not None and len(centroids) > 0:
    #    plt.hlines(0, np.amin(centroids[:, 0]), np.amax(centroids[:, 0]), color="k")

    if mask is not None:
        if len(centroids) > len(mask):
            mask = np.append(mask, np.zeros(len(centroids) - len(mask)))
        else:
            mask = mask[:len(centroids)]
        fast_vlines(centroids[mask.astype(bool)], mcolor, base, mfactor)

    if z != 0:
        isodist = create_isodist(centroids[np.argmax(centroids[:, 1]), 0], z, centroids)
        fast_vlines(isodist, zcolor, base, zfactor)

    if isodist is not None:
        fast_vlines(isodist, zcolor, base, zfactor)

    if centroids is not None:
        fast_vlines(centroids, color, base, factor)



def plot_pks(pks, data=None, centroids=None, scan=-1, show=False, labelz=False, title=None, ccolor="r", plotmass=True,
             zcolor=False, zcolormap="nipy_spectral", forcecolor=None, nocentroids=False, tickfont=12, labfont=14, labelpeaks=False):  # Added matched and actual parameters

    if plotmass:
        plt.subplot(121)
    if data is not None:
        plt.plot(data[:, 0], data[:, 1])

    if centroids is not None:
        cplot(centroids, color=ccolor)

    # Turn off top and right axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel("m/z", fontsize=labfont)
    plt.ylabel("Intensity", fontsize=labfont)
    plt.xticks(fontsize=tickfont)
    plt.yticks(fontsize=tickfont)

    if zcolor:
        cmap = plt.get_cmap(zcolormap)

    # Plotting peaks from pks
    if pks is not None:
        for p in pks.peaks:
            if scan == -1 or p.scan == scan:
                if forcecolor is not None:
                    color = forcecolor
                elif zcolor:
                    rescale = ((p.z - 1) ** 0.5) / (50 ** 0.5)
                    color = cmap(rescale)
                else:
                    color = p.color
                if labelpeaks:
                    intensity = 1000
                    if p.isodist is not None:
                        intensity = np.amax(p.isodist[:, 1])
                    plt.text(p.mz + 0.05, intensity, s=str(round(p.mz, 2)) + "," + str(p.z), fontsize=8, color=color)
                isodist = p.isodist
                if plotmass:
                    plt.subplot(121)
                cplot(isodist, color=color, factor=-1)
                if not nocentroids:
                    cplot(p.centroids, color=ccolor)
                if plotmass:
                    plt.subplot(122)
                    massdist = p.massdist
                    plt.gca().spines['top'].set_visible(False)
                    plt.gca().spines['right'].set_visible(False)
                    cplot(massdist, color=color)
                    mass = p.avgmass

    if title is not None:
        plt.suptitle(title)

    #plt.legend()  # Optional: Add a legend

    plt.connect('scroll_event', on_scroll)


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
    print("Time to plot 50", time.perf_counter() - starttime)
    plt.show()

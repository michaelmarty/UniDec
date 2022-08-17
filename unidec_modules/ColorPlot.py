import math

import numpy as np
import matplotlib.cm as cm
from scipy.stats import norm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar

from unidec_modules.PlottingWindow import PlottingWindow
from unidec_modules.unidectools import color_map_array

'''
# These functions are getting pulled from unidectools. You can uncomment this to avoid having to import them.

def make_alpha_cmap(rgb_tuple, alpha):
    """
    Make color map where RGB specified in tup
    :param rgb_tuple: Tuple of RGB vaalues [0,1]
    :param alpha: Maximum Alpha (Transparency) [0,1]
    :return: Color map dictionary with a constant color but varying transprency
    """
    cdict = {'red': ((0.0, rgb_tuple[0], rgb_tuple[0]),
                     (1.0, rgb_tuple[0], rgb_tuple[0])), 'green': ((0.0, rgb_tuple[1], rgb_tuple[1]),
                                                                   (1.0, rgb_tuple[1], rgb_tuple[1])),
             'blue': ((0.0, rgb_tuple[2], rgb_tuple[2]),
                      (1.0, rgb_tuple[2], rgb_tuple[2])), 'alpha': ((0.0, 0, 0),
                                                                    (1.0, alpha, alpha))}
    return cdict


def color_map_array(array, cmap, alpha):
    """
    For a specified array of values, map the intensity to a specified RGB color defined by cmap (output as topcm).
    For each color, create a color map where the color is constant but the transparency changes (output as cmarr).
    :param array: Values
    :param cmap: Color map
    :param alpha: Max Alpha value (transparency) [0,1]
    :return: cmarr, topcm (list of transparent color maps, list of RGB colors for array values)
    """
    rtab = array
    topcm = cm.ScalarMappable(cmap=cmap).to_rgba(rtab)[:, :3]
    cmarr = []
    for i in range(0, len(rtab)):
        cmarr.append(make_alpha_cmap(topcm[i], alpha))
    return cmarr, topcm

'''

__author__ = 'Michael.Marty'


class ColorPlot2D(PlottingWindow):
    """
    Method to perform a 3D plot by plotting the Z axis as 2D slices of a different color.
    Each 2D slice uses the alpha (trasparency) parameter to indicate intensity in square root scale.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize by inheriting from Plotting Window and setting the axes dimensions.
        :param args:
        :param kwargs:
        :return: ColorPlot2D object
        """
        PlottingWindow.__init__(self, *args, **kwargs)
        self._axes = [0.1, 0.1, 0.64, 0.8]

    def make_color_plot(self, mzgrid, mzax, dtax, ztab):
        """
        Creates the color plot
        :param mzgrid: mz x dt x z numpy array of intensity values
        :param mzax: m/z axis
        :param dtax: drift time axis
        :param ztab: charge axis
        :return: None
        """
        # Create Data Limits
        data_x_lim = (np.amin(mzax), np.amax(mzax))
        data_y_lim = (np.amin(dtax), np.amax(dtax))
        self.datalims = [data_x_lim[0], data_y_lim[0], data_x_lim[1], data_y_lim[1]]

        # Normalize Grid
        mzgrid = np.abs(mzgrid) / np.amax(mzgrid)

        # This is the magic for creating the color maps
        # Note: It skews the color map such that the maximum color change is when the charge state distribution is changing the most
        # You could make it linear, but a lot of the color change would happen at the boring fringes of the distribution
        zlen = len(ztab)
        ztot = np.sum(mzgrid, axis=(0, 1))
        zind = np.array(list(range(0, zlen)))
        ztab = ztab
        avg = np.average(zind, weights=ztot)
        std = math.sqrt(np.average(np.array(zind - avg) ** 2, weights=ztot))
        skew = norm.cdf(zind, loc=avg, scale=std * 1.75)
        topcmap = "gist_rainbow"
        cmarr, acolors = color_map_array(skew, topcmap, 1)

        # Create basics of the plot
        self.figure.clear()
        self.subplot1 = self.figure.add_axes(self._axes)

        # I think this clips things such that everything above 0.5 is a full strength of color
        normalization = cm.colors.Normalize(vmax=0.5, vmin=0)

        # Set up plot limits
        xdiff = mzax[1] - mzax[0]
        ydiff = dtax[1] - dtax[0]
        extent = (min(mzax) - 0.5 * xdiff, max(mzax) + 0.5 * xdiff, min(dtax) - 0.5 * ydiff, max(dtax) + 0.5 * ydiff)

        # Make the black background and title
        self.subplot1.imshow(np.transpose(np.ones_like(mzgrid[:, :, 0])), origin="lower", cmap="binary",
                             extent=extent, norm=normalization, aspect='auto')
        self.subplot1.set_title("Deconvolution in m/z vs. Arrival Time")

        # Loop through each charge state and make the color plot to layer on top of the black background
        for i in range(0, zlen):
            cm.register_cmap(name="newcmap", data=cmarr[i])
            grid_slice = np.sqrt(mzgrid[:, :, i])
            normalization = cm.colors.Normalize(vmax=0.5, vmin=0.00)
            self.subplot1.imshow(np.transpose(grid_slice), origin="lower", cmap="newcmap", extent=extent,
                                 norm=normalization,
                                 aspect='auto')

        # Labels and legends
        self.subplot1.set_xlabel("m/z (Th)")
        self.subplot1.set_ylabel("Arrival Time (ms)")
        cax = self.figure.add_axes([0.77, 0.1, 0.04, 0.8])
        cmap = colors.ListedColormap(acolors)
        normalization = colors.BoundaryNorm(ztab - 0.5, cmap.N)
        self.cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=normalization, orientation="vertical",
                                          ticks=ztab, drawedges="True")
        self.cbar.set_label("Charge")
        # Set tick marks to white
        for line in self.subplot1.xaxis.get_ticklines():
            line.set_color("white")
        for line in self.subplot1.yaxis.get_ticklines():
            line.set_color("white")

        # Set up zoom and repaint it
        self.setup_zoom([self.subplot1], 'box', data_lims=self.datalims)
        self.repaint()
        self.flag = True

    def make_compare_color_plot(self, data1, data2, test_kda=True, xlab="", ylab=""):
        # If data is coming in as 1D list in dat, reshape it

        # Test if we should plot kDa instead of Da
        if test_kda:
            self.kda_test(np.unique(data1[:, 0]))

        z1 = data1[:, 2]
        x1 = np.unique(data1[:, 0])/ self.kdnorm
        y1 = np.unique(data1[:, 1])
        z2 = data2[:, 2]
        x2 = np.unique(data2[:, 0])/ self.kdnorm
        y2 = np.unique(data2[:, 1])

        xlen1 = len(x1)
        ylen1 = len(y1)
        zgrid1 = np.reshape(z1, (xlen1, ylen1))
        xlen2 = len(x2)
        ylen2 = len(y2)
        zgrid2 = np.reshape(z2, (xlen2, ylen2))

        # Create Data Limits
        data_x_lim = (np.amin(np.concatenate([x1, x2])), np.amax(np.concatenate([x1, x2])))
        data_y_lim = (np.amin(np.concatenate([y1, y2])), np.amax(np.concatenate([y1, y2])))
        self.datalims = [data_x_lim[0], data_y_lim[0], data_x_lim[1], data_y_lim[1]]
        # Normalize Grids
        zgrid1 = np.abs(zgrid1) / np.amax(zgrid1)
        zgrid2 = np.abs(zgrid2) / np.amax(zgrid2)

        # This is the magic for creating the color maps
        # Note: It skews the color map such that the maximum color change is when the charge state distribution is changing the most
        # You could make it linear, but a lot of the color change would happen at the boring fringes of the distribution
        #zlen = len(ztab)
        #ztot = np.sum(mzgrid, axis=(0, 1))
        #zind = np.array(list(range(0, zlen)))

        #avg = np.average(zind, weights=ztot)
        #std = math.sqrt(np.average(np.array(zind - avg) ** 2, weights=ztot))
        #skew = norm.cdf(zind, loc=avg, scale=std * 1.75)
        topcmap = "seismic"
        cmarr, acolors = color_map_array([0,1], topcmap, 1)
        cm.register_cmap(name="newcmap1", data=cmarr[0])
        cm.register_cmap(name="newcmap2", data=cmarr[1])

        # Create basics of the plot
        self.figure.clear()
        self.subplot1 = self.figure.add_axes(self._axes)

        # I think this clips things such that everything above vmax is a full strength of color
        normalization = cm.colors.Normalize(vmax=1, vmin=0)

        # Set up plot limits
        xdiff = x1[1] - x1[0]
        ydiff = y1[1] - y1[0]
        extent = (data_x_lim[0] - 0.5 * xdiff, data_x_lim[1] + 0.5 * xdiff, data_y_lim[0] - 0.5 * ydiff, data_y_lim[1] + 0.5 * ydiff)
        extent2 = (np.amin(x2), np.amax(x2), np.amin(y2),np.amax(y2))
        extent1 = (np.amin(x1), np.amax(x1), np.amin(y1), np.amax(y1))

        # Make the background and title
        self.subplot1.imshow(np.transpose(np.zeros_like(zgrid1)), origin="lower", cmap="binary",
                             extent=extent, norm=normalization, aspect='auto')
        self.subplot1.set_title("Comparison")

        # Plot each comparison
        self.subplot1.imshow(np.transpose(zgrid1), origin="lower", cmap="newcmap1", extent=extent1,
                             norm=normalization,aspect='auto')
        self.subplot1.imshow(np.transpose(zgrid2), origin="lower", cmap="newcmap2", extent=extent2,
                             norm=normalization,aspect='auto')

        # Labels and legends
        self.subplot1.set_xlabel(self.xlabel)
        self.subplot1.set_ylabel(ylab)

        # Set up zoom and repaint it
        self.setup_zoom([self.subplot1], 'box', data_lims=self.datalims)
        self.repaint()
        self.flag = True


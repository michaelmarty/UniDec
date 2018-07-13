import math

import numpy as np
import matplotlib.cm as cm
from scipy.stats import norm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar

from unidec_modules.PlottingWindow import PlottingWindow
from unidec_modules.unidectools import color_map_array

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
        data_x_lim = (np.amin(mzax), np.amax(mzax))
        data_y_lim = (np.amin(dtax), np.amax(dtax))
        self.datalims = [data_x_lim[0], data_y_lim[0], data_x_lim[1], data_y_lim[1]]
        mzgrid = np.abs(mzgrid) / np.amax(mzgrid)

        zlen = len(ztab)
        ztot = np.sum(mzgrid, axis=(0, 1))
        zind = np.array(list(range(0, zlen)))
        ztab = ztab
        avg = np.average(zind, weights=ztot)
        std = math.sqrt(np.average(np.array(zind - avg) ** 2, weights=ztot))
        skew = norm.cdf(zind, loc=avg, scale=std * 1.75)
        topcmap = "gist_rainbow"
        cmarr, acolors = color_map_array(skew, topcmap, 1)

        self.figure.clear()
        self.subplot1 = self.figure.add_axes(self._axes)

        normalization = cm.colors.Normalize(vmax=0.5, vmin=0)

        xdiff = mzax[1] - mzax[0]
        ydiff = dtax[1] - dtax[0]
        extent = (min(mzax) - 0.5 * xdiff, max(mzax) + 0.5 * xdiff, min(dtax) - 0.5 * ydiff, max(dtax) + 0.5 * ydiff)

        self.subplot1.imshow(np.transpose(np.ones_like(mzgrid[:, :, 0])), origin="lower", cmap="binary",
                             extent=extent, norm=normalization, aspect='auto')
        self.subplot1.set_title("Deconvolution in m/z vs. Arrival Time")

        for i in range(0, zlen):
            cm.register_cmap(name="newcmap", data=cmarr[i])
            grid_slice = np.sqrt(mzgrid[:, :, i])
            normalization = cm.colors.Normalize(vmax=0.5, vmin=0.00)
            self.subplot1.imshow(np.transpose(grid_slice), origin="lower", cmap="newcmap", extent=extent,
                                 norm=normalization,
                                 aspect='auto')
        self.subplot1.set_xlabel("m/z (Th)")
        self.subplot1.set_ylabel("Arrival Time (ms)")
        cax = self.figure.add_axes([0.77, 0.1, 0.04, 0.8])
        cmap = colors.ListedColormap(acolors)
        normalization = colors.BoundaryNorm(ztab - 0.5, cmap.N)
        self.cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=normalization, orientation="vertical",
                                          ticks=ztab, drawedges="True")
        self.cbar.set_label("Charge")
        for line in self.subplot1.xaxis.get_ticklines():
            line.set_color("white")
        for line in self.subplot1.yaxis.get_ticklines():
            line.set_color("white")
        self.setup_zoom([self.subplot1], 'box', data_lims=self.datalims)
        self.repaint()
        self.flag = True

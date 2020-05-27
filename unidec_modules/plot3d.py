import numpy as np
from copy import deepcopy
from matplotlib.ticker import MaxNLocator

from unidec_modules.PlottingWindow import PlottingWindow
from unidec_modules.unidectools import interp_pos

__author__ = 'Michael.Marty'


class CubePlot(PlottingWindow):
    """
    Plotting class for cube plots. Performs isometric projections of three 2D plots onto three faces of a cube.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialize the plotting window. Specify axes. Specify that the plot does not have normal tick marks to repaint.
        :param args:
        :param kwargs:
        :return: CubePlot Object
        """
        PlottingWindow.__init__(self, *args, **kwargs)
        self._axes = [0.01, 0.01, 0.99, 0.99]
        self.normaltickes = False

    def cubeplot(self, xaxis, yaxis, zaxis, C, C2, C3, xlab="", ylab="", zlab="", cmap="jet"):
        """
        Make 2D isometric projection plot of three cube faces
        :param xaxis: X-axis values
        :param yaxis: Y-axis values
        :param zaxis: Z-axis values
        :param C: Intensity grid for X v. Y face
        :param C2: Intensity grid for X v. Z face
        :param C3: Intensity grid for Y v. Z face
        :param xlab: X-axis label
        :param ylab: Y-axis label
        :param zlab: Z-axis label
        :param cmap: Colormap for 2D plot
        :return: None
        """
        self.clear_plot("nopaint")
        self.subplot1 = self.figure.add_axes(self._axes)

        self.alpha = 35.264
        self.beta = 45
        nticks = 5
        numcountours = 50
        self.cmap = cmap
        self.set_tickcolor()


        self.make_isomatrices()

        self.xlen = len(xaxis)
        self.ylen = len(yaxis)
        self.zlen = len(zaxis)

        X, Y, X2, Y2, X3, Y3 = self.isogrids()

        self.subplot1.contourf(X, Y, C, numcountours, cmap=self.cmap)
        self.subplot1.contourf(X2, Y2, C2, numcountours, cmap=self.cmap)
        cax = self.subplot1.contourf(X3, Y3, C3, numcountours, cmap=self.cmap)

        self.isolines()
        self.subplot1.set_aspect(1)

        self.subplot1.set_ylim([np.amin(Y) - 0.1, np.amax(Y2) + 0.01])
        self.subplot1.set_xlim([np.amin(X) - 0.01, np.amax(X2) + 0.01])
        self.subplot1.set_xticks([])
        self.subplot1.set_yticks([])
        self.subplot1.axison = False

        self.xticlab = MaxNLocator(nbins=nticks).tick_values(np.amin(xaxis), np.amax(xaxis))
        self.yticlab = MaxNLocator(nbins=nticks).tick_values(np.amin(yaxis), np.amax(yaxis))
        self.zticlab = MaxNLocator(nbins=nticks).tick_values(np.amin(zaxis), np.amax(zaxis))

        self.xticlab = self.xticlab[
            np.logical_and(self.xticlab > np.amin(xaxis), self.xticlab < np.amax(xaxis))]
        self.yticlab = self.yticlab[
            np.logical_and(self.yticlab > np.amin(yaxis), self.yticlab < np.amax(yaxis))]
        self.zticlab = self.zticlab[
            np.logical_and(self.zticlab > np.amin(zaxis), self.zticlab < np.amax(zaxis))]

        self.xticloc = np.array([interp_pos(xaxis, i) for i in self.xticlab])
        self.yticloc = np.array([interp_pos(yaxis, i) for i in self.yticlab])
        self.zticloc = np.array([interp_pos(zaxis, i) for i in self.zticlab])

        self.isoticks(self.subplot1, self.tickcolor)
        self.isolabels(self.xticlab, self.yticlab, self.zticlab, self.subplot1, self.tickcolor)
        self.subplot1.text(0.5, 0.05, zlab, horizontalalignment="left")
        self.subplot1.text(-0.6, 0.05, xlab, horizontalalignment="right")
        self.subplot1.text(-0.9, 0.6, ylab, verticalalignment="bottom", rotation='vertical')

        self.cbar = self.figure.colorbar(cax, ax=self.subplot1, shrink=0.5, use_gridspec=True,
                                         ticks=[0, np.amax(C3) / 2, np.amax(C3)])
        self.cbar.ax.get_yaxis().set_tick_params(direction='out')
        self.cbar.ax.set_yticklabels(["0", "%", "100"])

        self.datalims = [np.amin(X) - 0.01, np.amin(Y) - 0.1, np.amax(X2) + 0.01, np.amax(Y2) + 0.01]
        self.setup_zoom([self.subplot1], 'box', data_lims=self.datalims)

        try:
            self.repaint()
        except MemoryError:
            print("Memory Error: Not updating 2D plot")
        self.flag = True

    def make_isomatrices(self):
        """
        Make projection matrix, self.pc, for calculating isometric projection
        :return: None
        """
        alpha = self.alpha * np.pi / 180.
        beta = self.beta * np.pi / 180.
        a = np.array([[1., 0., 0.], [0., np.cos(alpha), np.sin(alpha)], [0., -np.sin(alpha), np.cos(alpha)]])
        b = np.array([[np.cos(beta), 0., -np.sin(beta)], [0., 1., 0.], [np.sin(beta), 0., np.cos(beta)]])
        p = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 0.]])
        c = np.dot(a, b)
        self.pc = np.dot(p, c)

    def isometric_projection(self, x, y, z):
        """
        Give the isometric projection of a given (x,y,z) onto a plane set by self.alpha and self.beta
        :param x: x value (float)
        :param y: y value (float)
        :param z: z value (float)
        :return: o, isometric projection vector
        """
        vec = np.array([[x], [y], [z]])
        o = np.dot(self.pc, vec)
        return o

    def isogrids(self):
        """
        Sets jsometric projections of X, Y, and Z axis grids
        :return: None
        """
        xlen = self.xlen
        ylen = self.ylen
        zlen = self.zlen
        x_range = np.arange(0., xlen) / float(xlen - 1)
        y_range = np.arange(0., ylen) / float(ylen - 1)
        z_range = np.arange(0., zlen) / float(zlen - 1)

        iso = np.array(
            [[self.isometric_projection(0, y_range[j], x_range[xlen - i - 1]) for j in range(0, ylen)] for i in
             range(0, xlen)])
        xarr = iso[:, :, 0, 0]
        yarr = iso[:, :, 1, 0]

        iso = np.array(
            [[self.isometric_projection(z_range[j], 1, x_range[xlen - 1 - i]) for j in range(0, zlen)] for i in
             range(0, xlen)])
        xarr2 = iso[:, :, 0, 0]
        yarr2 = iso[:, :, 1, 0]

        iso = np.array(
            [[self.isometric_projection(z_range[j], y_range[i], 0) for j in range(0, zlen)] for i in range(0, ylen)])
        xarr3 = iso[:, :, 0, 0]
        yarr3 = iso[:, :, 1, 0]
        return xarr, yarr, xarr2, yarr2, xarr3, yarr3

    def isoticks(self, ax, col):
        """
        Draws isometric projection of tick marks.
        :param ax: Axis object
        :param col: Color of ticks
        :return: None
        """
        xlen = self.xlen
        ylen = self.ylen
        zlen = self.zlen
        xtics = deepcopy(self.xticloc)
        ytics = deepcopy(self.yticloc)
        ztics = deepcopy(self.zticloc)

        x_range = np.arange(0., xlen) / float(xlen - 1)
        y_range = np.arange(0., ylen) / float(ylen - 1)
        z_range = np.arange(0., zlen) / float(zlen - 1)

        xtics /= float(xlen - 1)
        ytics /= float(ylen - 1)
        ztics /= float(zlen - 1)

        xticlen = len(xtics)
        yticlen = len(ytics)
        zticlen = len(ztics)

        ticlen = 0.03
        cutmax = 0.95
        cutmin = 0.05

        for i in range(0, xticlen):
            if cutmin <= xtics[xticlen - i - 1] <= cutmax:
                iso = self.isometric_projection(0, y_range[0], xtics[xticlen - i - 1])
                iso2 = self.isometric_projection(0, y_range[0] + ticlen, xtics[xticlen - i - 1])
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(0, y_range[ylen - 1], xtics[xticlen - i - 1])
                iso2 = self.isometric_projection(0, y_range[ylen - 1] - ticlen, xtics[xticlen - i - 1])
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

                iso = self.isometric_projection(z_range[0], 1, xtics[xticlen - 1 - i])
                iso2 = self.isometric_projection(z_range[0] + ticlen, 1, xtics[xticlen - 1 - i])
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(z_range[zlen - 1], 1, xtics[xticlen - 1 - i])
                iso2 = self.isometric_projection(z_range[zlen - 1] - ticlen, 1, xtics[xticlen - 1 - i])
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

        for j in range(0, yticlen):
            if cutmin <= ytics[j] <= cutmax:
                iso = self.isometric_projection(0, ytics[j], x_range[0])
                iso2 = self.isometric_projection(0, ytics[j], x_range[0] + ticlen)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(0, ytics[j], x_range[xlen - i - 1])
                iso2 = self.isometric_projection(0, ytics[j], x_range[xlen - i - 1] - ticlen)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

                iso = self.isometric_projection(z_range[0], ytics[j], 0)
                iso2 = self.isometric_projection(z_range[0] + ticlen, ytics[j], 0)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(z_range[zlen - 1], ytics[j], 0)
                iso2 = self.isometric_projection(z_range[zlen - 1] - ticlen, ytics[j], 0)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

        for k in range(0, zticlen):
            if cutmin <= ztics[k] <= cutmax:
                iso = self.isometric_projection(ztics[k], 1, x_range[0])
                iso2 = self.isometric_projection(ztics[k], 1, x_range[0] + ticlen)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(ztics[k], 1, x_range[xlen - i - 1])
                iso2 = self.isometric_projection(ztics[k], 1, x_range[xlen - i - 1] - ticlen)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

                iso = self.isometric_projection(ztics[k], y_range[0], 0)
                iso2 = self.isometric_projection(ztics[k], y_range[0] + ticlen, 0)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)
                iso = self.isometric_projection(ztics[k], y_range[ylen - 1], 0)
                iso2 = self.isometric_projection(ztics[k], y_range[ylen - 1] - ticlen, 0)
                ax.plot([iso[0, 0], iso2[0, 0]], [iso[1, 0], iso2[1, 0]], color=col)

    def isolabels(self, xticlab, yticlab, zticlab, ax, col):
        """
        Draw isometric projection of axis labels
        :param xticlab: X-axis label (string)
        :param yticlab: Y-axis label (string)
        :param zticlab: Z-axis label (string)
        :param ax: Axis object
        :param col: Color (not presently used)
        :return: None
        """
        # TODO: Clean up some of the duplicate code between this function and isoticks
        xlen = self.xlen
        ylen = self.ylen
        zlen = self.zlen
        xtics = deepcopy(self.xticloc)
        ytics = deepcopy(self.yticloc)
        ztics = deepcopy(self.zticloc)

        x_range = np.arange(0., xlen) / float(xlen - 1)
        y_range = np.arange(0., ylen) / float(ylen - 1)

        xtics /= float(xlen - 1)
        ytics /= float(ylen - 1)
        ztics /= float(zlen - 1)

        xticlen = len(xtics)
        yticlen = len(ytics)
        zticlen = len(ztics)

        ticlen = 0.03
        cutmin = 0.05
        cutmax = 0.95

        for i in range(0, xticlen):
            if cutmin <= xtics[xticlen - i - 1] <= cutmax:
                iso2 = self.isometric_projection(0, y_range[0] - ticlen, xtics[xticlen - i - 1])
                lab = xticlab[i]
                if lab > 100:
                    lab = int(lab)
                ax.text(iso2[0, 0], iso2[1, 0], str(lab), horizontalalignment="right", verticalalignment="top")

        for j in range(0, yticlen):
            if cutmin <= ytics[j] <= cutmax and ytics[j] >= cutmin:
                iso2 = self.isometric_projection(0, ytics[j], x_range[xlen - xticlen - 1] + ticlen)
                lab = yticlab[j]
                if lab > 100:
                    lab = int(lab)
                ax.text(iso2[0, 0], iso2[1, 0], str(lab), horizontalalignment="right", verticalalignment="center",
                        rotation='vertical')

        for k in range(0, zticlen):
            if cutmin <= ztics[k] <= cutmax:
                iso2 = self.isometric_projection(ztics[k], y_range[0] - ticlen, 0)
                lab = zticlab[k]
                if lab > 100 or True:
                    lab = int(lab)
                ax.text(iso2[0, 0], iso2[1, 0], str(lab), horizontalalignment="left", verticalalignment="top")

    def isolines(self):
        """
        Draw isometric projection of lines at the edge of the cube.
        :return: None
        """
        lines = [
            [[1, 0, 0], [0, 0, 0]],
            [[1, 1, 0], [0, 1, 0]],
            [[0, 1, 0], [0, 0, 0]],
            [[1, 1, 0], [1, 0, 0]],
            [[0, 0, 0], [0, 0, 1]],
            [[0, 1, 0], [0, 1, 1]],
            [[0, 0, 1], [0, 1, 1]],
            [[0, 1, 1], [1, 1, 1]],
            [[1, 1, 0], [1, 1, 1]]
        ]
        for line in lines:
            pt1 = line[0]
            pt2 = line[1]
            iso = self.isometric_projection(pt1[0], pt1[1], pt1[2])
            iso2 = self.isometric_projection(pt2[0], pt2[1], pt2[2])
            self.subplot1.plot(np.array([iso[0, 0], iso2[0, 0]]), np.array([iso[1, 0], iso2[1, 0]]), color=self.tickcolor)

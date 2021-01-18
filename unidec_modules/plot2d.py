import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FixedLocator
import numpy as np
import scipy.ndimage.filters as filt
from matplotlib.image import NonUniformImage
from unidec_modules.PlottingWindow import PlottingWindow
from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


# TODO: The 100% line on the color bar doesn't come right up to the top. Fix it so it is square.
class Plot2d(PlottingWindow):
    """
    Plotting class for 2D contour plots
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize parameters for plotting window
        :param args: Arguments passed to PlottingWindow
        :param kwargs: Keywords passed to PlottingWindow
        :return: Plot2d object
        """
        PlottingWindow.__init__(self, *args, **kwargs)

    def contourplot(self, dat=None, config=None, xvals=None, yvals=None, zgrid=None, xlab='m/z (Th)', ylab="Charge",
                    title='', normflag=1, normrange=[0, 1], repaint=True, nticks=None, test_kda=False, discrete=None,
                    ticloc=None, ticlab=None):
        """
        Make 2D plot.

        Data can be added using two methods:
            1. If dat is specified, it will look for an N x 3 list of x,y,z values
            2. If xvals, yvals, and zgrid are filled, it will plot zgrid assuming its shape is (len(xvals),len(yvals))
        :param dat: N x 3 list in [x,y,z] format of data to be plotted
        :param config: UniDecConfig object
        :param xvals: x-axis values
        :param yvals: y-axis values
        :param zgrid: numpy array of shape (len(xvals),len(yvals))
        :param xlab: Label for x-axis
        :param ylab: Label for y-axis
        :param title: Plot title
        :param normflag: Set to 1 to normalize the plot range from min to max. If 0, will not normalize.
        :param normrange: Range to normalize to if normflag is 0.
        :param repaint: If true, will repaint the plot after it is done.
        :param nticks: Number of ticks on the x-axis. If None, will use default.
        :param test_kda: If true, will decide if the x-axis is better displayed in Da or kDa units.
        :return: None
        """
        # Clear Plot
        self.clear_plot('nopaint')
        # Set xlabel and ylabel
        self.xlabel = xlab
        self.ylabel = ylab

        # Get values from config
        if config is not None:
            speedplot = config.discreteplot
            publicationmode = config.publicationmode
            self.cmap = config.cmap
        else:
            speedplot = 0
            publicationmode = 0
            self.cmap = u"jet"

        if discrete is not None:
            speedplot = discrete
        # Set Tick colors
        self.set_tickcolor()

        # If data is coming in as 1D list in dat, reshape it
        if xvals is None or yvals is None or zgrid is None:
            zgrid = dat[:, 2]
            xvals = np.unique(dat[:, 0])
            yvals = np.unique(dat[:, 1])
        xlen = len(xvals)
        ylen = len(yvals)
        newgrid = np.reshape(zgrid, (xlen, ylen))

        # Save Data
        if dat is None:
            X2, Y2 = np.meshgrid(xvals, yvals, indexing="ij")
            X2 = np.ravel(X2)
            Y2 = np.ravel(Y2)
            Z2 = np.ravel(newgrid.transpose())
            dat = np.transpose([X2, Y2, Z2])
        self.data = dat

        # Test if we should plot kDa instead of Da
        if test_kda:
            self.kda_test(xvals)

        # Decide whether or not to normalize
        if normflag == 1:
            norm = cm.colors.Normalize(vmax=np.amax(newgrid), vmin=np.amin(newgrid))
        else:
            norm = cm.colors.Normalize(vmax=normrange[1], vmin=normrange[0])

        # Add axes
        self.subplot1 = self.figure.add_axes(self._axes)
        # Plot
        # speedplot=0
        if speedplot == 0:
            # Slow contour plot that interpolates grid
            b1 = newgrid > 0.01 * np.amax(newgrid)
            #If the data is sparse, use the tricontourf, otherwise, use the regular contourf
            if np.sum(b1)/len(newgrid.ravel()) < 0.1:
                try:
                    b1 = b1.astype(float)
                    b1 = filt.uniform_filter(b1, size=3) > 0
                    b1[0, 0] = True
                    b1[0, -1] = True
                    b1[-1, 0] = True
                    b1[-1, -1] = True
                    X2, Y2 = np.meshgrid(xvals, yvals, indexing="ij")
                    X2 = np.ravel(X2[b1])
                    Y2 = np.ravel(Y2[b1])
                    Z2 = np.ravel(newgrid[b1].transpose())
                    cax = self.subplot1.tricontourf(X2 / self.kdnorm, Y2, Z2, 100, cmap=self.cmap, norm=norm)
                except Exception as e:
                    print("Error with fast tricontourf plot", e)
                    cax = self.subplot1.contourf(xvals / self.kdnorm, yvals, np.transpose(newgrid), 100, cmap=self.cmap,
                                                 norm=norm)
            else:
                cax = self.subplot1.contourf(xvals / self.kdnorm, yvals, np.transpose(newgrid), 100, cmap=self.cmap,
                                             norm=norm)

            datalims = [np.amin(xvals) / self.kdnorm, np.amin(yvals), np.amax(xvals) / self.kdnorm, np.amax(yvals)]
        else:
            # Fast discrete plot using imshow
            try:
                xdiff = (xvals[1] - xvals[0]) / self.kdnorm
                ydiff = yvals[1] - yvals[0]
            except:
                xdiff = 1
                ydiff = 1
            extent = (np.amin(xvals) / self.kdnorm - 0.5 * xdiff, np.amax(xvals) / self.kdnorm + 0.5 * xdiff,
                      np.amin(yvals) - 0.5 * ydiff, np.amax(yvals) + 0.5 * ydiff)

            try:
                ax = self.subplot1
                im = NonUniformImage(ax, interpolation="nearest", extent=extent, cmap=self.cmap, norm=norm,)
                im.set_data(xvals / self.kdnorm, yvals, np.transpose(newgrid))
                ax.images.append(im)
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
                cax = im
            except Exception as e:
                print("Error in NonUniformImage:", e)
                cax = self.subplot1.imshow(np.transpose(newgrid), origin="lower", cmap=self.cmap, extent=extent,
                                           aspect='auto', norm=norm, interpolation='nearest')
            datalims = [extent[0], extent[2], extent[1], extent[3]]
            print(newgrid.shape)
        # Set X and Y axis labels
        self.subplot1.set_xlabel(self.xlabel)
        self.subplot1.set_ylabel(self.ylabel)
        # Set Title
        if publicationmode == 0:
            self.subplot1.set_title(title)
        # Set colorbar
        if normflag == 1:
            self.cbar = self.figure.colorbar(cax, ax=None, use_gridspec=True,
                                             ticks=[0, np.amax(newgrid) / 2, np.amax(newgrid)])
            self.cbar.ax.get_yaxis().set_tick_params(direction='out')
            self.cbar.ax.set_yticklabels(["0", "%", "100"])
        else:
            self.cbar = self.figure.colorbar(cax, ax=None, use_gridspec=True)
        # Change tick colors
        if nticks is not None:
            self.subplot1.xaxis.set_major_locator(MaxNLocator(nbins=nticks))
        if ticloc is not None and ticlab is not None:
            self.subplot1.xaxis.set_major_locator(FixedLocator(ticloc))
            self.subplot1.set_xticklabels(ticlab, rotation=90, fontsize=8)

        '''
        for line in self.subplot1.xaxis.get_ticklines():
            line.set_color(self.tickcolor)
        for line in self.subplot1.yaxis.get_ticklines():
            line.set_color(self.tickcolor)
        '''
        # Setup zoom and repaint
        self.setup_zoom([self.subplot1], 'box', data_lims=datalims)
        if repaint:
            self.repaint()
        self.flag = True

    def plot_native_z(self, offset, col, xvals, width=0, alpha=1, shape=0):
        """
        Plots a showing the native charge offsets.
        :param offset: Offset value
        :param col: Color
        :param xvals: x-axis (mass)
        :param width: Width of charge distribution.
        :param alpha: Transparency of plot overlay
        :param shape: Indicates whether the width defines a step function of a Gaussian distribution
        :return: None
        """
        x1, x2, y1, y2 = self.subplot1.axis()
        yvals = ud.simchargefit(xvals) + offset
        self.subplot1.plot(xvals, yvals, color=col)
        if width > 0 and shape == 0:
            self.subplot1.fill_between(xvals / self.kdnorm, yvals - width, yvals + width, color=col, alpha=alpha)
        elif width > 0 and shape == 1:
            zbuff = width * 3
            start = -zbuff
            end = zbuff + 1
            zrange = range(start, end)
            for i in range(1, len(zrange)):
                weight = np.exp(-(zrange[i]) ** 2 / (2. * width * width))
                self.subplot1.fill_between(xvals / self.kdnorm, yvals + zrange[i - 1] + 0.5, yvals + zrange[i] + 0.5,
                                           color=col, alpha=alpha * weight, linewidth=0.0)
        self.subplot1.axis([x1, x2, y1, y2])
        self.nativez.append([offset, col])
        self.repaint()

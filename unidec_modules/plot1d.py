from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from unidec_modules.isolated_packages.ZoomSpan import ZoomSpan

from unidec_modules.PlottingWindow import PlottingWindow
from matplotlib.collections import LineCollection
import matplotlib.colorbar as colorbar
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt


class Plot1d(PlottingWindow):
    """
    Class for 1D plots.
    """

    def __init__(self, *args, **kwargs):
        """
        Inherit from PlottingWindow
        :param args:
        :param kwargs:
        :return:
        """
        PlottingWindow.__init__(self, *args, **kwargs)
        self.mlist = []
        self.x1, self.x2 = None, None
        self.colors = []

    def plotrefreshtop(self, xvals, yvals, title="", xlabel="", ylabel="", label="", config=None, color="black",
                       marker=None, zoom="box", nopaint=False, test_kda=False, integerticks=False, **kwargs):
        """
        Create a new 1D plot.
        :param xvals: x values
        :param yvals: y values
        :param title: Plot title
        :param xlabel: x axis label
        :param ylabel: y axis label
        :param label: label for legend
        :param config: UniDecConfig object
        :param color: Plot color
        :param marker: Marker type
        :param zoom: zoom type ('box' or 'span')
        :param nopaint: Boolean, whether to repaint or not
        :param test_kda: Boolean, whether to attempt to plot kDa rather than Da
        :param integerticks: Boolean, whether to use inter tickmarks only
        :param kwargs: Keywords
        :return: None
        """
        self.clear_plot("nopaint")
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zoomtype = zoom
        if "nticks" in kwargs:
            nticks = kwargs["nticks"]
            del kwargs['nticks']
        else:
            nticks = None

        if test_kda:
            self.kda_test(xvals)
        else:
            self.kdnorm = 1

        pubflag = 0
        if config is not None:
            if config.publicationmode != 0:
                pubflag = 1

        if pubflag == 0:
            self.subplot1 = self.figure.add_axes(self._axes)
            self.subplot1.plot(np.array(xvals) / self.kdnorm, yvals, color=color, label=label, marker=marker, **kwargs)
            self.subplot1.set_ylabel(self.ylabel)
            self.subplot1.set_title(title)
        else:
            self.subplot1 = self.figure.add_axes(self._axes)
            self.subplot1.plot(np.array(xvals) / self.kdnorm, yvals, color=color, label=label, marker=marker, **kwargs)
            self.subplot1.spines['top'].set_visible(False)
            self.subplot1.spines['right'].set_visible(False)
            self.subplot1.get_xaxis().tick_bottom()
            self.subplot1.get_yaxis().tick_left()
            self.subplot1.get_yaxis().set_tick_params(direction='out')
            self.subplot1.get_xaxis().set_tick_params(direction='out')
            if config.peaknorm != 2:
                self.subplot1.get_yaxis().set_ticks([0, np.amax(yvals) / 2, np.amax(yvals)])
                self.subplot1.get_yaxis().set_ticklabels(["0", '%', "100"])
            else:
                self.subplot1.set_ylabel("Relative Intensity")

        if nticks is not None:
            self.subplot1.xaxis.set_major_locator(MaxNLocator(nbins=nticks))
        if integerticks:
            self.subplot1.xaxis.set_major_locator(MaxNLocator(integer=True))

        self.subplot1.set_xlabel(self.xlabel)
        self.subplot1.set_clip_on(True)

        self.setup_zoom([self.subplot1], self.zoomtype)
        if not nopaint:
            # self.setup_zoom([self.subplot1], self.zoomtype)
            self.repaint()
        self.flag = True
        self.mlist = []
        self.x1, self.x2 = None, None
        self.colors = []

    def plotadd(self, xvals, yvals, colval="black", newlabel="", nopaint=True, **kwargs):
        """
        Adds a plot.
        :param xvals: x values
        :param yvals: y values
        :param colval: Color
        :param newlabel: Label
        :param nopaint: Boolean, whether to repaint or not
        :return: None
        """
        self.subplot1.plot(np.array(xvals) / self.kdnorm, yvals, color=colval, label=newlabel, **kwargs)
        self.setup_zoom([self.subplot1], self.zoomtype)
        if not nopaint:
            # self.setup_zoom([self.subplot1], self.zoomtype)
            self.repaint()

    def multiplot(self, x, y, a, b):

        mpl.rcParams['ps.useafm'] = True
        mpl.rcParams['ps.fonttype'] = 42
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['lines.linewidth'] = 0.75
        mpl.rcParams['errorbar.capsize'] = 3
        mpl.rcParams['patch.force_edgecolor'] = True
        mpl.rcParams['patch.facecolor'] = 'b'
        mpl.rcParams['lines.markersize'] = 7
        mpl.rcParams['font.size'] = 8
        mpl.rcParams['font.sans-serif'] = "Arial"

        self.sp1 = self.figure.add_subplot(221)
        self.sp2 = self.figure.add_subplot(222)
        self.sp3 = self.figure.add_subplot(223)
        self.sp4 = self.figure.add_subplot(224)

        self.subplots = [self.sp1, self.sp2, self.sp3, self.sp4]

        for subplots in self.subplots:
            subplots.plot(x, y)
            subplots.spines['top'].set_visible(False)
            subplots.spines['right'].set_visible(False)
            subplots.get_yaxis().set_ticks([0, 0.5 * max(y), max(y)])
            subplots.get_yaxis().set_ticklabels(["0", '%', "100"])
            subplots.set_ylim(min(y), max(y))
            subplots.set_xlim(min(x), max(x))
            subplots.set_xlabel("m/z")

        self.axins1 = inset_axes(self.sp1, width="40%", height="40%",
                                 bbox_to_anchor=(.65, .65, 1.0, 1.0),
                                 bbox_transform=self.sp1.transAxes, loc=3)
        self.axins2 = inset_axes(self.sp2, width="40%", height="40%",
                                 bbox_to_anchor=(.65, .65, 1.0, 1.0),
                                 bbox_transform=self.sp2.transAxes, loc=3)
        self.axins3 = inset_axes(self.sp3, width="40%", height="40%",
                                 bbox_to_anchor=(.65, .65, 1.0, 1.0),
                                 bbox_transform=self.sp3.transAxes, loc=3)
        self.axins4 = inset_axes(self.sp4, width="40%", height="40%",
                                 bbox_to_anchor=(.65, .65, 1.0, 1.0),
                                 bbox_transform=self.sp4.transAxes, loc=3)
        self.insetaxes = [self.axins1, self.axins2, self.axins3, self.axins4]

        for insetaxes in self.insetaxes:
            insetaxes.plot(a / 1000, b)
            insetaxes.spines['top'].set_visible(False)
            insetaxes.spines['right'].set_visible(False)
            insetaxes.get_yaxis().set_ticks([0, 0.5 * max(b), max(b)])
            insetaxes.get_yaxis().set_ticklabels(["0", '%', "100"])
            insetaxes.set_ylim(min(b), max(b))
            insetaxes.set_xlim(min(a), max(a))
            insetaxes.set_xlabel("Mass (kDa)")

        self.setup_zoom(self.insetaxes, "box", link_axes=True)
        self.setup_zoom(self.subplots, "box", link_axes=True)

        self.repaint()

    def insetplot(self, x, y, a, b):
        self.subplot1 = self.figure.add_axes(self._axes)

        self.axins = inset_axes(self.subplot1, width="40%", height="40%",
                                bbox_to_anchor=(.65, .65, 1.0, 1.0),
                                bbox_transform=self.subplot1.transAxes, loc=3)

        mpl.rcParams['figure.dpi'] = 150
        mpl.rcParams['lines.linewidth'] = 0.75

        self.subplot1.plot(x, y)
        self.subplot1.get_yaxis().set_ticks([0, 0.5 * max(y), max(y)])
        self.subplot1.get_yaxis().set_ticklabels(["0", '%', "100"])
        self.subplot1.set_ylim(min(y), max(y))
        self.subplot1.set_xlim(min(x), max(x))
        self.subplot1.set_xlabel("m/z")
        self.subplot1.spines['top'].set_visible(False)
        self.subplot1.spines['right'].set_visible(False)

        self.axins.plot(a / 1000, b)
        self.axins.get_yaxis().set_ticks([0, 0.5 * max(b), max(b)])
        self.axins.get_yaxis().set_ticklabels(["0", '%', "100"])
        self.axins.set_ylim(min(b), max(b))
        self.axins.set_xlim(min(a) / 1000, max(a) / 1000)
        self.axins.set_xlabel("Mass (KDa)")
        self.axins.spines['top'].set_visible(False)
        self.axins.spines['right'].set_visible(False)

        self.setup_zoom([self.subplot1, self.axins], "box", link_axes=False)

        self.repaint()

    def errorbars(self, xvals, yvals, xerr=None, yerr=None, color="black", newlabel="", nopaint=True, **kwargs):
        self.subplot1.errorbar(np.array(xvals) / self.kdnorm, yvals, xerr=xerr, yerr=yerr, color=color, label=newlabel,
                               **kwargs)
        self.setup_zoom([self.subplot1], self.zoomtype, pad=0.02)
        if not nopaint:
            self.repaint()

    def addtext(self, txt, x, y, vlines=True, color="k", ymin=0, ymax=None, **kwargs):
        """
        Adds text and lines. Puts things that have been added in self.lines and self.text lists.
        :param txt: String of text to add
        :param x: x position for where to add
        :param y: y position for where to add
        :param vlines: Whether to add vertical lines to the text as well.
        :param color: Color of text and lines
        :param kwargs: Keywords
        If range=(a,b) is specified, adds a line from a to b and vertical lines at a and b.
        :return: None
        """
        if ymax is None:
            ymax = y
        text = self.subplot1.text(np.array(x) / self.kdnorm, y, txt, horizontalalignment="center",
                                  verticalalignment="top", color=color)
        self.text.append(text)
        if vlines:
            if "range" in kwargs:
                line_range = kwargs["range"]
                line = self.subplot1.vlines(line_range[0] / self.kdnorm, ymin, y * 0.6, color=color)
                self.lines.append(line)
                line = self.subplot1.vlines(line_range[1] / self.kdnorm, ymin, y * 0.6, color=color)
                self.lines.append(line)
                line = self.subplot1.hlines(y * 0.3, line_range[0] / self.kdnorm, line_range[1] / self.kdnorm,
                                            linestyles="dashed", color=color)
                self.lines.append(line)
                pass
            else:
                line = self.subplot1.vlines(np.array(x) / self.kdnorm, ymin, y - 0.05 * ymax, linestyles="dashed",
                                            color=color)
                self.lines.append(line)
        self.repaint()

    def textremove(self):
        """
        Remove text and lines previous placed in the self.text and self.lines containers
        :return: None
        """
        if len(self.text) > 0:
            for i in range(0, len(self.text)):
                self.text[i].remove()
                try:
                    self.lines[i].remove()
                except:
                    print(self.text[i])
        self.text = []
        self.lines = []
        self.repaint()

    def filledplot(self, x, y, color="black"):
        """
        Creates a plot filled between the y values and the y axis
        :param x: x values
        :param y: y values
        :param color: color
        :return: None
        """
        self.subplot1.fill_between(np.array(x) / self.kdnorm, y, y2=0, facecolor=color, alpha=0.75)

    def histogram(self, xarray, labels=None, xlab="", ylab="", title=""):
        """
        Creates a histogram plot.
        :param xarray: Array of values in (N x M) array
        :param labels: Labels for value
        :param xlab: x axis label
        :param ylab: y axis label
        :param title: Plot Title
        :return: None
        """
        self.clear_plot("nopaint")
        self.subplot1 = self.figure.add_axes(self._axes)

        ymax = 0
        xmax = 0
        for i in range(0, len(xarray)):
            if labels is not None:
                label = "KD" + str(labels[i])
            else:
                label = ""

            xvals = xarray[i]
            n, bins, patches = self.subplot1.hist(xvals, label=label, histtype="stepfilled", alpha=0.4, density=1)
            ytempmax = np.amax(n)
            xtempmax = np.amax(bins)
            if ytempmax > ymax:
                ymax = ytempmax
            if xtempmax > xmax:
                xmax = xtempmax

        self.subplot1.set_xlabel(xlab)
        self.subplot1.set_ylabel(ylab)
        self.subplot1.set_title(title)

        self.add_legend(location=2)
        self.setup_zoom([self.subplot1], 'box')
        self.repaint()

    def barplottop(self, xarr, yarr, peakval, colortab, xlabel="", ylabel="", title="", zoom="box", repaint=True):
        """
        Create a bar plot.
        :param xarr: x value array
        :param yarr: y value array
        :param peakval: Bar labels to be listed below bars
        :param colortab: List of colors for the various bars
        :param xlabel: Label for x axis
        :param ylabel: Label for y axis
        :param title: Plot title
        :param zoom: Type of zoom ('box' or 'span)
        :param repaint: Boolean, whether to repaint or not when done.
        :return: None
        """
        self.clear_plot("nopaint")
        self.zoomtype = zoom
        xticloc = np.array(xarr)
        peaklab = [str(p) for p in peakval]
        self.subplot1 = self.figure.add_axes(self._axes, xticks=xticloc)
        self.subplot1.set_xticklabels(peaklab, rotation=90, fontsize=8)
        self.subplot1.bar(xarr, yarr, color=colortab, label="Intensities", width=1)
        self.subplot1.set_xlabel(xlabel)
        self.subplot1.set_ylabel(ylabel)
        self.subplot1.set_title(title)
        self.subplot1.spines['top'].set_visible(False)
        self.subplot1.spines['right'].set_visible(False)
        self.subplot1.set_clip_on(False)
        self.setup_zoom([self.subplot1], self.zoomtype)
        self.flag = True

    # TODO make the axes work for negative and positive bars
    def barplottoperrors(self, xarr, yarr, peakval, colortab, xlabel="", ylabel="", title="", zoom="box", repaint=True,
                         xerr=0, yerr=0):
        """
        Create a bar plot.
        :param xarr: x value array
        :param yarr: y value array
        :param peakval: Bar labels to be listed below bars
        :param colortab: List of colors for the various bars
        :param xlabel: Label for x axis
        :param ylabel: Label for y axis
        :param title: Plot title
        :param zoom: Type of zoom ('box' or 'span)
        :param repaint: Boolean, whether to repaint or not when done.
        :param xerr: for error bars
        :param yerr: for error bars
        :return: None
        """
        self.clear_plot("nopaint")
        self.zoomtype = zoom
        xticloc = np.array(xarr)
        peaklab = [str(p) for p in peakval]
        self.subplot1 = self.figure.add_axes(self._axes, xticks=xticloc, ymargin=1)
        self.subplot1.set_xticklabels(peaklab, rotation=90, fontsize=8)
        self.subplot1.bar(xarr, yarr, color=colortab, label="Intensities", width=1, xerr=xerr, yerr=yerr)
        self.subplot1.set_xlabel(xlabel)
        self.subplot1.set_ylabel(ylabel)
        self.subplot1.set_title(title)
        # Adjust axes for error bars
        # Negative bars
        if np.amin(np.asarray(yarr)) < 0 and np.amax(np.asarray(yarr)) <= 0:
            minindex = 0
            for x in range(len(yarr)):
                if yarr[x] - yerr[x] < yarr[minindex] - yerr[minindex]:
                    minindex = x
            left = yarr[minindex] - yerr[minindex]
            right = np.amax(np.asarray(yarr))
        # Positive bars
        elif np.amax(np.asarray(yarr)) > 0 and np.amin(np.asarray(yarr)) >= 0:
            maxindex = 0
            for x in range(len(yarr)):
                if yarr[x] + yerr[x] > yarr[maxindex] + yerr[maxindex]:
                    maxindex = x
            left = np.amin(np.asarray(yarr))
            right = yarr[maxindex] + yerr[maxindex]
        # Negative and positive bars
        else:
            maxindex = 0
            minindex = 0
            for x in range(len(yarr)):
                if yarr[x] + yerr[x] > yarr[maxindex] + yerr[maxindex]:
                    maxindex = x
                if yarr[x] - yerr[x] < yarr[minindex] - yerr[minindex]:
                    minindex = x
            left = yarr[minindex] - yerr[minindex]
            right = yarr[maxindex] + yerr[maxindex]
        self.setup_zoom([self.subplot1], self.zoomtype)
        self.flag = True
        self.subplot1.set_ylim(left, right)
        self.subplot1.spines['top'].set_visible(False)
        self.subplot1.spines['right'].set_visible(False)
        self.subplot1.set_clip_on(False)
        if repaint:
            self.repaint()

    def colorplotMD(self, xvals, yvals, cvals, title="", xlabel="", ylabel="", clabel="Mass Defect", cmap="hsv",
                    config=None,
                    color="black", max=1,
                    marker=None, zoom="box", nopaint=False, test_kda=False, integerticks=False, **kwargs):
        self._axes = [0.11, 0.1, 0.64, 0.8]
        self.clear_plot("nopaint")
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zoomtype = zoom
        if "nticks" in kwargs:
            nticks = kwargs["nticks"]
            del kwargs['nticks']
        else:
            nticks = None

        if test_kda:
            self.kda_test(xvals)
            xvals = xvals / self.kdnorm

        pubflag = 0
        if config is not None:
            if config.publicationmode != 0:
                pubflag = 1

        self.subplot1 = self.figure.add_axes(self._axes)

        points = np.array([xvals, yvals]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # print segments
        t = cvals  # np.linspace(0, 1, len(xvals), endpoint=True)
        lc = LineCollection(segments, cmap=cmap)
        lc.set_array(t)
        lc.set_clim(0, max)
        # print(max, np.amax(t))
        self.subplot1.add_collection(lc)

        if pubflag == 0:
            self.subplot1.set_ylabel(self.ylabel)
            self.subplot1.set_title(title)
        else:

            self.subplot1.spines['top'].set_visible(False)
            self.subplot1.spines['right'].set_visible(False)
            self.subplot1.get_xaxis().tick_bottom()
            self.subplot1.get_yaxis().tick_left()
            self.subplot1.get_yaxis().set_tick_params(direction='out')
            self.subplot1.get_xaxis().set_tick_params(direction='out')
            if config.peaknorm != 2:
                self.subplot1.get_yaxis().set_ticks([0, np.amax(yvals) / 2, np.amax(yvals)])
                self.subplot1.get_yaxis().set_ticklabels(["0", '%', "100"])
            else:
                self.subplot1.set_ylabel("Relative Intensity")

        if nticks is not None:
            self.subplot1.xaxis.set_major_locator(MaxNLocator(nbins=nticks))
        if integerticks:
            self.subplot1.xaxis.set_major_locator(MaxNLocator(integer=True))

        self.subplot1.set_xlabel(self.xlabel)
        cax = self.figure.add_axes([0.77, 0.1, 0.04, 0.8])
        ticks = np.linspace(0., 1., 11, endpoint=True)

        cmap2 = cm.get_cmap(cmap)
        self.cbar = colorbar.ColorbarBase(cax, cmap=cmap2, orientation="vertical", ticks=ticks)
        ticks = ticks * max
        if max > 10:
            ticks = np.round(ticks).astype(np.int)
        else:
            ticks = np.round(ticks, 1)
        ticks = ticks.astype(str)
        self.cbar.set_ticklabels(ticks)
        self.cbar.set_label(clabel)

        self.setup_zoom([self.subplot1], self.zoomtype,
                        data_lims=[np.amin(xvals), np.amin(yvals), np.max(xvals), np.amax(yvals)])

        if not nopaint:
            self.repaint()
        self.flag = True
        self.mlist = []
        self.x1, self.x2 = None, None
        self.colors = []

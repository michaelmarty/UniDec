# base class for all plotting windows used in UniDec
# contains basic setup functionality
from unidec_modules.isolated_packages.ZoomCommon import *
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
import matplotlib.cm as cm
import numpy as np

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7

# rcParams['axes.linewidth']=1
# rcParams['font.size']=18
# matplotlib.rc('font', family='sans-serif')
# matplotlib.rc('font', serif='Helvetica')


class PlotBase(object):
    """
    Class for matplotlib plots
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize plot window parameters.

        Optional keywords:
        figsize: size of figure in inches

        :param args: Arguments
        :param kwargs: Keywords
        :return:
        """
        self.defaultfigsize = (6. * 0.9, 5. * 0.9)
        if "figsize" in kwargs:
            figsize = kwargs["figsize"]
            del kwargs["figsize"]
        else:
            figsize = self.defaultfigsize
        self.figsize = figsize

        if "axes" in kwargs:
            self._axes = kwargs["axes"]
            del kwargs["axes"]
        else:
            if figsize[0] < 5:
                self._axes = [0.2, 0.2, 0.7, 0.7]
            else:
                self._axes = [0.11, 0.11, 0.8, 0.8]

        self.figure = Figure(figsize=figsize)  # , dpi=
        self.subplot1 = None
        self.zoom = None
        self.subplot1 = None
        self.resize = 1
        self.flag = False
        self.kda = False
        self.kdnorm = 1.
        self.normalticks = True
        self.nativez = []
        self.text = []
        self.lines = []
        self.cbar = None
        self.datalims = None
        self.data = None
        self.cmap = None
        self.set_color()
        self.xlabel = ""
        self.ylabel = ""
        self.zoomtype = "box"
        self.tickcolor = "black"
        self.mouse_active = False

    def repaint(self, setupzoom=False):
        if setupzoom:
            self.setup_zoom()
        try:
            self.zoomout()
        except:
            pass

    def zoomout(self):
        set_clipping(self.subplot1)

    def setup_zoom(self, plots=None, zoom=None, data_lims=None, pad=0, groups=None):
        # Autoscale Y
        xmin, ymin, xmax, ymax = GetMaxes(self.subplot1)
        self.subplot1.set_ylim((ymin, ymax))
        self.subplot1.set_xlim((xmin, xmax))
        pass

    def get_blank_axis(self, scale=None):
        if scale is None:
            scale = self._axes
        self.clear_plot("nopaint")
        self.subplot1 = self.figure.add_axes(scale)
        return self.subplot1

    def on_save_fig(self, evt, path, **kwargs):
        """
        Save figure to path.
        :param evt: wx.Event (unused)
        :param path: Path to save figure to
        :param kwargs: keywords passed to save_figure
        :return: None
        """
        if path is not None:
            self.save_figure(path, **kwargs)

    def save_figure(self, path, **kwargs):
        """
        Saves Figure to path.
        :param path: Path to save figure at.
        :param kwargs: Keywords passed to matplotlib.figure.savefig (note only specific ones are passed)
        :return: None
        """
        if "transparent" in kwargs:
            t = kwargs["transparent"]
        else:
            t = True
        if "dpi" in kwargs:
            dpi = kwargs["dpi"]
        else:
            dpi = None
        self.figure.savefig(path, transparent=t, dpi=dpi)
        print("Saved Figure: ", path)

    def kda_test(self, xvals):
        """
        Test whether the axis should be normalized to convert mass units from Da to kDa.
        Will use kDa if: xvals[int(len(xvals) / 2)] > 100000 or xvals[len(xvals) - 1] > 1000000

        If kDa is used, self.kda=True and self.kdnorm=1000. Otherwise, self.kda=False and self.kdnorm=1.
        :param xvals: mass axis
        :return: None
        """
        try:
            if xvals[int(len(xvals) / 2)] > 20000 or xvals[len(xvals) - 1] > 150000:
                self.kdnorm = 1000.
                self.xlabel = "Mass (kDa)"
                self.kda = True
            else:
                self.xlabel = "Mass (Da)"
                self.kda = False
                self.kdnorm = 1.
        except (TypeError, ValueError):
            self.xlabel = "Mass (Da)"
            self.kdnorm = 1.
            self.kda = False

    def plotadddot(self, x, y, colval, markval, label=""):
        """
        Adds a scatter plot to the figure. May be one or more.
        :param x: x values
        :param y: y values
        :param colval: Color
        :param markval: Marker
        :param label: Label for Plot
        :return: None
        """
        self.subplot1.plot(np.array(x) / self.kdnorm, y, color=colval, marker=markval, linestyle='None', clip_on=True,
                           markeredgecolor="k", label=label)

    def addtext(self, txt, x, y, vlines=True, hlines=False, color="k", ymin=0, ymax=None, verticalalignment="top",
                xmin=0, xmax=None, nopaint=False, **kwargs):
        """
        Adds text and lines. Puts things that have been added in self.lines and self.text lists.
        :param txt: String of text to add
        :param x: x position for where to add
        :param y: y position for where to add
        :param vlines: Whether to add vertical lines to the text as well.
        :param hlines: Whether to add horizontal lines to the text as well.
        :param color: Color of text and lines
        :param ymin: Minimum value of vlines
        :param ymax: Maximum value of vlines
        :param verticalalignment: Vertical alignment of the text
        :param xmin: Minimum value of hlines
        :param xmax: Maximum value of hlines
        :param nopaint: Don't paint it
        :param kwargs: Keywords
        If range=(a,b) is specified, adds a line from a to b and vertical lines at a and b.
        :return: None
        """
        if ymax is None:
            ymax = y
        if xmax is None:
            xmax = x

        text = self.subplot1.text(np.array(x) / self.kdnorm, y, txt, horizontalalignment="center",
                                  verticalalignment=verticalalignment, color=color)
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

        if hlines:
            line = self.subplot1.hlines(y, xmin / self.kdnorm, xmax / self.kdnorm - 0.05 * xmax / self.kdnorm,
                                        linestyles="dashed",
                                        color=color)
            self.lines.append(line)

        if not nopaint:
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

    def clear_plot(self, *args):
        """
        Clear the plot and rest some of the parameters.
        :param args: Arguments
        :return:
        """
        self.figure.clear()
        self.flag = False
        self.nativez = []
        self.text = []
        self.lines = []
        self.kda = False
        self.kdnorm = 1.
        if "nopaint" not in args:
            self.repaint()

    def set_nticks(self, bins):
        """
        Set the number of ticks in the x-axis.
        :param bins: Number of ticks in the x-axis
        :return: None
        """
        if self.normalticks:
            self.subplot1.tick_params(axis="x", labelsize=12)
            self.subplot1.tick_params(axis="y", labelsize=12)
            self.subplot1.xaxis.set_major_locator(MaxNLocator(nbins=bins))
        self.repaint()

    def add_legend(self, location=1, anchor=None):
        """
        Adds a legend to the plot.
        :param location: Integer code for location
        :param anchor: BBox to anchor parameter
        :return: None
        """
        handles, labels = self.subplot1.get_legend_handles_labels()
        if anchor is None:
            anchor = (1, 1)
        if location == 1:
            self.subplot1.legend(handles, labels, loc=location, bbox_to_anchor=anchor)
        else:
            self.subplot1.legend(handles, labels, loc=location)
        self.repaint()

    def add_title(self, title=""):
        self.subplot1.set_title(title)
        self.repaint()

    def set_color(self, rgbtuple=None):
        """
        Sets background color
        :param rgbtuple: background color
        :return:
        """
        # Set figure and canvas colours to be the same
        if not rgbtuple:
            rgbtuple = [255., 255., 255.]
        col = [c / 255.0 for c in rgbtuple]
        self.figure.set_facecolor(col)
        self.figure.set_edgecolor(col)

    def set_tickcolor(self):
        """
        Sets tick colors based on the colormap set at self.cmap
        :return: None
        """
        if self.cmap[:2] == "b'":
            self.cmap = self.cmap[2:-1]
        try:
            self.cmap = str(self.cmap, encoding="utf-8")
        except:
            pass

        output = cm.ScalarMappable(norm=None, cmap=str(self.cmap)).to_rgba(0)

        if sum(output[:2]) > 0.9:
            self.tickcolor = u"black"
        else:
            self.tickcolor = u"white"

    def write_data(self, path):
        if self.data is not None:
            print("Saving Data to", path)
            print("Data Dimensions:", self.data.shape)
            np.savetxt(path, self.data)

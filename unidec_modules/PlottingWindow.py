from __future__ import unicode_literals
# base class for all plotting windows used in UniDec
# contains basic setup functionality

import wx
import tempfile, os
from matplotlib import interactive
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
# from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
# from matplotlib.backends.backend_wx import _load_bitmap
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
# import matplotlib
import matplotlib.cm as cm
import numpy as np

from unidec_modules.isolated_packages.ZoomSpan import ZoomSpan
from unidec_modules.isolated_packages.ZoomBox import ZoomBox
from unidec_modules.isolated_packages.NoZoomSpan import NoZoomSpan
from unidec_modules.isolated_packages import FileDialogs
from unidec_modules.miscwindows import DoubleInputDialog

# import matplotlib.style as mplstyle
# mplstyle.use('fast')

interactive(True)

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


class PlottingWindow(wx.Window):
    """
    Class for wx window with embedded matplotlib plots
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize plot window parameters.

        Optional keywords:
        figsize: size of figure in inches
        integrate: 0 or 1 of whether the plot should send the integrate pubsub when a right click is activated.
        smash: 0 or 1 of whether the plot should send the integrate pubsub when a right click is activated.

        :param args: Arguments
        :param kwargs: Keywords
        :return:
        """
        self.displaysize = wx.GetDisplaySize()
        self.defaultfigsize = (6. * 0.9, 5. * 0.9)
        if "figsize" in kwargs:
            figsize = kwargs["figsize"]
            del kwargs["figsize"]
        else:
            figsize = self.defaultfigsize

        if "axes" in kwargs:
            self._axes = kwargs["axes"]
            del kwargs["axes"]
        else:
            if figsize[0] < 5:
                self._axes = [0.2, 0.2, 0.7, 0.7]
            else:
                self._axes = [0.11, 0.11, 0.8, 0.8]
        self.figsize = figsize

        if "integrate" in kwargs:
            self.int = kwargs["integrate"]
            del kwargs["integrate"]
        else:
            self.int = 0
        if "smash" in kwargs:
            self.smash = kwargs["smash"]
            del kwargs["smash"]
        else:
            self.smash = 0

        wx.Window.__init__(self, *args, **kwargs)
        self.figure = Figure(figsize=figsize)  # , dpi=
        self.subplot1 = None
        self.zoom = None
        self.subplot1 = None
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.Bind(wx.EVT_SIZE, self.size_handler)
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
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('key_press_event', self.on_key)

    def on_release(self, event):
        """
        Function triggered on button release event from plot.
        Currently wired to trigger on_save_figure_dialog on middle button.
        :param event: wx.Event
        :return: None
        """
        if event.button == 1:
            if wx.GetKeyState(wx.WXK_ALT):
                try:
                    self.zoom.switch_label()
                except:
                    print("Could not switch on labels")
        if event.button == 2:
            if wx.GetKeyState(wx.WXK_CONTROL):
                dlg = DoubleInputDialog(self)
                dlg.initialize_interface("Matplotlib RC Parameters", "RC Param Name:", 'lines.markersize',
                                         "Value:", "6")
                dlg.ShowModal()
                rcname = dlg.value
                rcval = dlg.value2
                print(rcname, rcval)
                rcParams[rcname] = rcval
            elif wx.GetKeyState(wx.WXK_ALT):
                dlg = DoubleInputDialog(self)
                dlg.initialize_interface("Set Plot X Range", "Min:", '',
                                         "Max:", "")
                dlg.ShowModal()
                minval = dlg.value
                maxval = dlg.value2

                try:
                    minval = float(minval)
                    maxval = float(maxval)
                    self.zoom.set_manual(minval, maxval)
                    print("Manually Set Zoom:", minval, maxval)
                except:
                    print("Error converting string to float:", minval, maxval)
            elif wx.GetKeyState(wx.WXK_SHIFT):
                dlg = DoubleInputDialog(self)
                dlg.initialize_interface("Set Plot Y Range", "Min:", '',
                                         "Max:", "")
                dlg.ShowModal()
                minval = dlg.value
                maxval = dlg.value2

                try:
                    minval = float(minval)
                    maxval = float(maxval)
                    self.zoom.set_manual_y(minval, maxval)
                    print("Manually Set Zoom:", minval, maxval)
                except:
                    print("Error converting string to float:", minval, maxval)
            elif wx.GetKeyState(wx.WXK_SPACE):
                try:
                    self.zoom.switch_label()
                except:
                    print("Could not switch on labels")
            else:
                self.on_save_fig_dialog(event)

    def on_key(self, evt):
        # print("you pressed", evt.key)
        if evt.key == "ctrl+c":
            self.copy_to_clipboard()
        if evt.key == "ctrl+u":
            self.on_write_dialog(evt)

    def on_save_fig_dialog(self, evt):
        """
        Open a save figure dialog for specified plot.
        :param evt: wx.Event (unused)
        :return: None
        """
        path = FileDialogs.save_file_dialog()
        if path is not None:
            self.save_figure(path)

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
        :return: None
        """
        self.subplot1.plot(np.array(x) / self.kdnorm, y, color=colval, marker=markval, linestyle='None', clip_on=True
                           , markeredgecolor="k", label=label)

    def addtext(self, txt, x, y, vlines=True, hlines=False, color="k", ymin=0, ymax=None, verticalalignment="top",
                xmin=0, xmax=None, nopaint=False, **kwargs):
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

    def repaint(self, setupzoom=False):
        """
        Redraw and refresh the plot.
        :return: None
        """
        if setupzoom:
            self.setup_zoom([self.subplot1], self.zoomtype)
        self.canvas.draw()

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
        # self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))

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
        '''
        if self.cmap[-1] == "r":
            self.tickcolor = "black"
        else:
            self.tickcolor = "white"
        '''

    def size_handler(self, *args, **kwargs):
        """
        Resizes the plots
        :param args:
        :param kwargs:
        :return: None
        """
        if self.resize == 1:
            self.canvas.SetSize(self.GetSize())

    def copy_to_clipboard(self, *args, **kwargs):
        obj = tempfile.NamedTemporaryFile(delete=False)
        self.canvas.print_figure(obj, format="png", dpi=300)
        obj.close()
        img = wx.Image(obj.name)
        btm = wx.Bitmap(img)
        bobj = wx.BitmapDataObject(btm)
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(bobj)
            # wx.TheClipboard.SetData(wx.TextDataObject("Test"))
            wx.TheClipboard.Close()
            print("Image Copied")
        os.remove(obj.name)

    def on_write_dialog(self, evt):
        """
        Open a save figure dialog for specified plot.
        :param evt: wx.Event (unused)
        :return: None
        """
        if self.data is not None:
            path = FileDialogs.save_file_dialog()
            if path is not None:
                self.write_data(path)
        else:
            print("Data object empty. Plot likely unsupported for text file output.")

    def write_data(self, path):
        if self.data is not None:
            print("Saving Data to", path)
            print("Data Dimensions:", self.data.shape)
            np.savetxt(path, self.data)

    def setup_zoom(self, plots, zoom, data_lims=None, pad=0, groups=None):
        """
        Set up zoom on axes.
        :param plots: Axes objects to setup
        :param zoom: Type of zoom ('span' or 'box')
        :param data_lims: Optional manual description of the data limits (where to go when fully zoomed out)
        :param groups: Group of axes to link zoom together for
        :return: None
        """
        # setup for zoom box
        if groups is None:
            groups = [1, 2, 3, 4, 5, 6, 7, 8]
        if zoom == 'span':
            self.zoom = ZoomSpan(
                plots,
                None,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))
        if zoom == 'box':
            self.zoom = ZoomBox(
                plots,
                None,
                groups=groups,
                drawtype='box',
                useblit=True,
                button=1,  # so zoombox is left button
                onmove_callback=None,
                spancoords='data',
                rectprops=dict(alpha=0.2, facecolor='yellow'),
                data_lims=data_lims,
                integrate=self.int, smash=self.smash, pad=pad)
        if zoom == "fixed_span":
            self.zoom = NoZoomSpan(
                plots,
                None,
                minspan=0,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))

    '''
    def add_toolbar(self):
        self.toolbar = MyNavigationToolbar(self.canvas, True)
        self.toolbar.Realize()'''


"""
This technically could be used, but it didn't work well when I tried it.
class MyNavigationToolbar(NavigationToolbar):
    def __init__(self, canvas, cankill):
        NavigationToolbar.__init__(self, canvas)

        # for simplicity I'm going to reuse a bitmap from wx, you'll
        # probably want to add your own.
        tool = self.AddTool(wx.ID_ANY, 'Click me', _load_bitmap('back.png'),
                            'Activate custom contol')
        self.Bind(wx.EVT_TOOL, self._on_custom, id=tool.GetId())

    def _on_custom(self, evt):
        # add some text to the axes in a random location in axes (0,1)
        # coords) with a random color

        # get the axes
        ax = self.canvas.figure.axes[0]

        # generate a random location can color
        x, y = np.random.rand(2)
        rgb = np.random.rand(3)

        # add the text and draw
        ax.text(x, y, 'You clicked me',
                transform=ax.transAxes,
                color=rgb)
        self.canvas.draw()
        evt.Skip()"""

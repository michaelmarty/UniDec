# base class for all plotting windows used in unidec
# contains basic setup functionality

import wx
import tempfile
import os
from matplotlib import interactive
from matplotlib import rcParams
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

from unidec.modules.isolated_packages.ZoomSpan import ZoomSpan
from unidec.modules.isolated_packages.ZoomBox import ZoomBox
from unidec.modules.isolated_packages.NoZoomSpan import NoZoomSpan
from unidec.modules.isolated_packages import FileDialogs
from unidec.modules.miscwindows import DoubleInputDialog
from unidec.modules.PlotBase import PlotBase
from unidec.modules.plot1d import Plot1dBase
from unidec.modules.plot2d import Plot2dBase

# import matplotlib.style as mplstyle
# mplstyle.use('fast')

interactive(True)


class PlottingWindowBase(PlotBase, wx.Window):
    """
    Class for wx window with embedded matplotlib plots
    Inherits from PlotBase, which is the core functions stripped of everything interactive
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

        PlotBase.__init__(self, *args, **kwargs)
        self.displaysize = wx.GetDisplaySize()

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

        wx.Window.__init__(self, *args)
        self.Bind(wx.EVT_SIZE, self.size_handler)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('figure_enter_event', self.mouse_activate)
        self.canvas.mpl_connect('figure_leave_event', self.mouse_inactivate)

    def repaint(self, setupzoom=True):
        """
        Redraw and refresh the plot.
        :return: None
        """
        if setupzoom:
            try:
                self.setup_zoom([self.subplot1], self.zoomtype)
            except:
                pass
        try:
            self.zoomout()
        except:
            pass
        self.canvas.draw()

    def zoomout(self):
        self.zoom.zoomout()

    def mouse_activate(self, event):
        self.mouse_active = True

    def mouse_inactivate(self, event):
        self.mouse_active = False

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
        if event.button == 2 or (event.button == 1 and wx.GetKeyState(wx.WXK_ESCAPE)):
            if wx.GetKeyState(wx.WXK_CONTROL):
                dlg = DoubleInputDialog(self)
                dlg.initialize_interface("Matplotlib RC Parameters", "RC Param Name:", 'lines.markersize',
                                         "Value:", "6")
                dlg.ShowModal()
                rcname = dlg.value
                rcval = dlg.value2

                if rcname == "aspect":
                    self.set_aspect(rcval)
                    self.repaint(setupzoom=False)
                else:
                    print(rcname, rcval, "Replot to activate")
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

    def setup_zoom(self, plots, zoom, data_lims=None, pad=0, groups=None):
        """
        Set up zoom on axes.
        :param pad: Padding around plot
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
                zoombutton=3,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))


class Plot1d(PlottingWindowBase, Plot1dBase):
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
        Plot1dBase.__init__(self, *args, **kwargs)
        PlottingWindowBase.__init__(self, *args, **kwargs)


class Plot2d(PlottingWindowBase, Plot2dBase):
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
        Plot2dBase.__init__(self, *args, **kwargs)
        PlottingWindowBase.__init__(self, *args, **kwargs)
        self.is2d = True


class PlotAny(PlottingWindowBase, Plot1dBase, Plot2dBase):
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
        PlottingWindowBase.__init__(self, *args, **kwargs)
        Plot1dBase.__init__(self, *args, **kwargs)
        Plot2dBase.__init__(self, *args, **kwargs)


class NetworkFrame(PlottingWindowBase):
    def __init__(self, *args, **kwargs):
        PlottingWindowBase.__init__(self, *args, **kwargs)
        self.axes = self.figure.add_axes(self._axes)
        self.flag = True

    def clear_frame(self):
        self.figure.clear()
        self.axes = self.figure.add_axes(self._axes)
        self.repaint()

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
from pubsub import pub

# import matplotlib.style as mplstyle
# mplstyle.use('fast')

interactive(True)

# Create ScanSelected custom event
ScanSelectedEventType = wx.NewEventType()
EVT_SCANS_SELECTED = wx.PyEventBinder(ScanSelectedEventType, 1)


class ScanSelectedEvent(wx.PyCommandEvent):
    def __init__(self, evttype, id, smin=None, smax=None):
        wx.PyCommandEvent.__init__(self, evttype, id)
        self.smin = smin
        self.smax = smax


# Create mzlimits custom event
MZLimitsEventType = wx.NewEventType()
EVT_MZLIMITS = wx.PyEventBinder(MZLimitsEventType, 1)


class MZLimitsEvent(wx.PyCommandEvent):
    def __init__(self, evttype, id, minval=None, maxval=None):
        wx.PyCommandEvent.__init__(self, evttype, id)
        self.minval = minval
        self.maxval = maxval

# Create swoop drag custom event
SwoopDragEventType = wx.NewEventType()
EVT_SWOOP_DRAG = wx.PyEventBinder(SwoopDragEventType, 1)

class SwoopDragEvent(wx.PyCommandEvent):
    def __init__(self, evttype, id, sarray=None):
        wx.PyCommandEvent.__init__(self, evttype, id)
        self.sarray = sarray


class PlottingWindowBase(PlotBase, wx.Panel):
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
        self.EVT_SCANS_SELECTED = EVT_SCANS_SELECTED
        self.EVT_MZLIMITS = EVT_MZLIMITS
        self.EVT_SWOOP_DRAG = EVT_SWOOP_DRAG

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

        if "parent" in kwargs:
            self.parent = kwargs["parent"]
            del kwargs["parent"]
        else:
            self.parent = None

        wx.Window.__init__(self, *args)
        self.Bind(wx.EVT_SIZE, self.size_handler)
        # self.Bind(wx.EVT_IDLE, self.on_idle)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('figure_enter_event', self.mouse_activate)
        self.canvas.mpl_connect('figure_leave_event', self.mouse_inactivate)
        self.canvas.mpl_connect('draw_event', self.on_draw)

        # self.sizer = wx.BoxSizer(wx.VERTICAL)
        # self.sizer.Add(self.canvas, 1, wx.EXPAND)
        # self.SetSizer(self.sizer)
        # self.Show()
        # self.Layout()

    def repaint(self, setupzoom=True, resetzoom=False):
        """
        Redraw and refresh the plot.
        :param setupzoom: Boolean, whether to setup zooming
        :param resetzoom: Boolean, whether to reset zooming
        :return: None
        """
        if resetzoom:
            self.reset_zoom()

        if setupzoom:
            try:
                self.setup_zoom([self.subplot1], self.zoomtype)
            except:
                pass
        else:
            try:
                self.zoomout()
            except:
                pass
        if self.zoomvals is not None:
            self.set_zoomvals()
        self.canvas.draw()

    def zoomout(self):
        self.zoom.initialize()

    def mouse_activate(self, event):
        self.mouse_active = True

    def mouse_inactivate(self, event):
        self.mouse_active = False

    def on_draw(self, event):
        self.get_zoomvals()

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
        if event.button == 2 or (event.button == 1 and wx.GetKeyState(wx.WXK_DOWN)):
            # Middle Clicks
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

    def size_handler(self, event=None, *args, **kwargs):
        """
        Resizes the plots
        :param args:
        :param kwargs:
        :return: None
        """
        # print("Size Handler")
        # event.Skip()
        # wx.CallAfter(self.on_size, event)
        if self.resize == 1:
            self.canvas.SetSize(self.GetSize())

    def on_size(self, event):
        """
        Resizes the plots
        :param event: wx.Event
        :return: None
        """
        if self.resize == 1:
            self.resizeflag = True

    def on_idle(self, event):
        """
        Function triggered on idle event from plot.
        :param event: wx.Event
        :return: None
        """
        if self.resizeflag:
            # print(self.GetSize())
            # print(self.GetClientSize())
            if self.parent is not None:
                pass
                '''
                print("Client Size:", self.parent.GetClientSize())
                print("Parent Size:", self.parent.GetSize())
                psize = self.parent.GetClientSize()
                sizer = self.parent.GetSizer()
                cellsize = sizer.GetCellSize(0, 0)
                print("Cell Size:", cellsize)

                halfcolumn = [psize[0] / 2, psize[1] / 2]
                if cellsize[0] != 0:
                    newcellsize = (cellsize[0] * halfcolumn[0] / cellsize[0], cellsize[1] * halfcolumn[0] / cellsize[0])
                    print(newcellsize)
                    # self.canvas.resize(*newcellsize)
                    # self.canvas.SetSize(newcellsize)
                    # self.canvas.SetMinSize(newcellsize)
                    self.set_resize(newcellsize)'''

            # self.canvas.SetMinSize(self.GetSize())
            # self.SetMinSize(self.GetSize())

            # self.SetSize(self.GetSize())
            # self.canvas.SetClientSize(self.GetSize())
            # self.canvas.SetMinClientSize(self.GetSize())
            # print(self.GetClientSize())
            # self.canvas.resize(self.GetSize())
            # self.canvas.draw()
            # self.SetClientSize(self.GetSize())
            self.resizeflag = False

    def set_resize(self, newsize):
        self.SetSize(newsize)
        # self.canvas.Destroy()
        # self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.canvas.SetSize(newsize)
        self.figure.set_size_inches(float(newsize[0]) / self.figure.get_dpi(),
                                    float(newsize[1]) / self.figure.get_dpi())
        self.canvas.draw()
        self.Layout()
        self.Refresh()

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
                parent=self,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))
        if zoom == 'box':
            self.zoom = ZoomBox(
                plots,
                None,
                parent=self,
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
                parent=self,
                minspan=0,
                zoombutton=3,
                useblit=True,
                onmove_callback=None,
                rectprops=dict(alpha=0.2, facecolor='yellow'))

    def on_newxy(self, x, y):
        # print(x, y)
        pub.sendMessage('newxy', xpos=x, ypos=y)

    def on_scans_selected(self, smin, smax):
        # pub.sendMessage('scans_selected', min=smin, max=smax)
        event = ScanSelectedEvent(ScanSelectedEventType, self.GetId(), smin=smin, smax=smax)
        self.GetEventHandler().ProcessEvent(event)

    def on_left_click(self, x, y):
        pub.sendMessage('left_click', xpos=x, ypos=y)

    def on_right_click(self, event=None):
        if self.int == 1:
            pub.sendMessage('integrate')
        elif self.smash == 1:
            if event.dblclick:
                pub.sendMessage('smash')
            else:
                pub.sendMessage('mzlimits')
                event = ScanSelectedEvent(MZLimitsEventType, self.GetId())
                self.GetEventHandler().ProcessEvent(event)
        elif self.smash == 2:
            event = ScanSelectedEvent(MZLimitsEventType, self.GetId())
            self.GetEventHandler().ProcessEvent(event)
        else:
            event = ScanSelectedEvent(ScanSelectedEventType, self.GetId())
            self.GetEventHandler().ProcessEvent(event)

    def on_swoop_drag(self, sarray):
        event = SwoopDragEvent(SwoopDragEventType, self.GetId(), sarray=sarray)
        self.GetEventHandler().ProcessEvent(event)

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

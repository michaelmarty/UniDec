from unidec.modules.gui_elements.mainwindow_base import MainwindowBase
import wx
import wx.lib.scrolledpanel as scrolled
import os
from unidec.modules.plotting import PlottingWindow
from unidec.modules.gui_elements import peaklistsort
from unidec.modules.gui_elements import IsoDecControls
from unidec.modules.gui_elements import IsoDecMenu
from isodec.fragment_view import plot_fragment_matches
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
import math

class IsoDecView(MainwindowBase):
    def __init__(self, parent, title, config, iconfile=None, tabbed=None):
        super().__init__(parent, title, config, iconfile, tabbed)
        self.parent = parent

        self.menu = IsoDecMenu.main_menu(self, self.config, self.pres)
        self.SetMenuBar(self.menu.menuBar)

        self.setup_main_panel()

        keys = [["O", self.pres.on_open, self.menu.menuOpen],
                ["G", self.pres.on_paste_spectrum, self.menu.menupastespectrum],
                ["D", self.pres.on_dataprep_button, self.controls.dataprepbutton],
                ["N", self.pres.on_replot, self.controls.replotbutton],
                ["Q", self.on_exit, self.menu.menuExit],
                ["E", self.pres.on_run_all, self.controls.runallbutton]
                ]

        keys = keys + self.menu.menukeys
        self.setup_shortcuts(keys)
        self.import_config_to_gui()


    def setup_main_panel(self):
        # Create Status Bar
        statusbar_log_silencer = wx.LogNull()
        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, -6, -4, -4, -5, -3, -3])
        del statusbar_log_silencer
        # Sizers to develop layout
        # s1 = (min(self.displaysize[0], 1851), self.displaysize[1])
        # s2 = (550, self.displaysize[1])
        self.splitterwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D | wx.SP_BORDER)
        splitterwindow2 = wx.SplitterWindow(self.splitterwindow, -1, style=wx.SP_3D | wx.SP_BORDER)
        # self.splitterwindow.SetSashGravity(0)
        # splitterwindow2.SetSashGravity(0.5)
        panelp = wx.Panel(splitterwindow2, -1)
        panel = scrolled.ScrolledPanel(splitterwindow2, -1)  # wx.Panel(splitterwindow2, -1)
        splitterwindow2.SplitVertically(panelp, panel, sashPosition=-270)

        file_drop_target = MyFileDropTarget(self)
        self.splitterwindow.SetDropTarget(file_drop_target)
        # .................................
        #
        #    Layout the Plots
        #
        # ...................................

        # Scrolled panel view of plots

        # TODO: Line up plots on left hand side so that they share an m/z axis
        plotwindow = scrolled.ScrolledPanel(self.splitterwindow)
        self.splitterwindow.SplitVertically(plotwindow, splitterwindow2, sashPosition=-550)
        self.sizerplot = wx.GridBagSizer()
        figsize = self.config.figsize
        self.plot1 = PlottingWindow.Plot1d(plotwindow, smash=1, figsize=figsize, parent=plotwindow)
        self.plot2 = PlottingWindow.Plot1d(plotwindow, integrate=1, figsize=figsize, parent=plotwindow)
        self.fragment_panel = wx.Panel(plotwindow)
        self.fragment_figure = Figure(figsize=(12, 3))
        self.fragment_ax = self.fragment_figure.add_subplot(111)
        self.fragment_canvas = FigureCanvasWxAgg(self.fragment_panel, -1, self.fragment_figure)
        self.fragment_has_matches = False
        fragment_sizer = wx.BoxSizer(wx.VERTICAL)
        fragment_sizer.Add(self.fragment_canvas, 1, wx.EXPAND)
        self.fragment_panel.SetSizer(fragment_sizer)
        self.fragment_panel.SetMinSize((1200, 250))
        self.clear_fragment_plot()

        self.sizerplot.Add(self.plot1, (0, 0), span=(1, 1), flag=wx.EXPAND)
        self.sizerplot.Add(self.plot2, (0, 1), span=(1, 1), flag=wx.EXPAND)
        self.sizerplot.Add(self.fragment_panel, (1, 0), span=(1, 2), flag=wx.EXPAND)

        # plotwindow.SetScrollbars(1, 1,1,1)
        if self.system == "Linux":
            plotwindow.SetSizer(self.sizerplot)
            self.sizerplot.Fit(self)
        else:
            plotwindow.SetSizerAndFit(self.sizerplot)
        plotwindow.SetupScrolling()
        plotwindow.SetFocus()
        plotwindow.Bind(wx.EVT_SET_FOCUS, self.onFocus)
        self.plotpanel = plotwindow

        self.plots = [self.plot1, self.plot2]
        self.plotnames = ["Figure1", "Figure2"]


        # ...........................
        #
        #   Sizer for Peaks
        #
        # ...........................
        sizerpeaks = wx.BoxSizer(wx.VERTICAL)
        self.peakpanel = peaklistsort.PeakListCtrlPanel(panelp, size=(300, 600), isodec=True)
        self.bind_peakpanel()
        sizerpeaks.Add(self.peakpanel, 0, wx.EXPAND)
        panelp.SetSizer(sizerpeaks)
        sizerpeaks.Fit(self)

        # ..........................
        #
        # Setup Control Panel
        #
        # .............................
        sizercontrols = wx.BoxSizer(wx.VERTICAL)
        self.controls = IsoDecControls.MainControls(self, self.config, self.pres, panel)
        sizercontrols.Add(self.controls, 1, wx.EXPAND)
        panel.SetSizer(sizercontrols)
        sizercontrols.Fit(self)

        splitterwindow2.SetMinimumPaneSize(20)
        self.splitterwindow.SetMinimumPaneSize(20)
        # self.splitterwindow.SetMinSize((0,0))
        # splitterwindow2.SetMinSize((0,0))

        if self.system == "Linux":
            self.sizerplot.Fit(self.splitterwindow)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.splitterwindow, 1, wx.EXPAND)

        # Set everything up
        self.SetSizer(sizer)
        sizer.Fit(self)
        self.SetSize((self.GetSize().width, self.GetSize().height + 100))

        self.Layout()

        self.plotpanel.SetMinSize(wx.Size(-1, -1))
        self.plotpanel.Bind(wx.EVT_SIZE, self.resize_plots)

        self.splitterwindow.SetMinimumPaneSize(20)
        self.splitterwindow.SetSashGravity(0.99)

        splitterwindow2.SetMinimumPaneSize(20)
        splitterwindow2.SetSashGravity(0.5)

    def clear_fragment_plot(self):
        self.fragment_has_matches = False
        self.fragment_ax.clear()
        self.fragment_ax.set_axis_off()
        self.fragment_canvas.draw_idle()

    def show_fragment_matches(self, sequence, pks):
        residues_per_line = 70
        lines = math.ceil((len(pks.fragment_matches.index) + 1) / residues_per_line)
        ion_count = sum(column.endswith("_match") and pks.fragment_matches[column].notna().any()
                        for column in pks.fragment_matches)
        height = max(2.5, 0.35 + lines * (0.22 + 0.07 * ion_count))
        self.fragment_figure.set_size_inches(12, height)
        self.fragment_panel.SetMinSize((1200, int(height * self.fragment_figure.dpi)))
        plot_fragment_matches(self.fragment_ax, sequence, pks, residues_per_line=residues_per_line)
        self.fragment_has_matches = True
        self.fragment_figure.tight_layout(pad=0.3)
        self.fragment_canvas.draw_idle()
        self.sizerplot.Layout()
        self.plotpanel.SetupScrolling()

    def clear_all_plots(self, flag=0):
        super().clear_all_plots(flag)
        self.clear_fragment_plot()

    def save_all_figures(self, extension, extension2='', e=0, header=None, **kwargs):
        flags, files = super().save_all_figures(extension, extension2, e, header, **kwargs)
        if self.fragment_has_matches:
            path = f"{header or self.config.outfname}{extension2}_Figure3.{extension}"
            self.fragment_figure.savefig(path, **kwargs)
            flags.append(3)
            files.append([3, path])
        return flags, files


class MyFileDropTarget(wx.FileDropTarget):
    """"""

    def __init__(self, window):
        """Constructor"""
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, filenames):
        """
        When files are dropped, either open a single file or run in batch.
        """
        if len(filenames) == 1:
            # Open a single file
            path = filenames[0]
            directory, fname = os.path.split(path)
            if os.path.splitext(fname)[1] == ".raw" and os.path.isdir(path):
                print("Opening .raw file:", fname)
                self.window.pres.on_raw_open(0, path)
            elif fname[-9:] == "_conf.dat":
                print("Importing Configuration File:", path)
                self.window.pres.import_config(path)
            elif os.path.splitext(fname)[1] == ".zip":
                print("Loading State:", fname)
                self.window.pres.on_load_state(0, path)
            else:
                self.window.pres.on_open_file(fname, directory)
                # self.window.pres.on_auto() # Run the whole thing
        else:
            print("Error in file drop", filenames)
        return 0

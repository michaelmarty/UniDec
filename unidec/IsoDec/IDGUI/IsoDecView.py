from unidec.modules.gui_elements.mainwindow_base import MainwindowBase
import wx
import wx.lib.scrolledpanel as scrolled
import os
from modules.plotting import PlottingWindow
from unidec.modules.gui_elements import peaklistsort
from unidec.IsoDec.IDGUI import IsoDecControls
from unidec.IsoDec.IDGUI import IsoDecMenu

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
                ]

        keys = keys + self.menu.menukeys
        self.setup_shortcuts(keys)
        self.import_config_to_gui()


    def setup_main_panel(self):
        # Create Status Bar
        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, 300, 200, 200, 250, 150, 130])
        # Sizers to develop layout
        # s1 = (min(self.displaysize[0], 1851), self.displaysize[1])
        # s2 = (550, self.displaysize[1])
        self.splitterwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D | wx.SP_BORDER)
        splitterwindow2 = wx.SplitterWindow(self.splitterwindow, -1, style=wx.SP_3D | wx.SP_BORDER)
        self.splitterwindow.SetSashGravity(0)
        splitterwindow2.SetSashGravity(0.5)
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

        self.sizerplot.Add(self.plot1, (0, 0), span=(1, 1), flag=wx.EXPAND)
        self.sizerplot.Add(self.plot2, (0, 1), span=(1, 1), flag=wx.EXPAND)

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
        self.controls = IsoDecControls.main_controls(self, self.config, self.pres, panel, self.icon_path)
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
        sizer.Add(self.splitterwindow, 0, wx.EXPAND)

        # Set everything up
        self.SetSizer(sizer)
        sizer.Fit(self)

        self.Layout()


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
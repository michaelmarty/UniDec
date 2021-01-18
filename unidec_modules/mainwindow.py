import os

import numpy as np
import wx
import wx.lib.scrolledpanel as scrolled

from unidec_modules.gui_elements import ud_controls
from unidec_modules.gui_elements import ud_menu
from unidec_modules import ColorPlot, plot3d, plot1d, plot2d, miscwindows
from unidec_modules.gui_elements import peaklistsort
from unidec_modules.gui_elements.mainwindow_base import MainwindowBase

__author__ = 'Michael.Marty'


# noinspection PyAttributeOutsideInit,PyUnusedLocal,PyUnusedLocal
class Mainwindow(MainwindowBase):
    """
    Main UniDec GUI Window.
    """

    def __init__(self, parent, title, config, iconfile="logo.ico", tabbed=None):
        """
        initialize window and feed in links to presenter and config.

        :param parent: GUniDec Presenter -> self.pres
        :param title: Window title (string)
        :param config: UniDecConfig object ->self.config
        :return: None
        """
        MainwindowBase.__init__(self, parent, title, config, iconfile, tabbed)

        if tabbed is None:
            # If tabbed isn't specified, use the display size to decide what is best
            print("Display Size ", self.displaysize)
            self.tabbed = 0
            if self.displaysize[0] < 1400:
                self.tabbed = 1
            elif 1500 <= self.displaysize[0] < 1800:
                self.config.figsize = (5, 4)
            elif 1400 <= self.displaysize[0] < 1500:
                self.config.figsize = (4, 3)
            elif 1800 <= self.displaysize[0] < 2200:
                self.config.figsize = (6.5, 4.5)
            elif 2200 <= self.displaysize[0] < 3600:
                self.config.figsize = (9, 6)
            elif self.displaysize[0] >= 3600:
                self.config.figsize = (12, 6)
        else:
            self.tabbed = tabbed

        self.imflag = self.config.imflag

        self.twave = self.config.twaveflag > 0

        self.menu = ud_menu.main_menu(self, self.config, self.pres, self.tabbed)
        self.SetMenuBar(self.menu.menuBar)

        self.setup_main_panel()

        keys = [["O", self.pres.on_open, self.menu.menuOpen],
                ["G", self.pres.on_paste_spectrum, self.menu.menupastespectrum],
                ["I", self.pres.on_integrate, self.menu.menuintegrate],
                ["E", self.pres.on_auto, self.controls.autobutton],
                ["R", self.pres.on_unidec_button, self.controls.udbutton],
                ["D", self.pres.on_dataprep_button, self.controls.dataprepbutton],
                ["P", self.pres.on_pick_peaks, self.controls.ppbutton],
                ["K", self.pres.on_kendrick, self.menu.menukendrick],
                #["C", self.pres.on_plot_composite, self.controls.compositebutton],
                ["N", self.pres.on_replot, self.controls.replotbutton],
                ["F", self.pres.on_plot_offsets, self.menu.menuoffset],  # ["Z", self.pres.on_charge_plot],
                ["L", self.pres.on_load_everything, self.menu.menuLoadEverything],
                ["S", self.pres.on_save_state, self.menu.menuSaveState],
                ["B", self.pres.on_batch, self.menu.menuBatch],
                ["Q", self.on_exit, self.menu.menuExit],
                ["T", self.pres.on_mass_tools, self.menu.menuMassFile],
                ["M", self.pres.on_match, self.menu.menumatch],
                ["W", self.pres.on_auto_peak_width, self.menu.menuAutoWidth],
                ["Z", self.pres.on_undo, self.menu.menuundo],
                ["Y", self.pres.on_redo, self.menu.menuredo],
                ]

        keys = keys + self.menu.menukeys
        self.setup_shortcuts(keys)
        self.launch()
        pass

    # noinspection PyPep8,PyPep8
    def setup_main_panel(self):
        """
        Lays Out Main Panel. Binds some functions to presenter.
        :return: None
        """
        # Create Status Bar
        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, 300, 200, 200, 250, 150, 130])
        # Sizers to develop layout
        # s1 = (min(self.displaysize[0], 1851), self.displaysize[1])
        # s2 = (550, self.displaysize[1])
        splitterwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D | wx.SP_BORDER)
        splitterwindow2 = wx.SplitterWindow(splitterwindow, -1, style=wx.SP_3D | wx.SP_BORDER)
        panelp = wx.Panel(splitterwindow2, -1)
        panel = scrolled.ScrolledPanel(splitterwindow2, -1)  # wx.Panel(splitterwindow2, -1)
        splitterwindow2.SplitVertically(panelp, panel, sashPosition=-270 - self.config.imflag * 20)
        splitterwindow2.SetMinimumPaneSize(270)
        splitterwindow.SetMinimumPaneSize(250)
        # splitterwindow.SetMinSize((0,0))
        # splitterwindow2.SetMinSize((0,0))
        file_drop_target = MyFileDropTarget(self)
        splitterwindow.SetDropTarget(file_drop_target)
        # .................................
        #
        #    Layout the Plots
        #
        # ...................................

        # Tabbed view of plots
        if self.tabbed == 1:
            figsize = self.config.figsize
            plotwindow = wx.Notebook(splitterwindow)
            splitterwindow.SplitVertically(plotwindow, splitterwindow2, sashPosition=-550)
            tab1 = wx.Panel(plotwindow)
            tab2 = wx.Panel(plotwindow)
            tab3 = wx.Panel(plotwindow)
            tab4 = wx.Panel(plotwindow)
            tab5 = wx.Panel(plotwindow)
            tab6 = wx.Panel(plotwindow)

            self.plot1 = plot1d.Plot1d(tab1, smash=1, figsize=figsize)
            self.plot2 = plot1d.Plot1d(tab2, integrate=1, figsize=figsize)
            self.plot3 = plot2d.Plot2d(tab3, figsize=figsize)
            self.plot4 = plot1d.Plot1d(tab4, figsize=figsize)
            self.plot5 = plot2d.Plot2d(tab5, figsize=figsize)
            self.plot6 = plot1d.Plot1d(tab6, figsize=figsize)

            miscwindows.setup_tab_box(tab1, self.plot1)
            miscwindows.setup_tab_box(tab2, self.plot2)
            miscwindows.setup_tab_box(tab3, self.plot3)
            miscwindows.setup_tab_box(tab4, self.plot4)
            miscwindows.setup_tab_box(tab5, self.plot5)
            miscwindows.setup_tab_box(tab6, self.plot6)

            if self.config.imflag == 1:
                tab1im = wx.Panel(plotwindow)
                tab1fit = wx.Panel(plotwindow)
                tab2ccs = wx.Panel(plotwindow)
                tab3color = wx.Panel(plotwindow)
                tab5mccs = wx.Panel(plotwindow)
                tab5ccsz = wx.Panel(plotwindow)
                tab9 = wx.Panel(plotwindow)
                tab10 = wx.Panel(plotwindow)

                self.plot1im = plot2d.Plot2d(tab1im, figsize=figsize)
                self.plot1fit = plot2d.Plot2d(tab1fit, figsize=figsize)
                self.plot2ccs = plot1d.Plot1d(tab2ccs, figsize=figsize)
                self.plot5mccs = plot2d.Plot2d(tab5mccs, figsize=figsize)
                self.plot5ccsz = plot2d.Plot2d(tab5ccsz, figsize=figsize)
                self.plot3color = ColorPlot.ColorPlot2D(tab3color, figsize=figsize)
                self.plot9 = plot3d.CubePlot(tab9, figsize=figsize)
                self.plot10 = plot3d.CubePlot(tab10, figsize=figsize)

                miscwindows.setup_tab_box(tab1im, self.plot1im)
                miscwindows.setup_tab_box(tab1fit, self.plot1fit)
                miscwindows.setup_tab_box(tab2ccs, self.plot2ccs)
                miscwindows.setup_tab_box(tab3color, self.plot3color)
                miscwindows.setup_tab_box(tab5mccs, self.plot5mccs)
                miscwindows.setup_tab_box(tab5ccsz, self.plot5ccsz)
                miscwindows.setup_tab_box(tab9, self.plot9)
                miscwindows.setup_tab_box(tab10, self.plot10)

            plotwindow.AddPage(tab1, "MS Data v. Fit")
            if self.config.imflag == 1:
                plotwindow.AddPage(tab1im, "IM-MS Data")
                plotwindow.AddPage(tab1fit, "IM-MS Fit")
                plotwindow.AddPage(tab3color, "IM-MS Charges")
                plotwindow.AddPage(tab9, "m/z Cube")
            plotwindow.AddPage(tab3, "m/z Grid")
            plotwindow.AddPage(tab2, "Mass Distribution")
            if self.config.imflag == 1:
                plotwindow.AddPage(tab2ccs, "CCS Distribution")

            plotwindow.AddPage(tab4, "Individual Peaks")
            plotwindow.AddPage(tab5, "Mass vs. Charge")
            if self.config.imflag == 1:
                plotwindow.AddPage(tab5mccs, "Mass vs. CCS ")
                plotwindow.AddPage(tab5ccsz, "CCS vs. Charge")
                plotwindow.AddPage(tab10, "Mass Cube")
            plotwindow.AddPage(tab6, "Bar Chart")
        # Scrolled panel view of plots
        else:
            # TODO: Line up plots on left hand side so that they share an m/z axis
            plotwindow = scrolled.ScrolledPanel(splitterwindow)
            splitterwindow.SplitVertically(plotwindow, splitterwindow2, sashPosition=-550)
            sizerplot = wx.GridBagSizer()
            figsize = self.config.figsize
            self.plot1 = plot1d.Plot1d(plotwindow, smash=1, figsize=figsize)
            self.plot2 = plot1d.Plot1d(plotwindow, integrate=1, figsize=figsize)
            self.plot3 = plot2d.Plot2d(plotwindow, figsize=figsize)
            self.plot4 = plot1d.Plot1d(plotwindow, figsize=figsize)
            self.plot5 = plot2d.Plot2d(plotwindow, figsize=figsize)
            self.plot6 = plot1d.Plot1d(plotwindow, figsize=figsize)

            if self.config.imflag == 1:
                self.plot1im = plot2d.Plot2d(plotwindow, figsize=figsize)
                self.plot1fit = plot2d.Plot2d(plotwindow, figsize=figsize)
                self.plot2ccs = plot1d.Plot1d(plotwindow, figsize=figsize)
                self.plot5mccs = plot2d.Plot2d(plotwindow, figsize=figsize)
                self.plot5ccsz = plot2d.Plot2d(plotwindow, figsize=figsize)
                self.plot3color = ColorPlot.ColorPlot2D(plotwindow, figsize=figsize)
                self.plot9 = plot3d.CubePlot(plotwindow, figsize=figsize)
                self.plot10 = plot3d.CubePlot(plotwindow, figsize=figsize)

            if self.config.imflag == 0:
                sizerplot.Add(self.plot1, (0, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot2, (0, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot3, (1, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot4, (1, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot5, (2, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot6, (2, 1), span=(1, 1), flag=wx.EXPAND)
            else:
                sizerplot.Add(self.plot1, (0, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot1im, (0, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot1fit, (1, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot3color, (1, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot2, (2, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot2ccs, (3, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot3, (2, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot4, (4, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot5, (3, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot5mccs, (4, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot5ccsz, (5, 1), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot6, (5, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot9, (6, 0), span=(1, 1), flag=wx.EXPAND)
                sizerplot.Add(self.plot10, (6, 1), span=(1, 1), flag=wx.EXPAND)

            # plotwindow.SetScrollbars(1, 1,1,1)
            if self.system == "Linux":
                plotwindow.SetSizer(sizerplot)
                sizerplot.Fit(self)
            else:
                plotwindow.SetSizerAndFit(sizerplot)
            plotwindow.SetupScrolling()
            plotwindow.SetFocus()
            plotwindow.Bind(wx.EVT_SET_FOCUS, self.onFocus)
        self.plotpanel = plotwindow

        self.plots = [self.plot1, self.plot2, self.plot3, self.plot4, self.plot5, self.plot6]
        if self.config.imflag == 1:
            self.plots = self.plots + [self.plot1im, self.plot1fit, self.plot2ccs, self.plot5mccs, self.plot5ccsz,
                                       self.plot3color, self.plot9, self.plot10]
        self.plotnames = ["Figure1", "Figure2", "Figure3", "Figure4", "Figure5", "Figure6"]
        if self.config.imflag == 1:
            self.plotname = self.plotnames + ["Figure1im", "Figure1fit", "Figure2ccs", "Figure3color", "Figure5ccsz",
                                              "Figure5massccs", "mzCube", "massCube"]

        # ...........................
        #
        #   Sizer for Peaks
        #
        # ...........................
        sizerpeaks = wx.BoxSizer(wx.VERTICAL)
        self.peakpanel = peaklistsort.PeakListCtrlPanel(panelp)
        self.Bind(self.peakpanel.EVT_DELETE_SELECTION_2, self.pres.on_delete, self.peakpanel)
        self.Bind(self.peakpanel.EVT_CHARGE_STATE, self.pres.on_charge_states, self.peakpanel)
        self.Bind(self.peakpanel.EVT_DIFFERENCES, self.pres.on_differences, self.peakpanel)
        self.Bind(self.peakpanel.EVT_MASSES, self.pres.on_label_masses, self.peakpanel)
        sizerpeaks.Add(self.peakpanel, 0, wx.EXPAND)
        panelp.SetSizer(sizerpeaks)
        sizerpeaks.Fit(self)

        # ..........................
        #
        # Setup Control Panel
        #
        # .............................
        sizercontrols = wx.BoxSizer(wx.VERTICAL)
        self.controls = ud_controls.main_controls(self, self.config, self.pres, panel, self.icon_path)
        sizercontrols.Add(self.controls, 1, wx.EXPAND)
        panel.SetSizer(sizercontrols)
        sizercontrols.Fit(self)

        if self.system == "Linux" and self.tabbed != 1:
            sizerplot.Fit(splitterwindow)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(splitterwindow, 1, wx.EXPAND)

        # Set everything up
        self.SetSizer(sizer)
        sizer.Fit(self)
        pass


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

        elif len(filenames) > 1:
            # Batch process the files that were dropped
            if os.path.splitext(filenames[0])[1] == ".raw" and os.path.isdir(filenames[0]):
                print("Batch converting raw to txt")
                self.window.pres.on_batch_raw(0, filenames, clip=False)
            else:
                print("Running batch mode")
                self.window.pres.on_batch(batchfiles=filenames)
        else:
            print("Error in file drop", filenames)
        return 0

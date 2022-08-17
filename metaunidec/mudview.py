import os

import numpy as np
import wx
import wx.lib.scrolledpanel as scrolled
from matplotlib import rcParams

from unidec_modules import plot1d, plot2d, miscwindows
from unidec_modules.gui_elements import peaklistsort, mainwindow_base
from metaunidec.gui_elements.list_ctrls import ListCtrlPanel
from metaunidec.gui_elements.ud_cont_meta import main_controls
from metaunidec.gui_elements.ud_menu_meta import meta_menu

__author__ = 'Michael.Marty'

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42


# noinspection PyAttributeOutsideInit,PyUnusedLocal,PyUnusedLocal
class Mainwindow(mainwindow_base.MainwindowBase):
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
        mainwindow_base.MainwindowBase.__init__(self, parent, title, config, iconfile, tabbed)

        self.config.figsize = (5, 4)
        if tabbed is None:
            # If tabbed isn't specified, use the display size to decide what is best
            print("Display Size ", self.displaysize)
            self.tabbed = 0

            if self.displaysize[0] < 1600:
                self.tabbed = 1
            elif 1600 <= self.displaysize[0] < 1800:
                self.config.figsize = (5, 4)
            elif 1800 <= self.displaysize[0] < 2100:
                self.config.figsize = (5, 4)
            elif 2100 <= self.displaysize[0] < 3600:
                self.config.figsize = (8, 6)
            elif self.displaysize[0] >= 3600:
                self.config.figsize = (12, 6)

        else:
            self.tabbed = tabbed

        self.datachoices = {0: "Raw Data", 1: "Processed Data", 2: "Zero Charge Mass Spectrum"}
        self.extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position",
                               5: "Estimated Area"}
        self.extractlabels = {0: "Intensity", 1: "Intensity", 2: "Area", 3: "Mass", 4: "Mass", 5: "Area"}

        self.setupmainpanel()

        self.menu = meta_menu(self, self.config, self.pres)
        self.SetMenuBar(self.menu.menuBar)

        keys = [["E", self.pres.on_auto, self.controls.autobutton],
                ["G", self.pres.on_paste_spectrum, self.menu.menupastespectrum],
                ["R", self.pres.on_unidec_button, self.controls.udbutton],
                ["D", self.pres.on_dataprep_button, self.controls.dataprepbutton],
                ["O", self.pres.on_open, self.menu.openmenu],  # ["I", self.pres.on_integrate],
                ["P", self.pres.on_pick_peaks, self.controls.plotbutton],  # ["K", self.pres.on_plot_peaks],
                # ["C", self.pres.on_plot_composite, self.controls.compositebutton],
                ["N", self.pres.on_wizard, self.menu.wizardmenu],
                # ["F", self.pres.on_plot_offsets],  # ["Z", self.pres.on_charge_plot],
                # ["L", self.pres.on_load_state], ["S", self.pres.on_save_state],
                ["B", self.pres.on_batch_run, self.menu.menubatchrun],
                ["Q", self.on_exit, self.menu.menuExit],
                ["T", self.pres.on_mass_tools, self.menu.menuMassFile],
                ["M", self.pres.on_match, self.menu.menumatch],
                ["W", self.pres.on_auto_peak_width, self.menu.menuAutoWidth],
                ["Z", self.pres.on_undo, self.menu.menuundo],
                # ["Y", self.pres.on_redo, self.menu.menuredo],
                ["K", self.pres.on_kendrick, self.menu.menukendrick]
                ]
        self.setup_shortcuts(keys)

        self.launch()
        pass

    def setupmainpanel(self):
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel = wx.Panel(self)

        file_drop_target = MyFileDropTarget(self)
        panel.SetDropTarget(file_drop_target)

        self.ypanel = ListCtrlPanel(panel, self.pres, size=(300, 300))
        sizer.Add(self.ypanel, 0, wx.EXPAND)

        # Tabbed view of plots
        if self.tabbed == 1:
            figsize = (6, 5)
            plotwindow = wx.Notebook(panel)
            tab1 = wx.Panel(plotwindow)
            tab2 = wx.Panel(plotwindow)
            tab3 = wx.Panel(plotwindow)
            tab5 = wx.Panel(plotwindow)
            tab6 = wx.Panel(plotwindow)
            tab7 = wx.Panel(plotwindow)
            tab8 = wx.Panel(plotwindow)
            tab9 = wx.Panel(plotwindow)

            self.plot1 = plot1d.Plot1d(tab1, smash=1, figsize=figsize)
            self.plot2 = plot1d.Plot1d(tab2, integrate=1, figsize=figsize)
            self.plot3 = plot2d.Plot2d(tab3, figsize=figsize)
            self.plot5 = plot2d.Plot2d(tab5, figsize=figsize)
            self.plot6 = plot1d.Plot1d(tab6, figsize=figsize)
            self.plot7 = plot1d.Plot1d(tab7, figsize=figsize)
            self.plot8 = plot2d.Plot2d(tab8, figsize=figsize)
            self.plot9 = plot1d.Plot1d(tab9, figsize=figsize)

            miscwindows.setup_tab_box(tab1, self.plot1)
            miscwindows.setup_tab_box(tab2, self.plot2)
            miscwindows.setup_tab_box(tab3, self.plot3)
            miscwindows.setup_tab_box(tab5, self.plot5)
            miscwindows.setup_tab_box(tab6, self.plot6)
            miscwindows.setup_tab_box(tab7, self.plot7)
            miscwindows.setup_tab_box(tab8, self.plot8)
            miscwindows.setup_tab_box(tab9, self.plot9)

            plotwindow.AddPage(tab1, "MS Data")
            plotwindow.AddPage(tab9, "Charge Distributions")
            plotwindow.AddPage(tab2, "Mass Distribution")
            plotwindow.AddPage(tab7, "Extracts Line Plot")
            plotwindow.AddPage(tab8, "Extracts Grid Plot")
            plotwindow.AddPage(tab6, "Bar Chart")
            plotwindow.AddPage(tab3, "m/z Grid")
            plotwindow.AddPage(tab5, "Mass vs. Charge")
            sizer.Add(plotwindow, 1, wx.EXPAND)
        else:
            self.plotpanel = scrolled.ScrolledPanel(panel)
            sizerplot = wx.GridBagSizer()
            figsize = self.config.figsize
            self.plot1 = plot1d.Plot1d(self.plotpanel, smash=1, figsize=figsize)
            self.plot2 = plot1d.Plot1d(self.plotpanel, integrate=1, figsize=figsize)
            self.plot3 = plot2d.Plot2d(self.plotpanel, figsize=figsize)
            self.plot5 = plot2d.Plot2d(self.plotpanel, figsize=figsize)
            self.plot6 = plot1d.Plot1d(self.plotpanel, figsize=figsize)
            self.plot7 = plot1d.Plot1d(self.plotpanel, figsize=figsize)
            self.plot8 = plot2d.Plot2d(self.plotpanel, figsize=figsize)
            self.plot9 = plot1d.Plot1d(self.plotpanel, figsize=figsize)
            sizerplot.Add(self.plot1, (0, 0), span=(1, 1), flag=wx.EXPAND)
            sizerplot.Add(self.plot9, (0, 1), span=(1, 1), flag=wx.EXPAND)
            sizerplot.Add(self.plot2, (1, 0), span=(1, 2), flag=wx.EXPAND)
            sizerplot.Add(self.plot6, (3, 0), span=(1, 2), flag=wx.EXPAND)
            sizerplot.Add(self.plot7, (2, 0), span=(1, 1), flag=wx.EXPAND)
            sizerplot.Add(self.plot8, (2, 1), span=(1, 1), flag=wx.EXPAND)
            sizerplot.Add(self.plot3, (4, 0), span=(1, 1), flag=wx.EXPAND)
            sizerplot.Add(self.plot5, (4, 1), span=(1, 1), flag=wx.EXPAND)

            self.plotpanel.SetSizerAndFit(sizerplot)
            self.plotpanel.SetupScrolling()
            sizer.Add(self.plotpanel, 1, wx.EXPAND)

        self.plots = [self.plot1, self.plot2, self.plot6, self.plot7, self.plot8, self.plot9, self.plot3, self.plot5]
        self.plotnames = ["MetaFigure1", "MetaFigure2", "MetaFigure3", "MetaFigure4", "MetaFigure5", "MetaFigure6",
                          "MeatFigure7", "MetaFigure8"]

        self.peakpanel = peaklistsort.PeakListCtrlPanel(panel, meta=True)
        self.Bind(self.peakpanel.EVT_DELETE_SELECTION_2, self.pres.on_delete, self.peakpanel)
        self.Bind(self.peakpanel.EVT_CHARGE_STATE, self.pres.on_charge_states_mud, self.peakpanel)
        self.Bind(self.peakpanel.EVT_DIFFERENCES, self.pres.on_differences, self.peakpanel)
        self.Bind(self.peakpanel.EVT_MASSES, self.pres.on_label_masses, self.peakpanel)
        self.Bind(self.peakpanel.EVT_IMAGE, self.pres.make_image_plot, self.peakpanel)
        sizer.Add(self.peakpanel, 0, wx.EXPAND)

        self.controls = main_controls(self, self.config, self.pres, panel, self.icon_path)
        sizer.Add(self.controls, 0, wx.EXPAND)

        # Set everything up
        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, 600, 120, 0, 230, 250, 130])

    def export_gui_to_config(self):
        self.controls.export_gui_to_config()
        self.pres.export_vars()


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
        path = filenames[0]
        directory, fname = os.path.split(path)
        if os.path.splitext(fname)[1] == ".hdf5":
            print("Opening .hdf5 file:", fname)
            self.window.pres.open_file(path)
        elif os.path.splitext(fname)[1] == ".csv":
            self.window.pres.open_csv(path)
        else:
            self.window.pres.add_files(filenames)
        return 0

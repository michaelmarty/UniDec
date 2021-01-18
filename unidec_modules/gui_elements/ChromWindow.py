import os
import numpy as np
import wx
from metaunidec.gui_elements.list_ctrls import ListCtrlPanel
from metaunidec.gui_elements.ud_cont_meta import main_controls
from metaunidec.gui_elements.ud_menu_meta import meta_menu
from unidec_modules.gui_elements import peaklistsort, mainwindow_base
import wx.lib.scrolledpanel as scrolled
from pubsub import pub
from unidec_modules import plot1d, plot2d


class ChromWindow(mainwindow_base.MainwindowBase):
    def __init__(self, parent, title, config=None, iconfile="logo.ico", *args, **kwargs):
        mainwindow_base.MainwindowBase.__init__(self, parent, title, config, iconfile)
        # wx.Frame.__init__(self, None, title=title)  # ,size=(200,-1))
        self.pres = parent
        if config is None:
            self.config = self.pres.eng.config
        else:
            self.config = config
        self.icon_path = iconfile

        self.datachoices = {0: "Raw Data", 1: "Processed Data", 2: "Zero Charge Mass Spectrum"}
        self.extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position"}
        self.extractlabels = {0: "Intensity", 1: "Intensity", 2: "Area", 3: "Mass", 4: "Mass"}

        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, 200, 120, 200, 230, 250, 130])
        pub.subscribe(self.on_motion, 'newxy')
        pub.subscribe(self.pres.on_selection, 'scans_selected')

        self.menu = meta_menu(self, self.config, self.pres, type="Chrom")
        self.SetMenuBar(self.menu.menuBar)

        self.panel = wx.Panel(self)
        self.panel.SetDropTarget(ChromDropTarget(self))

        self.mainsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.leftsizer = wx.BoxSizer(wx.VERTICAL)

        labelfont = wx.Font(16, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)

        self.ctlsizer = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(self.panel, label="Chromatogram Parsing Tools", size=(300, 30))
        label.SetFont(labelfont)
        self.ctlsizer.Add(label, 0)

        self.manual_add_button = wx.Button(self.panel, label="Add From Manual Selection")
        self.Bind(wx.EVT_BUTTON, self.pres.on_manual_add, self.manual_add_button)
        self.ctlsizer.Add(self.manual_add_button)

        tsizer0 = wx.BoxSizer(wx.HORIZONTAL)
        self.cpeaksbutton = wx.Button(self.panel, label="Find Chrom. Peaks Near Width:")
        self.Bind(wx.EVT_BUTTON, self.pres.on_chrom_peaks, self.cpeaksbutton)
        tsizer0.Add(self.cpeaksbutton)

        self.ctlcpeaks_param1 = wx.TextCtrl(self.panel, value=str(self.config.chrom_peak_width), size=(50, 20))
        tsizer0.Add(self.ctlcpeaks_param1, 0, wx.ALIGN_CENTER_VERTICAL)
        tsizer0.Add(wx.StaticText(self.panel, label="min"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlsizer.Add(tsizer0)

        tsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.timepartbutton = wx.Button(self.panel, label="Partition in Time Steps of:")
        self.Bind(wx.EVT_BUTTON, self.pres.on_timepart, self.timepartbutton)
        tsizer1.Add(self.timepartbutton)

        self.ctlmin = wx.TextCtrl(self.panel, value=str(self.config.time_window), size=(50, 20))
        tsizer1.Add(self.ctlmin, 0, wx.ALIGN_CENTER_VERTICAL)
        tsizer1.Add(wx.StaticText(self.panel, label="min"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlsizer.Add(tsizer1)

        tsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.swbutton = wx.Button(self.panel, label="Sliding Window Width:")
        self.Bind(wx.EVT_BUTTON, self.pres.on_sliding_window, self.swbutton)
        tsizer2.Add(self.swbutton)

        self.ctlswwin = wx.TextCtrl(self.panel, value=str(self.config.sw_time_window), size=(50, 20))
        tsizer2.Add(self.ctlswwin, 0, wx.ALIGN_CENTER_VERTICAL)
        tsizer2.Add(wx.StaticText(self.panel, label="Offset:"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlswoffset = wx.TextCtrl(self.panel, value=str(int(self.config.sw_scan_offset)), size=(50, 20))
        tsizer2.Add(self.ctlswoffset, 0, wx.ALIGN_CENTER_VERTICAL)
        #tsizer2.Add(wx.StaticText(self.panel, label="min. Offset:"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlsizer.Add(tsizer2)

        self.clear_button = wx.Button(self.panel, label="Clear All Spectra")
        self.Bind(wx.EVT_BUTTON, self.pres.on_clear_spectra, self.clear_button)
        self.ctlsizer.Add(self.clear_button)

        label = wx.StaticText(self.panel, label="Parsed Spectra", size=(300, 30))
        label.SetFont(labelfont)
        self.ctlsizer.Add(label, 0)

        self.leftsizer.Add(self.ctlsizer)

        self.ypanel = ListCtrlPanel(self.panel, self.pres, size=(300, 300))
        self.leftsizer.Add(self.ypanel, 1, wx.EXPAND)

        self.ctlsizer2 = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(self.panel, label="UniDec of Manual Selection", size=(300, 30))
        label.SetFont(labelfont)
        self.ctlsizer2.Add(label, 0)

        self.run_ud_button = wx.Button(self.panel, label="Run UniDec On Selection")
        self.Bind(wx.EVT_BUTTON, self.pres.on_unidec_run, self.run_ud_button)
        self.ctlsizer2.Add(self.run_ud_button)

        self.pick_peaks_button_individual = wx.Button(self.panel, label="Pick Peaks On Selection")
        self.Bind(wx.EVT_BUTTON, self.pres.on_pick_peaks_individual, self.pick_peaks_button_individual)
        self.ctlsizer2.Add(self.pick_peaks_button_individual)

        label = wx.StaticText(self.panel, label="Peaks for Manual Selection", size=(300, 30))
        label.SetFont(labelfont)
        self.ctlsizer2.Add(label, 0)

        self.singlepeakpanel = peaklistsort.PeakListCtrlPanel(self.panel, meta=False, size=(300, 300))
        self.Bind(self.singlepeakpanel.EVT_DELETE_SELECTION_2, self.pres.on_single_delete, self.singlepeakpanel)
        self.Bind(self.singlepeakpanel.EVT_CHARGE_STATE, self.pres.on_single_charge_states, self.singlepeakpanel)
        self.Bind(self.singlepeakpanel.EVT_DIFFERENCES, self.pres.on_single_differences, self.singlepeakpanel)
        self.Bind(self.singlepeakpanel.EVT_MASSES, self.pres.on_single_label_masses, self.singlepeakpanel)
        self.ctlsizer2.Add(self.singlepeakpanel, 0, wx.EXPAND)

        self.leftsizer.Add(self.ctlsizer2, 0, wx.EXPAND)

        self.mainsizer.Add(self.leftsizer, 0, wx.EXPAND)

        plotwindow = scrolled.ScrolledPanel(self.panel)
        sizerplot = wx.GridBagSizer()

        figsize = (4.9, 3.5)
        self.plotc = plot1d.Plot1d(plotwindow, figsize=figsize)  # Chromatogram
        self.plotm = plot1d.Plot1d(plotwindow, figsize=figsize)  # Selection from chromatogram
        self.plot1 = plot1d.Plot1d(plotwindow, smash=1, figsize=figsize)  # MUD Plot 1 m/z cascade
        self.plot2 = plot1d.Plot1d(plotwindow, figsize=figsize)  # MUD Deconvolved Data
        self.plot7 = plot1d.Plot1d(plotwindow, figsize=figsize)  # MUD Extraction
        self.plot2s = plot1d.Plot1d(plotwindow, figsize=figsize)  # Selection mass
        self.plot3 = plot2d.Plot2d(plotwindow, figsize=figsize)  # MUD 2D m/z vs. time
        self.plot5 = plot2d.Plot2d(plotwindow, figsize=figsize)  # MUD 2D mass vs. time

        sizerplot.Add(self.plotc, (0, 0), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plotm, (1, 0), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot1, (0, 1), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot2, (1, 1), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot7, (2, 1), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot2s, (2, 0), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot3, (3, 0), span=(1, 1), flag=wx.EXPAND)
        sizerplot.Add(self.plot5, (3, 1), span=(1, 1), flag=wx.EXPAND)

        self.plots = [self.plotc, self.plotm, self.plot1, self.plot2, self.plot7, self.plot2s, self.plot3, self.plot5]
        self.plotnames = ["Chrom_TIC", "Chrom_mz_selected", "ChromFigure_mz", "ChromFigure_mass", "ChromFigure_XIC",
                          "Chrom_mass_selected","Chrom2Dmz", "Chrom2Dmass"]

        plotwindow.SetSizerAndFit(sizerplot)
        plotwindow.SetupScrolling()

        self.mainsizer.Add(plotwindow, 1, wx.EXPAND)

        self.sizer3 = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(self.panel, label=" Peaks for Parsed Spectra", size=(300, 30))
        label.SetFont(labelfont)
        self.sizer3.Add(label, 0)
        self.peakpanel = peaklistsort.PeakListCtrlPanel(self.panel, meta=True)
        self.Bind(self.peakpanel.EVT_DELETE_SELECTION_2, self.pres.on_delete, self.peakpanel)
        self.Bind(self.peakpanel.EVT_CHARGE_STATE, self.pres.on_charge_states_mud, self.peakpanel)
        self.Bind(self.peakpanel.EVT_DIFFERENCES, self.pres.on_differences, self.peakpanel)
        self.Bind(self.peakpanel.EVT_MASSES, self.pres.on_label_masses, self.peakpanel)
        self.sizer3.Add(self.peakpanel, 0, wx.EXPAND)
        self.mainsizer.Add(self.sizer3, 0, wx.EXPAND)

        self.controls = main_controls(self, self.config, self.pres, self.panel, self.icon_path)
        self.mainsizer.Add(self.controls, 0, wx.EXPAND)

        self.panel.SetSizer(self.mainsizer)
        self.mainsizer.Fit(self)

        keys = [["E", self.pres.on_auto, self.controls.autobutton],
                # ["G", self.pres.on_paste_spectrum, self.menu.menupastespectrum],
                ["R", self.pres.on_unidec_button, self.controls.udbutton],
                ["D", self.pres.on_dataprep_button, self.controls.dataprepbutton],
                ["O", self.pres.on_open, self.menu.openmenu],  # ["I", self.pres.on_integrate],
                ["P", self.pres.on_pick_peaks, self.controls.plotbutton],  # ["K", self.pres.on_plot_peaks],
                # ["C", self.pres.on_plot_composite, self.controls.compositebutton],
                # ["N", self.pres.on_wizard, self.menu.wizardmenu],
                # ["F", self.pres.on_plot_offsets],  # ["Z", self.pres.on_charge_plot],
                # ["L", self.pres.on_load_state], ["S", self.pres.on_save_state],
                # ["B", self.pres.on_batch_run, self.menu.menubatchrun],
                ["Q", self.on_exit, self.menu.menuExit],
                ["T", self.pres.on_mass_tools, self.menu.menuMassFile],
                ["M", self.pres.on_match, self.menu.menumatch],
                ["W", self.pres.on_auto_peak_width, self.menu.menuAutoWidth],
                # ["Z", self.pres.on_undo, self.menu.menuundo],
                # ["Y", self.pres.on_redo, self.menu.menuredo],
                ["K", self.pres.on_kendrick, self.menu.menukendrick]
                ]
        self.setup_shortcuts(keys)

        self.Centre()
        self.Show(True)

    def clear_plots(self, e=None):
        self.plot1.clear_plot()
        self.plot2.clear_plot()
        self.plot7.clear_plot()
        self.plot3.clear_plot()
        self.plot5.clear_plot()

    def export_gui_to_config(self):
        self.controls.export_gui_to_config()
        self.pres.export_vars()

    def on_motion(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=6)
        except:
            pass


class ChromDropTarget(wx.FileDropTarget):
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
            path = filenames[0]
            self.window.pres.open_file(path)
        elif len(filenames) > 1:
            for f in filenames:
                self.window.pres.open_file(f)
        else:
            print("Error in drag and drop.")
        return 0

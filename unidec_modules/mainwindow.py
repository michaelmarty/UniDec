import os
import wx
import numpy as np
from unidec_modules import ColorPlot, plot3d, plot1d, plot2d, peaklistsort, miscwindows
# import wx.lib.inspection
import unidec_modules.unidectools as ud
import wx.lib.scrolledpanel as scrolled
import platform
from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub

import wx.lib.agw.foldpanelbar as fpb

__author__ = 'Michael.Marty'


# noinspection PyAttributeOutsideInit,PyUnusedLocal,PyUnusedLocal
class Mainwindow(wx.Frame):
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
        wx.Frame.__init__(self, None, title=title)
        # Set presenter and config
        self.pres = parent
        self.config = config
        self.title = title

        self.version = self.pres.eng.version

        # Set Icon File
        if os.path.isfile(iconfile):
            favicon = wx.Icon(iconfile, wx.BITMAP_TYPE_ANY)
            wx.Frame.SetIcon(self, favicon)
            self.icon_path = os.path.abspath(iconfile)
        else:
            self.icon_path = None

        # Get a few tool bar icons
        tsize = (16, 16)
        try:
            self.open_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize)
            self.next_bmp = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD, wx.ART_TOOLBAR, tsize)
            self.report_bmp = wx.ArtProvider.GetBitmap(wx.ART_LIST_VIEW, wx.ART_TOOLBAR, tsize)
            self.A_bmp = wx.ArtProvider.GetBitmap(wx.ART_HELP_SETTINGS, wx.ART_TOOLBAR, tsize)
            try:
                self.ud_bmp = wx.BitmapFromImage(wx.Image(iconfile).Rescale(tsize[0], tsize[1]))
            except Exception, ex:
                self.ud_bmp = wx.ArtProvider.GetBitmap(wx.ART_HELP_SETTINGS, wx.ART_TOOLBAR, tsize)
                print ex
        except Exception, ex:
            self.open_bmp = None
            self.next_bmp = None
            self.report_bmp = None
            self.A_bmp = None
            self.ud_bmp = None
            print ex

        # Create Status Bar
        self.CreateStatusBar(7)
        self.SetStatusWidths([-1, 300, 200, 200, 250, 150, 130])

        # Get display size and intelligently reshape
        self.system = platform.system()
        self.displaysize = wx.GetDisplaySize()

        if tabbed is None:
            # If tabbed isn't specified, use the display size to decide what is best
            print "Display Size ", self.displaysize
            self.tabbed = 0
            if self.displaysize[0] < 1400:
                self.tabbed = 1
            elif 1500 <= self.displaysize[0] < 1800:
                self.config.figsize = (6, 5)
            elif 1400 <= self.displaysize[0] < 1500:
                self.config.figsize = (5, 4)
            elif self.displaysize[0]>=1800:
                self.config.figsize = (7, 5)
        else:
            self.tabbed = tabbed

        self.imflag = self.config.imflag

        self.twave = self.config.twaveflag > 0

        self.backgroundchoices = self.config.backgroundchoices

        pub.subscribe(self.on_motion, 'newxy')

        self.setup_menu()
        self.setup_main_panel()
        self.setup_shortcuts()

        self.Centre()
        self.Show(True)
        pass

    def setup_shortcuts(self):
        """
        Setup shortcuts in GUI. Binds key combinations to functions in presenter (self.pres)
        :return: None
        """
        keys = [["E", self.pres.on_auto], ["G", self.pres.on_paste_spectrum],
                ["R", self.pres.on_unidec_button], ["D", self.pres.on_dataprep_button],
                ["O", self.pres.on_open], ["I", self.pres.on_integrate],
                ["P", self.pres.on_pick_peaks], ["K", self.pres.on_plot_peaks],
                ["C", self.pres.on_plot_composite], ["N", self.pres.on_replot],
                ["F", self.pres.on_plot_offsets],  # ["Z", self.pres.on_charge_plot],
                ["L", self.pres.on_load_state], ["S", self.pres.on_save_state],
                ["B", self.pres.on_batch], ["Q", self.on_exit],
                ["T", self.pres.on_mass_tools], ["M", self.pres.on_match],
                ["W", self.pres.on_auto_peak_width],
                ["Z", self.pres.on_undo], ["Y", self.pres.on_redo]]
        ids = [wx.NewId() for a in keys]
        tab = []
        for i, k in enumerate(keys):
            self.Bind(wx.EVT_MENU, k[1], id=ids[i])
            tab.append((wx.ACCEL_CTRL, ord(k[0]), ids[i]))
        self.SetAcceleratorTable(wx.AcceleratorTable(tab))
        pass

    # noinspection PyPep8,PyPep8,PyPep8,PyPep8,PyPep8
    def setup_menu(self):
        """
        Sets menu and binds menu objects to functions in presenter and window
        :return: None
        """
        self.filemenu = wx.Menu()
        self.toolsmenu = wx.Menu()
        self.analysismenu = wx.Menu()
        self.advancedmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()

        # File Menu
        self.menuOpen = self.filemenu.Append(wx.ID_OPEN, "Open Text File\tCtrl+O", " Open a Text File in x y format")
        self.menuOpenRaw = self.filemenu.Append(wx.ID_ANY, "Open Raw File", " Open a Waters .Raw File")
        self.filemenu.AppendSeparator()

        self.menuLoadState = self.filemenu.Append(wx.ID_ANY, "Load State\tCtrl+L", "Load state from folder")
        self.menuSaveState = self.filemenu.Append(wx.ID_ANY, "Save State\tCtrl+S", "Save program state to fresh folder")
        self.filemenu.AppendSeparator()

        self.menupastespectrum = self.filemenu.Append(wx.ID_ANY, "Get Spectrum From Clipboard\tCtrl+G",
                                                      "Pastes the spectrum, formats, and loads")
        self.filemenu.AppendSeparator()

        self.menuLoad = self.filemenu.Append(wx.ID_ANY, "Load External Config File", "Load in a configuration file")
        self.menuLoadDefault = self.filemenu.Append(wx.ID_ANY, "Load Default Config File",
                                                    "Load in default configuration file")
        self.menuSaveDefault = self.filemenu.Append(wx.ID_ANY, "Save Default Config File",
                                                    "Save default configuration file")

        # Default Submenu
        self.defaultmenu = wx.Menu()
        self.menuDefault0 = self.defaultmenu.Append(999, "Low-resolution Native",
                                                    "General factory defaults for low-resolution native MS")
        self.menuDefault1 = self.defaultmenu.Append(1000, "High-resolution Native",
                                                    "Defaults for high-resolution data (Exactive EMR).")
        self.menuDefault2 = self.defaultmenu.Append(1001, "Zero-Charge Mass Deisotoping",
                                                    "Defaults for loading zero-charge mass spectrum as output")
        self.menuDefault3 = self.defaultmenu.Append(1002, "Isotopic Resolution",
                                                    "Defaults for isotopically resolved data.")
        self.menuDefault4 = self.defaultmenu.Append(1003, "Nanodiscs",
                                                    "Defaults for POPC Nanodiscs.")
        self.filemenu.AppendMenu(wx.ID_ANY, "Presets", self.defaultmenu)
        self.filemenu.AppendSeparator()

        self.menufigdialog = self.filemenu.Append(wx.ID_ANY, "Save Figure As",
                                                  "Dialog to define path and extension for figure saving")
        # Figure Submenu
        self.figmenu = wx.Menu()

        self.menuSaveFigure0 = self.figmenu.Append(wx.ID_ANY, "Save Figures as .pdf", "Save all figures to .pdf format")
        self.menuSaveFigure1s = self.figmenu.Append(wx.ID_ANY, "Save Figures as PDF thumbnail",
                                                    "Save all figures to small PDF format")
        self.menuSaveFigure1 = self.figmenu.Append(wx.ID_ANY, "Save Figures as .eps", "Save all figures to .eps format")
        self.menuSaveFigure2 = self.figmenu.Append(wx.ID_ANY, "Save Figures as .png", "Save all figures to .png format")
        if self.config.imflag == 0:
            self.menuSaveFigure4 = self.figmenu.Append(wx.ID_ANY, "Generate PDF Report", "Generate PDF Report")
        else:
            self.menuSaveFigure4 = self.figmenu.Append(wx.ID_ANY, "Save Figures as .pdf",
                                                       "Save all figures to .pdf format")
        self.filemenu.AppendMenu(wx.ID_ANY, 'Save Figure Presets', self.figmenu)
        self.filemenu.AppendSeparator()

        self.menuAbout = self.filemenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
        self.menuExit = self.filemenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", " Terminate the Program")

        # Setting Up the Tools Menu


        self.menuBatch = self.toolsmenu.Append(wx.ID_ANY, "Batch Process\tCtrl+B",
                                               "Run these conditions on multiple samples")
        self.menuBatch2 = self.toolsmenu.Append(wx.ID_ANY, "Batch Process - Independent Data Ranges",
                                                "Run these conditions with varying preset data ranges on multiple samples")
        self.toolsmenu.AppendSeparator()

        self.menuBatchRaw = self.toolsmenu.Append(wx.ID_ANY, "Simple Batch Process Raw To Txt",
                                                  "Batch convert .raw to .txt files with default parameters")
        self.menuImportWizard = self.toolsmenu.Append(wx.ID_ANY, "Raw to Txt Conversion Wizard",
                                                      "Batch convert .raw to .txt files with wizard")
        self.toolsmenu.AppendSeparator()

        self.menuAutoWidth = self.toolsmenu.Append(wx.ID_ANY, "Automatic Peak Width\tCtrl+W",
                                                   "Try to get peak width automatically.")
        self.menuWidth = self.toolsmenu.Append(wx.ID_ANY, "Peak Width Tool", "Help determine the peak width")
        self.toolsmenu.AppendSeparator()

        self.menuManualFile = self.toolsmenu.Append(wx.ID_ANY, "Manual Assignment",
                                                    "Manually set UniDec to preassign specific m/z values")
        self.toolsmenu.AppendSeparator()

        self.menuMassFile = self.toolsmenu.Append(wx.ID_ANY, "Oligomer and Mass Tools\tCtrl+T",
                                                  "Oligomer and Mass Tools")

        # Setting up Analysis Menu
        self.menuPlotZ = self.analysismenu.Append(wx.ID_ANY, "Native Charge/Mass Tools",
                                                  "Tool for exploring relationship between charge and native mass and extraction specific distributions")
        self.analysismenu.AppendSeparator()

        self.menucollect = self.analysismenu.Append(wx.ID_ANY, "Data Collector Utility",
                                                    "Collect results and extract values.  Experimental KD fitting.")
        self.analysismenu.AppendSeparator()

        self.menuExport = self.analysismenu.Append(wx.ID_ANY, "Export Peaks Parameters and Data",
                                                   "Export intensities of charge states, areas, average charge state, and other parameters for the peaks")
        self.menuFitNorm = self.analysismenu.Append(wx.ID_ANY, "Fit Peak Intensities",
                                                    "Fits masses and reports normalized and relative peak intensities")
        self.analysismenu.AppendSeparator()

        self.menukendrick = self.analysismenu.Append(wx.ID_ANY, "Kendrick Mass Tools", "Kendrick Mass Analysis")
        self.menu2Dgrid = self.analysismenu.Append(wx.ID_ANY, "2D Grid Analysis", "2D Grid Analysis")
        self.menuautocorr = self.analysismenu.Append(wx.ID_ANY, "Autocorrelation",
                                                     "Autocorrelation of Mass Distribution")
        self.analysismenu.AppendSeparator()

        self.menuintegrate = self.analysismenu.Append(wx.ID_ANY, "Integrate Peaks\tCtrl+I",
                                                      "Peak area with range set by Peak Detection Range or Integration Range")
        self.menumatch = self.analysismenu.Append(wx.ID_ANY, "Auto Match Peaks\tCtrl+M",
                                                  "Run \"Match to Mixed Oligomers\" in Oligomer and Mass Tools")
        self.analysismenu.AppendSeparator()

        self.menucom = self.analysismenu.Append(wx.ID_ANY, "Report Center of Mass",
                                                "Reports center of mass for the zoomed region of the zero-charge mass spectrum")
        self.analysismenu.AppendSeparator()

        self.menuchargeplot = self.analysismenu.Append(wx.ID_ANY, "Plot By Charge",
                                                       "Plots Mass Distributions as a Function of Charge")
        self.menuoffset = self.analysismenu.Append(wx.ID_ANY, "Plot Charge Offsets\tCtrl+F",
                                                   "Plots Mass vs. Charge Offset in 2D Plot")

        if self.config.imflag == 1:
            self.analysismenu.AppendSeparator()
            self.menuimtools = self.analysismenu.Append(wx.ID_ANY, "IM Parameters Tool",
                                                        "Tools for estimating IM parameters.")
            self.menuimtools2 = self.analysismenu.Append(wx.ID_ANY, "IM Extraction Tool",
                                                         "Tools for Extraction IM Results.")
            self.menunativeccs = self.analysismenu.Append(wx.ID_ANY, "Plot Predicted Native CCS",
                                                          "Plot Predicted Native CCS values on mass vs ccs plot.")
            self.Bind(wx.EVT_MENU, self.pres.on_im_tools, self.menuimtools)
            self.Bind(wx.EVT_MENU, self.pres.on_im_extract, self.menuimtools2)
            self.Bind(wx.EVT_MENU, self.pres.on_plot_nativeccs, self.menunativeccs)

        # Setting up the Advanced Menu
        self.menuReset = self.advancedmenu.Append(wx.ID_ANY, "Reset To Factory Default", "Reset Parameters to Default")
        self.advancedmenu.AppendSeparator()

        self.menuUnidecPath = self.advancedmenu.Append(wx.ID_FILE1, "UniDec File", "Find the UniDec executable file.")
        self.menuFileName = self.advancedmenu.Append(wx.ID_FILE2, "Rename Files",
                                                     "Rename the files into and out of UniDec.")
        self.menuOpenDir = self.advancedmenu.Append(wx.ID_ANY, "Open Saved File Directory",
                                                    "Opens the save directory in the file explorer")

        if self.config.imflag == 0:
            self.advancedmenu.AppendSeparator()
            self.advancedmenu.Append(401, "Best Mode", "Best UniDec Deconvolution Algorithm", wx.ITEM_RADIO)
            self.advancedmenu.Append(402, "Original Mode", "Deconvolution that is more aggressive, Marty et. Al. JASMS",
                                     wx.ITEM_RADIO)
            self.advancedmenu.Append(403, "Richardson-Lucy Mode", "Richardson-Lucy Deconvolution", wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.menu_401_403, id=401)
            self.Bind(wx.EVT_MENU, self.menu_401_403, id=402)
            self.Bind(wx.EVT_MENU, self.menu_401_403, id=403)

        self.advancedmenu.AppendSeparator()

        self.scalemenu = wx.Menu()

        self.scalemenu.Append(501, "Linear", "Normal linear intensity scale", wx.ITEM_RADIO)
        self.scalemenu.Append(502, "Logarithmic", "Logarithmic intensity scale",
                              wx.ITEM_RADIO)
        self.scalemenu.Append(503, "Square Root", "Square root intensity scale", wx.ITEM_RADIO)
        self.Bind(wx.EVT_MENU, self.menu_501_503, id=501)
        self.Bind(wx.EVT_MENU, self.menu_501_503, id=502)
        self.Bind(wx.EVT_MENU, self.menu_501_503, id=503)
        self.advancedmenu.AppendMenu(wx.ID_ANY, 'Intensity Scale', self.scalemenu)
        self.advancedmenu.AppendSeparator()

        if self.config.imflag == 0:
            self.menuflipmode = self.advancedmenu.Append(wx.ID_ANY, "Switch to Ion Mobility Mode",
                                                         "Switch interface to IM-MS Mode.")
        else:
            self.menuflipmode = self.advancedmenu.Append(wx.ID_ANY, "Switch to 1D Mass Spec Mode",
                                                         "Switch interface to MS Mode.")
        if self.tabbed == 0:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Tabbed Plots Mode",
                                                           "Put plots in individual tabs.")
        else:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Single Plot Window",
                                                           "Put plots in single large window.")

        # Experimental Menu
        self.menuundo = self.experimentalmenu.Append(wx.ID_ANY, "Undo Parameter Change\tCtrl+Z",
                                              "Go back to the previous set of parameters")
        self.menuredo = self.experimentalmenu.Append(wx.ID_ANY, "Redo Parameter Change\tCtrl+Y",
                                              "Go to the next set of parameters")
        self.experimentalmenu.AppendSeparator()
        self.Tweet = self.experimentalmenu.Append(wx.ID_ANY, "Twitter", "Twitter Extension")
        if self.config.imflag == 0:
            self.experimentalmenu.AppendSeparator()
            self.menuAdditionalParameters = self.experimentalmenu.Append(wx.ID_ANY, "Additional Parameters",
                                                                         "Adjust some experimental parameters")
            self.experimentalmenu.AppendSeparator()
            # self.menuMinimize = self.experimentalmenu.Append(wx.ID_ANY, "Minimize", "Minimize Peak List")
            # self.experimentalmenu.AppendSeparator()
            self.menuDeisotope = self.experimentalmenu.Append(wx.ID_ANY, "Load Zero-Charge Mass Spectrum",
                                                              "Load Zero-Charge Mass as input spectrum. Useful for deisotoping.")
            self.experimentalmenu.AppendSeparator()
            self.menuCrossValidate = self.experimentalmenu.Append(wx.ID_ANY, "Cross Validate",
                                                                  "Estimate errors through cross validation.")
            self.experimentalmenu.AppendSeparator()
            self.menucolor1d = self.experimentalmenu.Append(wx.ID_ANY, "Color Plots",
                                                            "Make a Different Colored 1D Plot")
            self.Bind(wx.EVT_MENU, self.pres.on_zerocharge_mass, self.menuDeisotope)
            self.Bind(wx.EVT_MENU, self.pres.on_cross_validate, self.menuCrossValidate)
            self.Bind(wx.EVT_MENU, self.pres.on_additional_parameters, self.menuAdditionalParameters)
            self.Bind(wx.EVT_MENU, self.pres.on_color_plot1d, self.menucolor1d)
            # self.Bind(wx.EVT_MENU, self.pres.on_minimize, self.menuMinimize)
        self.experimentalmenu.AppendSeparator()

        self.menusuperbatch = self.experimentalmenu.Append(wx.ID_ANY, "Speed Batch",
                                                           "Minimal Batch Run, only writing config and calling program.")
        self.experimentalmenu.AppendSeparator()

        self.menumassprocess = self.experimentalmenu.Append(wx.ID_ANY, "Process Zero-Charge Mass Spectrum",
                                                            "Apply smoothing, background subtraction, and intensity threshold to zero-charge mass spectrum")
        self.experimentalmenu.AppendSeparator()

        self.menuerrors = self.experimentalmenu.Append(wx.ID_ANY, "Get Errors")
        self.menucal = self.experimentalmenu.Append(wx.ID_ANY, "Apply Calibration")
        self.Bind(wx.EVT_MENU, self.pres.on_calibrate, self.menucal)
        self.experimentalmenu.AppendSeparator()
        self.menufft = self.experimentalmenu.Append(wx.ID_ANY, "FFT Window")
        self.Bind(wx.EVT_MENU, self.pres.on_fft_window, self.menufft)

        self.menugriddecon = self.experimentalmenu.Append(wx.ID_ANY, "Grid Deconvolution")
        self.Bind(wx.EVT_MENU, self.pres.on_grid_decon, self.menugriddecon)

        self.menuhdf5 = self.experimentalmenu.Append(wx.ID_ANY, "Write HDF5")
        self.Bind(wx.EVT_MENU, self.pres.on_write_hdf5, self.menuhdf5)

        # Setting Menu Bar
        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        self.menuBar.Append(self.analysismenu, "Analysis")
        self.menuBar.Append(self.advancedmenu, "Advanced")
        self.menuBar.Append(self.experimentalmenu, "Experimental")
        self.SetMenuBar(self.menuBar)

        # Set Events for Menu Bar

        # File Menu
        self.Bind(wx.EVT_MENU, self.pres.on_open, self.menuOpen)
        self.Bind(wx.EVT_MENU, self.pres.on_raw_open, self.menuOpenRaw)
        self.Bind(wx.EVT_MENU, self.pres.on_load_state, self.menuLoadState)
        self.Bind(wx.EVT_MENU, self.pres.on_save_state, self.menuSaveState)
        self.Bind(wx.EVT_MENU, self.pres.on_paste_spectrum, self.menupastespectrum)
        self.Bind(wx.EVT_MENU, self.pres.on_load_conf_file, self.menuLoad)
        self.Bind(wx.EVT_MENU, self.pres.on_save_default, self.menuSaveDefault)
        self.Bind(wx.EVT_MENU, self.pres.on_load_default, self.menuLoadDefault)
        self.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault0)
        self.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault1)
        self.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault2)
        self.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault3)
        self.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault4)
        self.Bind(wx.EVT_MENU, self.on_save_figure_dialog, self.menufigdialog)
        self.Bind(wx.EVT_MENU, self.on_save_figure_pdf, self.menuSaveFigure0)
        self.Bind(wx.EVT_MENU, self.on_save_figure_eps, self.menuSaveFigure1)
        self.Bind(wx.EVT_MENU, self.on_save_figure_small, self.menuSaveFigure1s)
        self.Bind(wx.EVT_MENU, self.on_save_figure_png, self.menuSaveFigure2)
        self.Bind(wx.EVT_MENU, self.pres.on_pdf_report, self.menuSaveFigure4)
        self.Bind(wx.EVT_MENU, self.on_about, self.menuAbout)
        self.Bind(wx.EVT_MENU, self.on_exit, self.menuExit)

        # Tools

        self.Bind(wx.EVT_MENU, self.pres.on_import_wizard, self.menuImportWizard)
        self.Bind(wx.EVT_MENU, self.pres.on_mass_tools, self.menuMassFile)
        self.Bind(wx.EVT_MENU, self.pres.on_batch, self.menuBatch)
        self.Bind(wx.EVT_MENU, self.pres.on_batch2, self.menuBatch2)
        self.Bind(wx.EVT_MENU, self.pres.on_batch_raw, self.menuBatchRaw)
        self.Bind(wx.EVT_MENU, self.pres.on_peak_width_tool, self.menuWidth)
        self.Bind(wx.EVT_MENU, self.pres.on_auto_peak_width, self.menuAutoWidth)
        self.Bind(wx.EVT_MENU, self.pres.on_manual, self.menuManualFile)
        self.Bind(wx.EVT_MENU, self.pres.on_undo, self.menuundo)
        self.Bind(wx.EVT_MENU, self.pres.on_redo, self.menuredo)

        # Analysis
        self.Bind(wx.EVT_MENU, self.pres.on_nativez_tools, self.menuPlotZ)
        self.Bind(wx.EVT_MENU, self.pres.on_data_collector, self.menucollect)
        self.Bind(wx.EVT_MENU, self.pres.on_export_params, self.menuExport)
        self.Bind(wx.EVT_MENU, self.pres.on_fit_masses, self.menuFitNorm)
        self.Bind(wx.EVT_MENU, self.pres.on_plot_offsets, self.menuoffset)
        self.Bind(wx.EVT_MENU, self.pres.on_charge_plot, self.menuchargeplot)
        self.Bind(wx.EVT_MENU, self.pres.on_center_of_mass, self.menucom)
        self.Bind(wx.EVT_MENU, self.pres.on_integrate, self.menuintegrate)
        self.Bind(wx.EVT_MENU, self.pres.on_2d_grid, self.menu2Dgrid)
        self.Bind(wx.EVT_MENU, self.pres.on_kendrick, self.menukendrick)
        self.Bind(wx.EVT_MENU, self.pres.on_autocorr_window, self.menuautocorr)
        self.Bind(wx.EVT_MENU, self.pres.on_match, self.menumatch)

        # Advanced
        self.Bind(wx.EVT_MENU, self.pres.on_reset, self.menuReset)
        self.Bind(wx.EVT_MENU, self.pres.on_unidec_path, self.menuUnidecPath)
        self.Bind(wx.EVT_MENU, self.pres.on_file_name, self.menuFileName)
        self.Bind(wx.EVT_MENU, self.on_open_dir, self.menuOpenDir)
        self.Bind(wx.EVT_MENU, self.pres.on_flip_tabbed, self.menufliptabbed)
        self.Bind(wx.EVT_MENU, self.pres.on_flip_mode, self.menuflipmode)

        # Experimental
        self.Bind(wx.EVT_MENU, self.pres.on_peak_errors, self.menuerrors)
        self.Bind(wx.EVT_MENU, self.pres.on_tweet, self.Tweet)
        self.Bind(wx.EVT_MENU, self.pres.on_mass_process, self.menumassprocess)
        self.Bind(wx.EVT_MENU, self.pres.on_super_batch, self.menusuperbatch)
        pass

    def on_defaults(self, e):
        """
        Resets the configuration to a default predefined in the unidecstructure file.
        :param e: Menu event
        :return: None
        """
        nid = e.GetId() % 100
        if nid == 0:
            self.config.default_high_res()
        elif nid == 1:
            self.config.default_zero_charge()
        elif nid == 2:
            self.config.default_isotopic_res()
        elif nid == 3:
            self.config.default_nanodisc()
        elif nid == 99:
            self.config.default_decon_params()
        self.pres.import_config(None)

    # noinspection PyPep8,PyPep8
    def setup_main_panel(self):
        """
        Lays Out Main Panel. Binds some functions to presenter.
        :return: None
        """
        # Sizers to develop layout
        # s1 = (min(self.displaysize[0], 1851), self.displaysize[1])
        # s2 = (550, self.displaysize[1])
        splitterwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D | wx.SP_BORDER)
        splitterwindow2 = wx.SplitterWindow(splitterwindow, -1, style=wx.SP_3D | wx.SP_BORDER)
        panelp = wx.Panel(splitterwindow2, -1)
        panel = wx.Panel(splitterwindow2, -1)
        splitterwindow2.SplitVertically(panelp, panel, sashPosition=-250 - self.config.imflag * 20)
        splitterwindow2.SetMinimumPaneSize(175)
        splitterwindow.SetMinimumPaneSize(175)
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
            figsize = (6, 5)
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
        sizerpeaks.Add(self.peakpanel, 0, wx.EXPAND)
        panelp.SetSizer(sizerpeaks)
        sizerpeaks.Fit(self)

        # ..........................
        #
        # Sizers for Controls
        #
        # .............................
        sizercontrol = wx.BoxSizer(wx.VERTICAL)

        # Small Toolbar
        buttonsizer = wx.BoxSizer(wx.HORIZONTAL)
        bsize = (50, 25)
        self.openbutton = wx.BitmapButton(panel, -1, self.open_bmp, size=bsize)
        self.procbutton = wx.BitmapButton(panel, -1, self.next_bmp, size=bsize)
        self.procbutton.SetBackgroundColour(wx.Colour(150, 150, 255))
        self.udbutton = wx.BitmapButton(panel, -1, self.ud_bmp, size=bsize)
        self.udbutton.SetBackgroundColour(wx.Colour(255, 255, 150))
        self.ppbutton = wx.BitmapButton(panel, -1, self.report_bmp, size=bsize)
        self.ppbutton.SetBackgroundColour(wx.Colour(255, 150, 150))
        self.autobutton = wx.Button(panel, -1, "All", size=bsize)  # self.A_bmp
        self.Bind(wx.EVT_BUTTON, self.pres.on_open, self.openbutton)
        self.Bind(wx.EVT_BUTTON, self.pres.on_dataprep_button, self.procbutton)
        self.Bind(wx.EVT_BUTTON, self.pres.on_unidec_button, self.udbutton)
        self.Bind(wx.EVT_BUTTON, self.pres.on_pick_peaks, self.ppbutton)
        self.Bind(wx.EVT_BUTTON, self.pres.on_auto, self.autobutton)
        buttons = [self.openbutton, self.procbutton, self.udbutton, self.ppbutton, self.autobutton]
        for button in buttons:
            buttonsizer.Add(button, 0, wx.EXPAND)
        sizercontrol.Add(buttonsizer, 0, wx.EXPAND)

        # Setting up main fold controls
        size1 = (75, -1)
        foldpanels = fpb.FoldPanelBar(panel, -1, agwStyle=fpb.FPB_VERTICAL)
        style1 = fpb.CaptionBarStyle()
        style1b = fpb.CaptionBarStyle()
        style1c = fpb.CaptionBarStyle()
        style2 = fpb.CaptionBarStyle()
        style2b = fpb.CaptionBarStyle()
        style3 = fpb.CaptionBarStyle()
        style3b = fpb.CaptionBarStyle()

        bright = 150
        bright2 = 200
        style1.SetFirstColour(wx.Colour(bright, bright, 255))
        style1b.SetFirstColour(wx.Colour(bright2, bright2, 255))
        style1c.SetFirstColour(wx.Colour(bright2, 255, bright2))
        style2.SetFirstColour(wx.Colour(255, 255, bright))
        style2b.SetFirstColour(wx.Colour(255, 255, bright2))
        style3.SetFirstColour(wx.Colour(255, bright, bright))
        style3b.SetFirstColour(wx.Colour(255, bright2, bright2))

        bright3 = 255
        bright4 = 255
        style1.SetSecondColour(wx.Colour(bright3, bright3, 255))
        style1b.SetSecondColour(wx.Colour(bright4, bright4, 255))
        style1c.SetSecondColour(wx.Colour(bright4, 255, bright4))
        style2.SetSecondColour(wx.Colour(255, 255, bright3))
        style2b.SetSecondColour(wx.Colour(255, 255, bright4))
        style3.SetSecondColour(wx.Colour(255, bright3, bright3))
        style3b.SetSecondColour(wx.Colour(255, bright4, bright4))

        # Panel to set the data prep controls
        foldpanel1 = foldpanels.AddFoldPanel(caption="Data Processing", collapsed=False, cbstyle=style1)
        panel1 = wx.Panel(foldpanel1, -1)
        sizercontrol1 = wx.GridBagSizer(wx.VERTICAL)

        self.ctlminmz = wx.TextCtrl(panel1, value="", size=(50, -1))
        self.ctlmaxmz = wx.TextCtrl(panel1, value="", size=(60, -1))
        mzrange = wx.BoxSizer(wx.HORIZONTAL)
        if self.config.imflag == 1:
            mzrange.Add(wx.StaticText(panel1, label="               "), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)
        mzrange.Add(wx.StaticText(panel1, label="m/z Range: "), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)
        mzrange.Add(self.ctlminmz)
        mzrange.Add(wx.StaticText(panel1, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        mzrange.Add(self.ctlmaxmz)
        mzrange.Add(wx.StaticText(panel1, label=" Th "), 0, wx.ALIGN_CENTER_VERTICAL)
        sizercontrol1.Add(mzrange, (0, 0), span=(1, 5))

        if self.config.imflag == 1:
            self.ctlmindt = wx.TextCtrl(panel1, value="", size=(50, -1))
            self.ctlmaxdt = wx.TextCtrl(panel1, value="", size=(60, -1))
            dtrange = wx.BoxSizer(wx.HORIZONTAL)
            dtrange.Add(wx.StaticText(panel1, label="Arrival Time Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
            dtrange.Add(self.ctlmindt)
            dtrange.Add(wx.StaticText(panel1, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
            dtrange.Add(self.ctlmaxdt)
            dtrange.Add(wx.StaticText(panel1, label=" ms "), 0, wx.ALIGN_CENTER_VERTICAL)
            sizercontrol1.Add(dtrange, (1, 0), span=(1, 5))

            self.ctlsmoothdt = wx.TextCtrl(panel1, value="", size=size1)
            self.ctlsubbuffdt = wx.TextCtrl(panel1, value="", size=size1)
            sizercontrol1.Add(self.ctlsubbuffdt, (3, 1))
            sizercontrol1.Add(self.ctlsmoothdt, (5, 1))

            sizercontrol1.Add(wx.StaticText(panel1, label=" m/z"), (2, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            sizercontrol1.Add(wx.StaticText(panel1, label=" m/z"), (4, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            sizercontrol1.Add(wx.StaticText(panel1, label=" A.T."), (3, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            sizercontrol1.Add(wx.StaticText(panel1, label=" A.T."), (5, 2), flag=wx.ALIGN_CENTER_VERTICAL)

            self.imflag = 1
        else:
            self.imflag = 0

        self.subtypectl = wx.Choice(panel1, -1, choices=self.backgroundchoices)
        self.dataprepbutton = wx.Button(panel1, -1, "Process Data")
        self.Bind(wx.EVT_BUTTON, self.pres.on_dataprep_button, self.dataprepbutton)
        self.ctlbuff = wx.TextCtrl(panel1, value="", size=size1)
        self.ctlsmooth = wx.TextCtrl(panel1, value="", size=size1)
        self.ctlbinsize = wx.TextCtrl(panel1, value="", size=size1)

        sizercontrol1.Add(self.subtypectl, (1 + self.config.imflag, 0))
        sizercontrol1.Add(self.ctlbuff, (1 + self.config.imflag, 1))
        sizercontrol1.Add(self.ctlsmooth, (2 + self.config.imflag * 2, 1))
        sizercontrol1.Add(wx.StaticText(panel1, label="Gaussian Smoothing: "), (2 + self.config.imflag * 2, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol1.Add(self.ctlbinsize, (3 + self.config.imflag * 3, 1))
        sizercontrol1.Add(wx.StaticText(panel1, label="Bin Every: "), (3 + self.config.imflag * 3, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol1.Add(self.dataprepbutton, (4 + self.config.imflag * 3, 0), span=(1, 2), flag=wx.EXPAND)

        panel1.SetSizer(sizercontrol1)
        sizercontrol1.Fit(panel1)
        foldpanels.AddFoldPanelWindow(foldpanel1, panel1, fpb.FPB_ALIGN_WIDTH)
        foldpanels.AddFoldPanelWindow(foldpanel1, wx.StaticText(foldpanel1, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Additional Data Processing Parameters
        foldpanel1b = foldpanels.AddFoldPanel(caption="Additional Data Processing Parameters", collapsed=True,
                                              cbstyle=style1b)
        panel1b = wx.Panel(foldpanel1b, -1)
        gbox1b = wx.GridBagSizer(wx.VERTICAL)

        if self.config.imflag == 1:
            self.ctlpusher = wx.TextCtrl(panel1b, value="", size=size1)
            gbox1b.Add(self.ctlpusher, (0, 1))
            gbox1b.Add(wx.StaticText(panel1b, label=u"Pusher Interval (\u03BCs)"), (0, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlintthresh = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlintthresh, (0 + self.config.imflag, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Intensity Threshold: "), (0 + self.config.imflag, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctladductmass = wx.TextCtrl(panel1b, value='', size=size1)
        gbox1b.Add(self.ctladductmass, (1 + self.config.imflag, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Adduct Mass (Da): "), (1 + self.config.imflag, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlaccelvolt = wx.TextCtrl(panel1b, value='', size=size1)
        gbox1b.Add(self.ctlaccelvolt, (2 + self.config.imflag, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Acceleration Voltage (kV): "), (2 + self.config.imflag, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        if self.config.imflag == 0:
            self.ctlbintype = wx.Choice(panel1b, -1, size=(240, 50),
                                        choices=["Linear m/z (Constant " + u'\N{GREEK CAPITAL LETTER DELTA}' + "m/z)",
                                                 "Linear resolution (Constant (m/z)/(" + u'\N{GREEK CAPITAL LETTER DELTA}' + "m/z))",
                                                 "Nonlinear", "Linear Interpolated", "Linear Resolution Interpolated"])
            gbox1b.Add(self.ctlbintype, (3 + self.config.imflag, 0), span=(1, 2))

        else:
            self.ctlconvertflag = wx.CheckBox(panel1b, label="Compress when converting to .txt")
            self.ctlconvertflag.SetValue(True)
            gbox1b.Add(self.ctlconvertflag, (3 + self.config.imflag, 0), span=(1, 2))

        panel1b.SetSizer(gbox1b)
        gbox1b.Fit(panel1b)
        foldpanels.AddFoldPanelWindow(foldpanel1b, panel1b, fpb.FPB_ALIGN_WIDTH)
        if self.config.imflag == 1:
            foldpanels.AddFoldPanelWindow(foldpanel1b, wx.StaticText(foldpanel1b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        if self.config.imflag == 1:
            foldpanel1c = foldpanels.AddFoldPanel(caption="Ion Mobility Parameters", collapsed=True, cbstyle=style1c)
            panel1c = wx.Panel(foldpanel1c, -1)
            gbox1c = wx.GridBagSizer(wx.VERTICAL)

            self.ctltwave = wx.RadioBox(panel1c, label="", choices=["Linear Cell", "Travelling Wave"])
            self.Bind(wx.EVT_RADIOBOX, self.pres.on_flip_twave, self.ctltwave)
            gbox1c.Add(self.ctltwave, (0, 0), span=(1, 2))

            self.twave = self.config.twaveflag > 0
            if not self.twave:
                self.ctltwave.SetSelection(0)

                self.ctlvolt = wx.TextCtrl(panel1c, value="", size=size1)
                gbox1c.Add(self.ctlvolt, (1, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Voltage (V): "), (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlpressure = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlpressure, (2, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Pressure (Torr): "), (2, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctltemp = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctltemp, (3, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label=u"Temperature (\u00B0C): "), (3, 0),
                           flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlgasmass = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlto = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlto, (5, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label=u"Dead Time (t\u2080 in ms): "), (5, 0),
                           flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctldriftlength = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctldriftlength, (6, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Drift Cell Length (m)"), (6, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            else:
                self.ctltwave.SetSelection(1)

                self.ctltcal1 = wx.TextCtrl(panel1c, value="", size=size1)
                gbox1c.Add(self.ctltcal1, (1, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Calibration Parameter 1: "), (1, 0),
                           flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctltcal2 = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctltcal2, (2, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Calibration Parameter 2: "), (2, 0),
                           flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctledc = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctledc, (3, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="EDC Parameter: "), (3, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlgasmass = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctltwavecaltype = wx.Choice(panel1c, -1, choices=self.config.twavedict.values())
                gbox1c.Add(self.ctltwavecaltype, (5, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Calibration Type: "), (5, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltwave.SetSelection(int(self.twave))
            panel1c.SetSizer(gbox1c)
            gbox1c.Fit(panel1c)

            foldpanels.AddFoldPanelWindow(foldpanel1c, panel1c, fpb.FPB_ALIGN_WIDTH)
            foldpanels.AddFoldPanelWindow(foldpanel1c, wx.StaticText(foldpanel1c, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for UniDec Parameters
        foldpanel2 = foldpanels.AddFoldPanel(caption="UniDec Parameters", collapsed=False, cbstyle=style2)
        panel2 = wx.Panel(foldpanel2, -1)
        sizercontrol2 = wx.GridBagSizer(wx.VERTICAL)

        self.ctlstartz = wx.TextCtrl(panel2, value="", size=(60, -1))
        self.ctlendz = wx.TextCtrl(panel2, value="", size=(60, -1))
        zrange = wx.BoxSizer(wx.HORIZONTAL)
        zrange.Add(wx.StaticText(panel2, label="Charge Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
        zrange.Add(self.ctlstartz)
        zrange.Add(wx.StaticText(panel2, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        zrange.Add(self.ctlendz)

        self.ctlmasslb = wx.TextCtrl(panel2, value="", size=(60, -1))
        self.ctlmassub = wx.TextCtrl(panel2, value="", size=(70, -1))
        massrange = wx.BoxSizer(wx.HORIZONTAL)
        massrange.Add(wx.StaticText(panel2, label="Mass Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
        massrange.Add(self.ctlmasslb)
        massrange.Add(wx.StaticText(panel2, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        massrange.Add(self.ctlmassub)
        massrange.Add(wx.StaticText(panel2, label=" Da  "), 0, wx.ALIGN_CENTER_VERTICAL)

        if self.config.imflag == 1:
            self.ctlccslb = wx.TextCtrl(panel2, value="", size=(60, -1))
            self.ctlccsub = wx.TextCtrl(panel2, value="", size=(70, -1))
            ccsrange = wx.BoxSizer(wx.HORIZONTAL)
            ccsrange.Add(wx.StaticText(panel2, label="CCS Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
            ccsrange.Add(self.ctlccslb)
            ccsrange.Add(wx.StaticText(panel2, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
            ccsrange.Add(self.ctlccsub)
            ccsrange.Add(wx.StaticText(panel2, label=u" \u212B\u00B2  "), 0, wx.ALIGN_CENTER_VERTICAL)
            sizercontrol2.Add(ccsrange, (2, 0), span=(1, 2))

            self.ctlccsbins = wx.TextCtrl(panel2, value="", size=size1)
            sizercontrol2.Add(self.ctlccsbins, (4, 1), span=(1, 2))
            sizercontrol2.Add(wx.StaticText(panel2, label=u"Sample CCS Every (\u212B\u00B2): "), (4, 0),
                              flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctldtsig = wx.TextCtrl(panel2, value="", size=size1)
            sizercontrol2.Add(self.ctldtsig, (6, 1), span=(1, 2))
            sizercontrol2.Add(wx.StaticText(panel2, label="Peak FWHM (ms): "), (6, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlmassbins = wx.TextCtrl(panel2, value="", size=size1)
        self.ctlmzsig = wx.TextCtrl(panel2, value="", size=size1)
        self.ctlpsfun = wx.RadioBox(panel2, label="Peak Shape Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        self.rununidec = wx.Button(panel2, -1, "Run UniDec")
        self.Bind(wx.EVT_BUTTON, self.pres.on_unidec_button, self.rununidec)

        sizercontrol2.Add(zrange, (0, 0), span=(1, 2))
        sizercontrol2.Add(massrange, (1, 0), span=(1, 2))
        sizercontrol2.Add(self.ctlmassbins, (2 + self.config.imflag, 1), span=(1, 2))
        sizercontrol2.Add(self.ctlmzsig, (3 + self.config.imflag * 2, 1), span=(1, 2))
        sizercontrol2.Add(wx.StaticText(panel2, label="Sample Mass Every (Da): "), (2 + self.config.imflag, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2.Add(wx.StaticText(panel2, label="Peak FWHM (Th): "), (3 + self.config.imflag * 2, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2.Add(self.ctlpsfun, (4 + self.config.imflag * 3, 0), span=(1, 2))
        sizercontrol2.Add(self.rununidec, (5 + self.config.imflag * 3, 0), span=(1, 2), flag=wx.EXPAND)

        panel2.SetSizer(sizercontrol2)
        sizercontrol2.Fit(panel2)
        foldpanels.AddFoldPanelWindow(foldpanel2, panel2, fpb.FPB_ALIGN_WIDTH)
        foldpanels.AddFoldPanelWindow(foldpanel2, wx.StaticText(foldpanel2, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Additional Restraints
        foldpanel2b = foldpanels.AddFoldPanel(caption="Additional Filters/Restraints", collapsed=True, cbstyle=style2b)
        panel2b = wx.Panel(foldpanel2b, -1)
        gbox2b = wx.GridBagSizer(wx.VERTICAL)

        self.ctlzzsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Charge Smooth Width: "), (0, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlzzsig, (0, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlmolig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Mass Difference (Da): "), (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlmolig, (1, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlmsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Mass Smooth Width: "), (2, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlmsig, (2, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        if self.config.imflag == 1:
            self.ctlcsig = wx.TextCtrl(panel2b, value="", size=size1)
            gbox2b.Add(wx.StaticText(panel2b, label="CCS Smooth Width: "), (3, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            gbox2b.Add(self.ctlcsig, (3, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlnumit = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label='Maximum # of Iterations: '), (3 + self.config.imflag, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlnumit, (3 + self.config.imflag, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlpoolflag = wx.RadioBox(panel2b, label="m/z to Mass Transformation",
                                       choices=["Integration", "Interpolation"])
        gbox2b.Add(self.ctlpoolflag, (4 + self.config.imflag, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)

        if self.config.imflag == 0:
            self.ctlisotopemode = wx.CheckBox(panel2b, label="Isotope Mode")
            gbox2b.Add(self.ctlisotopemode, (5 + self.config.imflag, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlmanualassign = wx.CheckBox(panel2b, label="Manual Mode")
        self.Bind(wx.EVT_CHECKBOX, self.on_check_manual, self.ctlmanualassign)
        gbox2b.Add(self.ctlmanualassign, (5 + self.config.imflag, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        mlsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmasslistflag = wx.CheckBox(panel2b, label="Mass List Window:")
        self.Bind(wx.EVT_CHECKBOX, self.on_mass_list, self.ctlmasslistflag)
        self.ctlmtabsig = wx.TextCtrl(panel2b, value="", size=(60, -1))
        mlsizer.Add(self.ctlmasslistflag, 0, wx.ALIGN_CENTER_VERTICAL)
        mlsizer.Add(self.ctlmtabsig, 0, wx.ALIGN_CENTER_VERTICAL)
        mlsizer.Add(wx.StaticText(panel2b, label=" Da "), 0, wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(mlsizer, (6 + self.config.imflag, 0), span=(1, 2))

        sb = wx.StaticBox(panel2b, label='Native Charge Offset Range')
        sbs = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        self.ctlminnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        self.ctlmaxnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        sbs.Add(self.ctlminnativez, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
        sbs.Add(wx.StaticText(panel2b, label=' to '), 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)
        sbs.Add(self.ctlmaxnativez, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
        gbox2b.Add(sbs, (7 + self.config.imflag, 0), span=(1, 2), flag=wx.EXPAND)

        if self.config.imflag == 1:
            sb2 = wx.StaticBox(panel2b, label='Native CCS Offset Range')
            sbs2 = wx.StaticBoxSizer(sb2, orient=wx.HORIZONTAL)
            self.ctlnativeccslb = wx.TextCtrl(panel2b, value='', size=(75, -1))
            self.ctlnativeccsub = wx.TextCtrl(panel2b, value='', size=(75, -1))
            sbs2.Add(self.ctlnativeccslb, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
            sbs2.Add(wx.StaticText(panel2b, label=' to '), 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)
            sbs2.Add(self.ctlnativeccsub, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
            sbs2.Add(wx.StaticText(panel2b, label=u" \u212B\u00B2 "), 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)
            gbox2b.Add(sbs2, (8 + self.config.imflag, 0), span=(1, 2), flag=wx.EXPAND)

        panel2b.SetSizer(gbox2b)
        gbox2b.Fit(panel2b)
        foldpanels.AddFoldPanelWindow(foldpanel2b, panel2b, fpb.FPB_ALIGN_WIDTH)
        foldpanels.AddFoldPanelWindow(foldpanel2b, wx.StaticText(foldpanel2b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Peak Selection and Plotting
        foldpanel3 = foldpanels.AddFoldPanel(caption="Peak Selection and Plotting", collapsed=False, cbstyle=style3)
        panel3 = wx.Panel(foldpanel3, -1)

        sizercontrol3 = wx.GridBagSizer(wx.VERTICAL)
        self.ctlwindow = wx.TextCtrl(panel3, value="", size=size1)
        self.ctlthresh = wx.TextCtrl(panel3, value="", size=size1)
        self.ctlnorm = wx.RadioBox(panel3, label="Peak Normalization", choices=["None", "Max", "Total"])

        self.plotbutton = wx.Button(panel3, -1, "Peak Detection")
        self.plotbutton2 = wx.Button(panel3, -1, "Plot Peaks")
        self.Bind(wx.EVT_BUTTON, self.pres.on_plot_peaks, self.plotbutton2)
        self.Bind(wx.EVT_BUTTON, self.pres.on_pick_peaks, self.plotbutton)

        sizercontrol3.Add(self.ctlwindow, (0, 1))
        sizercontrol3.Add(wx.StaticText(panel3, label="Peak Detection Range (Da): "), (0, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol3.Add(self.ctlthresh, (1, 1))
        sizercontrol3.Add(wx.StaticText(panel3, label="Peak Detection Threshold: "), (1, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol3.Add(self.ctlnorm, (2, 0), span=(1, 2), flag=wx.EXPAND)
        sizercontrol3.Add(self.plotbutton, (3, 0), span=(1, 1), flag=wx.EXPAND)
        sizercontrol3.Add(self.plotbutton2, (3, 1), span=(1, 1), flag=wx.EXPAND)

        panel3.SetSizer(sizercontrol3)
        sizercontrol3.Fit(panel3)
        foldpanels.AddFoldPanelWindow(foldpanel3, panel3, fpb.FPB_ALIGN_WIDTH)
        foldpanels.AddFoldPanelWindow(foldpanel3, wx.StaticText(foldpanel3, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Plotting Parameters
        foldpanel3b = foldpanels.AddFoldPanel(caption="Additional Plotting Parameters", collapsed=True, cbstyle=style3b)
        panel3b = wx.Panel(foldpanel3b, -1)

        gbox3b = wx.GridBagSizer(wx.VERTICAL)
        gbox3b.Add(wx.StaticText(panel3b, label='2D Color Map: '), (0, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctl2dcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)
        gbox3b.Add(wx.StaticText(panel3b, label='Peaks Color Map: '), (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlpeakcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)

        for mp in self.config.cmaps2:
            self.ctl2dcm.Append(mp)
        for mp in self.config.cmaps:
            self.ctlpeakcm.Append(mp)

        gbox3b.Add(self.ctl2dcm, (0, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(self.ctlpeakcm, (1, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctldiscrete = wx.CheckBox(panel3b, label="Discrete Plot")
        gbox3b.Add(self.ctldiscrete, (2, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlpublicationmode = wx.CheckBox(panel3b, label="Publication Mode")
        gbox3b.Add(self.ctlpublicationmode, (2, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlrawflag = wx.RadioBox(panel3b, label="", choices=["Reconvolved", "Raw Results"])
        gbox3b.Add(self.ctlrawflag, (3, 0), span=(1, 2), flag=wx.EXPAND)

        gbox3b.Add(wx.StaticText(panel3b, label="Marker Threshold: "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(wx.StaticText(panel3b, label="Species Separation: "), (5, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlthresh2 = wx.TextCtrl(panel3b, value="", size=size1)
        self.ctlsep = wx.TextCtrl(panel3b, value="", size=size1)
        gbox3b.Add(self.ctlthresh2, (4, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(self.ctlsep, (5, 1), flag=wx.ALIGN_CENTER_VERTICAL)

        sb2 = wx.StaticBox(panel3b, label='Integration Range')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.HORIZONTAL)
        self.ctlintlb = wx.TextCtrl(panel3b, value='', size=(75, -1))
        self.ctlintub = wx.TextCtrl(panel3b, value='', size=(75, -1))
        sbs2.Add(self.ctlintlb, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
        sbs2.Add(wx.StaticText(panel3b, label=' to '), 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)
        sbs2.Add(self.ctlintub, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, border=5)
        sbs2.Add(wx.StaticText(panel3b, label=' Da '), 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.EXPAND)
        gbox3b.Add(sbs2, (6, 0), span=(1, 2), flag=wx.EXPAND)

        self.replotbutton = wx.Button(panel3b, -1, "Replot")
        self.Bind(wx.EVT_BUTTON, self.pres.on_replot, self.replotbutton)
        gbox3b.Add(self.replotbutton, (7, 0), span=(1, 1), flag=wx.EXPAND)

        self.compositebutton = wx.Button(panel3b, -1, "Plot Composite")
        self.Bind(wx.EVT_BUTTON, self.pres.on_plot_composite, self.compositebutton)
        gbox3b.Add(self.compositebutton, (7, 1), span=(1, 1), flag=wx.EXPAND)

        if self.config.imflag == 1:
            self.cubeplotbutton = wx.Button(panel3b, -1, "Plot Cubes")
            self.Bind(wx.EVT_BUTTON, self.pres.make_cube_plot, self.cubeplotbutton)
            gbox3b.Add(self.cubeplotbutton, (8, 0), span=(1, 2), flag=wx.EXPAND)

        panel3b.SetSizer(gbox3b)
        gbox3b.Fit(panel3b)
        foldpanels.AddFoldPanelWindow(foldpanel3b, panel3b, fpb.FPB_ALIGN_WIDTH)
        foldpanels.AddFoldPanelWindow(foldpanel3b, wx.StaticText(foldpanel3b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        bright = 250
        foldpanel1.SetBackgroundColour(wx.Colour(bright, bright, 255))
        foldpanel1b.SetBackgroundColour(wx.Colour(bright, bright, 255))
        if self.config.imflag == 1:
            foldpanel1c.SetBackgroundColour(wx.Colour(bright, 255, bright))

        foldpanel2.SetBackgroundColour(wx.Colour(255, 255, bright))
        foldpanel2b.SetBackgroundColour(wx.Colour(255, 255, bright))

        foldpanel3.SetBackgroundColour(wx.Colour(255, bright, bright))
        foldpanel3b.SetBackgroundColour(wx.Colour(255, bright, bright))

        sizercontrol.SetMinSize((250 + self.config.imflag * 10, 0))

        # Add to top control sizer
        sizercontrol.Add(foldpanels, -1, wx.EXPAND)
        panel.SetSizer(sizercontrol)
        sizercontrol.Fit(self)

        if self.system == "Linux" and self.tabbed != 1:
            sizerplot.Fit(splitterwindow)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(splitterwindow, 1, wx.EXPAND)

        # Set everything up
        self.SetSizer(sizer)
        sizer.Fit(self)
        self.setup_tool_tips()
        pass

    # noinspection PyPep8
    def setup_tool_tips(self):
        """
        Sets Tool Tips for items on the Main Panel
        :return: None
        """
        self.ctlthresh2.SetToolTip(wx.ToolTip(
            "Set threshold for peaks to be plotted in m/z. Peak at given charge state must be greater than threshold * maximum m/z intensity."))
        self.ctlsep.SetToolTip(wx.ToolTip("Set distance between isolated peak m/z plots."))
        self.ctlwindow.SetToolTip(
            wx.ToolTip("Peak detection window. Peak must be maximum in a +/- window range in mass (Da)."))
        self.ctlthresh.SetToolTip(wx.ToolTip(
            "Peak detection threshold. Peak's intensity must be great than threshold times maximum mass intensity."))
        self.plotbutton.SetToolTip(wx.ToolTip("Pick peaks and plot. (Ctrl+P)"))
        self.plotbutton2.SetToolTip(wx.ToolTip("Plot individual peak species in m/z. (Ctrl+K)"))
        self.ctlmasslistflag.SetToolTip(wx.ToolTip(
            "Limit deconvolution to specific masses +/- some window.\nDefine in Tools>Oligomer and Mass Tools."))
        self.ctlmtabsig.SetToolTip(
            wx.ToolTip("Set window for mass limitations. Setting to 0 will force only listed masses."))
        self.ctlpsfun.SetToolTip(wx.ToolTip("Expected peak shape.\nSee Tools>Peak Width Tool for more tools."))
        self.rununidec.SetToolTip(wx.ToolTip("Write Configuration File, Run UniDec, and Plot Results. (Ctrl+R)"))
        self.ctlmzsig.SetToolTip(wx.ToolTip(
            "Expected peak FWHM in m/z (Th).\nFor nonlinear mode, minimum FWHM\nSee Tools>Peak Width Tool for more tools."))
        self.ctlzzsig.SetToolTip(wx.ToolTip(
            "Parameter for defining the width of the charge state smooth.\nUniDec will use a mean filter of width 2n+1 on log_e of the charge distribution"))
        self.ctlmassub.SetToolTip(wx.ToolTip(
            "Maximum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlmassbins.SetToolTip(wx.ToolTip("Sets the resolution of the zero-charge mass spectrum"))
        self.ctlmasslb.SetToolTip(wx.ToolTip(
            "Minimum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlstartz.SetToolTip(wx.ToolTip("Minimum allowed charge state in deconvolution."))
        self.ctlendz.SetToolTip(wx.ToolTip("Maximum allowed charge state in deconvolution."))
        self.ctlsmooth.SetToolTip(wx.ToolTip("Gaussian smooth sigma in units of data point number."))
        self.ctlbinsize.SetToolTip(wx.ToolTip(
            "Controls Linearization.\nConstant bin size (Th) for Linear m/z\nMinimum bin size (Th) for Linear Resolution\nNumber of data points compressed together for Nonlinear"))
        self.ctlintthresh.SetToolTip(
            wx.ToolTip("Set intensity threshold. Data points below threshold are excluded from deconvolution."))
        self.ctlbuff.SetToolTip(wx.ToolTip(
            "Background subtraction parameters.\nMinimum: 0=off 1=on\nLine: The first and last n data points will be averaged;\n a line between the two averages will be subtracted.\nCurved: Width of smoothed background"))
        self.ctlminmz.SetToolTip(wx.ToolTip("Set minimum m/z of data"))
        self.ctlmaxmz.SetToolTip(wx.ToolTip("Set maximum m/z of data"))
        self.dataprepbutton.SetToolTip(
            wx.ToolTip("Subtract, linearize, smooth, threshold, and write data to file. (Ctrl+D)"))
        self.ctladductmass.SetToolTip(wx.ToolTip("Mass of charge carrying adduct;\ntypically the mass of a proton"))
        self.ctlaccelvolt.SetToolTip(
            wx.ToolTip("QToF Acceleration Voltage: When set, will correct data for detector efficiency"))
        self.ctlmsig.SetToolTip(wx.ToolTip("Width of Mass Smooth Filter"))
        self.ctlmolig.SetToolTip(wx.ToolTip("Mass difference used for Mass Smooth Filter"))
        self.ctlminnativez.SetToolTip(wx.ToolTip("Minimum offset from a native charge state"))
        self.ctlmaxnativez.SetToolTip(wx.ToolTip("Maximum offset from a native charge state"))
        if self.config.imflag == 0:
            self.ctlisotopemode.SetToolTip(wx.ToolTip("Use isotopic distributions in deconvolution"))
            self.ctlmanualassign.SetToolTip(wx.ToolTip("Use manual assignments. See Tools>Manual Assignment"))
            self.ctlbintype.SetToolTip(wx.ToolTip(
                "Sets how to bin the data\nValue set by above with Bin Every\nLinear bins with linear m/z axis\nLinear Resolution bins with m/z axis that has a constant resolution\nNonlinear merges adjacent data points\nInterpolation uses the same axes but with interpolation instead of integration"))
        self.ctlnumit.SetToolTip(wx.ToolTip(
            "Maximum number of iterations. Note: Deconvolution will stop automically before this if it converges."))
        self.ctldiscrete.SetToolTip(wx.ToolTip("Set 2D plots to discrete rather than continuous"))
        self.ctlpublicationmode.SetToolTip(wx.ToolTip("Set plots to look good for publication rather than utility"))
        self.ctlrawflag.SetToolTip(
            wx.ToolTip(
                "Decide whether to outputs should be reconvolved with the peak shape or the raw deconvolution"))
        self.ctl2dcm.SetToolTip(wx.ToolTip("Set 2D plot color function"))
        self.ctlpeakcm.SetToolTip(wx.ToolTip("Set the color function for the peaks"))
        self.ctlintlb.SetToolTip(wx.ToolTip(
            "Controls range for integration.\nDefault is +/- Peak Detection Window.\nUses these boxes to manually set - and +"))
        self.ctlintub.SetToolTip(wx.ToolTip(
            "Controls range for integration.\nDefault is +/- Peak Detection Window.\nUses these boxes to manually set - and +"))
        self.ctlnorm.SetToolTip(wx.ToolTip(
            "Sets normalization of mass data.\nMaximum will normalize so that the maximum value is %100.\nTotal will normalize so that the sum of all peaks is %100"))
        self.replotbutton.SetToolTip(wx.ToolTip("Replot some of the plots. (Ctrl+N)"))
        self.compositebutton.SetToolTip(wx.ToolTip("Plot composite of simulated spectra from selected peaks. (Ctrl+C)"))
        self.openbutton.SetToolTip(wx.ToolTip("Open .txt or .mzML file (Ctrl+O)"))
        self.procbutton.SetToolTip(wx.ToolTip("Process Data (Ctrl+D)"))
        self.udbutton.SetToolTip(wx.ToolTip("Run UniDec (Ctrl+R)"))
        self.ppbutton.SetToolTip(wx.ToolTip("Pick Peaks (Ctrl+P)"))
        self.autobutton.SetToolTip(wx.ToolTip("Process Data, Run UniDec, Pick Peaks (Ctrl+E)"))
        if self.config.imflag == 1:
            self.cubeplotbutton.SetToolTip(wx.ToolTip("Plot the cubes!"))
            self.ctlmindt.SetToolTip(wx.ToolTip("Set minimum arrival time in ms"))
            self.ctlmaxdt.SetToolTip(wx.ToolTip("Set maximum arrival time in ms"))
            self.ctlsubbuffdt.SetToolTip(wx.ToolTip("Same as for m/z but in the arrival time dimension"))
            self.ctlsmoothdt.SetToolTip(wx.ToolTip("Same as for m/z but in the arrival time dimension"))
            self.ctlpusher.SetToolTip(
                wx.ToolTip(
                    "Pusher interval to convert ion mobility bin number into arrival time.\nAT=Pusher*Bin/1000"))
            self.ctltwave.SetToolTip(wx.ToolTip("Change between Linear Cell and Travelling Wave"))
            if not self.twave:
                self.ctlvolt.SetToolTip(wx.ToolTip("Linear voltage across drift cell."))
                self.ctlpressure.SetToolTip(wx.ToolTip("Pressure of drift cell in Torr."))
                self.ctltemp.SetToolTip(wx.ToolTip("Temperature of drift cell in Celsius."))
                self.ctlgasmass.SetToolTip(wx.ToolTip("Mass of gas in drift cell.\nDefault is for Helium."))
                self.ctlto.SetToolTip(wx.ToolTip("Set dead time in instrument.\nDrift Time=Arrival Time - Dead Time"))
                self.ctlto.SetToolTip(wx.ToolTip("Set Drift Cell Length in meters."))
            else:
                self.ctltcal1.SetToolTip(wx.ToolTip("T-Wave calibration parameter 1."))
                self.ctltcal2.SetToolTip(wx.ToolTip("T-Wave calibration parameter 2."))
                self.ctledc.SetToolTip(wx.ToolTip("T-Wave instrument EDC parameter."))
                self.ctlgasmass.SetToolTip(wx.ToolTip("Mass of gas in drift cell.\nDefault is for Nitrogen."))
            self.ctlccslb.SetToolTip(wx.ToolTip("Minimum allowed CCS."))
            self.ctlccsub.SetToolTip(wx.ToolTip("Maximum allowed CCS."))
            self.ctlccsbins.SetToolTip(wx.ToolTip("Sample CCS at this resolution."))
            self.ctldtsig.SetToolTip(wx.ToolTip("Peak width in the arrival time dimension."))
            self.ctlcsig.SetToolTip(wx.ToolTip("Width of CCS smooth.\nCCS difference is taken as the CCS bin size."))
            self.ctlnativeccslb.SetToolTip(wx.ToolTip("Sets lower bound on offset from native CCS."))
            self.ctlnativeccsub.SetToolTip(wx.ToolTip("Sets upper bound on offset from native CCS."))
        self.ctlpoolflag.SetToolTip(wx.ToolTip(
            "Sets type of conversion from m/z to mass.\nIntegration:\n\tEach m/z bin goes to the nearest mass bin\n\tBest for undersampled masses\nInterpolation:\n\tEach mass value interpolates its value in m/z space\n\tBest for oversampled mass data"))
        pass

    def import_config_to_gui(self):
        """
        Imports parameters from the config object to the GUI.
        :return: None
        """
        if self.config.batchflag == 0:

            if self.config.imflag == 1 and self.twave != (self.config.twaveflag > 0):
                self.pres.on_flip_twave(0)
            self.ctlmassbins.SetValue(str(self.config.massbins))
            self.ctlstartz.SetValue(str(self.config.startz))
            self.ctlendz.SetValue(str(self.config.endz))
            self.ctlzzsig.SetValue(str(self.config.zzsig))
            self.ctlmzsig.SetValue(str(self.config.mzsig))
            self.ctlpsfun.SetSelection(self.config.psfun)
            self.ctlnorm.SetSelection(self.config.peaknorm)
            self.ctlmasslb.SetValue(str(self.config.masslb))
            self.ctlmassub.SetValue(str(self.config.massub))
            self.ctlmasslistflag.SetValue(self.config.mfileflag)
            self.ctlmtabsig.SetValue(str(self.config.mtabsig))
            self.ctlbuff.SetValue(str(self.config.subbuff))
            self.subtypectl.SetSelection(int(self.config.subtype))
            self.ctlsmooth.SetValue(str(self.config.smooth))
            self.ctlbinsize.SetValue(str(self.config.mzbins))
            self.ctlwindow.SetValue(str(self.config.peakwindow))
            self.ctlthresh.SetValue(str(self.config.peakthresh))
            self.ctlthresh2.SetValue(str(self.config.peakplotthresh))
            self.ctlsep.SetValue(str(self.config.separation))
            self.ctlintthresh.SetValue(str(self.config.intthresh))
            self.ctladductmass.SetValue(str(self.config.adductmass))
            self.ctlaccelvolt.SetValue(str(self.config.detectoreffva))
            self.ctlmsig.SetValue(str(self.config.msig))
            self.ctlmolig.SetValue(str(self.config.molig))
            self.ctlnumit.SetValue(str(self.config.numit))
            self.ctlminnativez.SetValue(str(self.config.nativezlb))
            self.ctlmaxnativez.SetValue(str(self.config.nativezub))
            self.ctlpoolflag.SetSelection(self.config.poolflag)
            if self.config.imflag == 0:
                self.ctlmanualassign.SetValue(self.config.manualfileflag)
                self.ctlisotopemode.SetValue(self.config.isotopemode)
                self.ctlbintype.SetSelection(int(self.config.linflag))
            self.ctldiscrete.SetValue(self.config.discreteplot)
            self.ctlpublicationmode.SetValue(self.config.publicationmode)
            self.ctlrawflag.SetSelection(self.config.rawflag)

            try:
                self.ctl2dcm.SetSelection(self.config.cmaps2.index(self.config.cmap))
                self.ctlpeakcm.SetSelection(self.config.cmaps.index(self.config.peakcmap))
            except ValueError:
                print "Could not find the specified color map. Try upgrading to the latest version of matplotlib."
                import matplotlib
                print "Current version:", matplotlib.__version__
                # Revert to the defaults
                self.ctl2dcm.SetSelection(self.config.cmaps.index("spectral"))
                self.ctlpeakcm.SetSelection(self.config.cmaps.index("rainbow"))

            if self.config.imflag == 1:
                self.ctlpusher.SetValue(str(self.config.pusher))
                self.ctlsubbuffdt.SetValue(str(self.config.subbufdt))
                self.ctlsmoothdt.SetValue(str(self.config.smoothdt))
                self.ctlccslb.SetValue(str(self.config.ccslb))
                self.ctlccsub.SetValue(str(self.config.ccsub))
                self.ctlccsbins.SetValue(str(self.config.ccsbins))
                self.ctldtsig.SetValue(str(self.config.dtsig))
                self.ctlcsig.SetValue(str(self.config.csig))
                self.ctlnativeccslb.SetValue(str(self.config.nativeccslb))
                self.ctlnativeccsub.SetValue(str(self.config.nativeccsub))
                if not self.twave:
                    self.ctlvolt.SetValue(str(self.config.volt))
                    self.ctltemp.SetValue(str(self.config.temp))
                    self.ctlpressure.SetValue(str(self.config.pressure))
                    self.ctlgasmass.SetValue(str(self.config.gasmass))
                    self.ctlto.SetValue(str(self.config.to))
                    self.ctldriftlength.SetValue(str(self.config.driftlength))
                else:
                    self.ctltcal1.SetValue(str(self.config.tcal1))
                    self.ctltcal2.SetValue(str(self.config.tcal2))
                    self.ctledc.SetValue(str(self.config.edc))
                    self.ctlgasmass.SetValue(str(self.config.gasmass))
                    self.ctltwavecaltype.SetSelection(self.config.twavedict.keys().index(self.config.twaveflag))

            try:
                x = float(self.config.integratelb)
                y = float(self.config.integrateub)
                self.ctlintlb.SetValue(str(x))
                self.ctlintub.SetValue(str(y))
            except (ValueError, TypeError):
                self.ctlintlb.SetValue("")
                self.ctlintub.SetValue("")
            if self.config.imflag == 0:
                if self.config.aggressiveflag == 1:
                    self.advancedmenu.Check(id=402, check=True)
                elif self.config.aggressiveflag == 2:
                    self.advancedmenu.Check(id=403, check=True)
                else:
                    self.advancedmenu.Check(id=401, check=True)

            if self.config.msig > 0:
                self.SetStatusText(
                    "Oligomer Blur Mass: " + str(self.config.molig) + " Std Dev: " + str(self.config.msig),
                    number=4)
            else:
                self.SetStatusText(" ", number=4)
        # If the batchflag is not 1, it will import the data range as well
        if self.config.batchflag != 1:
            self.ctlminmz.SetValue(str(self.config.minmz))
            self.ctlmaxmz.SetValue(str(self.config.maxmz))
            if self.config.imflag == 1:
                self.ctlmindt.SetValue(str(self.config.mindt))
                self.ctlmaxdt.SetValue(str(self.config.maxdt))

    def export_gui_to_config(self):
        """
        Exports parameters from the GUI to the config object.
        :return: None
        """
        self.config.minmz = ud.string_to_value(self.ctlminmz.GetValue())
        self.config.maxmz = ud.string_to_value(self.ctlmaxmz.GetValue())
        self.config.smooth = ud.string_to_value(self.ctlsmooth.GetValue())
        self.config.mzbins = ud.string_to_value(self.ctlbinsize.GetValue())
        self.config.subbuff = ud.string_to_value(self.ctlbuff.GetValue())
        self.config.subtype = self.subtypectl.GetSelection()
        self.config.intthresh = ud.string_to_value(self.ctlintthresh.GetValue())
        self.config.massbins = ud.string_to_value(self.ctlmassbins.GetValue())
        self.config.endz = ud.string_to_int(self.ctlendz.GetValue())
        self.config.startz = ud.string_to_int(self.ctlstartz.GetValue())
        self.config.zzsig = ud.string_to_value(self.ctlzzsig.GetValue())
        self.config.mzsig = ud.string_to_value(self.ctlmzsig.GetValue())
        self.config.massub = ud.string_to_value(self.ctlmassub.GetValue())
        self.config.masslb = ud.string_to_value(self.ctlmasslb.GetValue())
        self.config.mtabsig = ud.string_to_value(self.ctlmtabsig.GetValue())
        self.config.psfun = self.ctlpsfun.GetSelection()
        self.config.peaknorm = self.ctlnorm.GetSelection()
        self.config.mfileflag = int(self.ctlmasslistflag.GetValue())
        self.config.peakwindow = ud.string_to_value(self.ctlwindow.GetValue())
        self.config.peakthresh = ud.string_to_value(self.ctlthresh.GetValue())
        self.config.peakplotthresh = ud.string_to_value(self.ctlthresh2.GetValue())
        self.config.separation = ud.string_to_value(self.ctlsep.GetValue())
        self.config.adductmass = ud.string_to_value(self.ctladductmass.GetValue())
        self.config.detectoreffva = ud.string_to_value(self.ctlaccelvolt.GetValue())
        self.config.msig = ud.string_to_value(self.ctlmsig.GetValue())
        self.config.molig = ud.string_to_value(self.ctlmolig.GetValue())
        self.config.numit = ud.string_to_int(self.ctlnumit.GetValue())
        self.config.nativezlb = ud.string_to_value(self.ctlminnativez.GetValue())
        self.config.nativezub = ud.string_to_value(self.ctlmaxnativez.GetValue())
        self.config.integratelb = ud.string_to_value(self.ctlintlb.GetValue())
        self.config.integrateub = ud.string_to_value(self.ctlintub.GetValue())
        if self.config.imflag == 0:
            self.config.isotopemode = int(self.ctlisotopemode.GetValue())
            self.config.manualfileflag = int(self.ctlmanualassign.GetValue())
            self.config.linflag = self.ctlbintype.GetSelection()
            if self.config.mzbins == 0:
                self.config.linflag = 2
                self.ctlbintype.SetSelection(int(self.config.linflag))
        self.config.discreteplot = int(self.ctldiscrete.GetValue())
        self.config.publicationmode = int(self.ctlpublicationmode.GetValue())
        self.config.rawflag = self.ctlrawflag.GetSelection()

        self.config.cmap = self.ctl2dcm.GetStringSelection().encode('ascii')
        self.config.peakcmap = self.ctlpeakcm.GetStringSelection().encode('ascii')
        self.config.poolflag = self.ctlpoolflag.GetSelection()

        if self.config.imflag == 1:
            self.config.pusher = ud.string_to_value(self.ctlpusher.GetValue())
            self.config.mindt = ud.string_to_value(self.ctlmindt.GetValue())
            self.config.maxdt = ud.string_to_value(self.ctlmaxdt.GetValue())
            self.config.smoothdt = ud.string_to_value(self.ctlsmoothdt.GetValue())
            self.config.subbufdt = ud.string_to_value(self.ctlsubbuffdt.GetValue())
            self.config.ccslb = ud.string_to_value(self.ctlccslb.GetValue())
            self.config.ccsub = ud.string_to_value(self.ctlccsub.GetValue())
            self.config.ccsbins = ud.string_to_value(self.ctlccsbins.GetValue())
            self.config.dtsig = ud.string_to_value(self.ctldtsig.GetValue())
            self.config.csig = ud.string_to_value(self.ctlcsig.GetValue())
            self.config.nativeccslb = ud.string_to_value(self.ctlnativeccslb.GetValue())
            self.config.nativeccsub = ud.string_to_value(self.ctlnativeccsub.GetValue())

            if not self.twave:
                self.config.volt = ud.string_to_value(self.ctlvolt.GetValue())
                self.config.temp = ud.string_to_value(self.ctltemp.GetValue())
                self.config.pressure = ud.string_to_value(self.ctlpressure.GetValue())
                self.config.gasmass = ud.string_to_value(self.ctlgasmass.GetValue())
                self.config.to = ud.string_to_value(self.ctlto.GetValue())
                self.config.driftlength = ud.string_to_value(self.ctldriftlength.GetValue())
            else:
                self.config.tcal1 = ud.string_to_value(self.ctltcal1.GetValue())
                self.config.tcal2 = ud.string_to_value(self.ctltcal2.GetValue())
                self.config.edc = ud.string_to_value(self.ctledc.GetValue())
                self.config.gasmass = ud.string_to_value(self.ctlgasmass.GetValue())
                self.config.twaveflag = self.config.twavedict.keys()[self.ctltwavecaltype.GetSelection()]

            if not self.config.mindt and not ud.isempty(self.pres.eng.data.rawdata3):
                self.config.mindt = np.amin(self.pres.eng.data.rawdata3[:, 1])
            if not self.config.maxdt and not ud.isempty(self.pres.eng.data.rawdata3):
                self.config.maxdt = np.amax(self.pres.eng.data.rawdata3[:, 1])

        if not self.config.minmz and not ud.isempty(self.pres.eng.data.rawdata):
            self.config.minmz = np.amin(self.pres.eng.data.rawdata[:, 0])
        if not self.config.maxmz and not ud.isempty(self.pres.eng.data.rawdata):
            self.config.maxmz = np.amax(self.pres.eng.data.rawdata[:, 0])
        pass

    def menu_401_403(self, event):
        """
        Menu function to adjust the UniDec core function (agreesiveflag).
        :param event: wx Event
        :return: None
        """
        event_id = event.GetId()
        if event_id == 401:
            self.config.aggressiveflag = 0
        if event_id == 402:
            self.config.aggressiveflag = 1
        if event_id == 403:
            self.config.aggressiveflag = 2
        print self.config.aggressiveflag

    def menu_501_503(self, event):
        """
        Menu function to adjust the intensity scale.
        :param event: wx Event
        :return: None
        """
        event_id = event.GetId()
        if event_id == 501:
            self.config.intscale = "Linear"
        if event_id == 502:
            self.config.intscale = "Logarithmic"
        if event_id == 503:
            self.config.intscale = "Square Root"
        print self.config.intscale

    def clear_all_plots(self, flag=0):
        """
        Clear All Plots
        :return: None
        """
        self.plot1.clear_plot()
        self.plot2.clear_plot()
        self.plot3.clear_plot()
        self.plot4.clear_plot()
        self.plot5.clear_plot()
        self.plot6.clear_plot()
        try:
            if self.config.imflag == 1:
                if flag == 1:
                    self.plot1im.clear_plot()
                self.plot1fit.clear_plot()
                self.plot2ccs.clear_plot()
                self.plot3color.clear_plot()
                self.plot5ccsz.clear_plot()
                self.plot5mccs.clear_plot()
                self.plot9.clear_plot()
                self.plot10.clear_plot()
        except AttributeError:
            pass

    def on_motion(self, xpos, ypos):
        """
        Triggered by pubsub from plot windows. Reports text in Status Bar.
        :param xpos: x position fed from event
        :param ypos: y position fed from event
        :return: None
        """
        if xpos is not None and ypos is not None:
            self.SetStatusText("x=%.2f y=%.2f" % (xpos, ypos), number=6)
        pass

    # .......................................................
    #
    #  The File Menu
    #
    # .......................................................

    def on_about(self, e):
        """
        Displays message about program
        :param e:
        :return:
        """
        dlg = wx.MessageDialog(self,
                               "UniDec GUI version " + self.version +
                               "\nPlease contact mtmarty@email.arizona.edu with any questions, bugs, or features to add.\n"
                               "The latest version may be found at unidec.chem.ox.ac.uk.\n"
                               "If used in publication, please cite Marty et Al. Anal. Chem. 2015, DOI: 10.1021/acs.analchem.5b00140 ",
                               "About UniDec", wx.OK | wx.CENTER)
        dlg.ShowModal()
        dlg.Destroy()

    def on_exit(self, e):
        """
        Exit the Program
        :param e: Dummy wx event
        :return: None
        """
        self.Close(True)

    # .......................................................
    #
    #  Saving Figures
    #
    # .......................................................

    def save_all_figures(self, extension, extension2='', e=0, header=None, **kwargs):
        """
        Save All of the Figures. Will name as header+extension2+_FigureX.+exetension
        :param extension: Figure type (pdf, eps, png). Anything accepted by matplotlib
        :param extension2: Additional text to include in the figure header.
        :param e: Dummy wx Event
        :param header: Option to add different header. Default of none yields self.outfname as the path header
        :param kwargs: Any keywards to pass to the matplotlib savefig command such as Transparent or DPI
        :return: figureflags, files (the figures that were successfully saved and the files that they were saved to)
        """
        self.SetStatusText("Saving Figures", number=5)
        figureflags = []
        files = []
        if header is None:
            header = self.config.outfname + extension2
        else:
            header += extension2

        name1 = header + "_Figure1." + extension
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1, **kwargs)
            figureflags.append(1)
            files.append([1, name1])
        name2 = header + "_Figure2." + extension
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2, **kwargs)
            figureflags.append(2)
            files.append([2, name2])
        name3 = header + "_Figure3." + extension
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3, **kwargs)
            figureflags.append(3)
            files.append([3, name3])
        name4 = header + "_Figure4." + extension
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name4, **kwargs)
            figureflags.append(4)
            files.append([4, name4])
        name5 = header + "_Figure5." + extension
        if self.plot5.flag:
            self.plot5.on_save_fig(e, name5, **kwargs)
            figureflags.append(5)
            files.append([5, name5])
        name6 = header + "_Figure6." + extension
        if self.plot6.flag:
            self.plot6.on_save_fig(e, name6, **kwargs)
            figureflags.append(6)
            files.append([6, name6])
        if self.config.imflag == 1:
            name = header + "_Figure1im." + extension
            if self.plot1im.flag:
                self.plot1im.on_save_fig(e, name, **kwargs)
                figureflags.append(7)
                files.append([7, name])
            name = header + "_Figure1fit." + extension
            if self.plot1fit.flag:
                self.plot1fit.on_save_fig(e, name, **kwargs)
                figureflags.append(8)
                files.append([8, name])
            name = header + "_Figure2ccs." + extension
            if self.plot2ccs.flag:
                self.plot2ccs.on_save_fig(e, name, **kwargs)
                figureflags.append(9)
                files.append([9, name])
            name = header + "_Figure3color." + extension
            if self.plot3color.flag:
                self.plot3color.on_save_fig(e, name, **kwargs)
                figureflags.append(10)
                files.append([10, name])
            name = header + "_Figure5ccsz." + extension
            if self.plot5ccsz.flag:
                self.plot5ccsz.on_save_fig(e, name, **kwargs)
                figureflags.append(11)
                files.append([11, name])
            name = header + "_Figure5massccs." + extension
            if self.plot5mccs.flag:
                self.plot5mccs.on_save_fig(e, name, **kwargs)
                figureflags.append(12)
                files.append([12, name])
            name = header + "_mzCube." + extension
            if self.plot9.flag:
                self.plot9.on_save_fig(e, name, **kwargs)
                figureflags.append(13)
                files.append([13, name])
            name = header + "_massCube." + extension
            if self.plot10.flag:
                self.plot10.on_save_fig(e, name, **kwargs)
                figureflags.append(14)
                files.append([14, name])
        return figureflags, files

    def on_save_figure_eps(self, e):
        """
        Save all figures as EPS
        :param e: Dummy wx event
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        figureflags, files = self.save_all_figures("eps")
        self.SetStatusText("Saved to .eps", number=5)
        pass

    def on_save_figure_png(self, e, **kwargs):
        """
        Save all figures as PNG
        :param e: Dummy wx event
        :param kwargs: keywards to pass to matplotlib savefig
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        flags, self.pngs = self.save_all_figures("png", **kwargs)
        self.SetStatusText("Saved to .png", number=5)
        pass

    def on_save_figure_pdf(self, e):
        """
        Saves all figures as PDF
        :param e: Dummy wx event
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        figureflags, files = self.save_all_figures("pdf")
        self.SetStatusText("Saved to .pdf", number=5)
        return figureflags, files

    def on_save_figure_dialog(self, e):
        """
        Open dialog box to set the parameters for figure type, size, and path to save.
        :param e: Dummy wx event
        :return: None
        """
        figsize = self.plot1.GetSize()
        dpi = wx.ScreenDC().GetPPI()
        defaultsize = np.array([figsize[0] / dpi[0], figsize[1] / dpi[1]])
        defaultrect = np.array([0.1, 0.1, 0.8, 0.8])

        dlg = miscwindows.SaveFigureDialog(self)
        dlg.initialize_interface(self.config)
        code = dlg.ShowModal()
        if code == 0:
            directory = dlg.directory
            header = dlg.header
            extension = dlg.extension
            transparent = dlg.transparent
            dpi = dlg.dpi
            try:
                dpi = int(dpi)
            except Exception, e:
                print e, dpi
                dpi = None

            path = os.path.join(directory, header)
            self.figsize = np.array(dlg.figsize)
            self.rect = np.array(dlg.rect)

            if not np.all(defaultsize == self.figsize) or not np.all(defaultrect == self.rect):
                plots = self.shrink_all_figures()
                self.SetStatusText("Saving Figures", number=5)
                figureflags, files = self.save_all_figures(extension, extension2="", header=path,
                                                           transparent=transparent, dpi=dpi)
                for plot in plots:
                    plot.resize = 1
                    plot.size_handler()
            else:
                self.save_all_figures(extension, extension2="", header=path, transparent=transparent, dpi=dpi)
            # print self.directory
            self.SetStatusText("Saved Figures", number=5)
            # TODO: Remember Values from previous

    def shrink_figure(self, plot):
        """
        Automatically shrinks the plot to a figure size in inches set in self.figsize.
        :param plot: Plot object to shrink
        :return: None
        """
        if plot.flag:
            dpi = wx.ScreenDC().GetPPI()
            figsize2 = (int(self.figsize[0] * dpi[0]), int(self.figsize[1] * dpi[1]))
            plot.resize = 0
            plot.canvas.SetSize(figsize2)
            plot.canvas.draw()
            plot.set_nticks(5)
            plot.subplot1.set_position(self.rect)
            if plot.cbar is not None:
                plot.cbar.ax.set_position([0.85, 0.2, 0.05, 0.7])
            plot.repaint()

    def shrink_all_figures(self):
        """
        Shrinks all figures to the size specified in self.figsize
        :return: A list of plot objects that we shrunk
        """
        plots = [self.plot1, self.plot2, self.plot3, self.plot4, self.plot5, self.plot6]
        if self.config.imflag == 1:
            plots = plots + [self.plot1im, self.plot1fit, self.plot2ccs, self.plot5mccs, self.plot5ccsz,
                             self.plot3color, self.plot9, self.plot10]
        for plot in plots:
            self.shrink_figure(plot)
        return plots

    def on_save_figure_small(self, e):
        """
        Preset to shrink figures to 4.5 in by 3 in and save as PDF.
        :param e: Dummy wx event
        :return: None
        """
        self.figsize = (4.5, 3.0)
        self.rect = [0.2, 0.2, 0.6, 0.7]
        plots = self.shrink_all_figures()
        self.SetStatusText("Saving Figures", number=5)
        figureflags, files = self.save_all_figures("pdf", extension2="_Thumb")
        self.SetStatusText("Saved to Thumbnails", number=5)
        for plot in plots:
            plot.resize = 1
            plot.size_handler()
        pass

    # .......................................................
    #
    #  The Main Panel
    #
    # .......................................................

    def on_check_manual(self, e):
        """
        Checks the configuration to see if values for manual mode are set. If they are not,
        it opens the window to set the manual assignments.
        :param e: Dummy wx event passed on.
        :return: None
        """
        self.config.manualfileflag = self.ctlmanualassign.GetValue()
        if len(self.config.manuallist) < 1:
            self.pres.on_manual(e)
            if len(self.config.manuallist) < 1:
                self.ctlmanualassign.SetValue(False)

    def on_mass_list(self, e):
        """
        Checks the configuration to see if values for the mass list are set. If they are not,
        it opens the window to set the mass list.
        :param e: Dummy wx event passed on.
        :return: None
        """
        self.config.mfileflag = self.ctlmasslistflag.GetValue()
        if len(self.config.masslist) < 1:
            self.pres.on_mass_tools(e)
            if len(self.config.masslist) < 1:
                self.ctlmasslistflag.SetValue(False)

    def on_open_dir(self, e):
        save_dir = os.getcwd()
        print "Opening directory:", save_dir
        try:
            os.system(self.config.opencommand + save_dir)
        except Exception, err:
            print "Error opening directory", err


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
            if os.path.splitext(fname)[1] == ".raw":
                print "Opening .raw file:", fname
                self.window.pres.on_raw_open(0, path)
            else:
                self.window.pres.on_open_file(fname, directory)
                # self.window.pres.on_auto() # Run the whole thing
        elif len(filenames) > 1:
            # Batch process the files that were dropped
            if os.path.splitext(filenames[0])[1] == ".raw":
                print "Batch converting raw to txt"
                self.window.pres.on_batch_raw(0, filenames, clip=False)
            else:
                print "Running batch mode"
                self.window.pres.on_batch(batchfiles=filenames)
        else:
            print "Error in file drop", filenames

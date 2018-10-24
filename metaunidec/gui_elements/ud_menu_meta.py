import wx


class meta_menu(wx.Menu):
    def __init__(self, parent, config, pres):
        super(wx.Menu, self).__init__()
        self.pres = pres
        self.config = config
        self.parent = parent

        self.filemenu = wx.Menu()
        self.toolsmenu = wx.Menu()
        self.analysismenu = wx.Menu()
        self.advancedmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()
        self.helpmenu = wx.Menu()

        # File Menu
        self.openmenu = self.filemenu.Append(wx.ID_ANY, "Open File\tCtrl+O", "Open HDF5 file")
        self.filemenu.AppendSeparator()
        self.wizardmenu = self.filemenu.Append(wx.ID_ANY, "Import Wizard\tCtrl+N", "Import Data Into HDF5 File")
        self.filemenu.AppendSeparator()
        self.filemenu2 = wx.Menu()
        self.newmenu = self.filemenu2.Append(wx.ID_ANY, "New File", "Create New Blank HDF5 file")
        self.addmenu = self.filemenu2.Append(wx.ID_ANY, "Add Data Files",
                                             "Add data from individual text, mzml, or Thermo RAW files to current HDF5 file")
        self.menupastespectrum = self.filemenu2.Append(wx.ID_ANY, "Add Data From Clipboard",
                                                      "Add copied data from clipboard to current HDF5 file")
        self.menuexportdata = self.filemenu2.Append(wx.ID_ANY, "Export HDF5 to Text",
                                                      "Export data in HDF5 to individual text files")
        self.filemenu.AppendSubMenu(self.filemenu2, "Manual File Operations")
        self.filemenu.AppendSeparator()
        self.filemenu3 = wx.Menu()
        self.importmenu = self.filemenu3.Append(wx.ID_ANY, "Auto Import Chromatogram By Time",
                                               "Import mzML or Thermo RAW Ramps to HDF5 files")
        self.importmenu2 = self.filemenu3.Append(wx.ID_ANY, "Auto Import Chromatogram By Scans",
                                               "Import mzML or Thermo RAW Ramps to HDF5 files")
        self.importmenu3 = self.filemenu3.Append(wx.ID_ANY, "Auto Import Multiple Chromatograms By Range of Times",
                                               "Import mzML or Thermo RAW Ramps to HDF5 files")
        self.importmenu4 = self.filemenu3.Append(wx.ID_ANY, "Auto Import Multiple Chromatograms By Range of Scans",
                                                "Import mzML or Thermo RAW Ramps to HDF5 files")
        self.filemenu.AppendSubMenu(self.filemenu3, "Automated Chromatogram Parsing")
        self.filemenu.AppendSeparator()

        # Default Submenu
        self.varmenu = wx.Menu()
        self.menuvar1name = self.varmenu.Append(wx.ID_ANY, "Rename Variable 1", "Rename Variable 1")
        self.menuvar2name = self.varmenu.Append(wx.ID_ANY, "Rename Variable 2", "Rename Variable 1")
        self.menuvar1import = self.varmenu.Append(wx.ID_ANY, "Import Variable Metadata", "Import Variable Metadata")
        self.menuvar1export = self.varmenu.Append(wx.ID_ANY, "Export Variable Metadata", "Export Variable Metadata")
        self.filemenu.AppendSubMenu(self.varmenu, "Variables")
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
        self.filemenu.AppendSubMenu(self.defaultmenu, "Presets")
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

        self.filemenu.AppendSubMenu(self.figmenu, 'Save Figure Presets')
        self.filemenu.AppendSeparator()

        self.menuAbout = self.filemenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
        self.menuExit = self.filemenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", " Terminate the Program")

        # Setting Up the Tools Menu
        self.menubatchocnfig = self.toolsmenu.Append(wx.ID_ANY, "Batch assign configs",
                                                    "Apply current config to a batch of HDF5 files.")
        self.menubatchrun = self.toolsmenu.Append(wx.ID_ANY, "Batch Run",
                                                    "Apply current config and run deconvolution for batch of HDF5 files")
        self.menubatchextract = self.toolsmenu.Append(wx.ID_ANY, "Batch Extract",
                                                  "Run extraction on batch of HDF5 files")
        self.menubatchcre = self.toolsmenu.Append(wx.ID_ANY, "Batch Assign/Run/Extract",
                                                      "Apply current config, run, and extract on batch of HDF5 files")

        self.toolsmenu.AppendSeparator()
        self.menuExport = self.toolsmenu.Append(wx.ID_ANY, "Export Peaks Parameters and Data",
                                                "Export intensities of charge states, areas, average charge state, and other parameters for the peaks")
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
        self.menumatch = self.toolsmenu.Append(wx.ID_ANY, "Auto Match Peaks\tCtrl+M",
                                               "Run \"Match to Mixed Oligomers\" in Oligomer and Mass Tools")

        # Setting up Analysis Menu
        self.animatemenu = wx.Menu()

        self.menuanimate1 = self.animatemenu.Append(wx.ID_ANY, "Animate Zero-Charge Spectra",
                                                    "Animate 1D plots of zero-charge spectra")
        self.menuanimate15 = self.animatemenu.Append(wx.ID_ANY, "Animate Annotated Zero-Charge Spectra",
                                                     "Animate 1D plots of zero-charge mass spectra with markers")
        self.menuanimate2 = self.animatemenu.Append(wx.ID_ANY, "Animate Mass Spectra",
                                                     "Animate 1D plots of mass spectra")
        self.menuanimate25 = self.animatemenu.Append(wx.ID_ANY, "Animate Annotated Mass Spectra",
                                                    "Animate 1D plots of mass spectra with markers")
        self.menuanimate3 = self.animatemenu.Append(wx.ID_ANY, "Animate Mass v. Charge Grid",
                                                     "Animate 2D plots of mass vs. charge")
        self.menuanimate4 = self.animatemenu.Append(wx.ID_ANY, "Animate m/z v. Charge Grid",
                                                     "Animate 1D plots of m/z vs charge")
        self.analysismenu.AppendSubMenu(self.animatemenu, "Animate")
        self.analysismenu.AppendSeparator()
        self.menukendrick = self.analysismenu.Append(wx.ID_ANY, "Kendrick Mass Tools\tCtrl+K", "Kendrick Mass Analysis")
        self.menu2Dgrid = self.analysismenu.Append(wx.ID_ANY, "2D Grid Analysis", "2D Grid Analysis")
        self.menuautocorr = self.analysismenu.Append(wx.ID_ANY, "Autocorrelation",
                                                     "Autocorrelation of Mass Distribution")

        self.menufft = self.analysismenu.Append(wx.ID_ANY, "FFT Window")
        self.analysismenu.AppendSeparator()
        self.menuexpfit = self.analysismenu.Append(wx.ID_ANY, "Exponential Decay Fit",
                                                "Fit all plots to exponential decays")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_exp_fit, self.menuexpfit)
        self.menulinfit = self.analysismenu.Append(wx.ID_ANY, "Linear Fit",
                                                "Fit all plots to line")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_lin_fit, self.menulinfit)
        self.menusigfit = self.analysismenu.Append(wx.ID_ANY, "Logistic Fit",
                                                "Fit all plots to logistic equation")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_sig_fit, self.menusigfit)

        self.analysismenu.AppendSeparator()
        self.menuUM = self.analysismenu.Append(wx.ID_ANY, "Ultra Meta Data Collector")

        # Setting up the Advanced Menu
        self.menuReset = self.advancedmenu.Append(wx.ID_ANY, "Reset To Factory Default", "Reset Parameters to Default")
        self.advancedmenu.AppendSeparator()
        self.menuOpenDir = self.advancedmenu.Append(wx.ID_ANY, "Open Saved File Directory",
                                                    "Opens the save directory in the file explorer")
        self.advancedmenu.AppendSeparator()
        self.repackDirectory = self.advancedmenu.Append(wx.ID_ANY, "Repack Directory",
                                                        "Repack a directory that contains HDF5 files")

        if self.config.imflag == 0:
            self.advancedmenu.AppendSeparator()
            self.advancedmenu.Append(401, "Best Mode", "Best UniDec Deconvolution Algorithm", wx.ITEM_RADIO)
            self.advancedmenu.Append(402, "Auto Baseline", "Experimental version with automatic baseline fitting",
                                     wx.ITEM_RADIO)
            self.advancedmenu.Append(403, "Auto Baseline Subtract",
                                     "Experimental Version with automatic baseline subtraction", wx.ITEM_RADIO)
            self.parent.Bind(wx.EVT_MENU, self.menu_401_403, id=401)
            self.parent.Bind(wx.EVT_MENU, self.menu_401_403, id=402)
            self.parent.Bind(wx.EVT_MENU, self.menu_401_403, id=403)
        self.advancedmenu.AppendSeparator()
        if self.parent.tabbed == 0:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Tabbed Plots Mode",
                                                           "Put plots in individual tabs.")
        else:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Single Plot Window",
                                                           "Put plots in single large window.")

        # Experimental Menu
        self.menuundo = self.experimentalmenu.Append(wx.ID_ANY, "Undo Parameter Change\tCtrl+Z",
                                                     "Go back to the previous set of parameters")
        #self.menuredo = self.experimentalmenu.Append(wx.ID_ANY, "Redo Parameter Change\tCtrl+Y",
        #                                             "Go to the next set of parameters")
        self.experimentalmenu.AppendSeparator()
        #self.analysismenu.AppendSeparator()
        self.menuDC = self.experimentalmenu.Append(wx.ID_ANY, "Data Collector and KD Fitting")
        self.experimentalmenu.AppendSeparator()
        self.autoformat = self.experimentalmenu.Append(wx.ID_ANY, "Auto Format Monomer/Dimer",
                                                       "Mark matched monomers and dimers automatically")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_autoformat, self.autoformat)

        #self.menuAdditionalParameters = self.experimentalmenu.Append(wx.ID_ANY, "Additional Parameters",
        #                                                             "Adjust some experimental parameters")
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_additional_parameters, self.menuAdditionalParameters)
        # self.experimentalmenu.AppendSeparator()

        # self.menucal = self.experimentalmenu.Append(wx.ID_ANY, "Apply Calibration")
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_calibrate, self.menucal)

        #Help Menu
        self.getstarted = self.helpmenu.Append(wx.ID_ANY, "Getting Started")
        self.listconts = self.helpmenu.Append(wx.ID_ANY, "The Spectra Table")
        self.plotwindow = self.helpmenu.Append(wx.ID_ANY, "Plot Window")
        self.peakwindow = self.helpmenu.Append(wx.ID_ANY, "Peak Window")
        self.contmenu = wx.Menu()
        self.dataprocess = self.contmenu.Append(wx.ID_ANY, "Data Processing")
        self.unidecparas = self.contmenu.Append(wx.ID_ANY, "UniDec Parameters")
        self.additionalfilters = self.contmenu.Append(wx.ID_ANY, "Additional Filters/Restraints")
        self.peakselection = self.contmenu.Append(wx.ID_ANY, "Peak Selection, Extraction, and Plotting")
        self.additionalplotting = self.contmenu.Append(wx.ID_ANY, "Additional Plotting Parameters")
        self.helpmenu.AppendSubMenu(self.contmenu, "UniDec Controls")
        self.additionaltoolsmenu = wx.Menu()
        self.autoimport = self.additionaltoolsmenu.Append(wx.ID_ANY, "Auto Import Chromatograms")
        self.presets = self.additionaltoolsmenu.Append(wx.ID_ANY, "Presets")
        self.additionaltoolsmenu.AppendSeparator()
        self.batch = self.additionaltoolsmenu.Append(wx.ID_ANY, "Batch")
        self.peakwidthtool = self.additionaltoolsmenu.Append(wx.ID_ANY, "Peak Width Tool")
        self.oligomer = self.additionaltoolsmenu.Append(wx.ID_ANY, "Oligomer and Mass Tools")
        self.automatch = self.additionaltoolsmenu.Append(wx.ID_ANY, "Auto Match Peaks")
        self.additionaltoolsmenu.AppendSeparator()
        self.animate = self.additionaltoolsmenu.Append(wx.ID_ANY, "Animate")
        self.autocorr = self.additionaltoolsmenu.Append(wx.ID_ANY, "Autocorrelation")
        self.fft = self.additionaltoolsmenu.Append(wx.ID_ANY, "FFT Window")
        self.additionaltoolsmenu.AppendSeparator()
        self.baseline = self.additionaltoolsmenu.Append(wx.ID_ANY, "Baseline")
        self.helpmenu.AppendSubMenu(self.additionaltoolsmenu, "Other Useful Tools")
        # Set Events for Menu Bar

        # File Menu
        self.parent.Bind(wx.EVT_MENU, self.pres.on_new_file, self.newmenu)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_wizard, self.wizardmenu)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_add_file, self.addmenu)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_open, self.openmenu)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_import_mzml, self.importmenu)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_import_mzml_scans, self.importmenu2)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_import_multiple_times, self.importmenu3)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_import_multiple_scans, self.importmenu4)

        self.parent.Bind(wx.EVT_MENU, self.pres.rename_var1, self.menuvar1name)
        self.parent.Bind(wx.EVT_MENU, self.pres.rename_var2, self.menuvar2name)
        self.parent.Bind(wx.EVT_MENU, self.pres.import_vars, self.menuvar1import)
        self.parent.Bind(wx.EVT_MENU, self.pres.export_vars_dialog, self.menuvar1export)

        self.parent.Bind(wx.EVT_MENU, self.pres.on_paste_spectrum, self.menupastespectrum)
        self.parent.Bind(wx.EVT_MENU, self.pres.eng.export_spectra, self.menuexportdata)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_conf_file, self.menuLoad)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_save_default, self.menuSaveDefault)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_default, self.menuLoadDefault)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault0)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault1)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault2)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault3)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault4)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_dialog, self.menufigdialog)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_pdf, self.menuSaveFigure0)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_eps, self.menuSaveFigure1)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_small, self.menuSaveFigure1s)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_png, self.menuSaveFigure2)
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_pdf_report, self.menuSaveFigure4)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_about, self.menuAbout)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_exit, self.menuExit)

        # Tools
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_config, self.menubatchocnfig)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_run, self.menubatchrun)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_extract, self.menubatchextract)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_cre, self.menubatchcre)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_width_tool, self.menuWidth)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_peak_width, self.menuAutoWidth)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_manual, self.menuManualFile)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_mass_tools, self.menuMassFile)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_match, self.menumatch)

        # Analysis
        self.parent.Bind(wx.EVT_MENU, self.pres.on_2d_grid, self.menu2Dgrid)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_kendrick, self.menukendrick)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_autocorr_window, self.menuautocorr)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_fft_window, self.menufft)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_data_collector, self.menuDC)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_ultra_meta, self.menuUM)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_mass, self.menuanimate1)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_mz, self.menuanimate2)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_annotated_mz, self.menuanimate25)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_annotated_mass, self.menuanimate15)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_2d_mass, self.menuanimate3)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_2d_mz, self.menuanimate4)

        # Advanced
        self.parent.Bind(wx.EVT_MENU, self.pres.on_reset, self.menuReset)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_open_dir, self.menuOpenDir)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_export_params, self.menuExport)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_flip_tabbed, self.menufliptabbed)
        self.parent.Bind(wx.EVT_MENU, self.pres.recursive_repack_hdf5, self.repackDirectory)

        # Experimental
        self.parent.Bind(wx.EVT_MENU, self.pres.on_undo, self.menuundo)
        #self.parent.Bind(wx.EVT_MENU, self.pres.on_redo, self.menuredo)

        # Help
        self.parent.Bind(wx.EVT_MENU, self.pres.on_getting_started, self.getstarted)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_import_window, self.listconts)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_plot_window, self.plotwindow)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_window, self.peakwindow)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_data_processing, self.dataprocess)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_unidec_parameters, self.unidecparas)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_additional_filters, self.additionalfilters)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_selection, self.peakselection)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_additional_plotting, self.additionalplotting)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_import_help, self.autoimport)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_presets_help, self.presets)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_help, self.batch)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_width_tool_help, self.peakwidthtool)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_oligomer_help, self.oligomer)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_match_help, self.automatch)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_help, self.animate)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_autocorr_help, self.autocorr)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_fft_help, self.fft)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_baseline_help, self.baseline)

        # Setting Menu Bar
        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        self.menuBar.Append(self.analysismenu, "Analysis")
        self.menuBar.Append(self.advancedmenu, "Advanced")
        self.menuBar.Append(self.experimentalmenu, "Experimental")
        self.menuBar.Append(self.helpmenu, "Help")
        # self.Append(self.filemenu, "&File")
        # self.Append(self.toolsmenu, "Tools")
        # self.Append(self.analysismenu, "Analysis")
        # self.Append(self.advancedmenu, "Advanced")
        # self.Append(self.experimentalmenu, "Experimental")
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
        print(self.config.aggressiveflag)

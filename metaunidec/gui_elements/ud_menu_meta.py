import wx
import os
import numpy as np
import unidec_modules.isolated_packages.preset_manager as pm

class meta_menu(wx.Menu):
    def __init__(self, parent, config, pres, type="Meta"):
        super(wx.Menu, self).__init__()
        self.pres = pres
        self.config = config
        self.parent = parent
        self.type = type

        self.filemenu = wx.Menu()
        self.toolsmenu = wx.Menu()
        self.analysismenu = wx.Menu()
        self.advancedmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()
        self.helpmenu = wx.Menu()
        self.menuOpenRecent = wx.Menu()

        # File Menu
        if self.type == "Meta":
            self.openmenu = self.filemenu.Append(wx.ID_ANY, "Open HDF5 File\tCtrl+O", "Open HDF5 file")
            self.filemenu.AppendSubMenu(self.menuOpenRecent, "Open Recent File")
            self.filemenu.AppendSeparator()
            self.wizardmenu = self.filemenu.Append(wx.ID_ANY, "Import Wizard\tCtrl+N", "Import Data Into HDF5 File")
            self.filemenu.AppendSeparator()
            self.parent.Bind(wx.EVT_MENU, self.pres.on_open, self.openmenu)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_wizard, self.wizardmenu)


            self.filemenu2 = wx.Menu()
            self.newmenu = self.filemenu2.Append(wx.ID_ANY, "New File", "Create New Blank HDF5 file")
            self.addmenu = self.filemenu2.Append(wx.ID_ANY, "Add Data Files",
                                                 "Add data from individual text, mzml, or Thermo RAW files to current HDF5 file")
            self.menupastespectrum = self.filemenu2.Append(wx.ID_ANY, "Add Data From Clipboard",
                                                          "Add copied data from clipboard to current HDF5 file")
            self.menuexportdata = self.filemenu2.Append(wx.ID_ANY, "Export HDF5 to Text",
                                                          "Export data in HDF5 to individual text files")
            self.filemenu.AppendSubMenu(self.filemenu2, "Manual File Operations")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_new_file, self.newmenu)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_add_file, self.addmenu)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_paste_spectrum, self.menupastespectrum)
            self.parent.Bind(wx.EVT_MENU, self.pres.eng.export_spectra, self.menuexportdata)

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
            self.parent.Bind(wx.EVT_MENU, self.pres.on_import_mzml, self.importmenu)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_import_mzml_scans, self.importmenu2)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_import_multiple_times, self.importmenu3)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_import_multiple_scans, self.importmenu4)
            self.filemenu.AppendSeparator()

        else:
            self.openmenu = self.filemenu.Append(wx.ID_ANY, "Open mzML or Thermo Raw File\tCtrl+O", "Open mzML or Thermo Raw File")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_open, self.openmenu)

            self.openmenudir = self.filemenu.Append(wx.ID_ANY, "Open Waters or Agilent File",
                                                 "Open Waters or Agilent File")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_open_dir, self.openmenudir)
            self.filemenu.AppendSeparator()
            self.menuopenhdf5 = self.filemenu.Append(wx.ID_ANY, "Open HDF5 File",
                                                 "Open HDF5 File from Previous Deconvolution")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_open_hdf5, self.menuopenhdf5)

            self.filemenu.AppendSubMenu(self.menuOpenRecent, "Open Recent File")
            self.filemenu.AppendSeparator()
            self.menuexportdata = self.filemenu.Append(wx.ID_ANY, "Export HDF5 to Text",
                                                        "Export data in HDF5 to individual text files")
            self.parent.Bind(wx.EVT_MENU, self.pres.eng.export_spectra, self.menuexportdata)
            self.filemenu.AppendSeparator()

        # Default Submenu
        self.varmenu = wx.Menu()
        self.menuvar1name = self.varmenu.Append(wx.ID_ANY, "Rename Variable 1", "Rename Variable 1")
        self.menuvar2name = self.varmenu.Append(wx.ID_ANY, "Rename Variable 2", "Rename Variable 1")
        self.menuvar1import = self.varmenu.Append(wx.ID_ANY, "Import Variable Metadata", "Import Variable Metadata")
        self.menuvar1export = self.varmenu.Append(wx.ID_ANY, "Export Variable Metadata", "Export Variable Metadata")
        self.filemenu.AppendSubMenu(self.varmenu, "Variables")
        self.parent.Bind(wx.EVT_MENU, self.pres.rename_var1, self.menuvar1name)
        self.parent.Bind(wx.EVT_MENU, self.pres.rename_var2, self.menuvar2name)
        self.parent.Bind(wx.EVT_MENU, self.pres.import_vars, self.menuvar1import)
        self.parent.Bind(wx.EVT_MENU, self.pres.export_vars_dialog, self.menuvar1export)
        self.filemenu.AppendSeparator()

        self.menuLoad = self.filemenu.Append(wx.ID_ANY, "Load External Config File", "Load in a configuration file")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_conf_file, self.menuLoad)
        self.menuLoadDefault = self.filemenu.Append(wx.ID_ANY, "Load Default Config File",
                                                    "Load in default configuration file")
        self.menuSaveDefault = self.filemenu.Append(wx.ID_ANY, "Save Default Config File",
                                                    "Save default configuration file")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_save_default, self.menuSaveDefault)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_default, self.menuLoadDefault)
        # Default Submenu
        self.defaultmenu = wx.Menu()
        self.menuDefault2 = self.defaultmenu.Append(1001, "UniDec Default",
                                                    "General Default. Similar to High-Resolution Native.")
        self.menuDefault0 = self.defaultmenu.Append(999, "Low-resolution Native",
                                                    "General defaults for low-resolution native MS")
        self.menuDefault1 = self.defaultmenu.Append(1000, "High-resolution Native",
                                                    "Defaults for high-resolution data (Exactive EMR).")

        self.menuDefault3 = self.defaultmenu.Append(1002, "Isotopic Resolution",
                                                    "Defaults for isotopically resolved data.")
        self.menuDefault4 = self.defaultmenu.Append(1003, "Nanodiscs",
                                                    "Defaults for POPC Nanodiscs.")
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault0)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault1)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault2)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault3)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault4)
        # Custom Presets
        self.custommenu, self.masterd = pm.make_preset_menu(self.config.presetdir, topi=3500)
        self.defaultmenu.AppendSubMenu(self.custommenu, "Custom")
        for i, path, item in self.masterd:
            self.parent.Bind(wx.EVT_MENU, self.on_custom_defaults, item)
        # ..............

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
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_dialog, self.menufigdialog)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_pdf, self.menuSaveFigure0)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_eps, self.menuSaveFigure1)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_small, self.menuSaveFigure1s)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_png, self.menuSaveFigure2)
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_pdf_report, self.menuSaveFigure4)

        # Example Data
        self.examplemenu, self.masterd2 = pm.make_preset_menu(self.config.exampledatadir, exclude_dir="_unidecfiles",
                                                              topi=2500, ext="hdf5")
        self.filemenu.AppendSubMenu(self.examplemenu, "Load Example Data")
        for i, path, item in self.masterd2:
            # print(i, path, item)
            self.parent.Bind(wx.EVT_MENU, self.on_example_data, item)
        self.filemenu.AppendSeparator()

        self.menuAbout = self.filemenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
        self.menuExit = self.filemenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", " Terminate the Program")
        self.parent.Bind(wx.EVT_MENU, self.parent.on_about, self.menuAbout)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_exit, self.menuExit)


        # Setting Up the Tools Menu
        if self.type == "Meta":
            self.menubatchocnfig = self.toolsmenu.Append(wx.ID_ANY, "Batch assign configs",
                                                        "Apply current config to a batch of HDF5 files.")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_config, self.menubatchocnfig)
            self.menubatchrun = self.toolsmenu.Append(wx.ID_ANY, "Batch Run",
                                                        "Apply current config and run deconvolution for batch of HDF5 files")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_run, self.menubatchrun)
            self.menubatchextract = self.toolsmenu.Append(wx.ID_ANY, "Batch Extract",
                                                      "Run extraction on batch of HDF5 files")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_extract, self.menubatchextract)
            self.menubatchcre = self.toolsmenu.Append(wx.ID_ANY, "Batch Assign/Run/Extract",
                                                          "Apply current config, run, and extract on batch of HDF5 files")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_cre, self.menubatchcre)

            self.toolsmenu.AppendSeparator()

        else:
            self.menubatchrun = self.toolsmenu.Append(wx.ID_ANY, "Batch Run Files (Thermo/mzML/etc.)",
                                                      "Apply current config and time windows and run deconvolution for batch of  files")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_chrom1, self.menubatchrun)

            self.menubatchrun2 = self.toolsmenu.Append(wx.ID_ANY, "Batch Run Directories (Waters/Agilent)",
                                                      "Apply current config and time windows and run deconvolution for batch of  files")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_batch_chrom_dirs, self.menubatchrun2)

        self.menuExport = self.toolsmenu.Append(wx.ID_ANY, "Export Peaks Parameters and Data",
                                                "Export intensities of charge states, areas, average charge state, and other parameters for the peaks")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_export_params, self.menuExport)

        self.toolsmenu.AppendSeparator()
        self.menuAutoWidth = self.toolsmenu.Append(wx.ID_ANY, "Automatic Peak Width\tCtrl+W",
                                                   "Try to get peak width automatically.")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_peak_width, self.menuAutoWidth)
        self.menuWidth = self.toolsmenu.Append(wx.ID_ANY, "Peak Width Tool", "Help determine the peak width")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_width_tool, self.menuWidth)
        self.toolsmenu.AppendSeparator()

        self.menuManualFile = self.toolsmenu.Append(wx.ID_ANY, "Manual Assignment",
                                                    "Manually set UniDec to preassign specific m/z values")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_manual, self.menuManualFile)
        self.toolsmenu.AppendSeparator()

        self.menuMassFile = self.toolsmenu.Append(wx.ID_ANY, "Oligomer and Mass Tools\tCtrl+T",
                                                  "Oligomer and Mass Tools")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_mass_tools, self.menuMassFile)
        self.menumatch = self.toolsmenu.Append(wx.ID_ANY, "Auto Match Peaks\tCtrl+M",
                                               "Run \"Match to Mixed Oligomers\" in Oligomer and Mass Tools")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_match, self.menumatch)


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
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_mass, self.menuanimate1)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_mz, self.menuanimate2)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_annotated_mz, self.menuanimate25)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_annotated_mass, self.menuanimate15)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_2d_mass, self.menuanimate3)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_animate_2d_mz, self.menuanimate4)
        self.analysismenu.AppendSubMenu(self.animatemenu, "Animate")
        self.analysismenu.AppendSeparator()
        self.menukendrick = self.analysismenu.Append(wx.ID_ANY, "Mass Defect Tools\tCtrl+K", "Mass Defect Analysis")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_kendrick, self.menukendrick)
        self.menu2Dgrid = self.analysismenu.Append(wx.ID_ANY, "2D Grid Analysis", "2D Grid Analysis")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_2d_grid, self.menu2Dgrid)
        self.menuautocorr = self.analysismenu.Append(wx.ID_ANY, "Autocorrelation",
                                                     "Autocorrelation of Mass Distribution")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_autocorr_window_mud, self.menuautocorr)

        self.menufft = self.analysismenu.Append(wx.ID_ANY, "FFT Window")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_fft_window, self.menufft)

        if self.type == "Meta":
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
            self.parent.Bind(wx.EVT_MENU, self.pres.on_ultra_meta, self.menuUM)

        # Setting up the Advanced Menu
        self.menuReset = self.advancedmenu.Append(wx.ID_ANY, "Reset To Factory Default", "Reset Parameters to Default")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_reset, self.menuReset)
        self.advancedmenu.AppendSeparator()
        self.menuOpenDir = self.advancedmenu.Append(wx.ID_ANY, "Open Saved File Directory",
                                                    "Opens the save directory in the file explorer")
        self.parent.Bind(wx.EVT_MENU, self.parent.on_open_dir, self.menuOpenDir)

        self.advancedmenu.AppendSeparator()
        self.repackDirectory = self.advancedmenu.Append(wx.ID_ANY, "Repack Directory",
                                                        "Repack a directory that contains HDF5 files")
        self.parent.Bind(wx.EVT_MENU, self.pres.recursive_repack_hdf5, self.repackDirectory)

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

        if self.type == "Meta":
            self.advancedmenu.AppendSeparator()
            if self.parent.tabbed == 0:
                self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Tabbed Plots Mode",
                                                               "Put plots in individual tabs.")
            else:
                self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Single Plot Window",
                                                               "Put plots in single large window.")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_flip_tabbed, self.menufliptabbed)

            # Experimental Menu
            self.menuundo = self.experimentalmenu.Append(wx.ID_ANY, "Undo Parameter Change\tCtrl+Z",
                                                         "Go back to the previous set of parameters")
            # Experimental
            self.parent.Bind(wx.EVT_MENU, self.pres.on_undo, self.menuundo)

            #self.menuredo = self.experimentalmenu.Append(wx.ID_ANY, "Redo Parameter Change\tCtrl+Y",
            #                                             "Go to the next set of parameters")
            # self.parent.Bind(wx.EVT_MENU, self.pres.on_redo, self.menuredo)

            self.experimentalmenu.AppendSeparator()
            #self.analysismenu.AppendSeparator()
            self.menuDC = self.experimentalmenu.Append(wx.ID_ANY, "Data Collector and KD Fitting")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_data_collector, self.menuDC)
            self.experimentalmenu.AppendSeparator()
            self.autoformat = self.experimentalmenu.Append(wx.ID_ANY, "Auto Format Monomer/Dimer",
                                                           "Mark matched monomers and dimers automatically")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_autoformat, self.autoformat)
            self.experimentalmenu.AppendSeparator()

        self.menuwaterfall = self.experimentalmenu.Append(wx.ID_ANY, "Waterfall Plot",
                                                       "Make 3D Waterfall plot with mass distributions")
        self.parent.Bind(wx.EVT_MENU, self.pres.make_waterfall_plots, self.menuwaterfall)

        self.experimentalmenu.AppendSeparator()
        self.menulinreg = self.experimentalmenu.Append(wx.ID_ANY, "Linear Regression")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_linreg, self.menulinreg)
        self.menusubdiv = self.experimentalmenu.Append(wx.ID_ANY, "Subtract and Divide")
        self.parent.Bind(wx.EVT_MENU, self.pres.sub_div, self.menusubdiv)
        #self.menuAdditionalParameters = self.experimentalmenu.Append(wx.ID_ANY, "Additional Parameters",
        #                                                             "Adjust some experimental parameters")
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_additional_parameters, self.menuAdditionalParameters)
        self.experimentalmenu.AppendSeparator()

        self.menuscanpeaks = self.experimentalmenu.Append(wx.ID_ANY, "Get Scan Peaks")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_pick_scanpeaks, self.menuscanpeaks)

        self.menufilterpeaks = self.experimentalmenu.Append(wx.ID_ANY, "Filter Peaks by Score")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_filter_peaks_MUD, self.menufilterpeaks)

        if self.type == "Meta":
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
        if self.type == "Meta":
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
            self.config.default_decon_params()
        elif nid == 2:
            self.config.default_isotopic_res()
        elif nid == 3:
            self.config.default_nanodisc()
        elif nid == 99:
            self.config.default_low_res()
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

    def on_custom_defaults(self, e):
        # print("Clicked", e)
        try:
            nid = e.GetId()
            ids = self.masterd[:, 0].astype(np.float)
            pos = np.argmin(np.abs(ids - nid))
            path = self.masterd[pos, 1]
            print("Opening Config:", path)
            #self.pres.eng.load_config(path)
            self.pres.import_config(path)
        except Exception as e:
            print(e)

    def on_example_data(self, e):
        # print("Clicked", e)
        try:
            nid = e.GetId()
            ids = self.masterd2[:, 0].astype(np.float)
            pos = np.argmin(np.abs(ids - nid))
            path = self.masterd2[pos, 1]
            print("Opening Path:", path)
            self.pres.open_file(path)
        except Exception as e:
            print(e)

    def update_recent(self, ext=".hdf5"):
        menu_list = self.menuOpenRecent.GetMenuItems()
        for i in range(len(menu_list) - 1, -1, -1):  # clear menu
            self.menuOpenRecent.DestroyItem(menu_list[i])

        max_items = 5  # can be changed to whatever
        added = 0
        for file_path in self.pres.recent_files:
            if added >= max_items:
                break
            filename = os.path.basename(file_path)
            if os.path.splitext(filename)[1] == ext:
                self.add_to_recent(file_path)
                added += 1

    def add_to_recent(self, file_path):
        new_item = self.menuOpenRecent.Append(wx.ID_ANY, os.path.basename(file_path))
        self.parent.Bind(wx.EVT_MENU, lambda e: self.pres.open_file(file_path), new_item)

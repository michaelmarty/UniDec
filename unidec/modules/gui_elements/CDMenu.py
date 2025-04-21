

import wx
import unidec.modules.isolated_packages.preset_manager as pm

import numpy as np
import os


class CDMenu(wx.Menu):
    # noinspection PyMissingConstructor
    def __init__(self, parent, config, pres, tabbed, htmode=False):
        super(wx.Menu, self).__init__()
        self.pres = pres
        self.config = config
        self.parent = parent
        self.tabbed = tabbed
        self.htmode = htmode

        self.filemenu = wx.Menu()
        self.toolsmenu = wx.Menu()
        self.analysismenu = wx.Menu()
        self.advancedmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()
        self.menuOpenRecent = wx.Menu()

        # File Menu
        self.menuOpen = self.filemenu.Append(wx.ID_OPEN, "Open File (Text, mzML, or Thermo RAW)\tCtrl+O",
                                             " Open a Text File in x y text, mzML, or Thermo RAW format")
        self.filemenu.AppendSubMenu(self.menuOpenRecent, "Open Recent File")
        self.filemenu.AppendSeparator()

        self.menuLoad = self.filemenu.Append(wx.ID_ANY, "Load External Config File", "Load in a configuration file")

        if self.htmode:
            # Load chrom file
            self.filemenu.AppendSeparator()
            self.menuLoadChrom = self.filemenu.Append(wx.ID_ANY, "Load Chromatogram File",
                                                      "Load in a chromatogram file")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_load_chroms, self.menuLoadChrom)
            self.filemenu.AppendSeparator()

        '''
        self.menuLoadEverything = self.filemenu.Append(wx.ID_ANY, "Load Prior State for Current File\tCtrl+L",
                                                       "Load past deconvolution results from the current opened file")
        self.menuLoadState = self.filemenu.Append(wx.ID_ANY, "Load Zip State File", "Load state from zip file")
        self.menuSaveState = self.filemenu.Append(wx.ID_ANY, "Save Zip State File\tCtrl+S",
                                                  "Save program state to fresh folder")
        self.filemenu.AppendSeparator()
        
        self.menuLoadDefault = self.filemenu.Append(wx.ID_ANY, "Load Default Config File",
                                                    "Load in default configuration file")
        self.menuSaveDefault = self.filemenu.Append(wx.ID_ANY, "Save Default Config File",
                                                    "Save default configuration file")
        
        # Default Submenu
        
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
        '''
        # Custom Presets
        self.defaultmenu, self.masterd = pm.make_preset_menu(self.config.presetdirCD)
        for i, path, item in self.masterd:
            # print(i, path, item)
            self.parent.Bind(wx.EVT_MENU, self.on_custom_defaults, item)
        # ..............

        self.filemenu.AppendSubMenu(self.defaultmenu, "Presets")
        self.filemenu.AppendSeparator()
        self.menuSaveFigureHTML = self.filemenu.Append(wx.ID_ANY, "Generate HTML Report\tCtrl+H",
                                                       "Generate HTML Report")
        self.parent.Bind(wx.EVT_MENU, self.parent.pres.on_gen_html_report, self.menuSaveFigureHTML)
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

        # Example Data

        self.examplemenu, self.masterd2 = pm.make_preset_menu(self.config.exampledatadirCD, exclude_dir="_unidecfiles",
                                                              topi=2500, exclude_ext="hdf5")

        keys = []
        for i, d in enumerate(self.masterd2):
            if i < 10:
                keys.append([str(i + 1), self.pres.on_ex, d[2]])
        self.menukeys = keys

        self.filemenu.AppendSubMenu(self.examplemenu, "Load Example Data")
        for i, path, item in self.masterd2:
            # print(i, path, item)
            self.parent.Bind(wx.EVT_MENU, self.on_example_data, item)
        # ..............
        self.filemenu.AppendSeparator()

        self.menuAbout = self.filemenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
        self.menuExit = self.filemenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", " Terminate the Program")

        # Setting Up the Tools Menu
        self.menuBatch = self.toolsmenu.Append(wx.ID_ANY, "Batch Process\tCtrl+B",
                                               "Run these conditions on multiple samples")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch, self.menuBatch)
        self.toolsmenu.AppendSeparator()

        self.menusmashwindow = self.toolsmenu.Append(wx.ID_ANY, "Select Noise Peaks")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_smash_window, self.menusmashwindow)
        self.toolsmenu.AppendSeparator()

        self.menuMassFile = self.toolsmenu.Append(wx.ID_ANY, "Oligomer and Mass Tools\tCtrl+T",
                                                  "Oligomer and Mass Tools")
        self.toolsmenu.AppendSeparator()
        self.menuAutoWidth = self.toolsmenu.Append(wx.ID_ANY, "Automatic Peak Width\tCtrl+W",
                                                   "Try to get peak width automatically.")
        self.menuWidth = self.toolsmenu.Append(wx.ID_ANY, "Peak Width Tool", "Help determine the peak width")
        self.toolsmenu.AppendSeparator()

        self.menucal = self.toolsmenu.Append(wx.ID_ANY, "Calibration Tool")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_calibrate, self.menucal)

        if self.htmode:
            self.toolsmenu.AppendSeparator()
            self.menuexportHT = self.toolsmenu.Append(wx.ID_ANY, "Export Chromatograms")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_export_arrays, self.menuexportHT)

            self.menuplotkernel = self.toolsmenu.Append(wx.ID_ANY, "Plot Kernel")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_plot_kernel, self.menuplotkernel)

            self.menuautocycle = self.toolsmenu.Append(wx.ID_ANY, "Print Optimal Cycle Index")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_set_ct, self.menuautocycle)

        self.toolsmenu.AppendSeparator()
        self.menustori = self.toolsmenu.Append(wx.ID_ANY, "Convert STORI Folder of CSVs")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_stori, self.menustori)

        # Tools
        self.parent.Bind(wx.EVT_MENU, self.pres.on_mass_tools, self.menuMassFile)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_peak_width_tool, self.menuWidth)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_auto_peak_width, self.menuAutoWidth)

        '''
        

        self.menuManualFile = self.toolsmenu.Append(wx.ID_ANY, "Manual Assignment",
                                                    "Manually set UniDec to preassign specific m/z values")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_manual, self.menuManualFile)
        self.toolsmenu.AppendSeparator()
        '''

        # Setting up Analysis Menu
        self.menuPlotZ = self.analysismenu.Append(wx.ID_ANY, "Native Charge/Mass Tools",
                                                  "Tool for exploring relationship between charge and native mass "
                                                  "and extraction specific distributions")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_nativez_tools, self.menuPlotZ)
        self.analysismenu.AppendSeparator()
        self.menuExport = self.analysismenu.Append(wx.ID_ANY, "Export Peaks Parameters and Data",
                                                   "Export intensities of charge states, areas, average charge state, "
                                                   "and other parameters for the peaks")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_export_params, self.menuExport)
        self.menuFitNorm = self.analysismenu.Append(wx.ID_ANY, "Fit Peak Intensities",
                                                    "Fits masses and reports normalized and relative peak intensities")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_fit_masses, self.menuFitNorm)
        self.analysismenu.AppendSeparator()

        self.menukendrick = self.analysismenu.Append(wx.ID_ANY, "Mass Defect Tools\tCtrl+K", "Mass Defect Analysis")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_kendrick, self.menukendrick)

        self.menu2Dgrid = self.analysismenu.Append(wx.ID_ANY, "2D Grid Analysis", "2D Grid Analysis")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_2d_grid, self.menu2Dgrid)
        self.menuautocorr = self.analysismenu.Append(wx.ID_ANY, "Autocorrelation",
                                                     "Autocorrelation of Mass Distribution")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_autocorr_window, self.menuautocorr)
        self.analysismenu.AppendSeparator()

        self.menuintegrate = self.analysismenu.Append(wx.ID_ANY, "Integrate Peaks\tCtrl+I",
                                                      "Peak area with range set by Peak Detection Range or "
                                                      "Integration Range")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_integrate, self.menuintegrate)
        self.menumatch = self.analysismenu.Append(wx.ID_ANY, "Auto Match Peaks\tCtrl+M",
                                                  "Run \"Match to Mixed Oligomers\" in Oligomer and Mass Tools")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_match, self.menumatch)
        self.menucom = self.analysismenu.Append(wx.ID_ANY, "Report Center of Mass",
                                                "Reports center of mass for the zoomed region of the zero-charge "
                                                "mass spectrum")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_center_of_mass, self.menucom)

        self.maxcharge = self.analysismenu.Append(wx.ID_ANY, "Label Max Charge States",
                                                  "Labels the maximum charge state in each distribution")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_label_max_charge_states, self.maxcharge)

        self.maxcharge = self.analysismenu.Append(wx.ID_ANY, "Label Average Charge States",
                                                  "Labels the average charge state in each distribution")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_label_avg_charge_states, self.maxcharge)

        self.analysismenu.AppendSeparator()
        self.menubiocalc = self.analysismenu.Append(wx.ID_ANY, "Protein/RNA/DNA Mass Calculator")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_biopolymer, self.menubiocalc)

        self.analysismenu.AppendSeparator()

        '''
        # Analysis
        
        self.parent.Bind(wx.EVT_MENU, self.pres.on_data_collector, self.menucollect)
        
        
        self.parent.Bind(wx.EVT_MENU, self.pres.on_plot_offsets, self.menuoffset)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_charge_plot, self.menuchargeplot)

        self.menucollect = self.analysismenu.Append(wx.ID_ANY, "Data Collector Utility",
                                                    "Collect results and extract values.  Experimental KD fitting.")
        self.analysismenu.AppendSeparator()

        

        self.menufft = self.analysismenu.Append(wx.ID_ANY, "FFT Window")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_fft_window, self.menufft)
        
        self.analysismenu.AppendSeparator()

        self.menuchargeplot = self.analysismenu.Append(wx.ID_ANY, "Plot By Charge",
                                                       "Plots Mass Distributions as a Function of Charge")
        self.menuoffset = self.analysismenu.Append(wx.ID_ANY, "Plot Charge Offsets",
                                                   "Plots Mass vs. Charge Offset in 2D Plot")

        if self.config.imflag == 1:
            self.analysismenu.AppendSeparator()
            self.menuimtools = self.analysismenu.Append(wx.ID_ANY, "IM Parameters Tool",
                                                        "Tools for estimating IM parameters.")
            self.menuimtools2 = self.analysismenu.Append(wx.ID_ANY, "IM Extraction Tool",
                                                         "Tools for Extraction IM Results.")
            self.parent.Bind(wx.EVT_MENU, self.pres.on_im_tools, self.menuimtools)
            self.parent.Bind(wx.EVT_MENU, self.pres.on_im_extract, self.menuimtools2)
        '''

        # Setting up the Advanced Menu
        self.menuReset = self.advancedmenu.Append(wx.ID_ANY, "Reset To Factory Default", "Reset Parameters to Default")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_reset, self.menuReset)
        self.advancedmenu.AppendSeparator()





        self.menuUnidecPath = self.advancedmenu.Append(wx.ID_FILE1, "UniDec File", "Find the UniDec executable file.")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_unidec_path, self.menuUnidecPath)
        # self.menuFileName = self.advancedmenu.Append(wx.ID_FILE2, "Rename Files",
        #                                             "Rename the files into and out of unidec.")
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_file_name, self.menuFileName)
        self.menuOpenDir = self.advancedmenu.Append(wx.ID_ANY, "Open Saved File Directory",
                                                    "Opens the save directory in the file explorer")



        self.parent.Bind(wx.EVT_MENU, self.parent.on_open_dir, self.menuOpenDir)

        self.advancedmenu.AppendSeparator()

        self.scalemenu = wx.Menu()

        self.scalemenu.Append(501, "Linear", "Normal linear intensity scale", wx.ITEM_RADIO)
        self.scalemenu.Append(502, "Logarithmic", "Logarithmic intensity scale",
                              wx.ITEM_RADIO)
        self.scalemenu.Append(503, "Square Root", "Square root intensity scale", wx.ITEM_RADIO)
        self.parent.Bind(wx.EVT_MENU, self.menu_501_503, id=501)
        self.parent.Bind(wx.EVT_MENU, self.menu_501_503, id=502)
        self.parent.Bind(wx.EVT_MENU, self.menu_501_503, id=503)

        self.advancedmenu.AppendSubMenu(self.scalemenu, 'Intensity Scale')
        self.advancedmenu.AppendSeparator()

        if self.tabbed == 0:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Tabbed Plots Mode",
                                                           "Put plots in individual tabs.")
        else:
            self.menufliptabbed = self.advancedmenu.Append(wx.ID_ANY, "Switch to Single Plot Window",
                                                           "Put plots in single large window.")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_flip_tabbed, self.menufliptabbed)
        self.advancedmenu.AppendSeparator()
        self.gpumenu = wx.Menu()
        self.gpumenu.Append(603, "Exe", "Uses Exe", wx.ITEM_RADIO)
        self.gpumenu.Append(601, "Py", "Uses Python", wx.ITEM_RADIO)
        # self.gpumenu.Append(602, "CuPy", "Turns on GPU Acceleration with Python", wx.ITEM_RADIO)

        self.parent.Bind(wx.EVT_MENU, self.menu_601_602, id=601)
        # self.parent.Bind(wx.EVT_MENU, self.menu_601_602, id=602)
        self.parent.Bind(wx.EVT_MENU, self.menu_601_602, id=603)
        self.advancedmenu.AppendSubMenu(self.gpumenu, 'EXE Mode')

        # Experimental Menu
        self.menuundo = self.experimentalmenu.Append(wx.ID_ANY, "Undo Parameter Change\tCtrl+Z",
                                                     "Go back to the previous set of parameters")
        self.menuredo = self.experimentalmenu.Append(wx.ID_ANY, "Redo Parameter Change\tCtrl+Y",
                                                     "Go to the next set of parameters")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_undo, self.menuundo)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_redo, self.menuredo)
        self.experimentalmenu.AppendSeparator()
        self.menuscore = self.experimentalmenu.Append(wx.ID_ANY, "Filter Peak Scores", "Filter Peak Scores")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_score2, self.menuscore)

        self.menuscore2 = self.experimentalmenu.Append(wx.ID_ANY, "Peak Scores Window", "Launch Peak Scores Window")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_score_window, self.menuscore2)

        self.menuscore3 = self.experimentalmenu.Append(wx.ID_ANY, "Label Peak Scores", "Label Peak Scores on Plot")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_score_label, self.menuscore3)

        self.experimentalmenu.AppendSeparator()
        self.menulinreg = self.experimentalmenu.Append(wx.ID_ANY, "Linear Regression")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_linreg, self.menulinreg)
        self.menusubdiv = self.experimentalmenu.Append(wx.ID_ANY, "Subtract and Divide")
        self.parent.Bind(wx.EVT_MENU, self.pres.sub_div_window, self.menusubdiv)

        self.experimentalmenu.AppendSeparator()
        self.menunoise = self.experimentalmenu.Append(wx.ID_ANY, "Show Noise Level")
        self.parent.Bind(wx.EVT_MENU, self.pres.plotnoise, self.menunoise)
        self.menukernel = self.experimentalmenu.Append(wx.ID_ANY, "Plot Deconvolution Kernel PSF")
        self.parent.Bind(wx.EVT_MENU, self.pres.plotkernel, self.menukernel)

        self.experimentalmenu.AppendSeparator()
        self.menucompare = self.experimentalmenu.Append(wx.ID_ANY, "Import Mass Data to Compare")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_compare, self.menucompare)
        self.menucompare2 = self.experimentalmenu.Append(wx.ID_ANY, "Plot Compare")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_compare2, self.menucompare2)
        self.menucompare3 = self.experimentalmenu.Append(wx.ID_ANY, "Compare Processed and Unprocessed")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_compare_proc, self.menucompare3)

        self.experimentalmenu.AppendSeparator()
        self.menurefresh = self.experimentalmenu.Append(wx.ID_ANY, "Refresh Acquiring Data")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_refresh, self.menurefresh)
        '''    
        self.menucolor1d = self.experimentalmenu.Append(wx.ID_ANY, "Color Plots",
                                                        "Make a Different Colored 1D Plot")

        self.parent.Bind(wx.EVT_MENU, self.pres.on_additional_parameters, self.menuAdditionalParameters)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_color_plot1d, self.menucolor1d)


        self.experimentalmenu.AppendSeparator()

        self.experimentalmenu.AppendSeparator()
        self.menuisotopes = self.experimentalmenu.Append(wx.ID_ANY, "Plot Averagine Isotope Distributions")
        self.parent.Bind(wx.EVT_MENU, self.pres.on_plot_isotope_distribution, self.menuisotopes)

        self.menupdi = self.experimentalmenu.Append(wx.ID_ANY, "Print Polydispersity Index")
        self.parent.Bind(wx.EVT_MENU, self.pres.eng.polydispersity_index, self.menupdi)

        self.experimentalmenu.AppendSeparator()

        self.menutheomass = self.experimentalmenu.Append(wx.ID_ANY, "Plot Theoretical Mass")
        self.parent.Bind(wx.EVT_MENU, self.pres.plot_theo_mass, self.menutheomass)

        # self.menucentroid = self.experimentalmenu.Append(wx.ID_ANY, "Get Centroid at FWHM")
        # self.parent.Bind(wx.EVT_MENU, self.pres.on_centroid, self.menucentroid)

        self.experimentalmenu.AppendSeparator()
        self.menuRegister = self.experimentalmenu.Append(wx.ID_ANY, "Fix Agilent Imports",
                                                         "Registers the Agilent Interfaces. "
                                                         "MUST RUN AS ADMINISTRATOR")
        self.parent.Bind(wx.EVT_MENU, self.pres.register, self.menuRegister)

        
        '''
        # Set Events for Menu Bar

        # File Menu
        self.parent.Bind(wx.EVT_MENU, self.pres.on_open, self.menuOpen)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_conf_file, self.menuLoad)
        '''
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_state, self.menuLoadState)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_everything, self.menuLoadEverything)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_save_state, self.menuSaveState)
        
        self.parent.Bind(wx.EVT_MENU, self.pres.on_save_default, self.menuSaveDefault)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_load_default, self.menuLoadDefault)'''
        '''
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault0)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault1)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault2)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault3)
        self.parent.Bind(wx.EVT_MENU, self.on_defaults, self.menuDefault4)'''
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_dialog, self.menufigdialog)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_pdf, self.menuSaveFigure0)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_eps, self.menuSaveFigure1)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_small, self.menuSaveFigure1s)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_save_figure_png, self.menuSaveFigure2)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_about, self.menuAbout)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_exit, self.menuExit)

        # Setting Menu Bar
        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        self.menuBar.Append(self.analysismenu, "Analysis")
        self.menuBar.Append(self.advancedmenu, "Advanced")
        self.menuBar.Append(self.experimentalmenu, "Experimental")
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
        print(self.config.intscale)
        self.pres.on_replot()

    def menu_601_602(self, event):
        """
        Menu function to adjust the intensity scale.
        :param event: wx Event
        :return: None
        """
        event_id = event.GetId()
        exemode = True
        if event_id == 601:
            # gpumode = False
            exemode = False
        # if event_id == 602:
        #    gpumode = True
        #    exemode = False
        if event_id == 603:
            # gpumode = False
            exemode = True
        print("Exe Mode:", exemode)
        self.pres.on_exe_mode(exemode)

    def on_custom_defaults(self, e):
        # print("Clicked", e)
        try:
            nid = e.GetId()
            ids = self.masterd[:, 0].astype(float)
            pos = np.argmin(np.abs(ids - nid))
            path = self.masterd[pos, 1]
            print("Opening Config:", path)
            # self.pres.eng.load_config(path)
            self.pres.import_config(path)
        except Exception as e:
            print(e)

    def on_example_data(self, e):
        # print("Clicked", e)
        try:
            nid = e.GetId()
            ids = self.masterd2[:, 0].astype(float)
            pos = np.argmin(np.abs(ids - nid))
            self.load_example_data(pos)
        except Exception as e:
            print(e)

    def load_example_data(self, pos):
        path = self.masterd2[pos, 1]
        dirname = os.path.dirname(path)
        file = os.path.split(path)[1]
        print("Opening Path:", path, file, dirname)
        self.pres.on_open_file(file, dirname)

    def update_recent(self):
        menu_list = self.menuOpenRecent.GetMenuItems()
        for i in range(len(menu_list) - 1, -1, -1):  # clear menu
            self.menuOpenRecent.DestroyItem(menu_list[i])

        max_items = 5  # can be changed to whatever
        added = 0
        for file_path in self.pres.recent_files:
            if added >= max_items:
                break
            filename = os.path.basename(file_path)
            if os.path.splitext(filename)[1] != ".hdf5":
                self.add_to_recent(file_path)
                added += 1

    def add_to_recent(self, file_path):
        # This needs to be separate from update_recent() for binding to work
        filename = os.path.basename(file_path)
        dirname = os.path.dirname(file_path)
        new_item = self.menuOpenRecent.Append(wx.ID_ANY, filename)
        self.parent.Bind(wx.EVT_MENU, lambda e: self.pres.on_open_file(filename, dirname), new_item)

    def on_reset(self, event):
        wx.MessageBox("Reset to Factory Default", "Info", wx.OK | wx.ICON_INFORMATION)


import wx
import wx.lib.agw.foldpanelbar as fpb
import unidec.tools as ud
import numpy as np
import wx.lib.scrolledpanel as scrolled
from unidec.modules.HTEng import HTseqDict


class main_controls(wx.Panel):  # scrolled.ScrolledPanel):
    # noinspection PyMissingConstructor
    def __init__(self, parent, config, pres, panel, iconfile, htmode=False):
        super(wx.Panel, self).__init__(panel)
        # super(scrolled.ScrolledPanel, self).__init__(parent=panel, style=wx.ALL | wx.EXPAND)
        # self.SetAutoLayout(1)
        # self.SetupScrolling(scroll_x=False)
        self.htmode = htmode
        self.parent = parent
        self.config = config
        self.pres = pres
        self.backgroundchoices = self.config.backgroundchoices
        self.psigsettings = [0, 1, 5, 10]
        self.betasettings = [0, 1, 5, 10]
        self.update_flag = True

        # Get a few tool bar icons
        tsize = (16, 16)
        try:
            self.open_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize)
            self.next_bmp = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD, wx.ART_TOOLBAR, tsize)
            self.report_bmp = wx.ArtProvider.GetBitmap(wx.ART_LIST_VIEW, wx.ART_TOOLBAR, tsize)
            self.A_bmp = wx.ArtProvider.GetBitmap(wx.ART_HELP_SETTINGS, wx.ART_TOOLBAR, tsize)
            try:
                self.ud_bmp = wx.Bitmap(wx.Image(iconfile).Rescale(tsize[0], tsize[1]))
            except Exception as ex:
                self.ud_bmp = wx.ArtProvider.GetBitmap(wx.ART_HELP_SETTINGS, wx.ART_TOOLBAR, tsize)
                print(ex)
        except Exception as ex:
            self.open_bmp = None
            self.next_bmp = None
            self.report_bmp = None
            self.A_bmp = None
            self.ud_bmp = None
            print(ex)

        # ..........................
        #
        # Sizers for Controls
        #
        # .............................
        sizercontrol = wx.BoxSizer(wx.VERTICAL)

        # Small Toolbar
        buttonsizer = wx.BoxSizer(wx.HORIZONTAL)
        bsize = (54, 25)
        self.openbutton = wx.BitmapButton(self, -1, self.open_bmp, size=bsize)
        self.procbutton = wx.BitmapButton(self, -1, self.next_bmp, size=bsize)
        self.procbutton.SetBackgroundColour(wx.Colour(150, 150, 255))
        self.udbutton = wx.BitmapButton(self, -1, self.ud_bmp, size=bsize)
        self.udbutton.SetBackgroundColour(wx.Colour(255, 255, 150))
        self.ppbutton = wx.BitmapButton(self, -1, self.report_bmp, size=bsize)
        self.ppbutton.SetBackgroundColour(wx.Colour(255, 150, 150))
        self.autobutton = wx.Button(self, -1, "All", size=bsize)  # self.A_bmp
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_open, self.openbutton)
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_dataprep_button, self.procbutton)
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_unidec_button, self.udbutton)
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_pick_peaks, self.ppbutton)
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_auto, self.autobutton)
        buttons = [self.openbutton, self.procbutton, self.udbutton, self.ppbutton, self.autobutton]
        for button in buttons:
            buttonsizer.Add(button, 0, wx.EXPAND)
        sizercontrol.Add(buttonsizer, 0, wx.EXPAND)

        # Setting up main fold controls
        self.scrolledpanel = scrolled.ScrolledPanel(self, style=wx.ALL | wx.EXPAND)
        self.scrolledpanel.SetupScrolling()
        size1 = (75, -1)
        self.foldpanels = fpb.FoldPanelBar(self.scrolledpanel, -1, size=(250, 3400), agwStyle=fpb.FPB_VERTICAL)
        style1 = fpb.CaptionBarStyle()
        style1b = fpb.CaptionBarStyle()
        style1c = fpb.CaptionBarStyle()
        style2 = fpb.CaptionBarStyle()
        style2b = fpb.CaptionBarStyle()
        style3 = fpb.CaptionBarStyle()
        style3b = fpb.CaptionBarStyle()
        styleht = fpb.CaptionBarStyle()

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

        styleht.SetSecondColour(wx.Colour(0, bright4, 0))

        self.imflag = 0

        # Panel to set the data prep controls
        foldpanel1 = self.foldpanels.AddFoldPanel(caption="Data Processing", collapsed=False, cbstyle=style1)
        panel1 = wx.Panel(foldpanel1, -1)
        sizercontrol1 = wx.GridBagSizer(wx.VERTICAL)

        self.ctlminmz = wx.TextCtrl(panel1, value="", size=(50, -1))
        self.ctlmaxmz = wx.TextCtrl(panel1, value="", size=(60, -1))
        mzrange = wx.BoxSizer(wx.HORIZONTAL)
        if self.config.imflag == 1:
            mzrange.Add(wx.StaticText(panel1, label="               "), 0, wx.ALIGN_CENTER_VERTICAL)
        mzrange.Add(wx.StaticText(panel1, label="m/z: "), 0, wx.ALIGN_CENTER_VERTICAL)
        mzrange.Add(self.ctlminmz)
        mzrange.Add(wx.StaticText(panel1, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        mzrange.Add(self.ctlmaxmz)
        mzrange.Add(wx.StaticText(panel1, label=" Th "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.fullbutton = wx.Button(panel1, -1, "Full", size=(40, 25))
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_full, self.fullbutton)
        mzrange.Add(self.fullbutton)
        i = 0
        sizercontrol1.Add(mzrange, (i, 0), span=(1, 5))
        i += 1

        self.ctlstartz = wx.TextCtrl(panel1, value="", size=(60, -1))
        self.ctlendz = wx.TextCtrl(panel1, value="", size=(60, -1))
        zrange = wx.BoxSizer(wx.HORIZONTAL)
        zrange.Add(wx.StaticText(panel1, label="Charge Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
        zrange.Add(self.ctlstartz)
        zrange.Add(wx.StaticText(panel1, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        zrange.Add(self.ctlendz)
        sizercontrol1.Add(zrange, (i, 0), span=(1, 2))
        i += 1

        self.ctlbinsize = wx.TextCtrl(panel1, value="", size=size1)
        sizercontrol1.Add(self.ctlbinsize, (i, 1))
        sizercontrol1.Add(wx.StaticText(panel1, label="m/z Bin Size: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.subtypectl = wx.Choice(panel1, -1, choices=["Slope: S/N per z", "Slope: Intensity per z"])
        self.subtypectl.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
        self.subtypectl.SetSelection(0)
        sizercontrol1.Add(self.subtypectl, (i, 0))
        self.ctlslope = wx.TextCtrl(panel1, value="", size=size1)
        sizercontrol1.Add(self.ctlslope, (i, 1))
        # sizercontrol1.Add(wx.StaticText(panel1, label="Slope (Intensity/z)"), (i, 0),
        #                  flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlmassbins = wx.TextCtrl(panel1, value="", size=size1)
        sizercontrol1.Add(self.ctlmassbins, (i, 1), span=(1, 2))
        sizercontrol1.Add(wx.StaticText(panel1, label="Mass Bin Size (Da): "), (i, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if htmode:
            # Control for scan compression
            self.ctlscancompress = wx.TextCtrl(panel1, value="", size=size1)
            sizercontrol1.Add(self.ctlscancompress, (i, 1), span=(1, 2))
            sizercontrol1.Add(wx.StaticText(panel1, label="Scan Compression: "), (i, 0),
                              flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctldatanorm = wx.CheckBox(panel1, label="Normalize Data")
        self.ctldatanorm.SetValue(True)
        sizercontrol1.Add(self.ctldatanorm, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.dataprepbutton = wx.Button(panel1, -1, "Process Data")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_dataprep_button, self.dataprepbutton)
        sizercontrol1.Add(self.dataprepbutton, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        panel1.SetSizer(sizercontrol1)
        sizercontrol1.Fit(panel1)
        self.foldpanels.AddFoldPanelWindow(foldpanel1, panel1, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel1, wx.StaticText(foldpanel1, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Additional Data Processing Parameters
        foldpanel1b = self.foldpanels.AddFoldPanel(caption="Additional Data Processing Parameters", collapsed=True,
                                                   cbstyle=style1b)
        panel1b = wx.Panel(foldpanel1b, -1)
        gbox1b = wx.GridBagSizer(wx.VERTICAL)
        i = 0

        self.ctlscanstart = wx.TextCtrl(panel1b, value="", size=(60, -1))
        self.ctlscanend = wx.TextCtrl(panel1b, value="", size=(60, -1))
        scanrange = wx.BoxSizer(wx.HORIZONTAL)
        scanrange.Add(wx.StaticText(panel1b, label="Scan Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
        scanrange.Add(self.ctlscanstart)
        scanrange.Add(wx.StaticText(panel1b, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        scanrange.Add(self.ctlscanend)
        gbox1b.Add(scanrange, (i, 0), span=(1, 2))
        i += 1

        self.ctlzbins = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlzbins, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Charge Bin Size: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlcentroid = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlcentroid, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Centroid Filtering Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlsmooth = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlsmooth, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Gaussian Smoothing: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1
        gbox1b.Add(wx.StaticText(panel1b, label=" m/z"), (i - 1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox1b.Add(wx.StaticText(panel1b, label=" Z"), (i, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlsmoothdt = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlsmoothdt, (i, 1))
        i += 1

        self.ctlintthresh = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlintthresh, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Intensity Threshold: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlprethresh = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlprethresh, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Import Threshold: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        # Bind change in this text box to update function
        self.parent.Bind(wx.EVT_TEXT, self.update_prethresh, self.ctlprethresh)
        i += 1

        self.ctldatareductionpercent = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctldatareductionpercent, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Data Reduction (%): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        gbox1b.Add(wx.StaticText(panel1b, label="Background Sub.:"), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlbuff = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlbuff, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label=" m/z"), (i, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1
        self.ctlbuff.SetToolTip(wx.ToolTip(
            "Background subtraction width of smoothed background in m/z data points"))

        self.ctlsubbuffdt = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlsubbuffdt, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label=" Z"), (i, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1
        self.ctlbuff.SetToolTip(wx.ToolTip(
            "Background subtraction width of smoothed background in Z data points"))

        self.ctladductmass = wx.TextCtrl(panel1b, value='', size=size1)
        gbox1b.Add(self.ctladductmass, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Adduct Mass (Da): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlnegmode = wx.CheckBox(panel1b, label="Negative Mode")
        gbox1b.Add(self.ctlnegmode, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.parent.Bind(wx.EVT_CHECKBOX, self.export_gui_to_config, self.ctlnegmode)
        i += 1

        self.ctlsmashflag = wx.CheckBox(panel1b, label="Remove Noise Peaks")
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_check_smash, self.ctlsmashflag)
        gbox1b.Add(self.ctlsmashflag, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlsmashflag.SetToolTip(wx.ToolTip("Remove Noise Peaks. See Tools>Select Noise Peaks"))
        i += 1

        # Check box for inv inj time
        self.ctlinjtime = wx.CheckBox(panel1b, label="DMT Inverse Injection Time")
        gbox1b.Add(self.ctlinjtime, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlinjtime.SetToolTip(wx.ToolTip("Use inverse injection time for histogram."))
        i += 1

        gbox1b.Add(wx.StaticText(panel1b, label=""), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        panel1b.SetSizer(gbox1b)
        gbox1b.Fit(panel1b)

        self.foldpanels.AddFoldPanelWindow(foldpanel1b, panel1b, fpb.FPB_ALIGN_WIDTH)
        if htmode:
            self.foldpaneldm = self.foldpanels.AddFoldPanel(caption="Demultiplexing Controls", collapsed=False,
                                                            cbstyle=styleht)
            paneldm = wx.Panel(self.foldpaneldm, -1)
            sizercontrolht1 = wx.GridBagSizer(wx.VERTICAL)
            i = 0

            # Button for Run TIC HT
            self.runticht = wx.Button(paneldm, -1, "Run TIC Demultiplex")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.on_run_tic_ht, self.runticht)
            sizercontrolht1.Add(self.runticht, (i, 0), span=(1, 2), flag=wx.EXPAND)
            self.runticht.SetToolTip(wx.ToolTip("Demultiplexing of TIC"))
            i += 1

            # Button for Run EIC HT
            self.runeicht = wx.Button(paneldm, -1, "Run EIC Demultiplex")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.on_run_eic_ht, self.runeicht)
            sizercontrolht1.Add(self.runeicht, (i, 0), span=(1, 2), flag=wx.EXPAND)
            self.runeicht.SetToolTip(wx.ToolTip("Demultiplexing of each EIC selected. "
                                                "To select, zoom on a region of the 2D plot and right click."))
            i += 1

            self.runallht = wx.Button(paneldm, -1, "Run All Demultiplex")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.on_run_all_ht, self.runallht)
            sizercontrolht1.Add(self.runallht, (i, 0), span=(1, 2), flag=wx.EXPAND)
            self.runallht.SetToolTip(wx.ToolTip("Run demultiplex on all data points. Used to create crazy cubes."))
            i += 1

            # Button to Mass Transformation
            self.masstransform = wx.Button(paneldm, -1, "Run All Mass Transform")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.run_all_mass_transform, self.masstransform)
            sizercontrolht1.Add(self.masstransform, (i, 0), span=(1, 2), flag=wx.EXPAND)
            self.masstransform.SetToolTip(wx.ToolTip("Run Mass Transformation on all data points. "
                                                     "Click button above first to do HT."))
            i += 1

            # Control for decmultiplex choice
            self.ctlmultiplexmode = wx.Choice(paneldm, -1, choices=self.config.demultiplexchoices)
            self.ctlmultiplexmode.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
            self.ctlmultiplexmode.SetSelection(0)
            sizercontrolht1.Add(self.ctlmultiplexmode, (i, 1), span=(1, 2))
            sizercontrolht1.Add(wx.StaticText(paneldm, label="Demultiplex Mode: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctlmultiplexmode.SetToolTip(wx.ToolTip("Choose the demultiplexing mode."))
            self.ctlmultiplexmode.Bind(wx.EVT_CHOICE, self.update_demultiplex_mode)
            i += 1

            # Text Control for Kernel Smoothing
            self.ctlkernelsmooth = wx.TextCtrl(paneldm, value="", size=size1)
            sizercontrolht1.Add(self.ctlkernelsmooth, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(paneldm, label="Smoothing: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text Control for Analysis Time
            self.ctlanalysistime = wx.TextCtrl(paneldm, value="", size=size1)
            sizercontrolht1.Add(self.ctlanalysistime, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(paneldm, label="Analysis Time: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text control for max scans
            self.ctlmaxscans = wx.TextCtrl(paneldm, value="", size=size1)
            sizercontrolht1.Add(self.ctlmaxscans, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(paneldm, label="Max Scans: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text control for timepad
            self.ctltimepad = wx.TextCtrl(paneldm, value="", size=size1)
            sizercontrolht1.Add(self.ctltimepad, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(paneldm, label="Time Padding: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            paneldm.SetSizer(sizercontrolht1)
            sizercontrolht1.Fit(paneldm)

            self.foldpanels.AddFoldPanelWindow(self.foldpaneldm, paneldm, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(self.foldpaneldm, wx.StaticText(self.foldpaneldm, -1, " "),
                                               fpb.FPB_ALIGN_WIDTH)

            self.foldpaneldmplot = self.foldpanels.AddFoldPanel(caption="Demultiplexing Plotting", collapsed=True,
                                                                cbstyle=styleht)
            paneldmplot = wx.Panel(self.foldpaneldmplot, -1)
            sizercontrolht1p = wx.GridBagSizer(wx.VERTICAL)
            i = 0

            # Drop down menu for X-axis
            self.ctlxaxis = wx.Choice(paneldmplot, -1, choices=["Time", "Scans"])
            self.ctlxaxis.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
            self.ctlxaxis.SetSelection(0)
            sizercontrolht1p.Add(self.ctlxaxis, (i, 1), span=(1, 1))
            sizercontrolht1p.Add(wx.StaticText(paneldmplot, label="X-Axis: "), (i, 0),
                                 flag=wx.ALIGN_CENTER_VERTICAL)
            # Bind to on replot
            self.ctlxaxis.Bind(wx.EVT_CHOICE, self.pres.on_replot_chrom)
            i += 1

            # Time Range for Demultiplexing
            self.ctltime1 = wx.TextCtrl(paneldmplot, value="", size=size1)
            self.ctltime2 = wx.TextCtrl(paneldmplot, value="", size=size1)
            trange = wx.BoxSizer(wx.HORIZONTAL)
            trange.Add(wx.StaticText(paneldmplot, label="DM Time Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
            trange.Add(self.ctltime1)
            trange.Add(wx.StaticText(paneldmplot, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
            trange.Add(self.ctltime2)
            sizercontrolht1p.Add(trange, (i, 0), span=(1, 2))
            i += 1

            # Check box to show legends
            self.ctllegends = wx.CheckBox(paneldmplot, label="")
            sizercontrolht1p.Add(self.ctllegends, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            sizercontrolht1p.Add(wx.StaticText(paneldmplot, label="Show Legends: "), (i, 0),
                                 flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctllegends.SetValue(True)
            # Bind to on replot
            self.ctllegends.Bind(wx.EVT_CHECKBOX, self.pres.on_replot_chrom)
            i += 1

            # Button to make 2d time vs. charge plot
            self.maketvsc = wx.Button(paneldmplot, -1, "Time-Z Plot")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_charge_time_2dplot, self.maketvsc)
            sizercontrolht1p.Add(self.maketvsc, (i, 0), span=(1, 1), flag=wx.EXPAND)

            # Button to make 2d time vs. mz plot
            self.maketvsm = wx.Button(paneldmplot, -1, "Time-m/z Plot")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_mz_time_2dplot, self.maketvsm)
            sizercontrolht1p.Add(self.maketvsm, (i, 1), span=(1, 1), flag=wx.EXPAND)
            i += 1

            # Button to make 2d time vs. mass plot
            self.makevtm = wx.Button(paneldmplot, -1, "Time-Mass Plot")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_mass_time_2dplot, self.makevtm)
            sizercontrolht1p.Add(self.makevtm, (i, 0), span=(1, 1), flag=wx.EXPAND)

            self.plot5button2 = wx.Button(paneldmplot, -1, "Mass-Charge Plot")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.makeplot5, self.plot5button2)
            sizercontrolht1p.Add(self.plot5button2, (i, 1), span=(1, 1), flag=wx.EXPAND)
            i += 1

            # Button for Make m/z Cube plots
            self.makemzcube = wx.Button(paneldmplot, -1, "m/z Cube")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_cube_plot, self.makemzcube)
            sizercontrolht1p.Add(self.makemzcube, (i, 0), span=(1, 1), flag=wx.EXPAND)

            # Button for make mass cube plot
            self.makemasscube = wx.Button(paneldmplot, -1, "Mass Cube")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_mass_cube_plot, self.makemasscube)
            sizercontrolht1p.Add(self.makemasscube, (i, 1), span=(1, 1), flag=wx.EXPAND)
            i += 1

            paneldmplot.SetSizer(sizercontrolht1p)
            sizercontrolht1p.Fit(paneldmplot)

            self.foldpanels.AddFoldPanelWindow(self.foldpaneldmplot, paneldmplot, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(self.foldpaneldmplot, wx.StaticText(self.foldpaneldmplot, -1, " "),
                                               fpb.FPB_ALIGN_WIDTH)

            self.foldpanelht = self.foldpanels.AddFoldPanel(caption="HT Parameters", collapsed=False,
                                                            cbstyle=styleht)
            panelht = wx.Panel(self.foldpanelht, -1)
            sizercontrolht1 = wx.GridBagSizer(wx.VERTICAL)
            i = 0

            # Drop down for HTseq
            self.ctlhtseq = wx.Choice(panelht, -1, choices=list(HTseqDict.keys()))
            self.ctlhtseq.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
            self.ctlhtseq.SetSelection(3)
            sizercontrolht1.Add(self.ctlhtseq, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelht, label="HT Sequence Bit: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text Control for HTtimeshift
            self.ctlhttimeshift = wx.TextCtrl(panelht, value="", size=size1)
            sizercontrolht1.Add(self.ctlhttimeshift, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelht, label="Time Shift: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text control for cycle index
            self.ctlcycleindex = wx.TextCtrl(panelht, value="", size=size1)
            sizercontrolht1.Add(self.ctlcycleindex, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelht, label="Cycle Length: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text control for HTmaskn
            self.ctlhtmaskn = wx.TextCtrl(panelht, value="", size=size1)
            sizercontrolht1.Add(self.ctlhtmaskn, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelht, label="mHT Iterations: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Text control for HTwin
            self.ctlhtwin = wx.TextCtrl(panelht, value="", size=size1)
            sizercontrolht1.Add(self.ctlhtwin, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelht, label="mHT Masks/Iter: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            panelht.SetSizer(sizercontrolht1)
            sizercontrolht1.Fit(panelht)

            self.foldpanels.AddFoldPanelWindow(self.foldpanelht, panelht, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(self.foldpanelht, wx.StaticText(self.foldpanelht, -1, " "),
                                               fpb.FPB_ALIGN_WIDTH)

            self.foldpanelft = self.foldpanels.AddFoldPanel(caption="FT Parameters", collapsed=True,
                                                            cbstyle=styleht)
            panelft = wx.Panel(self.foldpanelft, -1)
            sizercontrolht1 = wx.GridBagSizer(wx.VERTICAL)
            i = 0

            # Text Control for FTstart and FTend
            self.ctlftstart = wx.TextCtrl(panelft, value="", size=size1)
            sizercontrolht1.Add(self.ctlftstart, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelft, label="FT Start (Hz): "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctlftend = wx.TextCtrl(panelft, value="", size=size1)
            sizercontrolht1.Add(self.ctlftend, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelft, label="FT End (Hz): "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            # Create Two checkboxes for FTflatten and FTapodize
            self.ctlftflatten = wx.CheckBox(panelft, label="Flatten")
            sizercontrolht1.Add(self.ctlftflatten, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctlftflatten.SetValue(self.config.FTflatten)
            self.ctlftflatten.SetToolTip(wx.ToolTip("Flatten the FT data"))
            self.ctlftapodize = wx.CheckBox(panelft, label="Apodize")
            sizercontrolht1.Add(self.ctlftapodize, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctlftapodize.SetValue(self.config.FTapodize)
            self.ctlftapodize.SetToolTip(wx.ToolTip("Apodize the FT data"))
            i += 1

            # Add text input for post smoothing
            self.ctlftsmooth = wx.TextCtrl(panelft, value=str(self.config.FTsmooth), size=size1)
            sizercontrolht1.Add(self.ctlftsmooth, (i, 1), span=(1, 1))
            sizercontrolht1.Add(wx.StaticText(panelft, label="Post Smoothing: "), (i, 0),
                                flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            panelft.SetSizer(sizercontrolht1)
            sizercontrolht1.Fit(panelft)

            self.foldpanels.AddFoldPanelWindow(self.foldpanelft, panelft, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(self.foldpanelft, wx.StaticText(self.foldpanelft, -1, " "),
                                               fpb.FPB_ALIGN_WIDTH)

            self.foldpanelim = self.foldpanels.AddFoldPanel(caption="Ion Mobility Parameters", collapsed=False,
                                                            cbstyle=style1c)
            panel1c = wx.Panel(self.foldpanelim, -1)
            gbox1c = wx.GridBagSizer(wx.VERTICAL)
            i = 0

            # Create Run CCS Calc button
            self.runccsbutton = wx.Button(panel1c, -1, "Run All DT to CCS")
            self.Bind(wx.EVT_BUTTON, self.pres.on_run_ccs, self.runccsbutton)
            gbox1c.Add(self.runccsbutton, (i, 0), span=(1, 2), flag=wx.EXPAND)
            i += 1

            # Button to run EIC CCS
            self.runeicccs = wx.Button(panel1c, -1, "Run EIC CCS")
            self.Bind(wx.EVT_BUTTON, self.pres.on_run_eic_ccs, self.runeicccs)
            gbox1c.Add(self.runeicccs, (i, 0), span=(1, 2), flag=wx.EXPAND)
            i += 1

            self.ctlvolt = wx.TextCtrl(panel1c, value="", size=size1)
            gbox1c.Add(self.ctlvolt, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Voltage (V): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctlpressure = wx.TextCtrl(panel1c, value='', size=size1)
            gbox1c.Add(self.ctlpressure, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Pressure (Torr): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctltemp = wx.TextCtrl(panel1c, value='', size=size1)
            gbox1c.Add(self.ctltemp, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Temperature (\u00B0C): "), (i, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctlgasmass = wx.TextCtrl(panel1c, value='', size=size1)
            gbox1c.Add(self.ctlgasmass, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Gas Mass (Da): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctlto = wx.TextCtrl(panel1c, value='', size=size1)
            gbox1c.Add(self.ctlto, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Dead Time (t\u2080 in ms): "), (i, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctldriftlength = wx.TextCtrl(panel1c, value='', size=size1)
            gbox1c.Add(self.ctldriftlength, (i, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(panel1c, label="Drift Cell Length (m)"), (i, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctlccsbins = wx.TextCtrl(panel1c, value="", size=size1)
            gbox1c.Add(self.ctlccsbins, (i, 1), span=(1, 2))
            gbox1c.Add(wx.StaticText(panel1c, label="Sample CCS Every (\u212B\u00B2): "), (i, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

            panel1c.SetSizer(gbox1c)
            gbox1c.Fit(panel1c)

            self.foldpanels.AddFoldPanelWindow(self.foldpanelim, panel1c, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(self.foldpanelim, wx.StaticText(self.foldpanelim, -1, " "),
                                               fpb.FPB_ALIGN_WIDTH)

        # Panel for unidec Parameters
        foldpanel2 = self.foldpanels.AddFoldPanel(caption="UniDec Parameters", collapsed=False, cbstyle=style2)
        panel2 = wx.Panel(foldpanel2, -1)
        sizercontrol2 = wx.GridBagSizer(wx.VERTICAL)

        self.ctlmasslb = wx.TextCtrl(panel2, value="", size=(60, -1))
        self.ctlmassub = wx.TextCtrl(panel2, value="", size=(70, -1))
        massrange = wx.BoxSizer(wx.HORIZONTAL)
        massrange.Add(wx.StaticText(panel2, label="Mass Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
        massrange.Add(self.ctlmasslb)
        massrange.Add(wx.StaticText(panel2, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        massrange.Add(self.ctlmassub)
        massrange.Add(wx.StaticText(panel2, label=" Da  "), 0, wx.ALIGN_CENTER_VERTICAL)

        i = 0
        sizercontrol2.Add(massrange, (i, 0), span=(1, 2))
        i += 1

        self.ctlmzsig = wx.TextCtrl(panel2, value="", size=size1)
        sizercontrol2.Add(self.ctlmzsig, (i, 1), span=(1, 2))
        sizercontrol2.Add(wx.StaticText(panel2, label="m/z Spread FWHM (Th): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlcsig = wx.TextCtrl(panel2, value="", size=size1)
        sizercontrol2.Add(wx.StaticText(panel2, label="Charge Spread FWHM: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2.Add(self.ctlcsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.rununidec = wx.Button(panel2, -1, "Run UniDec")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_unidec_button, self.rununidec)
        sizercontrol2.Add(self.rununidec, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        panel2.SetSizer(sizercontrol2)
        sizercontrol2.Fit(panel2)
        self.foldpanels.AddFoldPanelWindow(foldpanel2, panel2, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2, wx.StaticText(foldpanel2, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Quick Access Panel
        foldpanel2a = self.foldpanels.AddFoldPanel(caption="Quick Controls", collapsed=False, cbstyle=style2)
        panel2a = wx.Panel(foldpanel2a, -1)
        sizercontrol2a = wx.GridBagSizer(wx.VERTICAL)
        i = 0

        self.ctlzsmoothcheck = wx.CheckBox(panel2a, label="Smooth Charge States Distributions", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_z_smooth, self.ctlzsmoothcheck)
        sizercontrol2a.Add(self.ctlzsmoothcheck, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlpselect = wx.RadioBox(panel2a, label="Smooth Nearby Points",
                                      choices=["None", "Some", "Lots", "Other"])
        self.parent.Bind(wx.EVT_RADIOBOX, self.on_p_select, self.ctlpselect)
        sizercontrol2a.Add(self.ctlpselect, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlbselect = wx.RadioBox(panel2a, label="Suppress Artifacts",
                                      choices=["None", "Some", "Lots", "Other"])
        self.parent.Bind(wx.EVT_RADIOBOX, self.on_b_select, self.ctlbselect)
        sizercontrol2a.Add(self.ctlbselect, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlmolig = wx.TextCtrl(panel2a, value="", size=size1)
        self.ctlmsmoothcheck = wx.CheckBox(panel2a, label="Mass Differences (Da): ", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_m_smooth, self.ctlmsmoothcheck)
        sizercontrol2a.Add(self.ctlmsmoothcheck, (i, 0), span=(1, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2a.Add(self.ctlmolig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        '''
        self.ctlpeakwidthcheck = wx.CheckBox(panel2a, label="Use Automatic m/z Peak Width", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_pw_check, self.ctlpeakwidthcheck)
        sizercontrol2a.Add(self.ctlpeakwidthcheck, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1        
        '''
        panel2a.SetSizer(sizercontrol2a)
        sizercontrol2a.Fit(panel2a)
        self.foldpanels.AddFoldPanelWindow(foldpanel2a, panel2a, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2a, wx.StaticText(foldpanel2a, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Additional Restraints
        foldpanel2b = self.foldpanels.AddFoldPanel(caption="Additional Deconvolution Parameters", collapsed=True,
                                                   cbstyle=style2b)
        panel2b = wx.Panel(foldpanel2b, -1)
        gbox2b = wx.GridBagSizer(wx.VERTICAL)

        i = 0

        self.ctlpsfun = wx.RadioBox(panel2b, label="m/z Peak Shape Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        gbox2b.Add(self.ctlpsfun, (i, 0), span=(1, 2))
        i += 1

        self.ctlpsfun2 = wx.RadioBox(panel2b, label="Charge Peak Shape Function",
                                     choices=["Gaussian", "Lorentzian", "Split L/G"])
        gbox2b.Add(self.ctlpsfun2, (i, 0), span=(1, 2))
        i += 1

        self.ctlnumit = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label='Maximum # of Iterations: '), (i, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlnumit, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlzzsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Charge Smooth Floor: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlzzsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlpsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Point Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlpsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlbeta = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Beta: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlbeta, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlmsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Mass Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlmsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlpoolflag = wx.RadioBox(panel2b, label="m/z to Mass Transformation",
                                       choices=["Integrate", "Interpolate"])
        gbox2b.Add(self.ctlpoolflag, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        sb = wx.StaticBox(panel2b, label='Native Charge Offset Range')
        sbs = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        self.ctlminnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        self.ctlmaxnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        sbs.Add(self.ctlminnativez, flag=wx.LEFT | wx.EXPAND, border=5)
        sbs.Add(wx.StaticText(panel2b, label=' to '), 0, wx.EXPAND)
        sbs.Add(self.ctlmaxnativez, flag=wx.LEFT | wx.EXPAND, border=5)
        gbox2b.Add(sbs, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        '''
        
        self.ctlmanualassign = wx.CheckBox(panel2b, label="Manual Mode")
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_check_manual, self.ctlmanualassign)
        gbox2b.Add(self.ctlmanualassign, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        mlsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmasslistflag = wx.CheckBox(panel2b, label="Mass List Window:")
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_mass_list, self.ctlmasslistflag)
        self.ctlmtabsig = wx.TextCtrl(panel2b, value="", size=(60, -1))
        mlsizer.Add(self.ctlmasslistflag, 0, wx.ALIGN_CENTER_VERTICAL)
        mlsizer.Add(self.ctlmtabsig, 0, wx.ALIGN_CENTER_VERTICAL)
        mlsizer.Add(wx.StaticText(panel2b, label=" Da "), 0, wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(mlsizer, (i, 0), span=(1, 2))
        i += 1

        
        '''
        panel2b.SetSizer(gbox2b)
        gbox2b.Fit(panel2b)
        self.foldpanels.AddFoldPanelWindow(foldpanel2b, panel2b, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2b, wx.StaticText(foldpanel2b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Peak Selection and Plotting
        foldpanel3 = self.foldpanels.AddFoldPanel(caption="Peak Selection and Plotting", collapsed=False,
                                                  cbstyle=style3)
        panel3 = wx.Panel(foldpanel3, -1)

        sizercontrol3 = wx.GridBagSizer(wx.VERTICAL)
        self.ctlwindow = wx.TextCtrl(panel3, value="", size=size1)
        self.ctlthresh = wx.TextCtrl(panel3, value="", size=size1)

        self.plotbutton = wx.Button(panel3, -1, "Peak Detection")
        self.plotbutton2 = wx.Button(panel3, -1, "Plot Peaks")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_plot_peaks, self.plotbutton2)
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_pick_peaks, self.plotbutton)

        sizercontrol3.Add(self.ctlwindow, (0, 1))
        sizercontrol3.Add(wx.StaticText(panel3, label="Peak Detection Range (Da): "), (0, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol3.Add(self.ctlthresh, (1, 1))
        sizercontrol3.Add(wx.StaticText(panel3, label="Peak Detection Threshold: "), (1, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)

        sizercontrol3.Add(self.plotbutton, (2, 0), span=(1, 1), flag=wx.EXPAND)
        sizercontrol3.Add(self.plotbutton2, (2, 1), span=(1, 1), flag=wx.EXPAND)

        panel3.SetSizer(sizercontrol3)
        sizercontrol3.Fit(panel3)
        self.foldpanels.AddFoldPanelWindow(foldpanel3, panel3, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel3, wx.StaticText(foldpanel3, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Plotting Parameters
        foldpanel3b = self.foldpanels.AddFoldPanel(caption="Additional Plotting Parameters", collapsed=True,
                                                   cbstyle=style3b)
        panel3b = wx.Panel(foldpanel3b, -1)

        gbox3b = wx.GridBagSizer(wx.VERTICAL)

        i = 0

        self.ctl2dcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)
        self.ctl2dcm.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)

        self.ctlpeakcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)
        self.ctlpeakcm.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)

        self.ctl2dcm.AppendItems(self.config.cmaps2)
        self.ctlpeakcm.AppendItems(self.config.cmaps)

        gbox3b.Add(self.ctl2dcm, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(wx.StaticText(panel3b, label='2D Color Map: '), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        gbox3b.Add(self.ctlpeakcm, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(wx.StaticText(panel3b, label='Peaks Color Map: '), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctldiscrete = wx.CheckBox(panel3b, label="Discrete Plot")
        gbox3b.Add(self.ctldiscrete, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlpublicationmode = wx.CheckBox(panel3b, label="Publication Mode")
        gbox3b.Add(self.ctlpublicationmode, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlrawflag = wx.RadioBox(panel3b, label="", choices=["Reconvolved/Profile", "Raw/Centroid"])
        gbox3b.Add(self.ctlrawflag, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.ctlnorm = wx.RadioBox(panel3b, label="Peak Normalization", choices=["None", "Max", "Total"])
        gbox3b.Add(self.ctlnorm, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.ctlthresh2 = wx.TextCtrl(panel3b, value="", size=size1)
        gbox3b.Add(wx.StaticText(panel3b, label="Marker Threshold: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(self.ctlthresh2, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlsep = wx.TextCtrl(panel3b, value="", size=size1)
        gbox3b.Add(wx.StaticText(panel3b, label="Species Separation: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(self.ctlsep, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        sb2 = wx.StaticBox(panel3b, label='Integration Range')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.HORIZONTAL)
        self.ctlintlb = wx.TextCtrl(panel3b, value='', size=(75, -1))
        self.ctlintub = wx.TextCtrl(panel3b, value='', size=(75, -1))
        sbs2.Add(self.ctlintlb, flag=wx.LEFT | wx.EXPAND, border=5)
        sbs2.Add(wx.StaticText(panel3b, label=' to '), 0, flag=wx.EXPAND)
        sbs2.Add(self.ctlintub, flag=wx.LEFT | wx.EXPAND, border=5)
        sbs2.Add(wx.StaticText(panel3b, label=' Da '), 0, flag=wx.EXPAND)
        gbox3b.Add(sbs2, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.replotbutton = wx.Button(panel3b, -1, "Replot")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_replot, self.replotbutton)
        gbox3b.Add(self.replotbutton, (i, 0), span=(1, 1), flag=wx.EXPAND)

        self.compositebutton = wx.Button(panel3b, -1, "Plot Composite")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_plot_composite, self.compositebutton)
        gbox3b.Add(self.compositebutton, (i, 1), span=(1, 1), flag=wx.EXPAND)
        i += 1

        self.plot5button = wx.Button(panel3b, -1, "Plot Mass v. Charge Grid")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.makeplot5, self.plot5button)
        gbox3b.Add(self.plot5button, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.plot5button2 = wx.Button(panel3b, -1, "Plot m/z v. Mass Grid (Fast)")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.make_mzmass_plot, self.plot5button2)
        gbox3b.Add(self.plot5button2, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.plot5button3 = wx.Button(panel3b, -1, "Plot m/z v. Mass Continuous (Slow)")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.make_mzmass_plot_continuous, self.plot5button3)
        gbox3b.Add(self.plot5button3, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        panel3b.SetSizer(gbox3b)
        gbox3b.Fit(panel3b)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, panel3b, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, wx.StaticText(foldpanel3b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        bright = 250
        foldpanel1.SetBackgroundColour(wx.Colour(bright, bright, 255))
        foldpanel1b.SetBackgroundColour(wx.Colour(bright, bright, 255))

        foldpanel2.SetBackgroundColour(wx.Colour(255, 255, bright))
        foldpanel2a.SetBackgroundColour(wx.Colour(255, 255, bright))
        foldpanel2b.SetBackgroundColour(wx.Colour(255, 255, bright))

        foldpanel3.SetBackgroundColour(wx.Colour(255, bright, bright))
        foldpanel3b.SetBackgroundColour(wx.Colour(255, bright, bright))

        sizercontrol.SetMinSize((250 + self.config.imflag * 10, 0))

        # Add to sizer and setup
        sizercontrolscrolled = wx.BoxSizer(wx.VERTICAL)
        sizercontrolscrolled.Add(self.foldpanels, 1, wx.EXPAND)
        self.scrolledpanel.SetSizer(sizercontrolscrolled)
        sizercontrolscrolled.Fit(self.scrolledpanel)

        # Add to top control sizer
        sizercontrol.Add(self.scrolledpanel, 1, wx.EXPAND)

        # Bottom Buttons
        if self.config.imflag == 0:
            buttonsizer2 = wx.BoxSizer(wx.HORIZONTAL)
            bsize = (45, 25)
            self.mainbutton = wx.Button(self, -1, "Main", size=bsize)
            self.expandbutton = wx.Button(self, -1, "All", size=bsize)
            self.collapsebutton = wx.Button(self, -1, "None", size=bsize)
            self.bluebutton = wx.Button(self, -1, "Blue", size=bsize)
            self.yellowbutton = wx.Button(self, -1, "Yellow", size=bsize)
            self.redbutton = wx.Button(self, -1, "Red", size=bsize)
            self.bluebutton.SetBackgroundColour(wx.Colour(150, 150, 255))
            self.yellowbutton.SetBackgroundColour(wx.Colour(255, 255, 150))
            self.redbutton.SetBackgroundColour(wx.Colour(255, 150, 150))
            self.parent.Bind(wx.EVT_BUTTON, self.on_expand_main, self.mainbutton)
            self.parent.Bind(wx.EVT_BUTTON, self.on_expand_all, self.expandbutton)
            self.parent.Bind(wx.EVT_BUTTON, self.on_collapse_all, self.collapsebutton)
            self.parent.Bind(wx.EVT_BUTTON, self.on_expand_blue, self.bluebutton)
            self.parent.Bind(wx.EVT_BUTTON, self.on_expand_yellow, self.yellowbutton)
            self.parent.Bind(wx.EVT_BUTTON, self.on_expand_red, self.redbutton)
            buttons = [self.mainbutton, self.expandbutton, self.collapsebutton, self.bluebutton, self.yellowbutton,
                       self.redbutton]
            for button in buttons:
                buttonsizer2.Add(button, 0, wx.EXPAND)
            sizercontrol.Add(buttonsizer2, 0, wx.EXPAND)

        # Set top sizer
        self.SetSizer(sizercontrol)
        sizercontrol.Fit(self)

        # Add remaining things
        self.setup_tool_tips()
        self.bind_changes()

    def import_config_to_gui(self):
        """
        Imports parameters from the config object to the GUI.
        :return: None
        """
        self.Freeze()
        self.update_flag = False
        if self.config.batchflag == 0:
            self.ctlmassbins.SetValue(str(self.config.massbins))

            self.ctlstartz.SetValue(str(self.config.startz))
            self.ctlendz.SetValue(str(self.config.endz))
            self.ctlslope.SetValue(str(self.config.CDslope))
            self.ctlmzsig.SetValue(str(self.config.mzsig))
            self.ctlnorm.SetSelection(int(self.config.peaknorm))
            self.ctlmasslb.SetValue(str(self.config.masslb))
            self.ctlmassub.SetValue(str(self.config.massub))
            self.ctlzzsig.SetValue(str(self.config.zzsig))
            self.ctlbuff.SetValue(str(self.config.subbuff))
            self.ctlsubbuffdt.SetValue(str(self.config.subbufdt))
            self.ctlminnativez.SetValue(str(self.config.nativezlb))
            self.ctlmaxnativez.SetValue(str(self.config.nativezub))
            self.ctlbeta.SetValue(str(self.config.beta))
            self.ctlpsig.SetValue(str(self.config.psig))
            try:
                self.ctlpoolflag.SetSelection(self.config.poolflag)
            except:
                self.ctlpoolflag.SetSelection(1)

            self.subtypectl.SetSelection(int(self.config.subtype))
            self.ctlpsfun.SetSelection(int(self.config.psfun))
            self.ctlpsfun2.SetSelection(int(self.config.psfunz))
            self.ctlmsig.SetValue(str(self.config.msig))
            self.ctlmolig.SetValue(str(self.config.molig))

            '''
            self.ctlmasslistflag.SetValue(self.config.mfileflag)
            self.ctlmtabsig.SetValue(str(self.config.mtabsig))
            self.ctlmanualassign.SetValue(self.config.manualfileflag)            
            '''
            self.ctlscanstart.SetValue(str(self.config.CDScanStart))
            self.ctlscanend.SetValue(str(self.config.CDScanEnd))

            self.ctlzbins.SetValue(str(self.config.CDzbins))
            self.ctlcentroid.SetValue(str(self.config.CDres))
            self.ctlbinsize.SetValue(str(self.config.mzbins))
            self.ctlsmooth.SetValue(str(self.config.smooth))
            self.ctlsmoothdt.SetValue(str(self.config.smoothdt))
            self.ctlwindow.SetValue(str(self.config.peakwindow))
            self.ctlthresh.SetValue(str(self.config.peakthresh))
            self.ctlthresh2.SetValue(str(self.config.peakplotthresh))
            self.ctlsep.SetValue(str(self.config.separation))
            self.ctlintthresh.SetValue(str(self.config.intthresh))
            self.ctlprethresh.SetValue(str(self.config.CDprethresh))
            self.ctladductmass.SetValue(str(self.config.adductmass))
            self.ctldatanorm.SetValue(int(self.config.datanorm))
            self.ctlnumit.SetValue(str(self.config.numit))
            self.ctldatareductionpercent.SetValue(str(self.config.reductionpercent))

            self.ctldiscrete.SetValue(self.config.discreteplot)
            self.ctlsmashflag.SetValue(self.config.smashflag)
            self.ctlinjtime.SetValue(self.config.CDiitflag)
            self.ctlpublicationmode.SetValue(self.config.publicationmode)
            self.ctlrawflag.SetSelection(self.config.rawflag)

            if self.htmode:
                self.ctlscancompress.SetValue(str(self.config.CDScanCompress))
                self.ctlkernelsmooth.SetValue(str(self.config.HTksmooth))
                self.ctlhtseq.SetStringSelection(str(self.config.htbit))
                self.ctlhttimeshift.SetValue(str(self.config.HTtimeshift))
                self.ctltimepad.SetValue(str(self.config.HTtimepad))
                self.ctlanalysistime.SetValue(str(self.config.HTanalysistime))
                self.ctlmaxscans.SetValue(str(self.config.HTmaxscans))
                self.ctlcycleindex.SetValue(str(self.config.HTcycleindex))
                self.ctlhtmaskn.SetValue(str(self.config.HTmaskn))
                self.ctlhtwin.SetValue(str(self.config.HTwin))
                self.ctlxaxis.SetStringSelection(self.config.HTxaxis)
                self.ctlmultiplexmode.SetStringSelection(self.config.demultiplexmode)
                self.ctlftstart.SetValue(str(self.config.FTstart))
                self.ctlftend.SetValue(str(self.config.FTend))
                self.ctlftflatten.SetValue(self.config.FTflatten)
                self.ctlftapodize.SetValue(self.config.FTapodize)
                self.ctlftsmooth.SetValue(str(self.config.FTsmooth))

                self.ctlvolt.SetValue(str(self.config.volt))
                self.ctltemp.SetValue(str(self.config.temp))
                self.ctlpressure.SetValue(str(self.config.pressure))
                self.ctlgasmass.SetValue(str(self.config.gasmass))
                self.ctlto.SetValue(str(self.config.to))
                self.ctldriftlength.SetValue(str(self.config.driftlength))
                self.ctlccsbins.SetValue(str(self.config.ccsbins))

                self.ctltime1.SetValue(str(self.config.HToutputlb))
                self.ctltime2.SetValue(str(self.config.HToutputub))

            try:
                if float(self.config.adductmass) < 0:
                    self.ctlnegmode.SetValue(1)
                else:
                    self.ctlnegmode.SetValue(0)
            except Exception as e:
                self.ctlnegmode.SetValue(0)
                self.config.adductmass = 1.007276467

            try:
                self.ctl2dcm.SetSelection(self.config.cmaps2.index(self.config.cmap))
                self.ctlpeakcm.SetSelection(self.config.cmaps.index(self.config.peakcmap))
            except ValueError:
                print("Could not find the specified color map. Try upgrading to the latest version of matplotlib.")
                import matplotlib
                print("Current version:", matplotlib.__version__)
                # Revert to the defaults
                self.ctl2dcm.SetSelection(self.config.cmaps.index("nipy_spectral"))
                self.ctlpeakcm.SetSelection(self.config.cmaps.index("rainbow"))

            self.ctlcsig.SetValue(str(self.config.csig))

            try:
                x = float(self.config.integratelb)
                y = float(self.config.integrateub)
                self.ctlintlb.SetValue(str(x))
                self.ctlintub.SetValue(str(y))
            except (ValueError, TypeError):
                self.ctlintlb.SetValue("")
                self.ctlintub.SetValue("")

            if self.config.msig > 0:
                self.parent.SetStatusText(
                    "Oligomer Blur Mass: " + str(self.config.molig) + " Std Dev: " + str(self.config.msig),
                    number=4)
            else:
                self.parent.SetStatusText(" ", number=4)

        # If the batchflag is not 1, it will import the data range as well
        if self.config.batchflag != 1:
            self.ctlminmz.SetValue(str(self.config.minmz))
            self.ctlmaxmz.SetValue(str(self.config.maxmz))

        self.update_flag = True
        try:
            self.update_quick_controls()
        except Exception as e:
            print("Error updating quick controls", e)

        self.Thaw()

        if self.htmode:
            self.update_demultiplex_mode()

    def export_gui_to_config(self, e=None):
        """
        Exports parameters from the GUI to the config object.
        :return: None
        """
        self.config.minmz = ud.string_to_value(self.ctlminmz.GetValue())
        self.config.maxmz = ud.string_to_value(self.ctlmaxmz.GetValue())
        self.config.CDzbins = ud.string_to_value(self.ctlzbins.GetValue())
        self.config.CDScanStart = ud.string_to_value(self.ctlscanstart.GetValue())
        self.config.CDScanEnd = ud.string_to_value(self.ctlscanend.GetValue())
        self.config.CDres = ud.string_to_value(self.ctlcentroid.GetValue())
        self.config.CDslope = ud.string_to_value(self.ctlslope.GetValue())
        self.config.subtype = self.subtypectl.GetSelection()
        self.config.mzbins = ud.string_to_value(self.ctlbinsize.GetValue())
        self.config.zzsig = ud.string_to_value(self.ctlzzsig.GetValue())
        self.config.smooth = ud.string_to_value(self.ctlsmooth.GetValue())
        self.config.smoothdt = ud.string_to_value(self.ctlsmoothdt.GetValue())
        self.config.subbuff = ud.string_to_value(self.ctlbuff.GetValue())
        self.config.subbufdt = ud.string_to_value(self.ctlsubbuffdt.GetValue())
        self.config.nativezlb = ud.string_to_value(self.ctlminnativez.GetValue())
        self.config.nativezub = ud.string_to_value(self.ctlmaxnativez.GetValue())
        self.config.psfun = self.ctlpsfun.GetSelection()
        self.config.psfunz = self.ctlpsfun2.GetSelection()
        self.config.msig = ud.string_to_value(self.ctlmsig.GetValue())
        self.config.molig = ud.string_to_value(self.ctlmolig.GetValue())
        '''          
        self.config.mtabsig = ud.string_to_value(self.ctlmtabsig.GetValue())
        self.config.mfileflag = int(self.ctlmasslistflag.GetValue())      
        self.config.manualfileflag = int(self.ctlmanualassign.GetValue())
        '''
        self.config.psig = ud.string_to_value(self.ctlpsig.GetValue())
        self.config.beta = ud.string_to_value(self.ctlbeta.GetValue())
        self.config.poolflag = self.ctlpoolflag.GetSelection()
        self.config.datanorm = int(self.ctldatanorm.GetValue())
        self.config.intthresh = ud.string_to_value(self.ctlintthresh.GetValue())
        self.config.CDprethresh = ud.string_to_value(self.ctlprethresh.GetValue())
        self.config.massbins = ud.string_to_value(self.ctlmassbins.GetValue())

        self.config.endz = ud.string_to_int(self.ctlendz.GetValue())
        self.config.startz = ud.string_to_int(self.ctlstartz.GetValue())

        self.config.mzsig = ud.string_to_value(self.ctlmzsig.GetValue())
        self.config.massub = ud.string_to_value(self.ctlmassub.GetValue())
        self.config.masslb = ud.string_to_value(self.ctlmasslb.GetValue())

        self.config.peaknorm = self.ctlnorm.GetSelection()

        self.config.peakwindow = ud.string_to_value(self.ctlwindow.GetValue())
        self.config.peakthresh = ud.string_to_value(self.ctlthresh.GetValue())
        self.config.peakplotthresh = ud.string_to_value(self.ctlthresh2.GetValue())
        self.config.separation = ud.string_to_value(self.ctlsep.GetValue())
        self.config.adductmass = ud.string_to_value(self.ctladductmass.GetValue())

        if self.htmode:
            self.config.CDScanCompress = ud.string_to_value(self.ctlscancompress.GetValue())
            self.config.HTksmooth = ud.string_to_value(self.ctlkernelsmooth.GetValue())
            #            self.config.htbit = int(self.ctlhtseq.GetStringSelection())
            self.config.htbit = self.ctlhtseq.GetStringSelection()
            self.config.htseq = HTseqDict[str(self.config.htbit)]
            self.config.HTxaxis = self.ctlxaxis.GetStringSelection()
            self.config.HTtimeshift = ud.string_to_value(self.ctlhttimeshift.GetValue())
            self.config.HTtimepad = ud.string_to_value(self.ctltimepad.GetValue())
            self.config.HTanalysistime = ud.string_to_value(self.ctlanalysistime.GetValue())
            self.config.HTmaxscans = ud.string_to_value(self.ctlmaxscans.GetValue())
            self.config.HTcycleindex = ud.string_to_value(self.ctlcycleindex.GetValue())
            self.config.HTmaskn = ud.string_to_value(self.ctlhtmaskn.GetValue())
            self.config.HTwin = ud.string_to_value(self.ctlhtwin.GetValue())
            self.config.demultiplexmode = self.ctlmultiplexmode.GetStringSelection()
            self.config.FTstart = ud.string_to_value(self.ctlftstart.GetValue())
            self.config.FTend = ud.string_to_value(self.ctlftend.GetValue())
            self.config.FTflatten = self.ctlftflatten.GetValue()
            self.config.FTapodize = self.ctlftapodize.GetValue()
            self.config.FTsmooth = ud.string_to_value(self.ctlftsmooth.GetValue())
            self.config.showlegends = self.ctllegends.GetValue()

            self.config.volt = ud.string_to_value(self.ctlvolt.GetValue())
            self.config.temp = ud.string_to_value(self.ctltemp.GetValue())
            self.config.pressure = ud.string_to_value(self.ctlpressure.GetValue())
            self.config.gasmass = ud.string_to_value(self.ctlgasmass.GetValue())
            self.config.to = ud.string_to_value(self.ctlto.GetValue())
            self.config.driftlength = ud.string_to_value(self.ctldriftlength.GetValue())
            self.config.ccsbins = ud.string_to_value(self.ctlccsbins.GetValue())

            self.config.HToutputlb = ud.string_to_value(self.ctltime1.GetValue())
            self.config.HToutputub = ud.string_to_value(self.ctltime2.GetValue())

        self.config.numit = ud.string_to_int(self.ctlnumit.GetValue())

        self.config.integratelb = ud.string_to_value(self.ctlintlb.GetValue())
        self.config.integrateub = ud.string_to_value(self.ctlintub.GetValue())
        self.config.reductionpercent = ud.string_to_value(self.ctldatareductionpercent.GetValue())

        self.config.smashflag = int(self.ctlsmashflag.GetValue())
        self.config.CDiitflag = int(self.ctlinjtime.GetValue())
        self.config.discreteplot = int(self.ctldiscrete.GetValue())
        self.config.publicationmode = int(self.ctlpublicationmode.GetValue())
        self.config.rawflag = self.ctlrawflag.GetSelection()

        try:
            self.config.adductmass = float(self.config.adductmass)
        except Exception:
            self.config.adductmass = 1.007276467

        if self.ctlnegmode.GetValue() == 1:
            # print("Negative Ion Mode")
            if self.config.adductmass > 0:
                self.config.adductmass *= -1
                self.ctladductmass.SetValue(str(self.config.adductmass))
        else:
            # print("Positive Ion Mode")
            if self.config.adductmass < 0:
                self.config.adductmass *= -1
                self.ctladductmass.SetValue(str(self.config.adductmass))

        self.config.cmap = str(self.ctl2dcm.GetStringSelection())
        self.config.peakcmap = str(self.ctlpeakcm.GetStringSelection())

        self.config.csig = ud.string_to_value(self.ctlcsig.GetValue())

        if not self.config.minmz and not ud.isempty(self.pres.eng.data.rawdata):
            self.config.minmz = np.amin(self.pres.eng.data.rawdata[:, 0])
        if not self.config.maxmz and not ud.isempty(self.pres.eng.data.rawdata):
            self.config.maxmz = np.amax(self.pres.eng.data.rawdata[:, 0])
        pass

        if self.config.msig > 0:
            self.parent.SetStatusText(
                "Oligomer Blur Mass: " + str(self.config.molig) + " Std Dev: " + str(self.config.msig),
                number=4)
        else:
            self.parent.SetStatusText(" ", number=4)
        # noinspection PyPep8

    def setup_tool_tips(self):
        """
        Sets Tool Tips for items on the Main Panel
        :return: None
        """
        self.ctlthresh2.SetToolTip(wx.ToolTip(
            "Set threshold for peaks to be plotted in m/z. "
            "Peak at given charge state must be greater than threshold * maximum m/z intensity."))
        self.ctlsep.SetToolTip(wx.ToolTip("Set distance between isolated peak m/z plots."))
        self.ctlwindow.SetToolTip(
            wx.ToolTip("Peak detection window. Peak must be maximum in a +/- window range in mass (Da)."))
        self.ctlthresh.SetToolTip(wx.ToolTip(
            "Peak detection threshold. Peak's intensity must be great than threshold times maximum mass intensity."))
        self.plotbutton.SetToolTip(wx.ToolTip("Pick peaks and plot. (Ctrl+P)"))
        self.plotbutton2.SetToolTip(wx.ToolTip("Plot individual peak species in m/z. (Ctrl+K)"))

        self.rununidec.SetToolTip(wx.ToolTip("Write Configuration File, Run UniDec, and Plot Results. (Ctrl+R)"))
        self.ctlmzsig.SetToolTip(wx.ToolTip(
            "Expected peak FWHM in m/z (Th).\nFor nonlinear mode, minimum FWHM"
            "\nSee Tools>Peak Width Tool for more tools."))

        self.ctlmassub.SetToolTip(wx.ToolTip(
            "Maximum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlmassbins.SetToolTip(wx.ToolTip("Sets the resolution of the zero-charge mass spectrum"))
        self.ctlmasslb.SetToolTip(wx.ToolTip(
            "Minimum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlstartz.SetToolTip(wx.ToolTip("Minimum allowed charge state in deconvolution."))
        self.ctlendz.SetToolTip(wx.ToolTip("Maximum allowed charge state in deconvolution."))
        self.ctlzbins.SetToolTip(wx.ToolTip("Charge Bin Size for Histogram."))
        self.ctlscanstart.SetToolTip(wx.ToolTip("Minimum scan number to include in deconvolution."))
        self.ctlscanend.SetToolTip(wx.ToolTip("Maximum scan number to include in deconvolution."))
        self.ctlcentroid.SetToolTip(wx.ToolTip("Width for Centroid Filtering in units of m/z."))
        self.ctlbinsize.SetToolTip(wx.ToolTip(
            "Controls Linearization.\nConstant bin size (Th) for Linear m/z"
            "\nMinimum bin size (Th) for Linear Resolution"
            "\nNumber of data points compressed together for Nonlinear"))
        self.ctlsmooth.SetToolTip(wx.ToolTip("Gaussian smooth sigma in units of data point number."))
        self.ctlsmoothdt.SetToolTip(wx.ToolTip("Gaussian smooth sigma in units of charge number."))
        self.ctlintthresh.SetToolTip(
            wx.ToolTip("Set intensity threshold. Data points below threshold are excluded from deconvolution."))
        self.ctlprethresh.SetToolTip(
            wx.ToolTip("Set intensity threshold applied during data import."
                       " Data points below threshold are excluded from import. "
                       "If you adjust this, you will need to reimport the data."))
        self.ctldatareductionpercent.SetToolTip(
            wx.ToolTip(
                "Reduces the amount of data by removing everything below a threshold."
                "\nSets the threshold to fit the percentage of data to remove."))
        self.ctlminmz.SetToolTip(wx.ToolTip("Set minimum m/z of data"))
        self.ctlmaxmz.SetToolTip(wx.ToolTip("Set maximum m/z of data"))
        self.dataprepbutton.SetToolTip(
            wx.ToolTip("Subtract, linearize, smooth, threshold, and write data to file. (Ctrl+D)"))
        self.ctladductmass.SetToolTip(wx.ToolTip("Mass of charge carrying adduct;\ntypically the mass of a proton"))
        self.ctlminnativez.SetToolTip(wx.ToolTip("Minimum offset from a native charge state"))
        self.ctlmaxnativez.SetToolTip(wx.ToolTip("Maximum offset from a native charge state"))
        self.ctlbselect.SetToolTip(wx.ToolTip(
            "Select whether to suppress deconvolution artifacts"))
        # TODO: Set several of these that are gone. Edit the others.
        '''
        self.ctlmasslistflag.SetToolTip(wx.ToolTip(
            "Limit deconvolution to specific masses +/- some window.\nDefine in Tools>Oligomer and Mass Tools."))
        self.ctlmtabsig.SetToolTip(
            wx.ToolTip("Set window for mass limitations. Setting to 0 will force only listed masses."))
        self.ctlpsfun.SetToolTip(wx.ToolTip("Expected peak shape.\nSee Tools>Peak Width Tool for more tools."))
        self.ctlzzsig.SetToolTip(wx.ToolTip(
            "Parameter for defining the width of the charge state smooth."
            "\nUniDec will use a mean filter of width 2n+1 on log_e of the charge distribution")) 
        self.ctlmsig.SetToolTip(wx.ToolTip("Width of Mass Smooth Filter"))
        self.ctlmolig.SetToolTip(wx.ToolTip("Mass difference used for Mass Smooth Filter"))
        
        self.ctlisotopemode.SetToolTip(wx.ToolTip(
            "Use isotopic distributions in deconvolution.\nOutput either monoisotopic or average masses"))
        self.ctlmanualassign.SetToolTip(wx.ToolTip("Use manual assignments. See Tools>Manual Assignment"))

        self.ctlpsig.SetToolTip(wx.ToolTip(
            "Parameter for defining the width of the data point smooth."
            "\nUniDec will weight +/- n data points to have the same charge state."))
        self.ctlbeta.SetToolTip(wx.ToolTip(
            "Parameter for defining the degree of Softmax distribution applied to the charge state vectors."
            "\n0 will shut it off."))
        self.ctlpselect.SetToolTip(wx.ToolTip(
            "Select whether to smooth nearby data points to have similar charge assignments"))
        
        self.ctlpoolflag.SetToolTip(wx.ToolTip(
            "Sets type of conversion from m/z to mass.\nIntegration:\n\tEach m/z bin goes to the nearest mass bin"
            "\n\tBest for undersampled masses\nInterpolation:\n\tEach mass value interpolates its value in m/z space"
            "\n\tBest for oversampled mass data"
            "\nSmart:\n\tDetermines interpolation vs. integration automatically\n\tfrom the data and blends the two."
            "\n\tUse Smart in most cases (default). \n\tNot available for IM-MS data."))
        self.ctlmsmoothcheck.SetToolTip(
            wx.ToolTip("Select whether to assume a smooth distribution spaced by a known mass difference"))
        self.ctlpeakwidthcheck.SetToolTip(
            wx.ToolTip("Select whether to incorporated the automatic peak width in deconvolution"))    
        '''
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
            "Controls range for integration.\nDefault is +/- Peak Detection Window."
            "\nUses these boxes to manually set - and +"))
        self.ctlintub.SetToolTip(wx.ToolTip(
            "Controls range for integration.\nDefault is +/- Peak Detection Window."
            "\nUses these boxes to manually set - and +"))
        self.ctlnorm.SetToolTip(wx.ToolTip(
            "Sets normalization of mass data.\nMaximum will normalize so that the maximum value is %100."
            "\nTotal will normalize so that the sum of all peaks is %100"))
        self.replotbutton.SetToolTip(wx.ToolTip("Replot some of the plots. (Ctrl+N)"))
        self.compositebutton.SetToolTip(
            wx.ToolTip("Plot composite of simulated spectra from selected peaks. (Ctrl+C)"))
        self.openbutton.SetToolTip(wx.ToolTip("Open .txt or .mzML file (Ctrl+O)"))
        self.procbutton.SetToolTip(wx.ToolTip("Process Data (Ctrl+D)"))
        self.udbutton.SetToolTip(wx.ToolTip("Run UniDec (Ctrl+R)"))
        self.ppbutton.SetToolTip(wx.ToolTip("Pick Peaks (Ctrl+P)"))
        self.autobutton.SetToolTip(wx.ToolTip("Process Data, Run UniDec, Pick Peaks (Ctrl+E)"))
        self.ctlslope.SetToolTip(wx.ToolTip("Slope for conversion of intensity to charge.\nz=I/slope"))
        self.ctlcsig.SetToolTip(wx.ToolTip("Width of Charge Distribution.\nFWHM of the charge distribution."))

        self.ctlzsmoothcheck.SetToolTip(
            wx.ToolTip("Select whether to assume a smooth charge state distribution"))

        if self.htmode:
            self.ctlscancompress.SetToolTip(wx.ToolTip(
                "Average each n scans together. This reduces the number of scans by a factor of n."))

            self.ctlkernelsmooth.SetToolTip(wx.ToolTip(
                "Smooth the data with a Gaussian kernel of this many data points."))
            self.ctlhtseq.SetToolTip(wx.ToolTip(
                "The Hadamard bit depth. Put bit{n} in file name to set automatically."))
            self.ctlhttimeshift.SetToolTip(wx.ToolTip(
                "Shift the HT start time by this amount to account for a delay at the start of the run. "
                "Units of retention time."))

            self.ctltimepad.SetToolTip(wx.ToolTip(
                "Pad the time axis by this amount to account for a delay at the end of the run. "
                "Units of retention time. Set cyc{n} and zp{m} in file name to set automatically as n and m."))
            self.ctlanalysistime.SetToolTip(wx.ToolTip(
                "Total analysis time. Must be set manually for DMT files."))
            self.ctlmaxscans.SetToolTip(wx.ToolTip(
                "Maximum number of scans to use in the analysis. Set to -1 to automatically determine from data"))
            self.ctlxaxis.SetToolTip(wx.ToolTip("Select the x-axis unit for plotting."))
            self.ctlcycleindex.SetToolTip(wx.ToolTip(
                "The number of scans per cycle. Determined automatically if -1. "
                "See tools menu or autocorrelation right click on TIC to get hints on setting this differently."))

        pass

    # .......................................................
    #
    #  The Main Panel
    #
    # .......................................................

    '''
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
                self.ctlmasslistflag.SetValue(False)'''

    def update_prethresh(self, e):
        self.config.CDprethresh = ud.string_to_value(self.ctlprethresh.GetValue())
        try:
            self.config.CDprethresh = float(self.config.CDprethresh)
        except Exception:
            self.config.CDprethresh = -1

    def on_z_smooth(self, e):
        value = self.ctlzsmoothcheck.Get3StateValue()
        if value == 1:
            self.ctlzzsig.SetValue("1e-6")
            self.ctlzbins.SetValue("1")
        elif value == 0:
            self.ctlzzsig.SetValue("0")
        self.export_gui_to_config()

    def on_m_smooth(self, e):
        value = self.ctlmsmoothcheck.Get3StateValue()
        if value == 1:
            self.ctlmsig.SetValue("1e-6")
        elif value == 0:
            self.ctlmsig.SetValue("0")
        self.export_gui_to_config()

    '''
    def on_pw_check(self, e):
        value = self.ctlpeakwidthcheck.Get3StateValue()
        if value == 1:
            self.ctlmzsig.SetValue(str(self.config.automzsig))
            self.ctlpsfun.SetSelection(self.config.autopsfun)
        elif value == 0:
            self.ctlmzsig.SetValue("0")
        self.export_gui_to_config()'''

    def on_check_smash(self, e):
        """
        Checks the configuration to see if values for noise smashing are set. If they are not,
        it opens the window to set the smash assignments.
        :param e: Dummy wx event passed on.
        :return: None
        """
        self.config.smashflag = self.ctlsmashflag.GetValue()
        if len(self.config.smashlist) < 1:
            self.pres.on_smash_window(e)
            if len(self.config.smashlist) < 1:
                self.ctlsmashflag.SetValue(False)

    def on_p_select(self, e):
        value = self.ctlpselect.GetSelection()
        self.ctlpsig.SetValue(str(self.psigsettings[value]))
        self.export_gui_to_config()

    def on_b_select(self, e):
        value = self.ctlbselect.GetSelection()
        self.ctlbeta.SetValue(str(self.betasettings[value]))
        self.export_gui_to_config()

    def on_collapse_all(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            self.foldpanels.Collapse(fp)
        pass

    def on_expand_all(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            self.foldpanels.Expand(fp)
        pass

    def expand_list(self, list):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in list:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)

    def on_expand_blue(self, e=None):
        self.expand_list([0, 1])

    def on_expand_yellow(self, e=None):
        if self.htmode:
            yellow = [9, 7, 8]
        else:
            yellow = [2, 3, 4]
        self.expand_list(yellow)

    def on_expand_red(self, e=None):
        if self.htmode:
            red = [10, 11]
        else:
            red = [5, 6]
        self.expand_list(red)

    def on_expand_main(self, e=None):
        if self.htmode:
            self.expand_list([0, 2, 8, 7, 10])
        else:
            self.expand_list([0, 2, 3, 5])

    def bind_changes(self, e=None):
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlzzsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlbeta)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlpsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlmsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlmzsig)
        '''
        
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlbuff)
        self.parent.Bind(wx.EVT_RADIOBOX, self.update_quick_controls, self.ctlpsfun)
        
        '''

    def on_mousewheel(self, e):
        # print("Wee")
        pass

    def update_demultiplex_mode(self, e=None):
        demultiplexmode = self.ctlmultiplexmode.GetStringSelection()
        if "HT" in demultiplexmode:
            self.foldpanels.Collapse(self.foldpanelft)
            self.foldpanels.Expand(self.foldpanelht)
            self.foldpanels.Collapse(self.foldpanelim)
        elif "FT" in demultiplexmode:
            self.foldpanels.Collapse(self.foldpanelht)
            self.foldpanels.Expand(self.foldpanelft)
        else:
            print("Unknown demultiplex mode:", demultiplexmode)

    def update_quick_controls(self, e=None):
        if self.update_flag:

            # Z Box
            try:
                value = float(self.ctlzzsig.GetValue())
            except:
                value = -1
            if value == 0:
                self.ctlzsmoothcheck.Set3StateValue(0)
            elif value == 1e-6:
                self.ctlzsmoothcheck.Set3StateValue(1)
            else:
                self.ctlzsmoothcheck.Set3StateValue(2)

            # beta box
            try:
                value = float(self.ctlbeta.GetValue())
            except:
                value = -1
            if value == self.betasettings[0]:
                self.ctlbselect.SetSelection(0)
            elif value == self.betasettings[1]:
                self.ctlbselect.SetSelection(1)
            elif value == self.betasettings[2]:
                self.ctlbselect.SetSelection(2)
            else:
                self.ctlbselect.SetSelection(3)
                self.betasettings[3] = value
                pass

            # point box
            try:
                value = float(self.ctlpsig.GetValue())
            except:
                value = -1
            if value == self.psigsettings[0]:
                self.ctlpselect.SetSelection(0)
            elif value == self.psigsettings[1]:
                self.ctlpselect.SetSelection(1)
            elif value == self.psigsettings[2]:
                self.ctlpselect.SetSelection(2)
            else:
                self.ctlpselect.SetSelection(3)
                self.psigsettings[3] = value
                pass

            # M box
            try:
                value = float(self.ctlmsig.GetValue())
            except:
                value = -1
            if value == 0:
                self.ctlmsmoothcheck.Set3StateValue(0)
            elif value == 1e-6:
                self.ctlmsmoothcheck.Set3StateValue(1)
            else:
                self.ctlmsmoothcheck.Set3StateValue(2)

            '''
            

            # PW box
            try:
                value = float(self.ctlmzsig.GetValue())
            except:
                value = -1
            psval = self.ctlpsfun.GetSelection()
            # try:
            #    fwhm, psfun, mid = ud.auto_peak_width(self.pres.eng.data.data2)
            # except:
            #    fwhm = -1
            #    psfun = -1
            if value == 0:
                self.ctlpeakwidthcheck.Set3StateValue(0)
            elif value == self.config.automzsig and psval == self.config.autopsfun:
                self.ctlpeakwidthcheck.Set3StateValue(1)
            else:
                self.ctlpeakwidthcheck.Set3StateValue(2)
                pass

                
            '''

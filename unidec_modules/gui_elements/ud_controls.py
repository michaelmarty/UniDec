import wx
import wx.lib.agw.foldpanelbar as fpb
import os
import unidec_modules.unidectools as ud
import numpy as np
import time
import wx.lib.scrolledpanel as scrolled


class main_controls(wx.Panel):#scrolled.ScrolledPanel):
    def __init__(self, parent, config, pres, panel, iconfile):
        super(wx.Panel, self).__init__(panel)
        #super(scrolled.ScrolledPanel, self).__init__(parent=panel, style=wx.ALL | wx.EXPAND)
        # self.SetAutoLayout(1)
        #self.SetupScrolling(scroll_x=False)

        self.parent = parent
        self.config = config
        self.pres = pres
        self.backgroundchoices = self.config.backgroundchoices
        self.psigsettings = [0, 1, 10, 100]
        self.betasettings = [0, 50, 500, 1000]
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
        self.scrolledpanel=scrolled.ScrolledPanel(self, style=wx.ALL | wx.EXPAND)
        self.scrolledpanel.SetupScrolling()
        size1 = (75, -1)
        self.foldpanels = fpb.FoldPanelBar(self.scrolledpanel, -1, size=(250, 1900), agwStyle=fpb.FPB_VERTICAL)
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

        if self.config.imflag == 1:
            self.ctlmindt = wx.TextCtrl(panel1, value="", size=(50, -1))
            self.ctlmaxdt = wx.TextCtrl(panel1, value="", size=(60, -1))
            dtrange = wx.BoxSizer(wx.HORIZONTAL)
            dtrange.Add(wx.StaticText(panel1, label="Arrival Time Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
            dtrange.Add(self.ctlmindt)
            dtrange.Add(wx.StaticText(panel1, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
            dtrange.Add(self.ctlmaxdt)
            dtrange.Add(wx.StaticText(panel1, label=" ms "), 0, wx.ALIGN_CENTER_VERTICAL)
            sizercontrol1.Add(dtrange, (i, 0), span=(1, 5))
            i += 1

            self.imflag = 1
        else:
            self.imflag = 0

        self.ctlbackcheck = wx.CheckBox(panel1, label="Use Background Subtraction", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_backcheck, self.ctlbackcheck)
        sizercontrol1.Add(self.ctlbackcheck, (i, 0), span=(1, 1), flag=wx.ALIGN_CENTER_VERTICAL)
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

        self.ctlsmooth = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlsmooth, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Gaussian Smoothing: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1
        if self.imflag == 1:
            gbox1b.Add(wx.StaticText(panel1b, label=" m/z"), (i - 1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            gbox1b.Add(wx.StaticText(panel1b, label=" A.T."), (i, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctlsmoothdt = wx.TextCtrl(panel1b, value="", size=size1)
            gbox1b.Add(self.ctlsmoothdt, (i, 1))
            i += 1

            self.ctlpusher = wx.TextCtrl(panel1b, value="", size=size1)
            gbox1b.Add(self.ctlpusher, (i, 1))
            gbox1b.Add(wx.StaticText(panel1b, label="Pusher Interval (\u03BCs)"), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.subtypectl = wx.Choice(panel1b, -1, choices=self.backgroundchoices)
        self.ctlbuff = wx.TextCtrl(panel1b, value="", size=size1)
        self.subtypectl.SetSelection(2)
        gbox1b.Add(self.subtypectl, (i, 0))
        gbox1b.Add(self.ctlbuff, (i, 1))
        i += 1
        if self.imflag == 1:
            gbox1b.Add(wx.StaticText(panel1b, label=" m/z"), (i - 1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            self.ctlsubbuffdt = wx.TextCtrl(panel1b, value="", size=size1)
            gbox1b.Add(self.ctlsubbuffdt, (i, 1))
            gbox1b.Add(wx.StaticText(panel1b, label=" A.T."), (i, 2), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctlaccelvolt = wx.TextCtrl(panel1b, value='', size=size1)
        gbox1b.Add(self.ctlaccelvolt, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Acceleration Voltage (kV): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlintthresh = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlintthresh, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Intensity Threshold: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 0:
            self.ctldatareductionpercent = wx.TextCtrl(panel1b, value="", size=size1)
            gbox1b.Add(self.ctldatareductionpercent, (i, 1), span=(1, 1))
            gbox1b.Add(wx.StaticText(panel1b, label="Data Reduction (%): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

            self.ctldatanorm = wx.CheckBox(panel1b, label="Normalize Data")
            self.ctldatanorm.SetValue(True)
            gbox1b.Add(self.ctldatanorm, (i, 0), span=(1, 2))
            i += 1

        self.ctlbinsize = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlbinsize, (i, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Bin Every: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 0:
            self.ctlbintype = wx.Choice(panel1b, -1, size=(240, 50),
                                        choices=["Linear m/z (Constant " + '\N{GREEK CAPITAL LETTER DELTA}' + "m/z)",
                                                 "Linear resolution (Constant (m/z)/(" + '\N{GREEK CAPITAL LETTER DELTA}' + "m/z))",
                                                 "Nonlinear", "Linear Interpolated", "Linear Resolution Interpolated"])
            gbox1b.Add(self.ctlbintype, (i, 0), span=(1, 2))
            i += 1
        else:
            self.ctlconvertflag = wx.CheckBox(panel1b, label="Compress when converting to .txt")
            self.ctlconvertflag.SetValue(True)
            gbox1b.Add(self.ctlconvertflag, (i, 0), span=(1, 2))
            i += 1

        panel1b.SetSizer(gbox1b)
        gbox1b.Fit(panel1b)

        self.foldpanels.AddFoldPanelWindow(foldpanel1b, panel1b, fpb.FPB_ALIGN_WIDTH)

        if self.config.imflag == 1:
            self.foldpanels.AddFoldPanelWindow(foldpanel1b, wx.StaticText(foldpanel1b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        if self.config.imflag == 1:
            foldpanel1c = self.foldpanels.AddFoldPanel(caption="Ion Mobility Parameters", collapsed=False,
                                                       cbstyle=style1c)
            panel1c = wx.Panel(foldpanel1c, -1)
            gbox1c = wx.GridBagSizer(wx.VERTICAL)

            self.ctltwave = wx.RadioBox(panel1c, label="", choices=["Linear Cell", "Travelling Wave"])
            self.parent.Bind(wx.EVT_RADIOBOX, self.pres.on_flip_twave, self.ctltwave)
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
                gbox1c.Add(wx.StaticText(panel1c, label="Temperature (\u00B0C): "), (3, 0),
                           flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlgasmass = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

                self.ctlto = wx.TextCtrl(panel1c, value='', size=size1)
                gbox1c.Add(self.ctlto, (5, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Dead Time (t\u2080 in ms): "), (5, 0),
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

                self.ctltwavecaltype = wx.Choice(panel1c, -1, choices=list(self.config.twavedict.values()))
                gbox1c.Add(self.ctltwavecaltype, (5, 1), span=(1, 1))
                gbox1c.Add(wx.StaticText(panel1c, label="Calibration Type: "), (5, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltwave.SetSelection(int(self.twave))
            panel1c.SetSizer(gbox1c)
            gbox1c.Fit(panel1c)

            self.foldpanels.AddFoldPanelWindow(foldpanel1c, panel1c, fpb.FPB_ALIGN_WIDTH)
            self.foldpanels.AddFoldPanelWindow(foldpanel1c, wx.StaticText(foldpanel1c, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for UniDec Parameters
        foldpanel2 = self.foldpanels.AddFoldPanel(caption="UniDec Parameters", collapsed=False, cbstyle=style2)
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

        i = 0
        sizercontrol2.Add(zrange, (i, 0), span=(1, 2))
        i += 1
        sizercontrol2.Add(massrange, (i, 0), span=(1, 2))
        i += 1

        if self.config.imflag == 1:
            self.ctlccslb = wx.TextCtrl(panel2, value="", size=(60, -1))
            self.ctlccsub = wx.TextCtrl(panel2, value="", size=(70, -1))
            ccsrange = wx.BoxSizer(wx.HORIZONTAL)
            ccsrange.Add(wx.StaticText(panel2, label="CCS Range: "), 0, wx.ALIGN_CENTER_VERTICAL)
            ccsrange.Add(self.ctlccslb)
            ccsrange.Add(wx.StaticText(panel2, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
            ccsrange.Add(self.ctlccsub)
            ccsrange.Add(wx.StaticText(panel2, label=" \u212B\u00B2  "), 0, wx.ALIGN_CENTER_VERTICAL)
            sizercontrol2.Add(ccsrange, (i, 0), span=(1, 2))
            i += 1

            self.ctlccsbins = wx.TextCtrl(panel2, value="", size=size1)
            sizercontrol2.Add(self.ctlccsbins, (i, 1), span=(1, 2))
            sizercontrol2.Add(wx.StaticText(panel2, label="Sample CCS Every (\u212B\u00B2): "), (i, 0),
                              flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctlmassbins = wx.TextCtrl(panel2, value="", size=size1)
        sizercontrol2.Add(self.ctlmassbins, (i, 1), span=(1, 2))
        sizercontrol2.Add(wx.StaticText(panel2, label="Sample Mass Every (Da): "), (i, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
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

        self.ctlpeakwidthcheck = wx.CheckBox(panel2a, label="Use Automatic m/z Peak Width", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_pw_check, self.ctlpeakwidthcheck)
        sizercontrol2a.Add(self.ctlpeakwidthcheck, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 0:
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
        # sizercontrol2.Add(wx.StaticText(panel2, label="Mass Difference (Da): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2a.Add(self.ctlmsmoothcheck, (i, 0), span=(1, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        sizercontrol2a.Add(self.ctlmolig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        panel2a.SetSizer(sizercontrol2a)
        sizercontrol2a.Fit(panel2a)
        self.foldpanels.AddFoldPanelWindow(foldpanel2a, panel2a, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2a, wx.StaticText(foldpanel2a, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Additional Restraints
        foldpanel2b = self.foldpanels.AddFoldPanel(caption="Additional Deconvolution Parameters", collapsed=True,
                                                   cbstyle=style2b)
        panel2b = wx.Panel(foldpanel2b, -1)
        gbox2b = wx.GridBagSizer(wx.VERTICAL)

        self.ctlmzsig = wx.TextCtrl(panel2b, value="", size=size1)
        self.ctlpsfun = wx.RadioBox(panel2b, label="Peak Shape Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        i = 0
        gbox2b.Add(self.ctlmzsig, (i, 1), span=(1, 2))
        gbox2b.Add(wx.StaticText(panel2b, label="Peak FWHM (Th): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1
        gbox2b.Add(self.ctlpsfun, (i, 0), span=(1, 2))
        i += 1

        if self.config.imflag == 1:
            self.ctldtsig = wx.TextCtrl(panel2b, value="", size=size1)
            gbox2b.Add(self.ctldtsig, (i, 1), span=(1, 2))
            gbox2b.Add(wx.StaticText(panel2b, label="Peak FWHM (ms): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        if self.config.imflag == 0:
            self.ctlbeta = wx.TextCtrl(panel2b, value="", size=size1)
            gbox2b.Add(wx.StaticText(panel2b, label="Beta: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            gbox2b.Add(self.ctlbeta, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctlzzsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Charge Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlzzsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 0:
            self.ctlpsig = wx.TextCtrl(panel2b, value="", size=size1)
            gbox2b.Add(wx.StaticText(panel2b, label="Point Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            gbox2b.Add(self.ctlpsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctlmsig = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label="Mass Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlmsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 1:
            self.ctlcsig = wx.TextCtrl(panel2b, value="", size=size1)
            gbox2b.Add(wx.StaticText(panel2b, label="CCS Smooth Width: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
            gbox2b.Add(self.ctlcsig, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
            i += 1

        self.ctlnumit = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label='Maximum # of Iterations: '), (i, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlnumit, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlpoolflag = wx.RadioBox(panel2b, label="m/z to Mass Transformation",
                                       choices=["Integrate", "Interpolate", "Smart"])
        gbox2b.Add(self.ctlpoolflag, (i, 0), span=(1, 2), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctladductmass = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(self.ctladductmass, (i, 1), span=(1, 1))
        gbox2b.Add(wx.StaticText(panel2b, label="Adduct Mass (Da): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlnegmode = wx.CheckBox(panel2b, label="Negative Mode")
        gbox2b.Add(self.ctlnegmode, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        self.parent.Bind(wx.EVT_CHECKBOX, self.export_gui_to_config, self.ctlnegmode)

        self.ctlmanualassign = wx.CheckBox(panel2b, label="Manual Mode")
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_check_manual, self.ctlmanualassign)
        gbox2b.Add(self.ctlmanualassign, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        if self.config.imflag == 0:
            # self.ctlisotopemode = wx.CheckBox(panel2b, label="Isotope Mode")
            self.ctlisotopemode = wx.Choice(panel2b, -1, size=(100, -1), choices=self.config.isotopechoices)
            gbox2b.Add(self.ctlisotopemode, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlorbimode = wx.CheckBox(panel2b, label="Charge Scaling")
            gbox2b.Add(self.ctlorbimode, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
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

        sb = wx.StaticBox(panel2b, label='Native Charge Offset Range')
        sbs = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        self.ctlminnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        self.ctlmaxnativez = wx.TextCtrl(panel2b, value='', size=(75, -1))
        sbs.Add(self.ctlminnativez, flag=wx.LEFT | wx.EXPAND, border=5)
        sbs.Add(wx.StaticText(panel2b, label=' to '), 0, wx.EXPAND)
        sbs.Add(self.ctlmaxnativez, flag=wx.LEFT | wx.EXPAND, border=5)
        gbox2b.Add(sbs, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        if self.config.imflag == 1:
            sb2 = wx.StaticBox(panel2b, label='Native CCS Offset Range')
            sbs2 = wx.StaticBoxSizer(sb2, orient=wx.HORIZONTAL)
            self.ctlnativeccslb = wx.TextCtrl(panel2b, value='', size=(75, -1))
            self.ctlnativeccsub = wx.TextCtrl(panel2b, value='', size=(75, -1))
            sbs2.Add(self.ctlnativeccslb, flag=wx.LEFT | wx.EXPAND, border=5)
            sbs2.Add(wx.StaticText(panel2b, label=' to '), 0, wx.EXPAND)
            sbs2.Add(self.ctlnativeccsub, flag=wx.LEFT | wx.EXPAND, border=5)
            sbs2.Add(wx.StaticText(panel2b, label=" \u212B\u00B2 "), 0, wx.EXPAND)
            gbox2b.Add(sbs2, (i, 0), span=(1, 2), flag=wx.EXPAND)
            i += 1  # Check

        ddsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.doubledecflag = wx.CheckBox(panel2b, label="Double Deconvolution")
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_double_dec, self.doubledecflag)
        self.doubledecbutton = wx.Button(panel2b, label="Open kernel file")
        self.parent.Bind(wx.EVT_BUTTON, self.on_open_kernel, self.doubledecbutton)
        ddsizer.Add(self.doubledecflag, 0, wx.ALIGN_CENTER_VERTICAL)
        ddsizer.Add(self.doubledecbutton, 0, wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(ddsizer, (i, 0), span=(1, 2))

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
        self.ctlnorm = wx.RadioBox(panel3b, label="Peak Normalization", choices=["None", "Max", "Total"])
        gbox3b.Add(self.ctlnorm, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        self.ctl2dcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)

        self.ctlpeakcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)

        # for mp in self.config.cmaps2:
        #    self.ctl2dcm.Append(mp)
        # for mp in self.config.cmaps:
        #    self.ctlpeakcm.Append(mp)
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

        if self.config.imflag == 1:
            self.cubeplotbutton = wx.Button(panel3b, -1, "Plot Cubes")
            self.parent.Bind(wx.EVT_BUTTON, self.pres.make_cube_plot, self.cubeplotbutton)
            gbox3b.Add(self.cubeplotbutton, (i, 0), span=(1, 2), flag=wx.EXPAND)
            i += 1

        panel3b.SetSizer(gbox3b)
        gbox3b.Fit(panel3b)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, panel3b, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, wx.StaticText(foldpanel3b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        bright = 250
        foldpanel1.SetBackgroundColour(wx.Colour(bright, bright, 255))
        foldpanel1b.SetBackgroundColour(wx.Colour(bright, bright, 255))
        if self.config.imflag == 1:
            foldpanel1c.SetBackgroundColour(wx.Colour(bright, 255, bright))

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
        sizercontrol.Add(self.scrolledpanel,1, wx.EXPAND)

        #Bottom Buttons
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
            buttons = [self.mainbutton, self.expandbutton, self.collapsebutton, self.bluebutton, self.yellowbutton, self.redbutton]
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

        tstart = time.perf_counter()
        self.Freeze()
        self.update_flag = False
        if self.config.batchflag == 0:

            if self.config.imflag == 1 and self.twave != (self.config.twaveflag > 0):
                self.pres.on_flip_twave(0)
            self.ctlmassbins.SetValue(str(self.config.massbins))
            self.ctlstartz.SetValue(str(self.config.startz))
            self.ctlendz.SetValue(str(self.config.endz))
            self.ctlzzsig.SetValue(str(self.config.zzsig))
            self.ctlmzsig.SetValue(str(self.config.mzsig))
            self.ctlpsfun.SetSelection(int(self.config.psfun))
            self.ctlnorm.SetSelection(int(self.config.peaknorm))
            self.ctlmasslb.SetValue(str(self.config.masslb))
            self.ctlmassub.SetValue(str(self.config.massub))
            self.ctlmasslistflag.SetValue(self.config.mfileflag)
            self.doubledecflag.SetValue(self.config.doubledec)  # DoubleDec config import
            if self.config.kernel != "":
                kernel_name = os.path.splitext(os.path.basename(self.config.kernel))[0]
                self.doubledecbutton.SetLabel(kernel_name)
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
                self.ctldatareductionpercent.SetValue(str(self.config.reductionpercent))
                self.ctlmanualassign.SetValue(self.config.manualfileflag)
                self.ctlisotopemode.SetSelection(self.config.isotopemode)
                self.ctlorbimode.SetValue(int(self.config.orbimode))
                self.ctldatanorm.SetValue(int(self.config.datanorm))
                self.ctlbintype.SetSelection(int(self.config.linflag))
                self.ctlpsig.SetValue(str(self.config.psig))
                self.ctlbeta.SetValue(str(self.config.beta))
            self.ctldiscrete.SetValue(self.config.discreteplot)
            self.ctlpublicationmode.SetValue(self.config.publicationmode)
            self.ctlrawflag.SetSelection(self.config.rawflag)

            if self.config.adductmass < 0:
                self.ctlnegmode.SetValue(1)
            else:
                self.ctlnegmode.SetValue(0)

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
                    self.ctltwavecaltype.SetSelection(list(self.config.twavedict.keys()).index(self.config.twaveflag))

            try:
                x = float(self.config.integratelb)
                y = float(self.config.integrateub)
                self.ctlintlb.SetValue(str(x))
                self.ctlintub.SetValue(str(y))
            except (ValueError, TypeError):
                self.ctlintlb.SetValue("")
                self.ctlintub.SetValue("")
            if self.config.imflag == 0:
                try:
                    if self.config.aggressiveflag == 1:
                        self.parent.menu.advancedmenu.Check(id=402, check=True)
                    elif self.config.aggressiveflag == 2:
                        self.parent.menu.advancedmenu.Check(id=403, check=True)
                    else:
                        self.parent.menu.advancedmenu.Check(id=401, check=True)
                except:
                    print("No Menu Found")

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
            if self.config.imflag == 1:
                self.ctlmindt.SetValue(str(self.config.mindt))
                self.ctlmaxdt.SetValue(str(self.config.maxdt))
        # print("4: %.2gs" % (time.perf_counter() - tstart))
        self.update_flag = True
        try:
            self.update_quick_controls()
        except Exception as e:
            print("Error updating quick controls", e)
        # print("5: %.2gs" % (time.perf_counter() - tstart))
        self.Thaw()

    def export_gui_to_config(self, e=None):
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
        self.config.doubledec = self.doubledecflag.GetValue()  # DoubleDec config export
        # self.config.kernel =
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
            self.config.reductionpercent = ud.string_to_value(self.ctldatareductionpercent.GetValue())
            self.config.isotopemode = int(self.ctlisotopemode.GetSelection())
            self.config.orbimode = int(self.ctlorbimode.GetValue())
            self.config.datanorm = int(self.ctldatanorm.GetValue())
            self.config.psig = ud.string_to_value(self.ctlpsig.GetValue())
            self.config.beta = ud.string_to_value(self.ctlbeta.GetValue())
            self.config.manualfileflag = int(self.ctlmanualassign.GetValue())
            self.config.linflag = self.ctlbintype.GetSelection()
            if self.config.mzbins == 0:
                self.config.linflag = 2
                self.ctlbintype.SetSelection(int(self.config.linflag))
        self.config.discreteplot = int(self.ctldiscrete.GetValue())
        self.config.publicationmode = int(self.ctlpublicationmode.GetValue())
        self.config.rawflag = self.ctlrawflag.GetSelection()

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
                self.config.twaveflag = list(self.config.twavedict.keys())[self.ctltwavecaltype.GetSelection()]

            if not self.config.mindt and not ud.isempty(self.pres.eng.data.rawdata3):
                self.config.mindt = np.amin(self.pres.eng.data.rawdata3[:, 1])
            if not self.config.maxdt and not ud.isempty(self.pres.eng.data.rawdata3):
                self.config.maxdt = np.amax(self.pres.eng.data.rawdata3[:, 1])

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
        self.ctlmasslistflag.SetToolTip(wx.ToolTip(
            "Limit deconvolution to specific masses +/- some window.\nDefine in Tools>Oligomer and Mass Tools."))
        self.ctlmtabsig.SetToolTip(
            wx.ToolTip("Set window for mass limitations. Setting to 0 will force only listed masses."))
        self.ctlpsfun.SetToolTip(wx.ToolTip("Expected peak shape.\nSee Tools>Peak Width Tool for more tools."))
        self.rununidec.SetToolTip(wx.ToolTip("Write Configuration File, Run UniDec, and Plot Results. (Ctrl+R)"))
        self.ctlmzsig.SetToolTip(wx.ToolTip(
            "Expected peak FWHM in m/z (Th).\nFor nonlinear mode, minimum FWHM"
            "\nSee Tools>Peak Width Tool for more tools."))
        self.ctlzzsig.SetToolTip(wx.ToolTip(
            "Parameter for defining the width of the charge state smooth."
            "\nUniDec will use a mean filter of width 2n+1 on log_e of the charge distribution"))
        self.ctlmassub.SetToolTip(wx.ToolTip(
            "Maximum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlmassbins.SetToolTip(wx.ToolTip("Sets the resolution of the zero-charge mass spectrum"))
        self.ctlmasslb.SetToolTip(wx.ToolTip(
            "Minimum allowed mass in deconvolution.\nTip: A negative value will force the axis to the absolute value."))
        self.ctlstartz.SetToolTip(wx.ToolTip("Minimum allowed charge state in deconvolution."))
        self.ctlendz.SetToolTip(wx.ToolTip("Maximum allowed charge state in deconvolution."))
        self.ctlsmooth.SetToolTip(wx.ToolTip("Gaussian smooth sigma in units of data point number."))
        self.ctlbinsize.SetToolTip(wx.ToolTip(
            "Controls Linearization.\nConstant bin size (Th) for Linear m/z"
            "\nMinimum bin size (Th) for Linear Resolution"
            "\nNumber of data points compressed together for Nonlinear"))
        self.ctlintthresh.SetToolTip(
            wx.ToolTip("Set intensity threshold. Data points below threshold are excluded from deconvolution."))
        self.ctlbuff.SetToolTip(wx.ToolTip(
            "Background subtraction parameters.\nMinimum: 0=off 1=on"
            "\nLine: The first and last n data points will be averaged;"
            "\n a line between the two averages will be subtracted.\nCurved: Width of smoothed background"))
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
            self.ctldatanorm.SetToolTip(wx.ToolTip("Normalize Data and Results"))
            self.ctldatareductionpercent.SetToolTip(
                wx.ToolTip(
                    "Reduces the amount of data by removing everything below a threshold.\nSets the threshold to fit the percentage of data to remove."))
            self.ctlisotopemode.SetToolTip(wx.ToolTip(
                "Use isotopic distributions in deconvolution.\nOutput either monoisotopic or average masses"))
            self.ctlorbimode.SetToolTip(wx.ToolTip("Scale intensities by 1/z to compensate for induced charge effects"))
            self.ctlmanualassign.SetToolTip(wx.ToolTip("Use manual assignments. See Tools>Manual Assignment"))
            self.ctlbintype.SetToolTip(wx.ToolTip(
                "Sets how to bin the data\nValue set by above with Bin Every\nLinear bins with linear m/z axis"
                "\nLinear Resolution bins with m/z axis that has a constant resolution"
                "\nNonlinear merges adjacent data points"
                "\nInterpolation uses the same axes but with interpolation instead of integration"))
            self.ctlpsig.SetToolTip(wx.ToolTip(
                "Parameter for defining the width of the data point smooth."
                "\nUniDec will weight +/- n data points to have the same charge state."))
            self.ctlbeta.SetToolTip(wx.ToolTip(
                "Parameter for defining the degree of Softmax distribution applied to the charge state vectors."
                "\n0 will shut it off."))
            self.ctlpselect.SetToolTip(wx.ToolTip(
                "Select whether to smooth nearby data points to have similar charge assignments"))
            self.ctlbselect.SetToolTip(wx.ToolTip(
                "Select whether to suppress deconvolution artifacts"))
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
                self.ctlto.SetToolTip(
                    wx.ToolTip("Set dead time in instrument.\nDrift Time=Arrival Time - Dead Time"))
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
            self.ctlcsig.SetToolTip(
                wx.ToolTip("Width of CCS smooth.\nCCS difference is taken as the CCS bin size."))
            self.ctlnativeccslb.SetToolTip(wx.ToolTip("Sets lower bound on offset from native CCS."))
            self.ctlnativeccsub.SetToolTip(wx.ToolTip("Sets upper bound on offset from native CCS."))
        self.ctlpoolflag.SetToolTip(wx.ToolTip(
            "Sets type of conversion from m/z to mass.\nIntegration:\n\tEach m/z bin goes to the nearest mass bin"
            "\n\tBest for undersampled masses\nInterpolation:\n\tEach mass value interpolates its value in m/z space"
            "\n\tBest for oversampled mass data"
            "\nSmart:\n\tDetermines interpolation vs. integration automatically\n\tfrom the data and blends the two."
            "\n\tUse Smart in most cases (default). \n\tNot available for IM-MS data."))
        self.ctlmsmoothcheck.SetToolTip(
            wx.ToolTip("Select whether to assume a smooth distribution spaced by a known mass difference"))
        self.ctlzsmoothcheck.SetToolTip(
            wx.ToolTip("Select whether to assume a smooth charge state distribution"))
        self.ctlpeakwidthcheck.SetToolTip(
            wx.ToolTip("Select whether to incorporated the automatic peak width in deconvolution"))
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

    def on_double_dec(self, e):
        """
        Sets the DoubleDec flag in the config to the appropriate value
        :param e: Dummy wx event. Unused.
        :return: None
        """
        self.config.doubledec = self.doubledecflag.GetValue()
        print("DoubleDec mode is", self.config.doubledec)

    def on_open_kernel(self, e):
        """
        Opens a kernel file and saves the path in the config. The button text is
        also updated to match. Note that the user is supposed to select the raw
        m/z file--the corresponding mass file will be found if it exists. This way,
        possible file naming issues are avoided.
        :param e: Triggering event. Unused.
        :return: None
        """
        kernel_path = self.pres.on_open_kernel()  # The file the user chose
        kernel_path2 = kernel_path  # The _mass.txt file
        kernel_name = os.path.basename(kernel_path)

        if kernel_name.split('_')[-1] != "mass.txt":
            # If the filename doesn't end in "_mass.txt", it's not a mass file
            is_mass = False
        else:
            # Otherwise, it might be a mass file. Check if linear
            kernel_dat = np.loadtxt(kernel_path)
            kernel_diff = np.diff(kernel_dat[:, 0])
            is_mass = np.all(kernel_diff == kernel_diff[0])
        if not is_mass:  # If not a mass file, find the mass file
            bare_name = os.path.splitext(kernel_name)[0]
            kernel_path2 = os.path.dirname(kernel_path) + "\\" + bare_name + "_unidecfiles\\" + \
                           bare_name + "_mass.txt"

        try:
            with open(kernel_path2, "r") as f:
                self.doubledecbutton.SetLabel(os.path.splitext(os.path.basename(kernel_path2))[0])
                self.config.kernel = kernel_path2
        except IOError:
            print("Please deconvolve the m/z file [" + kernel_name + "] with UniDec first.")
        print(self.config.kernel)

    def on_z_smooth(self, e):
        value = self.ctlzsmoothcheck.Get3StateValue()
        if value == 1:
            self.ctlzzsig.SetValue("1")
        elif value == 0:
            self.ctlzzsig.SetValue("0")
        self.export_gui_to_config()

    def on_backcheck(self, e):
        value = self.ctlbackcheck.Get3StateValue()
        if value == 1:
            self.ctlbuff.SetValue("100")
            self.subtypectl.SetSelection(2)
        elif value == 0:
            self.ctlbuff.SetValue("0")
        self.export_gui_to_config()
        self.update_quick_controls()

    def on_m_smooth(self, e):
        value = self.ctlmsmoothcheck.Get3StateValue()
        if value == 1:
            self.ctlmsig.SetValue("1")
        elif value == 0:
            self.ctlmsig.SetValue("0")
        self.export_gui_to_config()

    def on_pw_check(self, e):
        value = self.ctlpeakwidthcheck.Get3StateValue()
        if value == 1:
            self.ctlmzsig.SetValue(str(self.config.automzsig))
            self.ctlpsfun.SetSelection(self.config.autopsfun)
        elif value == 0:
            self.ctlmzsig.SetValue("0")
        self.export_gui_to_config()

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

    def on_expand_blue(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in [0,1]:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)
        pass

    def on_expand_yellow(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in [2,3,4]:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)
        pass

    def on_expand_red(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in [5, 6]:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)
        pass

    def on_expand_main(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in [0, 2, 3, 5]:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)
        pass

    def bind_changes(self, e=None):
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlzzsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlmsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlmzsig)
        self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlbuff)
        self.parent.Bind(wx.EVT_CHOICE, self.update_quick_controls, self.subtypectl)
        self.parent.Bind(wx.EVT_RADIOBOX, self.update_quick_controls, self.ctlpsfun)
        if self.config.imflag == 0:
            self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlbeta)
            self.parent.Bind(wx.EVT_TEXT, self.update_quick_controls, self.ctlpsig)

    def update_quick_controls(self, e=None):
        if self.update_flag:
            # Background Subtraction
            try:
                value = float(self.ctlbuff.GetValue())
            except:
                value = -1
            selection = self.subtypectl.GetSelection()
            if value == 100 and selection == 2:
                self.ctlbackcheck.Set3StateValue(1)
            elif value == 0:
                self.ctlbackcheck.Set3StateValue(0)
            else:
                self.ctlbackcheck.Set3StateValue(2)

            # Z Box
            try:
                value = float(self.ctlzzsig.GetValue())
            except:
                value = -1
            if value == 0:
                self.ctlzsmoothcheck.Set3StateValue(0)
            elif value == 1:
                self.ctlzsmoothcheck.Set3StateValue(1)
            else:
                self.ctlzsmoothcheck.Set3StateValue(2)

            # M box
            try:
                value = float(self.ctlmsig.GetValue())
            except:
                value = -1
            if value == 0:
                self.ctlmsmoothcheck.Set3StateValue(0)
            elif value == 1:
                self.ctlmsmoothcheck.Set3StateValue(1)
            else:
                self.ctlmsmoothcheck.Set3StateValue(2)

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

            if self.config.imflag == 0:
                # beta box
                try:
                    value = float(self.ctlbeta.GetValue())
                except:
                    valu = -1
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

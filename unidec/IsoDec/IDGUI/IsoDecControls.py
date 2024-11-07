import wx
import wx.lib.agw.foldpanelbar as fpb
import os
import unidec.tools as ud
import numpy as np
import time
import wx.lib.scrolledpanel as scrolled

class main_controls(wx.Panel):  # scrolled.ScrolledPanel):
    # noinspection PyMissingConstructor
    def __init__(self, parent, config, pres, panel, iconfile):
        super(wx.Panel, self).__init__(panel)
        self.parent = parent
        self.config = config
        self.pres = pres
        self.backgroundchoices = self.config.backgroundchoices
        self.psigsettings = [0, 1, 10, 100]
        self.betasettings = [0, 50, 500, 1000]
        self.update_flag = True

        # ..........................
        #
        # Sizers for Controls
        #
        # .............................
        sizercontrol = wx.BoxSizer(wx.VERTICAL)

        # Setting up main fold controls
        self.scrolledpanel = scrolled.ScrolledPanel(self, style=wx.ALL | wx.EXPAND)
        self.scrolledpanel.SetupScrolling()
        size1 = (75, -1)
        self.foldpanels = fpb.FoldPanelBar(self.scrolledpanel, -1, size=(250, 800), agwStyle=fpb.FPB_VERTICAL)
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

        self.ctlbackcheck = wx.CheckBox(panel1, label="Use Background Subtraction", style=wx.CHK_3STATE)
        self.parent.Bind(wx.EVT_CHECKBOX, self.on_backcheck, self.ctlbackcheck)
        sizercontrol1.Add(self.ctlbackcheck, (i, 0), span=(1, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlcentroided = wx.CheckBox(panel1, label="Data is Centroided")
        sizercontrol1.Add(self.ctlcentroided, (i, 0), span=(1, 1), flag=wx.ALIGN_CENTER_VERTICAL)
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

        self.subtypectl = wx.Choice(panel1b, -1, choices=self.backgroundchoices)
        self.subtypectl.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
        self.ctlbuff = wx.TextCtrl(panel1b, value="", size=size1)
        self.subtypectl.SetSelection(2)
        gbox1b.Add(self.subtypectl, (i, 0))
        gbox1b.Add(self.ctlbuff, (i, 1))
        i += 1

        self.ctlintthresh = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctlintthresh, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Intensity Threshold: "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctldatareductionpercent = wx.TextCtrl(panel1b, value="", size=size1)
        gbox1b.Add(self.ctldatareductionpercent, (i, 1), span=(1, 1))
        gbox1b.Add(wx.StaticText(panel1b, label="Data Reduction (%): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctldatanorm = wx.CheckBox(panel1b, label="Normalize Data")
        self.ctldatanorm.SetValue(True)
        gbox1b.Add(self.ctldatanorm, (i, 0), span=(1, 2))
        i += 1

        panel1b.SetSizer(gbox1b)
        gbox1b.Fit(panel1b)

        self.foldpanels.AddFoldPanelWindow(foldpanel1b, panel1b, fpb.FPB_ALIGN_WIDTH)
        i=0

        # Panel for unidec Parameters
        foldpanel2 = self.foldpanels.AddFoldPanel(caption="UniDec Parameters", collapsed=False, cbstyle=style2)
        panel2 = wx.Panel(foldpanel2, -1)
        sizercontrol2 = wx.GridBagSizer(wx.VERTICAL)


        self.ctlmassbins = wx.TextCtrl(panel2, value="", size=size1)
        sizercontrol2.Add(self.ctlmassbins, (i, 1), span=(1, 2))
        sizercontrol2.Add(wx.StaticText(panel2, label="Sample Mass Every (Da): "), (i, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.rununidec = wx.Button(panel2, -1, "Run IsoDec")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_unidec_button, self.rununidec)
        sizercontrol2.Add(self.rununidec, (i, 0), span=(1, 2), flag=wx.EXPAND)
        i += 1

        panel2.SetSizer(sizercontrol2)
        sizercontrol2.Fit(panel2)
        self.foldpanels.AddFoldPanelWindow(foldpanel2, panel2, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2, wx.StaticText(foldpanel2, -1, " "), fpb.FPB_ALIGN_WIDTH)


        # Panel for Additional Restraints
        foldpanel2b = self.foldpanels.AddFoldPanel(caption="Additional Deconvolution Parameters", collapsed=True,
                                                   cbstyle=style2b)
        panel2b = wx.Panel(foldpanel2b, -1)
        gbox2b = wx.GridBagSizer(wx.VERTICAL)

        i=0
        self.ctlwindow = wx.TextCtrl(panel2b, value="", size=size1)
        self.ctlthresh = wx.TextCtrl(panel2b, value="", size=size1)
        gbox2b.Add(self.ctlwindow, (i, 1))
        gbox2b.Add(wx.StaticText(panel2b, label="Peak Detection Range (Da): "), (i, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        i+=1
        gbox2b.Add(self.ctlthresh, (i, 1))
        gbox2b.Add(wx.StaticText(panel2b, label="Peak Detection Threshold: "), (i, 0),
                          flag=wx.ALIGN_CENTER_VERTICAL)
        i+=1
        self.ctlnumit = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(wx.StaticText(panel2b, label='Knockdown Rounds: '), (i, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        gbox2b.Add(self.ctlnumit, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctladductmass = wx.TextCtrl(panel2b, value='', size=size1)
        gbox2b.Add(self.ctladductmass, (i, 1), span=(1, 1))
        gbox2b.Add(wx.StaticText(panel2b, label="Adduct Mass (Da): "), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlnegmode = wx.CheckBox(panel2b, label="Negative Mode")
        gbox2b.Add(self.ctlnegmode, (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.parent.Bind(wx.EVT_CHECKBOX, self.export_gui_to_config, self.ctlnegmode)

        # self.ctlisotopemode = wx.CheckBox(panel2b, label="Isotope Mode")
        self.ctlisotopemode = wx.Choice(panel2b, -1, size=(100, -1), choices=self.config.isotopechoices)
        self.ctlisotopemode.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
        gbox2b.Add(self.ctlisotopemode, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)


        panel2b.SetSizer(gbox2b)
        gbox2b.Fit(panel2b)
        self.foldpanels.AddFoldPanelWindow(foldpanel2b, panel2b, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel2b, wx.StaticText(foldpanel2b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        # Panel for Peak Selection and Plotting
        foldpanel3 = self.foldpanels.AddFoldPanel(caption="Peak Selection and Plotting", collapsed=False,
                                                  cbstyle=style3)
        panel3 = wx.Panel(foldpanel3, -1)

        sizercontrol3 = wx.GridBagSizer(wx.VERTICAL)


        self.plotbutton2 = wx.Button(panel3, -1, "Plot Peaks")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_plot_peaks, self.plotbutton2)
        sizercontrol3.Add(self.plotbutton2, (0, 1), span=(1, 1), flag=wx.EXPAND)

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

        self.ctlpeakcm = wx.ComboBox(panel3b, wx.ID_ANY, style=wx.CB_READONLY)
        self.ctlpeakcm.Bind(wx.EVT_MOUSEWHEEL, self.on_mousewheel)
        self.ctlpeakcm.AppendItems(self.config.cmaps)

        gbox3b.Add(self.ctlpeakcm, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
        gbox3b.Add(wx.StaticText(panel3b, label='Peaks Color Map: '), (i, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        i += 1

        self.ctlpublicationmode = wx.CheckBox(panel3b, label="Publication Mode")
        gbox3b.Add(self.ctlpublicationmode, (i, 1), flag=wx.ALIGN_CENTER_VERTICAL)
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

        self.replotbutton = wx.Button(panel3b, -1, "Replot")
        self.parent.Bind(wx.EVT_BUTTON, self.pres.on_replot, self.replotbutton)
        gbox3b.Add(self.replotbutton, (i, 0), span=(1, 1), flag=wx.EXPAND)

        panel3b.SetSizer(gbox3b)
        gbox3b.Fit(panel3b)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, panel3b, fpb.FPB_ALIGN_WIDTH)
        self.foldpanels.AddFoldPanelWindow(foldpanel3b, wx.StaticText(foldpanel3b, -1, " "), fpb.FPB_ALIGN_WIDTH)

        bright = 250
        foldpanel1.SetBackgroundColour(wx.Colour(bright, bright, 255))
        foldpanel1b.SetBackgroundColour(wx.Colour(bright, bright, 255))
        foldpanel2.SetBackgroundColour(wx.Colour(255, 255, bright))
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

        # Set top sizer
        self.SetSizer(sizercontrol)
        sizercontrol.Fit(self)

        # Add remaining things
        self.setup_tool_tips()


    def import_config_to_gui(self):
        """
        Imports parameters from the config object to the GUI.
        :return: None
        """

        tstart = time.perf_counter()
        self.Freeze()
        self.update_flag = False
        if self.config.batchflag == 0:
            self.ctlcentroided.SetValue(self.config.centroided)
            self.ctlmassbins.SetValue(str(self.config.massbins))
            self.ctlnorm.SetSelection(int(self.config.peaknorm))
            self.ctlbuff.SetValue(str(self.config.subbuff))
            self.subtypectl.SetSelection(int(self.config.subtype))
            self.ctlwindow.SetValue(str(self.config.peakwindow))
            self.ctlthresh.SetValue(str(self.config.peakthresh))
            self.ctlthresh2.SetValue(str(self.config.peakplotthresh))
            self.ctlsep.SetValue(str(self.config.separation))
            self.ctlintthresh.SetValue(str(self.config.intthresh))
            self.ctladductmass.SetValue(str(self.config.adductmass))
            self.ctlnumit.SetValue(str(self.config.numit))
            self.ctldatareductionpercent.SetValue(str(self.config.reductionpercent))
            self.ctlisotopemode.SetSelection(self.config.isotopemode)
            self.ctldatanorm.SetValue(int(self.config.datanorm))
            self.ctlpublicationmode.SetValue(self.config.publicationmode)

            if self.config.adductmass < 0:
                self.ctlnegmode.SetValue(1)
            else:
                self.ctlnegmode.SetValue(0)

            try:
                self.ctlpeakcm.SetSelection(self.config.cmaps.index(self.config.peakcmap))
            except ValueError:
                print("Could not find the specified color map. Try upgrading to the latest version of matplotlib.")
                import matplotlib
                print("Current version:", matplotlib.__version__)
                # Revert to the defaults
                self.ctlpeakcm.SetSelection(self.config.cmaps.index("rainbow"))

        # If the batchflag is not 1, it will import the data range as well
        if self.config.batchflag != 1:
            self.ctlminmz.SetValue(str(self.config.minmz))
            self.ctlmaxmz.SetValue(str(self.config.maxmz))

        self.Thaw()

    def export_gui_to_config(self, e=None):
        """
        Exports parameters from the GUI to the config object.
        :return: None
        """
        self.config.centroided = self.ctlcentroided.GetValue()
        self.config.minmz = ud.string_to_value(self.ctlminmz.GetValue())
        self.config.maxmz = ud.string_to_value(self.ctlmaxmz.GetValue())
        self.config.subbuff = ud.string_to_value(self.ctlbuff.GetValue())
        self.config.subtype = self.subtypectl.GetSelection()
        self.config.intthresh = ud.string_to_value(self.ctlintthresh.GetValue())
        self.config.massbins = ud.string_to_value(self.ctlmassbins.GetValue())
        self.config.peaknorm = self.ctlnorm.GetSelection()
        self.config.peakwindow = ud.string_to_value(self.ctlwindow.GetValue())
        self.config.peakthresh = ud.string_to_value(self.ctlthresh.GetValue())
        self.config.peakplotthresh = ud.string_to_value(self.ctlthresh2.GetValue())
        self.config.separation = ud.string_to_value(self.ctlsep.GetValue())
        self.config.adductmass = ud.string_to_value(self.ctladductmass.GetValue())
        self.config.numit = ud.string_to_int(self.ctlnumit.GetValue())
        self.config.reductionpercent = ud.string_to_value(self.ctldatareductionpercent.GetValue())
        self.config.isotopemode = int(self.ctlisotopemode.GetSelection())
        self.config.datanorm = int(self.ctldatanorm.GetValue())
        self.config.publicationmode = int(self.ctlpublicationmode.GetValue())

        try:
            test = float(self.config.adductmass)
        except:
            self.config.adductmass=0

        if self.ctlnegmode.GetValue() == 1:
            # print("Negative Ion Mode")
            if self.config.adductmass > 0:
                self.config.adductmass *= -1
                self.ctladductmass.SetValue(str(self.config.adductmass))
        else:
            # print("Positive Ion Mode")\
            if self.config.adductmass < 0:
                self.config.adductmass *= -1
                self.ctladductmass.SetValue(str(self.config.adductmass))

        self.config.peakcmap = str(self.ctlpeakcm.GetStringSelection())


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
        self.plotbutton2.SetToolTip(wx.ToolTip("Plot individual peak species in m/z. (Ctrl+K)"))
        self.rununidec.SetToolTip(wx.ToolTip("Write Configuration File, Run UniDec, and Plot Results. (Ctrl+R)"))
        self.ctlmassbins.SetToolTip(wx.ToolTip("Sets the resolution of the zero-charge mass spectrum"))
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

        self.ctldatanorm.SetToolTip(wx.ToolTip("Normalize Data and Results"))
        self.ctldatareductionpercent.SetToolTip(
            wx.ToolTip(
                "Reduces the amount of data by removing everything below a threshold.\nSets the threshold to fit the percentage of data to remove."))
        self.ctlisotopemode.SetToolTip(wx.ToolTip(
            "Use isotopic distributions in deconvolution.\nOutput either monoisotopic or average masses"))
        self.ctlnumit.SetToolTip(wx.ToolTip(
            "Maximum number of iterations. Note: Deconvolution will stop automically before this if it converges."))

        self.ctlpublicationmode.SetToolTip(wx.ToolTip("Set plots to look good for publication rather than utility"))

        self.ctlpeakcm.SetToolTip(wx.ToolTip("Set the color function for the peaks"))

        self.ctlnorm.SetToolTip(wx.ToolTip(
            "Sets normalization of mass data.\nMaximum will normalize so that the maximum value is %100."
            "\nTotal will normalize so that the sum of all peaks is %100"))
        self.replotbutton.SetToolTip(wx.ToolTip("Replot some of the plots. (Ctrl+N)"))
        pass

    # .......................................................
    #
    #  The Main Panel
    #
    # .......................................................

    def on_backcheck(self, e):
        value = self.ctlbackcheck.Get3StateValue()
        if value == 1:
            self.ctlbuff.SetValue("100")
            self.subtypectl.SetSelection(2)
        elif value == 0:
            self.ctlbuff.SetValue("0")
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
            if i in [0, 1]:
                self.foldpanels.Expand(fp)
            else:
                self.foldpanels.Collapse(fp)
        pass

    def on_expand_yellow(self, e=None):
        num = self.foldpanels.GetCount()
        for i in range(0, num):
            fp = self.foldpanels.GetFoldPanel(i)
            if i in [2, 3, 4]:
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

    def on_mousewheel(self, e):
        pass



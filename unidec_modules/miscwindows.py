import os
import wx

from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


def setup_tab_box(tab, plot):
    box1 = wx.BoxSizer(wx.VERTICAL)
    box1.Add(plot, 0, wx.EXPAND)
    tab.SetSizerAndFit(box1)


class SingleInputDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        # self.SetSize((400, 200))

    def InitUI(self, title="", message="", defaultvalue=""):
        self.SetTitle(title)
        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox = wx.TextCtrl(self.pnl, value=defaultvalue)
        self.hbox.Add(wx.StaticText(self.pnl, label=message), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox.Add(self.inputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.pnl.SetSizer(self.hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        hboxend.Add(okButton)

        self.vbox.Add(self.pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)

    def OnClose(self, e):
        self.value = self.inputbox.GetValue()
        self.Destroy()
        self.EndModal(0)


class AdditionalParameters(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((400, 200))
        self.SetTitle("Set Experimental Parameters")

    def InitUI(self, config):
        self.config = config

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Set Experimental Parameters')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox5 = wx.TextCtrl(self.pnl, value=str(self.config.inflate))
        self.hbox5.Add(wx.StaticText(self.pnl, label='Peak Inflation Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox5.Add(self.inputbox5, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox5.Add(wx.StaticText(self.pnl, label=' x Peak Shape '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox5, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox6 = wx.TextCtrl(self.pnl, value=str(self.config.damp))
        self.hbox6.Add(wx.StaticText(self.pnl, label='Damping Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox6.Add(self.inputbox6, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox6.Add(wx.StaticText(self.pnl, label=' '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox6, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox6b = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox6b = wx.ComboBox(self.pnl, wx.ID_ANY, style=wx.CB_READONLY,
                                      choices=["None", "Harmonics", "Satellites", "Both", "Total", "Ultimate",
                                               "Isotope"])
        self.inputbox6b.SetSelection(self.config.suppression)
        self.hbox6b.Add(wx.StaticText(self.pnl, label='Suppression Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox6b.Add(self.inputbox6b, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox6b.Add(wx.StaticText(self.pnl, label=' Be careful...'), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox6b, 1, wx.ALIGN_CENTER_VERTICAL)

        self.pnl.SetSizer(self.sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

    def OnClose(self, e):
        self.config.suppression = ud.string_to_int(self.inputbox6b.GetSelection())
        self.config.inflate = ud.string_to_value(self.inputbox5.GetValue())
        self.config.damp = ud.string_to_value(self.inputbox6.GetValue())
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)


class SaveFigureDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((500, 300))
        self.SetTitle("Save Figures")

    def InitUI(self, config):
        self.config = config
        self.directory = None
        self.header = None
        self.extension = None
        self.transparent = None
        self.rect = None
        self.figsize = None

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Save Settings')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(wx.StaticText(self.pnl, label="Directory:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.dirinput = wx.TextCtrl(self.pnl, value=os.getcwd(), size=(300, 20))
        self.hbox.Add(self.dirinput, 0, wx.EXPAND)
        self.dirbutton = wx.Button(self.pnl, label="...", size=(20, 20))
        self.hbox.Add(self.dirbutton, 1, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_choose_dir, self.dirbutton)
        self.sbs.Add(self.hbox, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.headerbox = wx.TextCtrl(self.pnl, value=str(self.config.outfname), size=(300, 20))
        self.hbox5.Add(wx.StaticText(self.pnl, label='File Name Header: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox5.Add(self.headerbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox5, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.extbox = wx.TextCtrl(self.pnl, value=str("png"))
        self.hbox6.Add(wx.StaticText(self.pnl, label='Extension: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox6.Add(self.extbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.tbox = wx.CheckBox(self.pnl, label=str("Transparent"))
        self.hbox6.Add(self.tbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox6, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        self.widebox = wx.TextCtrl(self.pnl, value=str(6), size=(40, 20))
        self.tallbox = wx.TextCtrl(self.pnl, value=str(5), size=(40, 20))
        self.hbox7.Add(wx.StaticText(self.pnl, label='Figure Size: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox7.Add(self.widebox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox7.Add(wx.StaticText(self.pnl, label=' wide by '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox7.Add(self.tallbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox7.Add(wx.StaticText(self.pnl, label=' tall (inches) '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox7, 1, wx.ALIGN_CENTER_VERTICAL)

        self.hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        size = (40, 20)
        self.b1 = wx.TextCtrl(self.pnl, value=str(0.1), size=size)
        self.b2 = wx.TextCtrl(self.pnl, value=str(0.1), size=size)
        self.b3 = wx.TextCtrl(self.pnl, value=str(0.8), size=size)
        self.b4 = wx.TextCtrl(self.pnl, value=str(0.8), size=size)
        self.hbox8.Add(wx.StaticText(self.pnl, label='Plot Position: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox8.Add(self.b1, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox8.Add(self.b2, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox8.Add(self.b3, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox8.Add(self.b4, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox8.Add(wx.StaticText(self.pnl, label=' Left Bottom Width Height (%) '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox8, 1, wx.ALIGN_CENTER_VERTICAL)

        self.pnl.SetSizer(self.sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

    def OnClose(self, e):
        self.directory = self.dirinput.GetValue()
        self.header = self.headerbox.GetValue()
        self.extension = self.extbox.GetValue()
        self.transparent = self.tbox.GetValue()
        self.rect = [float(self.b1.GetValue()), float(self.b2.GetValue()), float(self.b3.GetValue()),
                     float(self.b4.GetValue())]
        self.figsize = [float(self.widebox.GetValue()), float(self.tallbox.GetValue())]
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def on_choose_dir(self, e):
        if self.directory is not None:
            dlg = wx.DirDialog(None, "Choose Top Directory", "",
                               wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)  # ,defaultPath=self.directory)
        else:
            dlg = wx.DirDialog(None, "Choose Top Directory", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.directory = dlg.GetPath()
            self.dirinput.SetValue(self.directory)
            # print self.directory
        dlg.Destroy()


class FileNameDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((250, 300))
        self.SetTitle("Change File Names")

    def InitUI(self, config):
        self.config = config
        self.defaultin = self.config.infname
        self.defaultout = self.config.outfname
        self.defaultconf = self.config.confname
        self.defaultmassfile = self.config.mfile
        self.defaulttruncfile = self.config.manualfile
        self.defaultofile = self.config.ofile
        self.defaultmatchfile = self.config.matchfile
        self.defaultpeaksfile = self.config.peaksfile
        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Rename Files Into and Out of UniDec')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox = wx.TextCtrl(self.pnl, value=self.defaultin)
        self.hbox1.Add(wx.StaticText(self.pnl, label='Input File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox1.Add(self.inputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox1, 0, wx.EXPAND)

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.confbox = wx.TextCtrl(self.pnl, value=self.defaultconf)
        self.hbox2.Add(wx.StaticText(self.pnl, label='Configuration File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.confbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox2, 0, wx.EXPAND)

        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.outputbox = wx.TextCtrl(self.pnl, value=self.defaultout)
        self.hbox3.Add(wx.StaticText(self.pnl, label='Output File Headers: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox3.Add(self.outputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox3, 0, wx.EXPAND)

        self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.massbox = wx.TextCtrl(self.pnl, value=self.defaultmassfile)
        self.hbox4.Add(wx.StaticText(self.pnl, label='Mass List File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox4.Add(self.massbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox4, 0, wx.EXPAND)

        self.hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.truncbox = wx.TextCtrl(self.pnl, value=self.defaulttruncfile)
        self.hbox5.Add(wx.StaticText(self.pnl, label='Manual File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox5.Add(self.truncbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox5, 0, wx.EXPAND)

        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.obox = wx.TextCtrl(self.pnl, value=self.defaultofile)
        self.hbox6.Add(wx.StaticText(self.pnl, label='Oligomer File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox6.Add(self.obox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox6, 0, wx.EXPAND)

        self.hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        self.matchbox = wx.TextCtrl(self.pnl, value=self.defaultmatchfile)
        self.hbox7.Add(wx.StaticText(self.pnl, label='Match File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox7.Add(self.matchbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox7, 0, wx.EXPAND)

        self.hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.peakbox = wx.TextCtrl(self.pnl, value=self.defaultpeaksfile)
        self.hbox8.Add(wx.StaticText(self.pnl, label='Peak File Name: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox8.Add(self.peakbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.sbs.Add(self.hbox8, 0, wx.EXPAND)

        self.pnl.SetSizer(self.sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

        self.CenterOnParent()

    def OnClose(self, e):
        self.inboxval = self.inputbox.GetValue()
        self.outboxval = self.outputbox.GetValue()
        self.confboxval = self.confbox.GetValue()
        self.truncval = self.truncbox.GetValue()
        self.massfileval = self.massbox.GetValue()
        self.ofileval = self.obox.GetValue()
        self.matchfileval = self.matchbox.GetValue()
        self.peaksfileval = self.peakbox.GetValue()
        if self.inboxval != "":
            self.config.infname = self.inboxval
            print self.inboxval
        if self.outboxval != "":
            self.config.outfname = self.outboxval
            print self.outboxval
        if self.confboxval != "":
            self.config.confname = self.confboxval
            print self.confboxval
        if self.massfileval != "":
            self.config.mfile = self.massfileval
            print self.massfileval
        if self.truncval != "":
            self.config.manualfile = self.truncval
            print self.truncval
        if self.peaksfileval != "":
            self.config.peaksfile = self.peaksfileval
            print self.peaksfileval
        if self.matchfileval != "":
            self.config.matchfile = self.matchfileval
            print self.matchfileval
        if self.ofileval != "":
            self.config.ofile = self.ofileval
            print self.ofileval
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

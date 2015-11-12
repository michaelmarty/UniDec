import os
import wx
from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


def setup_tab_box(tab, plot):
    """
    Shortcut function for adding an object (plot) to a box sizer and fitting it to a notebook tab.
    :param tab: wx.Notebook tab
    :param plot: Object to add to box sizer
    :return: None
    """
    box1 = wx.BoxSizer(wx.VERTICAL)
    box1.Add(plot, 0, wx.EXPAND)
    tab.SetSizerAndFit(box1)


class SingleInputDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Simple dialog for a single textctrl input.
        :param args: Passed to wx.Dialog.
        :param kwargs: Passed to wx.Dialog.
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        # self.SetSize((400, 200))
        self.value = None
        self.inputbox = None

    def initialize_interface(self, title="", message="", defaultvalue=""):
        """
        Create a simple dialog window interface.
        :param title: Title for the frame.
        :param message: Message inside the window.
        :param defaultvalue: String default value for the window.
        :return: None
        """
        self.SetTitle(title)
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox = wx.TextCtrl(pnl, value=defaultvalue)
        hbox.Add(wx.StaticText(pnl, label=message), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox.Add(self.inputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        pnl.SetSizer(hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        hboxend.Add(okbutton)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        self.SetSizer(vbox)

        okbutton.Bind(wx.EVT_BUTTON, self.on_close)

    def on_close(self, e):
        """
        Close the window and set self.value to the value returned by the inputbox.
        :param e: Unused event
        :return: self.value
        """
        self.value = self.inputbox.GetValue()
        self.Destroy()
        self.EndModal(0)
        return self.value


class AdditionalParameters(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create a dialog for setting some obscure additional parameters.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((400, 200))
        self.SetTitle("Set Experimental Parameters")
        self.config = None
        self.inputbox5 = None
        self.inputbox6 = None
        self.inputbox6b = None

    def initialize_interface(self, config):
        """
        Initialize the GUI.
        :param config: UniDecConfig object
        :return: None
        """
        self.config = config

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Set Experimental Parameters')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox5 = wx.TextCtrl(pnl, value=str(self.config.inflate))
        hbox5.Add(wx.StaticText(pnl, label='Peak Inflation Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox5.Add(self.inputbox5, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox5.Add(wx.StaticText(pnl, label=' x Peak Shape '), 0, wx.ALIGN_CENTER_VERTICAL)
        sbs.Add(hbox5, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox6 = wx.TextCtrl(pnl, value=str(self.config.damp))
        hbox6.Add(wx.StaticText(pnl, label='Damping Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox6.Add(self.inputbox6, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox6.Add(wx.StaticText(pnl, label=' '), 0, wx.ALIGN_CENTER_VERTICAL)
        sbs.Add(hbox6, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox6b = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox6b = wx.ComboBox(pnl, wx.ID_ANY, style=wx.CB_READONLY,
                                      choices=["None", "Harmonics", "Satellites", "Both", "Total", "Ultimate",
                                               "Isotope"])
        self.inputbox6b.SetSelection(self.config.suppression)
        hbox6b.Add(wx.StaticText(pnl, label='Suppression Parameter: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox6b.Add(self.inputbox6b, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox6b.Add(wx.StaticText(pnl, label=' Be careful...'), 0, wx.ALIGN_CENTER_VERTICAL)
        sbs.Add(hbox6b, 1, wx.ALIGN_CENTER_VERTICAL)

        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)

        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)

    def on_close(self, e):
        """
        Update parameters in self.config from the GUI. Close the window.
        :param e: Unused Event
        :return: None
        """
        self.config.suppression = ud.string_to_int(self.inputbox6b.GetSelection())
        self.config.inflate = ud.string_to_value(self.inputbox5.GetValue())
        self.config.damp = ud.string_to_value(self.inputbox6.GetValue())
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but don't update the parameters in self.config.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)


class SaveFigureDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Window for setting save figure parameters.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((500, 300))
        self.SetTitle("Save Figures")
        self.config = None
        self.headerbox = None
        self.dirinput = None
        self.extbox = None
        self.widebox = None
        self.tallbox = None
        self.tbox = None
        self.b1 = None
        self.b2 = None
        self.b3 = None
        self.b4 = None
        # Key parameters
        self.directory = None
        self.header = None
        self.extension = None
        self.transparent = None
        self.rect = None
        self.figsize = None

    def initialize_interface(self, config):
        """
        Initialize the GUI.
        :param config: UniDecConfig
        :return: 
        """
        self.config = config

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Save Settings')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(pnl, label="Directory:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.dirinput = wx.TextCtrl(pnl, value=os.getcwd(), size=(300, 20))
        hbox.Add(self.dirinput, 0, wx.EXPAND)
        dirbutton = wx.Button(pnl, label="...", size=(20, 20))
        hbox.Add(dirbutton, 1, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_choose_dir, dirbutton)
        sbs.Add(hbox, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.headerbox = wx.TextCtrl(pnl, value=str(self.config.outfname), size=(300, 20))
        hbox5.Add(wx.StaticText(pnl, label='File Name Header: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox5.Add(self.headerbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox5, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.extbox = wx.TextCtrl(pnl, value=str("png"))
        hbox6.Add(wx.StaticText(pnl, label='Extension: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox6.Add(self.extbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.tbox = wx.CheckBox(pnl, label=str("Transparent"))
        hbox6.Add(self.tbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox6, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox6b = wx.BoxSizer(wx.HORIZONTAL)
        self.widebox = wx.TextCtrl(pnl, value=str(6), size=(40, 20))
        self.tallbox = wx.TextCtrl(pnl, value=str(5), size=(40, 20))
        hbox6b.Add(wx.StaticText(pnl, label='Figure Size: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox6b.Add(self.widebox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox6b.Add(wx.StaticText(pnl, label=' wide by '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox6b.Add(self.tallbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox6b.Add(wx.StaticText(pnl, label=' tall (inches) '), 0, wx.ALIGN_CENTER_VERTICAL)
        sbs.Add(hbox6b, 1, wx.ALIGN_CENTER_VERTICAL)

        hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        size = (40, 20)
        self.b1 = wx.TextCtrl(pnl, value=str(0.1), size=size)
        self.b2 = wx.TextCtrl(pnl, value=str(0.1), size=size)
        self.b3 = wx.TextCtrl(pnl, value=str(0.8), size=size)
        self.b4 = wx.TextCtrl(pnl, value=str(0.8), size=size)
        hbox8.Add(wx.StaticText(pnl, label='Plot Position: '), 0, wx.ALIGN_CENTER_VERTICAL)
        hbox8.Add(self.b1, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox8.Add(self.b2, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox8.Add(self.b3, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox8.Add(self.b4, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox8.Add(wx.StaticText(pnl, label=' Left Bottom Width Height (%) '), 0, wx.ALIGN_CENTER_VERTICAL)
        sbs.Add(hbox8, 1, wx.ALIGN_CENTER_VERTICAL)

        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1,
                 flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend,
                 flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)

        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)

    def on_close(self, e):
        """
        Close the window and save the parameters to the class.
        :param e: Unused event
        :return: None
        """
        # TODO: Put all of this into self.config
        self.directory = self.dirinput.GetValue()
        self.header = self.headerbox.GetValue()
        self.extension = self.extbox.GetValue()
        self.transparent = self.tbox.GetValue()
        self.rect = [float(self.b1.GetValue()), float(self.b2.GetValue()), float(self.b3.GetValue()),
                     float(self.b4.GetValue())]
        self.figsize = [float(self.widebox.GetValue()), float(self.tallbox.GetValue())]
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Destroy the window without setting anything.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def on_choose_dir(self, e):
        """
        Open a dialog to choose the directory.
        :param e: Unused event
        :return: None
        """
        dlg = wx.DirDialog(None, "Choose Top Directory", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            directory = dlg.GetPath()
            self.dirinput.SetValue(directory)
            # print self.directory
        dlg.Destroy()


class FileNameDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Window for renaming default file names.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((250, 300))
        self.SetTitle("Change File Names")
        self.config = None
        self.inputbox = None
        self.confbox = None
        self.outputbox = None
        self.massbox = None
        self.manualbox = None
        self.obox = None
        self.matchbox = None
        self.peakbox = None

    def initialize_interface(self, config):
        """
        Initialize the interface.
        :param config: UniDecConfig object
        :return: None
        """
        self.config = config

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Rename Files Into and Out of UniDec')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox = wx.TextCtrl(pnl, value=self.config.infname)
        hbox1.Add(wx.StaticText(pnl, label='Input File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(self.inputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox1, 0, wx.EXPAND)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.confbox = wx.TextCtrl(pnl, value=self.config.confname)
        hbox2.Add(wx.StaticText(pnl, label='Configuration File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox2.Add(self.confbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox2, 0, wx.EXPAND)

        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.outputbox = wx.TextCtrl(pnl, value=self.config.outfname)
        hbox3.Add(wx.StaticText(pnl, label='Output File Headers: '), wx.ALIGN_CENTER_VERTICAL)
        hbox3.Add(self.outputbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox3, 0, wx.EXPAND)

        hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.massbox = wx.TextCtrl(pnl, value=self.config.mfile)
        hbox4.Add(wx.StaticText(pnl, label='Mass List File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox4.Add(self.massbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox4, 0, wx.EXPAND)

        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.manualbox = wx.TextCtrl(pnl, value=self.config.manualfile)
        hbox5.Add(wx.StaticText(pnl, label='Manual File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox5.Add(self.manualbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox5, 0, wx.EXPAND)

        hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.obox = wx.TextCtrl(pnl, value=self.config.ofile)
        hbox6.Add(wx.StaticText(pnl, label='Oligomer File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox6.Add(self.obox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox6, 0, wx.EXPAND)

        hbox6b = wx.BoxSizer(wx.HORIZONTAL)
        self.matchbox = wx.TextCtrl(pnl, value=self.config.matchfile)
        hbox6b.Add(wx.StaticText(pnl, label='Match File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox6b.Add(self.matchbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox6b, 0, wx.EXPAND)

        hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.peakbox = wx.TextCtrl(pnl, value=self.config.peaksfile)
        hbox8.Add(wx.StaticText(pnl, label='Peak File Name: '), wx.ALIGN_CENTER_VERTICAL)
        hbox8.Add(self.peakbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        sbs.Add(hbox8, 0, wx.EXPAND)

        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)

        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)

        self.CenterOnParent()

    def on_close(self, e):
        """
        Get the new filenames from the window.
        If each is not blank, update self.config.
        Close the window.
        :param e: Unused event
        :return: None
        """
        inboxval = self.inputbox.GetValue()
        outboxval = self.outputbox.GetValue()
        confboxval = self.confbox.GetValue()
        manualval = self.manualbox.GetValue()
        massfileval = self.massbox.GetValue()
        ofileval = self.obox.GetValue()
        matchfileval = self.matchbox.GetValue()
        peaksfileval = self.peakbox.GetValue()
        if inboxval != "":
            self.config.infname = inboxval
            print inboxval
        if outboxval != "":
            self.config.outfname = outboxval
            print outboxval
        if confboxval != "":
            self.config.confname = confboxval
            print confboxval
        if massfileval != "":
            self.config.mfile = massfileval
            print massfileval
        if manualval != "":
            self.config.manualfile = manualval
            print manualval
        if peaksfileval != "":
            self.config.peaksfile = peaksfileval
            print peaksfileval
        if matchfileval != "":
            self.config.matchfile = matchfileval
            print matchfileval
        if ofileval != "":
            self.config.ofile = ofileval
            print ofileval
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but don't update any of the parameters.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

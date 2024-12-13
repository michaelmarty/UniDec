import wx
import wx.lib.mixins.listctrl as listmix
from modules.plotting import PlottingWindow
from unidec.modules.CDEng import UniDecCD
from unidec.modules.isolated_packages import FileDialogs
import numpy as np

__author__ = 'Michael.Marty'

"""
This file contains a windows for fitting CD-MS data to determine the intensity/charge or S/N/charge slope.
"""


class CDListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_val, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create listctrl with 4 columns.
        :param parent: Passed to wx.ListCtrl
        :param id_val: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_val, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "File Name")
        self.InsertColumn(1, "Mass (Da)")
        self.InsertColumn(2, "Min Charge")
        self.InsertColumn(3, "Max Charge")
        self.SetColumnWidth(0, 700)
        self.SetColumnWidth(1, 100)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 100)

        self.popupID1 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_delete, id=self.popupID1)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self)

    def clear_list(self):
        """
        Clear list.
        :return: None
        """
        self.DeleteAllItems()

    def add_line(self, val="File Name"):
        """
        Add a new line to the list.
        :param val: Value for the first column. Default is 0. Default for second column is 0.
        :return: None
        """
        index = self.InsertItem(10000, str(val))
        #self.SetItem(index, 1, str(0))

    def populate(self, data, colors=None):
        """
        Add data from array or nested list to the listctrl.
        :param data: List or array of data values.
        :param colors: Background colors list
        :return: None
        """
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertItem(i, str(data[i][0]))
            self.SetItem(index, 1, str(data[i][1]))
            self.SetItem(index, 2, str(data[i][2]))
            self.SetItem(index, 3, str(data[i][3]))

    def get_list(self):
        """
        Return the list of values in the listctrl.
        :return: Nested list of listctrl output
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [str(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       float(self.GetItem(i, col=2).GetText()), float(self.GetItem(i, col=3).GetText())]
            list_out.append(sublist)
        return list_out

    def on_right_click(self, event):
        """
        Create a right click menu.
        :param event: Unused event
        :return: None
        """
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_delete(self, event):
        """
        Delete the selected item.
        :param event: Unused event.
        :return: None
        """
        item = self.GetFirstSelected()
        num = self.GetSelectedItemCount()
        selection = [item]
        for i in range(1, num):
            item = self.GetNextSelected(item)
            selection.append(item)
        for i in range(0, num):
            self.DeleteItem(selection[num - i - 1])


class CDListCtrlPanel(wx.Panel):
    def __init__(self, parent, size=(1000, 300)):
        """
        Creates the panel for the IMListCtrl
        :param parent: Parent window or panel
        :param size: Size in pixels foor list control
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        tid = wx.NewIdRef()
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = CDListCtrl(self, tid, size=size, style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)


class CDCalDialog(wx.Frame):
    def __init__(self, parent, *args, **kwargs):
        """
        Initialize the dialog window.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        # super(IMTools, self).__init__(*args, **kwargs)
        wx.Frame.__init__(self, parent, *args, **kwargs) #style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER,
        self.SetSize((1225, 850))
        self.SetTitle("CDMS Calibration Tool")
        self.plot = None
        self.masspanel = None
        self.config = None

    def initialize_interface(self, config):
        """
        Initialize the parameters, set up the panels, and display the interface.
        :param config: UniDecConfig object
        :return: None
        """
        self.config = config
        self.setup_panel()
        self.CenterOnParent()

    def setup_panel(self):
        """
        Make/remake the main panel.
        :return: None
        """
        self.pnl = wx.Panel(self)
        self.pnl2 = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        sb = wx.StaticBox(self.pnl, label='CDMS Calibration Plot')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        plotsizer = wx.BoxSizer(wx.HORIZONTAL)
        figsize = (6, 4)
        self.plot = PlottingWindow.Plot2d(self.pnl, figsize=figsize)
        self.plot2 = PlottingWindow.Plot2d(self.pnl, figsize=figsize)
        plotsizer.Add(self.plot, 1, wx.EXPAND)
        plotsizer.Add(self.plot2, 1, wx.EXPAND)

        sbs.Add(plotsizer, 1, wx.EXPAND)

        ctlsizer = wx.BoxSizer(wx.HORIZONTAL)


        vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.masspanel = CDListCtrlPanel(self.pnl)
        csvbutton = wx.Button(self.pnl, label="Load from CSV")
        self.Bind(wx.EVT_BUTTON, self.on_load_csv, csvbutton)
        vbox2.Add(csvbutton, 0, wx.EXPAND)

        addbutton = wx.Button(self.pnl, label="Add Line")
        self.Bind(wx.EVT_BUTTON, self.on_add, addbutton)
        vbox2.Add(addbutton, 0, wx.EXPAND)

        vbox2.Add(self.masspanel, 0, wx.EXPAND)

        plotbutton = wx.Button(self.pnl, label="Plot Species")
        self.Bind(wx.EVT_BUTTON, self.on_plot, plotbutton)
        vbox2.Add(plotbutton, 0, wx.EXPAND)

        ctlsizer.Add(vbox2, 1, wx.EXPAND)



        sb2 = wx.StaticBox(self.pnl, label='Fit Parameters')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)
        gbox1c = wx.GridBagSizer(wx.VERTICAL)
        size1 = (75, -1)

        # Readout
        self.ctlfit = wx.TextCtrl(self.pnl, value='', size=size1)
        gbox1c.Add(self.ctlfit, (0, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(self.pnl, label="Intensity Slope: "), (0, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlfit2 = wx.TextCtrl(self.pnl, value='', size=size1)
        gbox1c.Add(self.ctlfit2, (1, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(self.pnl, label="S/N Slope: "), (1, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)

        self.ctlmzwindow = wx.TextCtrl(self.pnl, value='', size=size1)
        gbox1c.Add(self.ctlmzwindow, (2, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(self.pnl, label="m/z Window: "), (2, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlmzwindow.SetValue("25")

        self.ctlnoise = wx.TextCtrl(self.pnl, value='', size=size1)
        gbox1c.Add(self.ctlnoise, (3, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(self.pnl, label="S/N Minimum: "), (3, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlnoise.SetValue("0")

        self.ctlsdmult = wx.TextCtrl(self.pnl, value='', size=size1)
        gbox1c.Add(self.ctlsdmult, (4, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(self.pnl, label="Std Dev Window: "), (4, 0),
                   flag=wx.ALIGN_CENTER_VERTICAL)
        self.ctlsdmult.SetValue("2")

        sbs2.Add(gbox1c, 0, wx.EXPAND)

        ctlsizer.Add(sbs2, 0, wx.EXPAND)

        sbs.Add(ctlsizer, 0, wx.EXPAND)
        vbox.Add(sbs, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)

        self.pnl.SetSizer(vbox)
        vbox.Fit(self.pnl)

        self.Centre()
        self.Show(True)
        self.Raise()

    def on_add(self, e):
        """
        Add lines to the masspanel.
        :param e: Unused event
        :return: None
        """
        self.masspanel.list.add_line()

    def on_load_csv(self, e=None):
        filename = FileDialogs.open_file_dialog("Open Calibration CSV File", file_types="*.csv*")
        print("Opening CSV File: ", filename)
        if filename is not None:
            self.load_csv(filename)

    def load_csv(self, path):
        caldat = np.genfromtxt(path, delimiter=",", dtype=np.str)[1:]
        self.masspanel.list.populate(caldat)

    def on_plot(self, e=None):
        """
        Gathers all of the parameters, fits, and plots the data.
        :param e: wx.Event
        :return: None
        """
        caldat = np.array(self.masspanel.list.get_list())
        try:
            mzwindow = float(self.ctlmzwindow.GetValue())
            print("Using MZ Window: ", mzwindow)
        except:
            mzwindow = 25
            print("Could Not Convert to Float:", self.ctlmzwindow.GetValue())
            print("Using Default for MZ window:", mzwindow)

        try:
            noisemult = float(self.ctlnoise.GetValue())
            print("Using Noise Threshold: ", noisemult)
        except:
            noisemult = 0
            print("Could Not Convert to Float:", self.ctlnoise.GetValue())
            print("Using Default for Noise Threshold:", noisemult)

        try:
            sdmult = float(self.ctlsdmult.GetValue())
            print("Using Std Dev Window: ", sdmult)
        except:
            sdmult = 0
            print("Could Not Convert to Float:", self.ctlsdmult.GetValue())
            print("Using Default for Std Dev Window:", sdmult)

        eng = UniDecCD()
        extracts = []
        snextracts = []
        for i, l in enumerate(caldat):
            f = l[0]
            eng.open_file(f)
            m = float(l[1])
            minz = float(l[2])
            maxz = float(l[3])
            ext, snext = eng.extract_intensities(m, minz, maxz, window=mzwindow, noise_mult=noisemult, sdmult=sdmult)

            if i == 0:
                extracts = ext
            else:
                extracts = np.concatenate((extracts, ext), axis=0)

            if i == 0:
                snextracts = snext
            else:
                snextracts = np.concatenate((snextracts, snext), axis=0)

            print("Length", len(ext))

        fits, fitdat = eng.get_fit(extracts)
        self.ctlfit.SetValue(str(fits[0]))
        snfits, snfitdat = eng.get_fit(snextracts)
        self.ctlfit2.SetValue(str(snfits[0]))
        self.plot.clear_plot()
        bins = [int(np.amax(fitdat[:, 0]) - np.amin(fitdat[:, 0])), 50]
        self.plot.hist2d(extracts[:, 0], extracts[:, 1], bins=bins, config=self.config, xlab="Known Charge",
                         ylab="Measured Intensity", title="Intensity vs. Charge")
        self.plot.subplot1.plot(fitdat[:, 0], fitdat[:, 1], color="red")
        self.plot.repaint()


        self.plot2.clear_plot()
        self.plot2.hist2d(snextracts[:, 0], snextracts[:, 1], bins=bins, config=self.config, xlab="Known Charge",
                         ylab="Measured S/N Ratio", title="S/N vs. Charge")
        self.plot2.subplot1.plot(snfitdat[:, 0], snfitdat[:, 1], color="red")
        self.plot2.repaint()


# Main App Execution
if __name__ == "__main__":
    #calpath = "Z:\Group Share\Marius Kostelic\CD-MS\Replicates\Calibration.csv"
    calpath = "C:\\Data\\CDMS\\AqpZ_STORI\\AqpZ_STORI\\caltest.csv"
    app = wx.App(False)
    frame = CDCalDialog(None)
    from unidec.modules.unidecstructure import UniDecConfig
    config = UniDecConfig()
    frame.initialize_interface(config)
    frame.load_csv(calpath)
    frame.on_plot(None)
    app.MainLoop()
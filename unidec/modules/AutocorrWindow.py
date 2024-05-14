import wx.lib.mixins.listctrl as listmix
import wx

from unidec.modules import PlottingWindow, peakstructure
import unidec.tools as ud

__author__ = 'Michael.Marty'

'''
Window for viewing autocorrleation results
'''


class CorrListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    """
    Class for the list control of peak values
    """

    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a two column list control
        :param parent: Parent window or panel
        :param id_value: ID passed to wx.ListCtrl
        :param pos: pos passed to wx.ListCtrl
        :param size: size passed to wx.ListCtrl
        :param style: styl passed to wx.ListCtrl
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(1, "Mass Diff (Da)")
        self.SetColumnWidth(1, width=100)  # , wx.LIST_AUTOSIZE)
        self.InsertColumn(0, "")
        self.SetColumnWidth(0, width=25)  # , wx.LIST_AUTOSIZE)

    def populate(self, pks):
        """
        Add Peaks class object to the panel. Adds p.textmarker to the first column and p.mass to the second.
        The background color is set to p.color.
        :param pks: Peaks object
        :return: None
        """
        self.DeleteAllItems()
        for i in range(0, pks.plen):
            p = pks.peaks[i]
            self.InsertItem(i, p.textmarker)
            self.SetItem(i, 1, str(p.mass))
            color = wx.Colour(int(round(p.color[0] * 255)), int(round(p.color[1] * 255)), int(round(p.color[2] * 255)),
                              alpha=255)
            self.SetItemBackgroundColour(i, col=color)

    def clear_list(self):
        """
        Clears the list.
        :return: None
        """
        self.DeleteAllItems()


class CorrListCtrlPanel(wx.Panel):
    """
    Panel for the ListCtrl
    """

    def __init__(self, parent):
        """
        Creates the panel.
        :param parent: Parent panel or window
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = CorrListCtrl(self, wx.NewIdRef(), size=(200, 550), style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class AutocorrWindow(wx.Dialog):
    """
    Dialog window for autocorrelation.
    """

    def __init__(self, *args, **kwargs):
        """
        Creates a dialog.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((800, 600))
        self.SetTitle("Autocorrelation Plots")
        self.config = None
        self.massdat = []
        self.plot1 = None
        self.listpanel = None
        self.window=None

    def initalize_dialog(self, config, massdat, window=None):
        """
        Initilizes dialog and the layout. Calls run_autocorr.
        :param config: UniDecConfig object
        :param massdat: Array of mass distribution (N x 2)
        :return: None
        """
        self.config = config
        self.massdat = massdat
        self.window=window

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        sb = wx.StaticBox(panel, label='Autocorrelation')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot1 = PlottingWindow.Plot1d(panel)
        sbs.Add(self.plot1)
        hbox.Add(sbs)
        self.listpanel = CorrListCtrlPanel(panel)
        hbox.Add(self.listpanel)
        panel.SetSizer(hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        # closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        # hboxend.Add(closebutton, flag=wx.LEFT, border=5)
        vbox.Add(panel, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        self.SetSizer(vbox)
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        # closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        self.CenterOnParent()
        self.run_autocorr(0)

    def run_autocorr(self, e):
        """
        Performs autocorrelation using ud.autocorr from unidectools.
        Plots the results and adds peaks to listctrl.
        :param e: Unused event
        :return: None
        """
        if self.window is not None:
            corr, peaks = ud.autocorr(self.massdat, config=None, window=self.window)
        else:
            corr, peaks = ud.autocorr(self.massdat, self.config)
        self.plot1.plotrefreshtop(corr[:, 0], corr[:, 1], "Autocorrelation", "Difference", "", "", self.config)
        pks2 = peakstructure.Peaks()
        pks2.add_peaks(peaks, self.config.massbins)
        pks2.default_params()
        if pks2.plen > 0:
            for p in pks2.peaks:
                if p.ignore == 0:
                    self.plot1.plotadddot(p.mass, p.height, p.color, p.marker)
            self.plot1.repaint()
        print(peaks)
        self.listpanel.list.populate(pks2)

    def on_close(self, e):
        """
        Close dialog.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Cancel the dialog.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(0)

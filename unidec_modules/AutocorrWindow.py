import wx.lib.mixins.listctrl as listmix
import wx

from unidec_modules import plot1d, peakstructure
import unidec_modules.unidectools as ud

__author__ = 'Michael.Marty'


class CorrListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(1, "Mass Diff (Da)")
        self.SetColumnWidth(1, width=100)  # , wx.LIST_AUTOSIZE)
        self.InsertColumn(0, "")
        self.SetColumnWidth(0, width=25)  # , wx.LIST_AUTOSIZE)

    def populate(self, pks):
        self.DeleteAllItems()
        for i in range(0, pks.plen):
            p = pks.peaks[i]
            self.InsertStringItem(i, p.textmarker)
            self.SetStringItem(i, 1, str(p.mass))
            color = wx.Colour(round(p.color[0] * 255), round(p.color[1] * 255), round(p.color[2] * 255), alpha=255)
            self.SetItemBackgroundColour(i, col=color)

    def clear_list(self):
        self.DeleteAllItems()


class CorrListCtrlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = CorrListCtrl(self, wx.NewId(), size=(200, 550), style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class AutocorrWindow(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((800, 600))
        self.SetTitle("Autocorrelation Plots")
        self.config = None
        self.massdat = []
        self.plot1 = None
        self.listpanel = None
        self.corr = []
        self.peaks = []
        self.pks2 = peakstructure.Peaks()

    def initalize_dialog(self, config, massdat):
        self.config = config
        self.massdat = massdat

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        sb = wx.StaticBox(panel, label='Autocorrelation')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot1 = plot1d.Plot1d(panel)
        sbs.Add(self.plot1)
        hbox.Add(sbs)
        self.listpanel = CorrListCtrlPanel(panel)
        hbox.Add(self.listpanel)
        panel.SetSizer(hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)
        vbox.Add(panel, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        self.SetSizer(vbox)
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        self.CenterOnParent()
        self.run_autocorr(0)

    def run_autocorr(self, e):
        self.corr, self.peaks = ud.autocorr(self.massdat, self.config)
        self.plot1.plotrefreshtop(self.corr[:, 0], self.corr[:, 1], "Autocorrelation", "Mass Difference", "", "",
                                  self.config)
        self.pks2 = peakstructure.Peaks()
        self.pks2.add_peaks(self.peaks)
        self.pks2.default_params()
        if self.pks2.plen > 0:
            for p in self.pks2.peaks:
                if p.ignore == 0:
                    self.plot1.plotadddot(p.mass, p.height, p.color, p.marker)
            self.plot1.repaint()
        print self.peaks
        self.listpanel.list.populate(self.pks2)

    def on_close(self, e):
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        self.Destroy()
        self.EndModal(0)

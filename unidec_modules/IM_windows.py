__author__ = 'Michael.Marty'

import sys
import os
from copy import deepcopy

import wx
import wx.lib.mixins.listctrl  as  listmix
import matplotlib.cm as cm

from unidec_modules.IM_functions import *
from unidec_modules import plot1d, plot2d


class TestListCtrl4(wx.ListCtrl,
                    listmix.ListCtrlAutoWidthMixin,
                    listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Mass (Da)")
        self.InsertColumn(1, u"CCS (\u212B\u00B2)")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 100)

    def Clear(self):
        self.DeleteAllItems()

    def AddLine(self, val=0):
        index = self.InsertStringItem(sys.maxint, str(val))
        self.SetStringItem(index, 1, str(0))

    def Populate(self, data, colors=None):
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertStringItem(sys.maxint, str(data[i][0]))
            self.SetStringItem(index, 1, str(data[i][1]))
            if colors is not None:
                color = colors[i]
                col = wx.Colour(round(color[0] * 255), round(color[1] * 255), round(color[2] * 255), alpha=255)
                self.SetItemBackgroundColour(index, col=col)

    def GetList(self):
        count = self.GetItemCount()
        list = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText())]
            list.append(sublist)
        return list


class TestListCtrlPanel4(wx.Panel):
    def __init__(self, parent, size=(200, 200)):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        tID = wx.NewId()
        sizer = wx.BoxSizer(wx.VERTICAL)
        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)
        self.list = TestListCtrl4(self, tID, size=size, style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")


class IMTools(wx.Dialog):
    def __init__(self, *args, **kwargs):
        super(IMTools, self).__init__(*args, **kwargs)
        self.SetSize((600, 800))
        self.SetTitle("Ion Mobility Tools")

    def InitUI(self, data3, config):
        self.defaultconfig = config
        self.config = deepcopy(config)
        self.data3 = data3
        self.pnl = wx.Panel(self)
        self.pnl2 = wx.Panel(self)
        self.SetupPanel()

        self.CenterOnParent()
        self.flag = 0

    def SetupPanel(self):
        for child in self.pnl.GetChildren():
            child.Destroy()
        for child in self.pnl2.GetChildren():
            child.Destroy()
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.sb = wx.StaticBox(self.pnl, label='Ion Mobility Parameters Tool')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.plot = plot2d.Plot2d(self.pnl)
        self.plot.contourplot(self.data3, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="IM-MS Data")
        self.sbs.Add(self.plot, 1, wx.EXPAND)

        self.ctlsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sb2 = wx.StaticBox(self.pnl, label='Instrumental Parameters')
        self.sbs2 = wx.StaticBoxSizer(self.sb2, orient=wx.VERTICAL)
        self.gbox1c = wx.GridBagSizer(wx.VERTICAL)
        size1 = (75, -1)
        self.ctltwave = wx.RadioBox(self.pnl, label="", choices=["Linear Cell", "Travelling Wave"])
        self.Bind(wx.EVT_RADIOBOX, self.OnFlipTWave, self.ctltwave)
        self.gbox1c.Add(self.ctltwave, (0, 0), span=(1, 5))

        if self.config.twaveflag == 0:
            self.ctlvolt = wx.TextCtrl(self.pnl, value="", size=size1)
            self.gbox1c.Add(self.ctlvolt, (1, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Voltage (V): "), (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlpressure = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctlpressure, (2, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Pressure (Torr): "), (2, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltemp = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctltemp, (3, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label=u"Temperature (\u00B0C): "), (3, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlgasmass = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlto = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctlto, (5, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label=u"Dead Time (t\u2080 in ms): "), (5, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctldriftlength = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctldriftlength, (6, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Drift Cell Length (m)"), (6, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)

            self.twave = 0

        else:

            self.ctltcal1 = wx.TextCtrl(self.pnl, value="", size=size1)
            self.gbox1c.Add(self.ctltcal1, (1, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Calibration Parameter 1: "), (1, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltcal2 = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctltcal2, (2, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Calibration Parameter 2: "), (2, 0),
                            flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctledc = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctledc, (3, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="EDC Parameter: "), (3, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlgasmass = wx.TextCtrl(self.pnl, value='', size=size1)
            self.gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
            self.gbox1c.Add(wx.StaticText(self.pnl, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.twave = 1

        self.sbs2.Add(self.gbox1c, 0, wx.EXPAND)
        self.ctlsizer.Add(self.sbs2, 0, wx.EXPAND)

        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.masspanel = TestListCtrlPanel4(self.pnl)
        self.addbutton = wx.Button(self.pnl, label="Add Species")
        self.Bind(wx.EVT_BUTTON, self.OnAdd, self.addbutton)
        self.vbox2.Add(self.addbutton, 0, wx.EXPAND)
        self.vbox2.Add(self.masspanel, 0, wx.EXPAND)
        self.plotbutton = wx.Button(self.pnl, label="Plot Species")
        self.Bind(wx.EVT_BUTTON, self.OnPlot, self.plotbutton)
        self.vbox2.Add(self.plotbutton, 0, wx.EXPAND)

        self.ctlsizer.Add(self.vbox2, -1, wx.EXPAND)

        self.sbs.Add(self.ctlsizer, 0, wx.EXPAND)
        self.pnl.SetSizer(self.sbs)

        self.hboxend = wx.BoxSizer(wx.HORIZONTAL)
        self.okButton = wx.Button(self.pnl2, label='Ok')
        self.closeButton = wx.Button(self.pnl2, label='Cancel')
        self.okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        self.closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)
        self.hboxend.Add(self.okButton)
        self.hboxend.Add(self.closeButton, flag=wx.LEFT, border=5)
        self.pnl2.SetSizer(self.hboxend)
        # self.hboxend.Fit(self.pnl2)

        self.vbox.Add(self.pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(self.pnl2, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        self.loadtogui(0)

    def OnClose(self, e):
        self.getfromgui(e)
        self.defaultconfig = self.config
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def loadtogui(self, e):
        self.ctltwave.SetSelection(self.config.twaveflag)
        if self.config.twaveflag == 0:
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

    def getfromgui(self, e):

        if self.config.twaveflag == 0:
            self.config.volt = ud.string_to_value(self.ctlvolt.GetValue())
            self.config.temp = ud.string_to_value(self.ctltemp.GetValue())
            self.config.pressure = ud.string_to_value(self.ctlpressure.GetValue())
            self.config.gasmass = ud.string_to_value(self.ctlgasmass.GetValue())
            self.config.to = ud.string_to_value(self.ctlto.GetValue())
            self.config.driftlength = ud.string_to_value(self.ctldriftlength.GetValue())
        elif self.config.twaveflag == 1:
            self.config.tcal1 = ud.string_to_value(self.ctltcal1.GetValue())
            self.config.tcal2 = ud.string_to_value(self.ctltcal2.GetValue())
            self.config.edc = ud.string_to_value(self.ctledc.GetValue())
            self.config.gasmass = ud.string_to_value(self.ctlgasmass.GetValue())

    def OnAdd(self, e):
        self.masspanel.list.AddLine()

    def OnPlot(self, e):
        self.getfromgui(e)
        outs = np.array(self.masspanel.list.GetList())
        ztab = np.arange(float(self.config.startz), float(self.config.endz + 1))
        self.plot.clear_plot()
        self.plot.contourplot(self.data3, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="IM-MS Data")
        for i in xrange(0, len(outs)):
            mass = outs[i, 0]
            ccs = outs[i, 1]
            mzvals = (mass + ztab * self.config.adductmass) / ztab
            if self.config.twaveflag == 0:
                dts = np.array([calc_linear_dt(mass, z, ccs, self.config) for z in ztab])
            else:
                dts = np.array([calc_twave_dt(mass, z, ccs, self.config) for z in ztab])
            dtdat = np.unique(self.data3[:, 1])
            maxdt = np.amax(dtdat)
            mindt = np.amin(dtdat)
            clipped = np.clip(dts, mindt, maxdt)
            print np.array([ztab, dts])
            if np.amin(clipped) == maxdt or np.amax(clipped) == mindt:
                dts = clipped
            self.plot.subplot1.plot(mzvals, dts, marker="o")
        self.plot.repaint()

    def OnFlipTWave(self, e):
        self.config.twaveflag = self.ctltwave.GetSelection()
        if self.config.twaveflag == 0:
            self.config.gasmass = 4.002602
            print "Using Linear Cell"
        elif self.config.twaveflag == 1:
            self.config.gasmass = 28.0134
            print "Using Travelling Wave"
        self.SetupPanel()
        self.ctltwave.SetSelection(self.config.twaveflag)


class IMToolExtract(wx.Dialog):
    def __init__(self, *args, **kwargs):
        super(IMToolExtract, self).__init__(*args, **kwargs)
        self.SetSize((800, 800))
        self.SetTitle("Ion Mobility Extraction Tools")

    def InitUI(self, massdat, ccsdat, mccsgrid, config, pks):
        self.config = config
        self.massdat = massdat
        self.ccsdat = ccsdat
        self.totalgrid = mccsgrid
        self.pks = pks
        self.ztab = np.arange(float(self.config.startz), float(self.config.endz + 1))
        self.zstrings = [str(int(z)) for z in self.ztab]
        self.zstrings.append("All")

        self.pnl = wx.Panel(self)

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.sb = wx.StaticBox(self.pnl, label='Ion Mobility Extraction Tool')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.HORIZONTAL)
        figsize = (6, 5)
        self.plot1 = plot2d.Plot2d(self.pnl, figsize=figsize)
        self.plot2 = plot1d.Plot1d(self.pnl, figsize=figsize)
        self.plotsizer = wx.BoxSizer(wx.VERTICAL)
        self.plotsizer.Add(self.plot1, 0, wx.EXPAND)
        self.plotsizer.Add(self.plot2, 0, wx.EXPAND)
        self.sbs.Add(self.plotsizer, 0, wx.EXPAND)

        self.sb2 = wx.StaticBox(self.pnl, label='Extraction Parameters')
        self.sbs2 = wx.StaticBoxSizer(self.sb2, orient=wx.VERTICAL)
        self.gbox1c = wx.GridBagSizer(wx.VERTICAL)
        size1 = (75, -1)

        self.ctlzout = wx.ComboBox(self.pnl, value="", size=size1, choices=self.zstrings)
        self.gbox1c.Add(self.ctlzout, (0, 1), span=(1, 1))
        self.gbox1c.Add(wx.StaticText(self.pnl, label="Charge State: "), (0, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        self.sbs2.Add(self.gbox1c, 0, wx.EXPAND)

        self.masspanel = TestListCtrlPanel4(self.pnl, size=(200, 700))
        self.addbutton = wx.Button(self.pnl, label="Add Species")
        self.plotbutton = wx.Button(self.pnl, label="Plot Species")
        self.sbs2.Add(self.plotbutton, 0, wx.EXPAND)
        self.sbs2.Add(self.masspanel, 0, wx.EXPAND)
        self.sbs2.Add(self.addbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.OnAdd, self.addbutton)
        self.Bind(wx.EVT_BUTTON, self.OnPlot, self.plotbutton)

        self.sbs.Add(self.sbs2, 0, wx.EXPAND)
        self.pnl.SetSizer(self.sbs)

        self.hboxend = wx.BoxSizer(wx.HORIZONTAL)
        self.okButton = wx.Button(self, label='Ok')
        self.closeButton = wx.Button(self, label='Cancel')
        self.okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        self.closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)
        self.hboxend.Add(self.okButton)
        self.hboxend.Add(self.closeButton, flag=wx.LEFT, border=5)
        self.vbox.Add(self.pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(self.hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

        self.CenterOnParent()
        self.loadpeaks(0)
        self.OnPlot(0)
        self.ctlzout.SetSelection(len(self.zstrings) - 1)

    def OnClose(self, e):
        self.config.zout = self.zout
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def loadpeaks(self, e):
        for p in self.pks.peaks:
            self.masspanel.list.AddLine(val=p.mass)

    def getfromgui(self, e):
        try:
            self.zout = int(self.ctlzout.GetStringSelection())
            print "Charge State: ", self.zout
        except ValueError:
            self.zout = 0

    def OnAdd(self, e):
        self.masspanel.list.AddLine()

    def OnPlot(self, e):
        self.getfromgui(e)
        outs = np.array(self.masspanel.list.GetList())
        fname = self.config.outfname + "_zout_" + str(self.zout) + ".bin"
        if os.path.isfile(fname):
            self.zoutgrid = np.fromfile(fname, dtype=float)
        else:
            self.zoutgrid = self.totalgrid

        self.plot1.contourplot(xvals=self.massdat[:, 0], yvals=self.ccsdat[:, 0], zgrid=self.zoutgrid,
                               config=self.config, ylab="CCS (${\AA}$$^2$)", title="Mass vs. CCS", test_kda=True)
        self.zoutgrid = np.reshape(self.zoutgrid, (len(self.massdat), len(self.ccsdat)))
        self.ccsproj = np.sum(self.zoutgrid, axis=0)
        self.plot2.plotrefreshtop(self.ccsdat[:, 0], self.ccsproj / np.amax(self.ccsproj), "CCS Extract",
                                  "CCS (${\AA}$$^2$)", "Normalized Intensity", "", self.config)

        colormap = cm.get_cmap(self.config.peakcmap, len(outs))
        xcolors = colormap(np.arange(len(outs)))
        for i, l in enumerate(outs):
            mass = l[0]
            ccsext = self.zoutgrid[ud.nearest(self.massdat[:, 0], mass)]
            ccsmax = self.ccsdat[np.argmax(ccsext), 0]
            l[1] = ccsmax
            self.plot2.plotadd(self.ccsdat[:, 0], ccsext / np.amax(ccsext), xcolors[i], '')
        self.plot2.repaint()
        self.masspanel.list.Populate(outs, colors=xcolors)

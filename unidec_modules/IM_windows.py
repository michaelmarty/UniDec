import sys
import os
from copy import deepcopy
import wx
import wx.lib.mixins.listctrl as listmix
import matplotlib.cm as cm
from unidec_modules.IM_functions import *
from unidec_modules import plot1d, plot2d

__author__ = 'Michael.Marty'

"""
This file contains two windows for manipulating IM-MS data.
IMTools is used for plotting predicted m/z and arrival time for a given set of parameters.
IMToolExtract is for extracting the CCS distribution for specific slices of the IM data.
"""


class IMListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_val, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create listctrl with 2 columns.
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
        self.InsertColumn(0, "Mass (Da)")
        self.InsertColumn(1, "CCS (\u212B\u00B2)")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 100)

        self.popupID1 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_delete, id=self.popupID1)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self)

    def clear_list(self):
        """
        Clear list.
        :return: None
        """
        self.DeleteAllItems()

    def add_line(self, val=0):
        """
        Add a new line to the list.
        :param val: Value for the first column. Default is 0. Default for second column is 0.
        :return: None
        """
        index = self.InsertItem(10000, str(val))
        self.SetItem(index, 1, str(0))

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
            if colors is not None:
                color = colors[i]
                col = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)), int(round(color[2] * 255)),
                                alpha=255)
                self.SetItemBackgroundColour(index, col=col)

    def get_list(self):
        """
        Return the list of values in the listctrl.
        :return: Nested list of listctrl output
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText())]
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


class IMListCtrlPanel(wx.Panel):
    def __init__(self, parent, size=(200, 200)):
        """
        Creates the panel for the IMListCtrl
        :param parent: Parent window or panel
        :param size: Size in pixels foor list control
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        tid = wx.NewIdRef()
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = IMListCtrl(self, tid, size=size, style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class IMTools(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Initialize the dialog window.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        # super(IMTools, self).__init__(*args, **kwargs)
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((600, 800))
        self.SetTitle("Ion Mobility Tools")
        self.pnl = None
        self.pnl2 = None
        self.plot = None
        self.ctltwave = None
        self.masspanel = None
        self.twave = False
        self.flag = 0
        self.data3 = []
        self.config = None
        self.defaultconfig = None

    def initialize_interface(self, data3, config):
        """
        Initialize the parameters, setup the panels, and display the interface.
        :param data3: IM-MS raw or processed data
        :param config: UniDecConfig object
        :return: None
        """
        self.config = config
        self.data3 = data3

        self.pnl = wx.Panel(self)
        self.pnl2 = wx.Panel(self)
        self.setup_panel()
        self.CenterOnParent()

    def setup_panel(self):
        """
        Make/remake the main panel.
        :return: None
        """
        # TODO: Inherit this from mainwindow somehow so controls are the same.
        for child in self.pnl.GetChildren():
            child.Destroy()
        for child in self.pnl2.GetChildren():
            child.Destroy()
        vbox = wx.BoxSizer(wx.VERTICAL)
        sb = wx.StaticBox(self.pnl, label='Ion Mobility Parameters Tool')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot = plot2d.Plot2d(self.pnl)
        self.plot.contourplot(self.data3, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="IM-MS Data")
        sbs.Add(self.plot, 1, wx.EXPAND)

        ctlsizer = wx.BoxSizer(wx.HORIZONTAL)
        sb2 = wx.StaticBox(self.pnl, label='Instrumental Parameters')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)
        gbox1c = wx.GridBagSizer(wx.VERTICAL)
        size1 = (75, -1)
        self.ctltwave = wx.RadioBox(self.pnl, label="", choices=["Linear Cell", "Travelling Wave"])
        self.Bind(wx.EVT_RADIOBOX, self.on_flip_twave, self.ctltwave)
        gbox1c.Add(self.ctltwave, (0, 0), span=(1, 5))

        self.twave = self.config.twaveflag > 0
        if not self.twave:
            self.ctltwave.SetSelection(0)
            # Linear Mode controls
            self.ctlvolt = wx.TextCtrl(self.pnl, value="", size=size1)
            gbox1c.Add(self.ctlvolt, (1, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Voltage (V): "), (1, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlpressure = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctlpressure, (2, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Pressure (Torr): "), (2, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltemp = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctltemp, (3, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Temperature (\u00B0C): "), (3, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlgasmass = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlto = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctlto, (5, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Dead Time (t\u2080 in ms): "), (5, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctldriftlength = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctldriftlength, (6, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Drift Cell Length (m)"), (6, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

        else:
            self.ctltwave.SetSelection(1)
            # T-Wave Controls
            self.ctltcal1 = wx.TextCtrl(self.pnl, value="", size=size1)
            gbox1c.Add(self.ctltcal1, (1, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Calibration Parameter 1: "), (1, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltcal2 = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctltcal2, (2, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Calibration Parameter 2: "), (2, 0),
                       flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctledc = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctledc, (3, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="EDC Parameter: "), (3, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctlgasmass = wx.TextCtrl(self.pnl, value='', size=size1)
            gbox1c.Add(self.ctlgasmass, (4, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Gas Mass (Da): "), (4, 0), flag=wx.ALIGN_CENTER_VERTICAL)

            self.ctltwavecaltype = wx.Choice(self.pnl, -1, choices=list(self.config.twavedict.values()))
            gbox1c.Add(self.ctltwavecaltype, (5, 1), span=(1, 1))
            gbox1c.Add(wx.StaticText(self.pnl, label="Calibration Type: "), (5, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        sbs2.Add(gbox1c, 0, wx.EXPAND)
        ctlsizer.Add(sbs2, 0, wx.EXPAND)

        vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.masspanel = IMListCtrlPanel(self.pnl)
        addbutton = wx.Button(self.pnl, label="Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add, addbutton)
        vbox2.Add(addbutton, 0, wx.EXPAND)
        vbox2.Add(self.masspanel, 0, wx.EXPAND)
        plotbutton = wx.Button(self.pnl, label="Plot Species")
        self.Bind(wx.EVT_BUTTON, self.on_plot, plotbutton)
        vbox2.Add(plotbutton, 0, wx.EXPAND)

        ctlsizer.Add(vbox2, -1, wx.EXPAND)

        sbs.Add(ctlsizer, 0, wx.EXPAND)
        self.pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self.pnl2, label='Ok')
        closebutton = wx.Button(self.pnl2, label='Cancel')
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)
        # TODO: There is a strange bug whereby the closebutton is not drawn when the window is flipped...
        self.pnl2.SetSizer(hboxend)

        vbox.Add(self.pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(self.pnl2, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        vbox.Fit(self)

        self.load_to_gui(0)

    def on_close(self, e):
        """
        Close the window and set self.defaultconfig to the new self.config.
        :param e: Unused event
        :return: None
        """
        self.get_from_gui(e)
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but do not update self.defaultconfig.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def load_to_gui(self, e):
        """
        Load from self.config into the GUI control boxes.
        :param e: Unused event
        :return: None
        """
        if not self.twave:
            self.ctltwave.SetSelection(0)
            self.ctlvolt.SetValue(str(self.config.volt))
            self.ctltemp.SetValue(str(self.config.temp))
            self.ctlpressure.SetValue(str(self.config.pressure))
            self.ctlgasmass.SetValue(str(self.config.gasmass))
            self.ctlto.SetValue(str(self.config.to))
            self.ctldriftlength.SetValue(str(self.config.driftlength))
        else:
            self.ctltwave.SetSelection(1)
            self.ctltcal1.SetValue(str(self.config.tcal1))
            self.ctltcal2.SetValue(str(self.config.tcal2))
            self.ctledc.SetValue(str(self.config.edc))
            self.ctlgasmass.SetValue(str(self.config.gasmass))
            self.ctltwavecaltype.SetSelection(list(self.config.twavedict.keys()).index(self.config.twaveflag))

    def get_from_gui(self, e):
        """
        Load from GUI to self.config
        :param e: Unused event
        :return: None
        """
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

    def on_add(self, e):
        """
        Add lines to the masspanel.
        :param e: Unused event
        :return: None
        """
        self.masspanel.list.add_line()

    def on_plot(self, e):
        """
        Gathers all of the parameters and simulates the m/z and dt values. Plots the results.
        :param e: wx.Event
        :return: None
        """
        self.get_from_gui(e)
        outs = np.array(self.masspanel.list.get_list())
        ztab = np.arange(float(self.config.startz), float(self.config.endz + 1))
        self.plot.clear_plot()
        self.plot.contourplot(self.data3, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="IM-MS Data")
        for i in range(0, len(outs)):
            mass = outs[i, 0]
            ccs = outs[i, 1]
            mzvals = (mass + ztab * self.config.adductmass) / ztab
            if self.config.twaveflag == 0:
                dts = np.array([calc_linear_dt(mass, z, ccs, self.config) for z in ztab])
            elif self.config.twaveflag == 1:
                dts = np.array([calc_twave_dt_log(mass, z, ccs, self.config) for z in ztab])
            elif self.config.twaveflag == 2:
                dts = np.array([calc_twave_dt_linear(mass, z, ccs, self.config) for z in ztab])
            elif self.config.twaveflag == 3:
                dts = np.array([calc_twave_dt_power(mass, z, ccs, self.config) for z in ztab])
            else:
                print("Error: twaveflat value not supported. Value was:", self.config.twaveflag)
            dtdat = np.unique(self.data3[:, 1])
            maxdt = np.amax(dtdat)
            mindt = np.amin(dtdat)
            clipped = np.clip(dts, mindt, maxdt)
            print(mass, np.array([ztab, dts]))
            if np.amin(clipped) == maxdt or np.amax(clipped) == mindt:
                dts = clipped
            self.plot.subplot1.plot(mzvals, dts, marker="o")
        self.plot.repaint()

    def on_flip_twave(self, e):
        """
        Flip configuration from linear to t-wave. Changes the gas from Ni to He or visa versa. Remakes the GUI.
        :param e: Unused event
        :return: None
        """
        self.config.twaveflag = self.ctltwave.GetSelection()
        if self.config.twaveflag == 0:
            self.config.gasmass = 4.002602
            print("Using Linear Cell")
        elif self.config.twaveflag > 0:
            self.config.gasmass = 28.0134
            print("Using Travelling Wave")
        self.setup_panel()
        self.ctltwave.SetSelection(int(self.twave))


class IMToolExtract(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Creates dialog window for performing extraction of IM data slices.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        # super(IMToolExtract, self).__init__(*args, **kwargs)
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((800, 800))
        self.SetTitle("Ion Mobility Extraction Tools")
        self.plot1 = None
        self.plot2 = None
        self.masspanel = None
        self.ctlzout = None
        self.config = None
        self.massdat = []
        self.ccsdat = []
        self.totalgrid = []
        self.pks = []
        self.ztab = []
        self.zout = 0

    def initialize_interface(self, massdat, ccsdat, mccsgrid, config, pks):
        """
        Initilizes the interface and plots the intial results.
        :param massdat: Mass distribution array (N x 2)
        :param ccsdat: CCS distribution array (M x 2)
        :param mccsgrid: Array of intensity values for corresponding mass and CCS values (N x M) array
        :param config: UniDecConfig object
        :param pks: Peaks object
        :return: None
        """
        self.config = config
        self.massdat = massdat
        self.ccsdat = ccsdat
        self.totalgrid = mccsgrid
        self.pks = pks
        self.ztab = np.arange(float(self.config.startz), float(self.config.endz + 1))
        zstrings = [str(int(z)) for z in self.ztab]
        zstrings.append("All")

        pnl = wx.Panel(self)

        vbox = wx.BoxSizer(wx.VERTICAL)
        sb = wx.StaticBox(pnl, label='Ion Mobility Extraction Tool')
        sbs = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        figsize = (6, 5)
        self.plot1 = plot2d.Plot2d(pnl, figsize=figsize)
        self.plot2 = plot1d.Plot1d(pnl, figsize=figsize)
        plotsizer = wx.BoxSizer(wx.VERTICAL)
        plotsizer.Add(self.plot1, 0, wx.EXPAND)
        plotsizer.Add(self.plot2, 0, wx.EXPAND)
        sbs.Add(plotsizer, 0, wx.EXPAND)

        sb2 = wx.StaticBox(pnl, label='Extraction Parameters')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)
        gbox1c = wx.GridBagSizer(wx.VERTICAL)
        size1 = (75, -1)

        self.ctlzout = wx.ComboBox(pnl, value="", size=size1, choices=zstrings)
        gbox1c.Add(self.ctlzout, (0, 1), span=(1, 1))
        gbox1c.Add(wx.StaticText(pnl, label="Charge State: "), (0, 0), flag=wx.ALIGN_CENTER_VERTICAL)

        sbs2.Add(gbox1c, 0, wx.EXPAND)

        self.masspanel = IMListCtrlPanel(pnl, size=(200, 700))
        addbutton = wx.Button(pnl, label="Add Species")
        plotbutton = wx.Button(pnl, label="Plot Species")
        sbs2.Add(plotbutton, 0, wx.EXPAND)
        sbs2.Add(self.masspanel, 0, wx.EXPAND)
        sbs2.Add(addbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_add, addbutton)
        self.Bind(wx.EVT_BUTTON, self.on_plot, plotbutton)

        sbs.Add(sbs2, 0, wx.EXPAND)
        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)
        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        vbox.Fit(self)

        self.CenterOnParent()
        self.loadpeaks(0)
        self.on_plot(0)
        self.ctlzout.SetSelection(len(zstrings) - 1)

    def on_close(self, e):
        """
        Close the dialog and set self.config.zout to self.zout.
        :param e: Unused event
        :return: None
        """
        self.config.zout = self.zout
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the dialog but do not update self.config.zout.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def loadpeaks(self, e):
        """
        Load masses from self.pks.peaks into the masspanel.
        :param e: Unused event
        :return: None
        """
        for p in self.pks.peaks:
            self.masspanel.list.add_line(val=p.mass)

    def get_from_gui(self, e):
        """
        Load from GUI to self.zout. If nothing is set, self.zout = 0.
        :param e: Unused event
        :return: None
        """
        try:
            self.zout = int(self.ctlzout.GetStringSelection())
            print("Charge State: ", self.zout)
        except ValueError:
            self.zout = 0

    def on_add(self, e):
        """
        Add a blank line to the masspanel.
        :param e: Unused event
        :return: None
        """
        self.masspanel.list.add_line()

    def on_plot(self, e):
        """
        First, it updates from the GUI.
        Second, if _zout_str(self.zout).bin file exists, it will be imported.
            This allows slices of specific charge states to be extracted. Otherwise, the default is all charge states.
        Third, the values are extracted and plotted.
        :param e: Unused event
        :return: None
        """
        # TODO: Split some of this off from the window into independent functions in IM_functions.py.
        # Get parameters
        self.get_from_gui(e)
        outs = np.array(self.masspanel.list.get_list())
        fname = self.config.outfname + "_zout_" + str(self.zout) + ".bin"
        if os.path.isfile(fname):
            zoutgrid = np.fromfile(fname, dtype=self.config.dtype)
        else:
            zoutgrid = self.totalgrid
        # 2D plot
        self.plot1.contourplot(xvals=self.massdat[:, 0], yvals=self.ccsdat[:, 0], zgrid=zoutgrid,
                               config=self.config, ylab="CCS (${\AA}$$^2$)", title="Mass vs. CCS", test_kda=True)
        # 1D CCS projection
        zoutgrid = np.reshape(zoutgrid, (len(self.massdat), len(self.ccsdat)))
        ccsproj = np.sum(zoutgrid, axis=0)
        self.plot2.plotrefreshtop(self.ccsdat[:, 0], ccsproj / np.amax(ccsproj), "CCS Extract",
                                  "CCS (${\AA}$$^2$)", "Normalized Intensity", "", self.config)
        colormap = cm.get_cmap(self.config.peakcmap, len(outs))
        xcolors = colormap(np.arange(len(outs)))
        for i, l in enumerate(outs):
            mass = l[0]
            ccsext = zoutgrid[ud.nearest(self.massdat[:, 0], mass)]
            ccsmax = self.ccsdat[np.argmax(ccsext), 0]
            l[1] = ccsmax
            self.plot2.plotadd(self.ccsdat[:, 0], ccsext / np.amax(ccsext), xcolors[i], '')
        self.plot2.repaint()
        # Update colors on mass list
        self.masspanel.list.populate(outs, colors=xcolors)

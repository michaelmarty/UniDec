import os
import numpy as np
import wx
from pubsub import pub
import wx.lib.mixins.listctrl as listmix
import unidec_modules.unidectools as ud
from unidec_modules import plot1d
import unidec
from unidec_modules.gui_elements import peaklistsort
from unidec_modules.unidec_presbase import UniDecPres
from unidec_modules.gui_elements.mainwindow_base import MainwindowBase
import multiprocessing


class UDListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
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
        self.InsertColumn(0, "File Name")
        self.InsertColumn(1, "Other")
        self.SetColumnWidth(0, 200)
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
        print(data)
        for i in range(0, len(data)):
            index = self.InsertItem(i, str(data[i][0]))
            self.SetItemData(index, i)
            try:
                self.SetItem(index, 1, str(data[i][1]))
            except:
                self.SetItem(index, 1, "")
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


class UDListCtrlPanel(wx.Panel):
    def __init__(self, parent, size=(300, 800)):
        """
        Creates the panel for the IMListCtrl
        :param parent: Parent window or panel
        :param size: Size in pixels foor list control
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        tid = wx.NewIdRef()
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = UDListCtrl(self, tid, size=size, style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class UniListWindow(MainwindowBase):
    def __init__(self, parent, title, config, iconfile="logo.ico", tabbed=None):
        """
        initialize window and feed in links to presenter and config.

        :param parent: GUniDec Presenter -> self.pres
        :param title: Window title (string)
        :param config: UniDecConfig object ->self.config
        :return: None
        """
        MainwindowBase.__init__(self, parent, title, config, iconfile, tabbed)
        self.CreateStatusBar(4)
        self.SetStatusWidths([-1, 150, 150, 150])
        pub.subscribe(self.on_motion, 'newxy')
        # pub.subscribe(self.on_selection, 'scans_selected')

        self.filemenu = wx.Menu()
        self.menuOpen = self.filemenu.Append(wx.ID_SAVE, "Open File Set", "Open File Set")
        # self.Bind(wx.EVT_MENU, self.on_open, self.menuOpen)

        self.menuSave = self.filemenu.Append(wx.ID_SAVE, "Save File Set", "Save File Set")
        # self.Bind(wx.EVT_MENU, self.on_open, self.menuOpen)

        self.menuAdd = self.filemenu.Append(wx.ID_SAVE, "Add Files To Set", "Add Files To Set")
        self.Bind(wx.EVT_MENU, self.pres.on_add, self.menuAdd)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.SetMenuBar(self.menuBar)

        self.panel = wx.Panel(self)
        # self.panel.SetDropTarget(ChromDropTarget(self))
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.inputsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.filelist = UDListCtrlPanel(self.panel)
        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_click, self.filelist.list)

        self.inputsizer.Add(self.filelist)

        self.plotsizer = wx.BoxSizer(wx.VERTICAL)
        self.plot4 = plot1d.Plot1d(self.panel)
        self.plot2 = plot1d.Plot1d(self.panel)
        self.insetplot = plot1d.Plot1d(self.panel)
        self.multiplot = plot1d.Plot1d(self.panel)

        self.plotsizer.Add(self.multiplot)
        self.plotsizer.Add(self.insetplot)
        self.plotsizer.Add(self.plot4)
        self.plotsizer.Add(self.plot2)

        self.inputsizer.Add(self.plotsizer)

        self.peakpanel = peaklistsort.PeakListCtrlPanel(self.panel)
        self.Bind(self.peakpanel.EVT_DELETE_SELECTION_2, self.pres.on_replot, self.peakpanel)
        self.Bind(self.peakpanel.EVT_CHARGE_STATE, self.pres.on_charge_states, self.peakpanel)
        self.Bind(self.peakpanel.EVT_DIFFERENCES, self.pres.on_differences, self.peakpanel)
        self.Bind(self.peakpanel.EVT_MASSES, self.pres.on_label_masses, self.peakpanel)
        self.inputsizer.Add(self.peakpanel, 0, wx.EXPAND)

        self.sizer.Add(self.inputsizer, 1, wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)

        self.Centre()
        self.Show(True)

    def populate(self, listdat):
        self.filelist.list.populate(listdat)

    def on_motion(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=3)
        except:
            pass

    def on_click(self, e):
        index = e.GetData()
        self.pres.update_gui(index)


class UniListApp(UniDecPres):
    def __init__(self, *args, **kwargs):
        """
        Initialize App
        :param args:
        :param kwargs:
        :return: UniDecApp object
        """
        UniDecPres.__init__(self, *args, **kwargs)
        self.init(*args, **kwargs)

    def init(self, *args, **kwargs):
        """
        Initialize Engine and View.
        :param args:
        :param kwargs:
        :return:
        """
        self.eng = unidec.UniDec()
        self.view = UniListWindow(self, "UniDec List Viewer", self.eng.config)

        # pub.subscribe(self.on_integrate, 'integrate')
        # pub.subscribe(self.on_smash, 'smash')
        # pub.subscribe(self.on_get_mzlimits, 'mzlimits')
        # pub.subscribe(self.on_left_click, 'left_click')

        if False:
            p1 = "C:\\Users\\margo\\Desktop\\Marty\\wanglabsample_protease_3.txt"
            p2 = "C:\\Users\\margo\\Desktop\\Marty\\wanglabsample_protease_3.txt"
            p3 = "C:\\Users\\margo\\Desktop\\Marty\\wanglabsample_protease_3.txt"

            paths = [[p1], [p2], [p3]]

            self.populate_list(paths)

    def on_add(self, e=0):
        """
        Open a dialog to add a file. Triggers self.add_file()
        :param e: Event
        :return: None
        """
        dlg = wx.FileDialog(self.view, "Choose a data file in x y list, mzML, or Thermo Raw format", '', "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            path = os.path.join(dirname, filename)
            print(path)
            self.add_file(path)
        dlg.Destroy()

    def add_file(self, path):
        """
        Adds path to the list and imports its engine.
        :param path: Path of data to add to the current list
        :return: None
        """
        self.listdat.append([path])
        self.view.populate(self.listdat)
        self.paths = [l[0] for l in self.listdat]
        self.engs.append(self.single_eng(path))

    def multiplot(self):
        mzdata = self.eng.data.data2
        massdata = self.eng.data.massdat
        self.view.multiplot.multiplot(mzdata[:, 0], mzdata[:, 1], massdata[:, 0], massdata[:, 1])

    def insetplot(self):
        massdata = self.eng.data.massdat
        mzdata = self.eng.data.data2
        self.view.insetplot.insetplot(mzdata[:, 0], mzdata[:, 1], massdata[:, 0], massdata[:, 1])

    def populate_list(self, listdat):
        """
        Add listdat items to the left hand list and get the engines for each. Populates the self.listdat and self.paths parameters.
        :param listdat: List of paths and associated data in N x M list.
        :return: None
        """
        self.listdat = listdat
        self.view.populate(listdat)
        self.paths = [l[0] for l in listdat]
        self.get_engs(self.paths)

    def get_engs(self, paths):
        """
        Creates a UniDecEngine object for each path and imports the data into that engine.
        :param paths: List of paths
        :return: None
        """
        self.engs = []
        for p in paths:
            eng = self.single_eng(p)
            self.engs.append(eng)
        return self.engs

    def single_eng(self, path):
        eng = unidec.UniDec()
        eng.open_file(path)
        eng.unidec_imports(everything=True, efficiency=True)
        eng.pick_peaks()
        # eng.dscore()
        return eng

    def update_gui(self, index):
        self.eng = self.engs[index]
        self.view.peakpanel.add_data(self.eng.pks)
        self.on_replot()

    def on_replot(self, e=None):
        """
        Makes Mass Plots (plot 2) and m/z Plots (plot 4)
        :param e: Event (unused)
        :return: None
        """

        self.multiplot()
        self.insetplot()
        self.makeplot4(0)
        self.makeplot2(0)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = UniListApp()
    app.start()

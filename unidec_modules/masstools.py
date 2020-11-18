import sys
import string
import wx.lib.mixins.listctrl as listmix
import wx
import numpy as np
from scipy.interpolate import interp1d

from unidec_modules.AutocorrWindow import AutocorrWindow
from unidec_modules.isolated_packages import FileDialogs, biopolymer_calculator
import unidec_modules.unidectools as ud

'''
Module for window defining the oligomers, the expected masses, and for matching peaks to the defined oligomers.
'''


# TODO: Inherit a subclass with basic list control functions like delete


class MassListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0, coltitle="Mass (Da)"):
        """
        Create the mass list ctrl with one column.
        :param parent: Passed to wx.ListCtrl
        :param id_value: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, coltitle)
        self.SetColumnWidth(0, width=190)  # , wx.LIST_AUTOSIZE)

        self.popupID1 = wx.NewIdRef()
        self.popupID3 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_masslist_delete, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_biopolymer, id=self.popupID3)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click_masslist, self)

    def populate(self, listctrldata):
        """
        Populate the the listctrl from a list of mass values
        :param listctrldata: mass value list
        :return: None
        """
        self.DeleteAllItems()
        for i in range(0, len(listctrldata)):
            index = self.InsertItem(i, str(listctrldata[i]))
            self.SetItemData(index, i)

    def clear(self):
        """
        Clear the list.
        :return: None
        """
        self.DeleteAllItems()

    def add_line(self):
        """
        Add a new line, default of 0.
        :return: None
        """
        self.InsertItem(10000, str(0))

    def get_list(self):
        """
        Return the list as a list of mass values
        :return: List of mass values
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            list_out.append(float(self.GetItemText(i)))
        return list_out

    def on_right_click_masslist(self, event):
        """
        Create a right click menu.
        :param event: Unused event
        :return: None
        """
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID3, "Calculate Mass from Sequence")
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_masslist_delete(self, event):
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

    def on_biopolymer(self, e=0):
        bp = biopolymer_calculator.BiopolymerFrame(self)
        result = bp.ShowModal()
        if result == 0:
            print("Imported Mass from Sequence:", bp.mass)
            calc_mass = 0
            try:
                calc_mass = float(bp.mass)
                if calc_mass > 0:
                    item = self.GetFirstSelected()
                    self.SetItem(item, 0, str(calc_mass))
            except:
                pass
        bp.Destroy()


class MassListCtrlPanel(wx.Panel):
    def __init__(self, parent, columntitle="Mass (Da)", size=(210, 380)):
        """
        ListCtrlPanel for the Mass list
        :param parent: Parent panel or window
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = MassListCrtl(self, wx.NewIdRef(), coltitle=columntitle, size=size, style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class OligomerListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a 5 column list ctrl.
        :param parent: Passed to wx.ListCrtl
        :param id_value: Passed to wx.ListCrtl
        :param pos: Passed to wx.ListCrtl
        :param size: Passed to wx.ListCrtl
        :param style: Passed to wx.ListCrtl
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Base Offset (Da)")
        self.InsertColumn(1, "Monomer Mass (Da)")
        self.InsertColumn(2, "Min # of Oligomers")
        self.InsertColumn(3, "Max # of Oligomers")
        self.InsertColumn(4, "Name")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 125)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 100)
        self.SetColumnWidth(4, 60)
        self.index = 0

        self.popupID2 = wx.NewIdRef()
        self.popupID3 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_oligo_delete, id=self.popupID2)
        self.Bind(wx.EVT_MENU, self.on_biopolymer, id=self.popupID3)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click_oligolist, self)

    def clear(self):
        """
        Clear the list.
        :return: None
        """
        self.DeleteAllItems()
        self.index = 0

    def add_line(self, a=0, b=0, c=0, d=1, e=None):
        """
        Add a blank line to the list.
        :return: None
        """
        index = self.InsertItem(10000, str(a))
        if e is None:
            e = string.ascii_uppercase[self.index]
        self.SetItem(index, 1, str(b))
        self.SetItem(index, 2, str(c))
        self.SetItem(index, 3, str(d))
        self.SetItem(index, 4, e)
        self.index += 1

    def populate(self, data):
        """
        Populate from an N x 4 or N x 5 array into the listctrl.
        :param data: Array
        :return: None
        """
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data)):
            index = self.InsertItem(i, str(data[i][0]))
            try:
                self.SetItem(index, 1, str(data[i][1]))
                self.SetItem(index, 2, str(data[i][2]))
                self.SetItem(index, 3, str(data[i][3]))
            except (ValueError, IndexError):
                self.SetItem(index, 1, "")
                self.SetItem(index, 2, "")
                self.SetItem(index, 3, "")
            try:
                self.SetItem(index, 4, str(data[i][4]))
            except (ValueError, IndexError):
                self.SetItem(index, 4, "")
                # self.SetItemData(index, i)

    def get_list(self):
        """
        Return the values on the apanel as a nested list.
        :return: List of values in N x 5 format.
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       int(self.GetItem(i, col=2).GetText()), int(self.GetItem(i, col=3).GetText()),
                       self.GetItem(i, col=4).GetText()]
            list_out.append(sublist)
        return list_out

    def on_right_click_oligolist(self, event):
        """
        Open a right click menu
        :param event: Unused event
        :return: None
        """
        menu = wx.Menu()
        menu.Append(self.popupID3, "Calculate Mass from Sequence")
        menu.Append(self.popupID2, "Delete")
        self.PopupMenu(menu)
        menu.Destroy()

    def on_oligo_delete(self, event):
        """
        Delete the selected item.
        :param event: Unused event
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

    def on_biopolymer(self, e=0):
        bp = biopolymer_calculator.BiopolymerFrame(self)
        result = bp.ShowModal()
        if result == 0:
            print("Imported Mass from Sequence:", bp.mass)
            calc_mass = 0
            try:
                calc_mass = float(bp.mass)
                if calc_mass > 0:
                    item = self.GetFirstSelected()
                    self.SetItem(item, 1, str(calc_mass))
            except:
                pass
        bp.Destroy()


class OligomerListCrtlPanel(wx.Panel):
    def __init__(self, parent):
        """
        ListCtrlPanel for the Oligomer list
        :param parent: Parent panel or window
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = OligomerListCtrl(self, wx.NewIdRef(), size=(500, 200), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class MatchListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a four column uneditable list control.
        :param parent: Parent panel or window
        :param id_value: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style=wx.LC_REPORT)
        #wx.ListCtrl.__init__(self, pos=wx.DefaultPosition, size=size,
                                     #style=wx.LC_REPORT | wx.BORDER_SUNKEN)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Peak Mass (Da)")
        self.InsertColumn(1, "Match")
        self.InsertColumn(2, "Error")
        self.InsertColumn(3, "Name")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 100)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 100)
        self.index = 0

    def clear(self):
        """
        Clear the list
        :return: None
        """
        self.DeleteAllItems()
        self.index = 0

    def populate(self, data1, data2, data3, data4):
        """
        Populate the list of matches from four arrays
        :param data1: The first column, the measured mass
        :param data2: The second column, the simulated mass
        :param data3: The third column, the error between measured and simulated.
        :param data4: The fourth column, the match name.
        :return: None
        """
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data1)):
            index = self.InsertItem(i, str(data1[i]))
            self.SetItem(index, 1, str(data2[i]))
            self.SetItem(index, 2, str(data3[i]))
            self.SetItem(index, 3, str(data4[i]))
            # self.SetItemData(index, i)

    def get_list(self):
        """
        :return: List items in a nested list
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       float(self.GetItem(i, col=2).GetText()), self.GetItem(i, col=3).GetText()]
            list_out.append(sublist)
        return list_out


class MatchListCrtlPanel(wx.Panel):
    def __init__(self, parent):
        """
        ListCtrlPanel for the Match list
        :param parent: Parent panel or window
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = MatchListCrtl(self, wx.NewIdRef(), size=(500, 200), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class DummyPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        self.parent = parent


class CommonMasses(wx.ListCtrl,  # listmix.ListCtrlAutoWidthMixin,
                   listmix.TextEditMixin, listmix.ColumnSorterMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a four column uneditable list control.
        :param parent: Parent panel or window
        :param id_value: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        self.parent = parent
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        # listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        listmix.ColumnSorterMixin.__init__(self, 3)
        self.InsertColumn(0, "Name")
        self.InsertColumn(1, "Mass")
        self.InsertColumn(2, "Type")
        self.SetColumnWidth(0, 200)
        self.SetColumnWidth(1, 100)
        self.SetColumnWidth(2, 100)
        self.index = 0

        self.data = []

        self.popupID1 = wx.NewIdRef()
        self.popupID2 = wx.NewIdRef()
        self.popupID3 = wx.NewIdRef()
        self.popupID4 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_transfer, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_delete, id=self.popupID2)
        self.Bind(wx.EVT_MENU, self.clear, id=self.popupID3)
        self.Bind(wx.EVT_MENU, self.repopulate, id=self.popupID4)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self)

    def on_right_click(self, event):
        """
        Open a right click menu
        :param event: Unused event
        :return: None
        """
        menu = wx.Menu()
        menu.Append(self.popupID1, "Add to Oligomer Builder")
        menu.Append(self.popupID2, "Delete")
        menu.Append(self.popupID3, "Delete All")
        menu.Append(self.popupID4, "Repopulate")
        self.PopupMenu(menu)
        menu.Destroy()

    def on_transfer(self, e):
        self.parent.parent.parent.on_common_to_oligo(e)

    def clear(self, e=None):
        """
        Clear the list
        :return: None
        """
        self.DeleteAllItems()
        self.index = 0
        self.data = []
        self.itemDataMap = self.data

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

    def populate(self, data1, data2, data3):
        """
        Populate the list of matches from four arrays
        :param data1: The first column, the measured mass
        :param data2: The second column, the simulated mass
        :param data3: The third column, the error between measured and simulated.

        :param data4: The fourth column, the match name.
        :return: None
        """
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data1)):
            index = self.InsertItem(i, str(data1[i]))
            self.SetItem(index, 1, str(data2[i]))
            self.SetItem(index, 2, str(data3[i]))

    def get_list(self):
        """
        :return: List items in a nested list
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [str(self.GetItem(i, col=0).GetText()), str(self.GetItem(i, col=1).GetText()),
                       str(self.GetItem(i, col=2).GetText())]
            list_out.append(sublist)
        return np.array(list_out)

    def get_selection(self):
        i = self.GetFirstSelected()
        num = self.GetSelectedItemCount()
        list_out = []
        sublist = [str(self.GetItem(i, col=0).GetText()), str(self.GetItem(i, col=1).GetText()),
                   str(self.GetItem(i, col=2).GetText())]
        list_out.append(sublist)
        for j in range(1, num):
            i = self.GetNextSelected(i)
            sublist = [str(self.GetItem(i, col=0).GetText()), str(self.GetItem(i, col=1).GetText()),
                       str(self.GetItem(i, col=2).GetText())]
            list_out.append(sublist)
        return np.array(list_out)

    def add_line(self):
        """
        Add a blank line to the list.
        :return: None
        """
        index = self.InsertItem(10000, string.ascii_uppercase[self.index])
        self.SetItem(index, 1, str(0))
        self.SetItem(index, 2, "User")
        self.index += 1

    def repopulate(self, e):
        self.populate(self.parent.parent.parent.commonmasses[:, 0], self.parent.parent.parent.commonmasses[:, 1],
                      self.parent.parent.parent.commonmasses[:, 2])
        pass

    def GetListCtrl(self):
        # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        count = self.GetItemCount()
        self.data = []
        for i in range(0, count):
            self.SetItemData(i, i)
            sublist = [str(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       str(self.GetItem(i, col=2).GetText())]
            self.data.append(sublist)
        self.itemDataMap = self.data
        return self


class CommonMassesPanel(wx.Panel):
    def __init__(self, parent):
        """
        ListCtrlPanel for the Match list
        :param parent: Parent panel or window
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.parent = parent
        self.list = CommonMasses(self, wx.NewIdRef(), size=(400, 500), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class MassSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Initialize the dialog
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((1200, 700))
        self.SetTitle("Mass and Oligomer Tools")

        self.pks = None
        self.config = None
        self.massdat = None
        self.commonmasses = None

        self.masslistbox = None
        self.oligomerlistbox = None
        self.matchlistbox = None

    def init_dialog(self, config, pks, massdat=None):
        """
        Creates the dialog window and displays it.
        :param config: UniDecConfig object
        :param pks: Peaks object
        :param massdat: Mass distribution
        :return: None
        """
        # massbins = config.massbins
        self.massdat = massdat
        self.config = config
        self.pks = pks

        panel = DummyPanel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        sb = wx.StaticBox(panel, label='Set the Mass List for Limited UniDec')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        peakpop = wx.Button(panel, label="Populate from Peak List")
        self.Bind(wx.EVT_BUTTON, self.pop_from_peaks, peakpop)

        importbutton = wx.Button(panel, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import_masses, importbutton)

        oligopop = wx.Button(panel, label="Populate from Isolated Oligomers")
        self.Bind(wx.EVT_BUTTON, self.pop_oligo_iso, oligopop)

        oligopop2 = wx.Button(panel, label="Populate from All Possible Oligomers")
        self.Bind(wx.EVT_BUTTON, self.pop_oligo_all, oligopop2)

        clearbutt = wx.Button(panel, label="Clear List")
        self.Bind(wx.EVT_BUTTON, self.on_clear_masslist, clearbutt)

        addbutton = wx.Button(panel, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_mass, addbutton)

        simbutton = wx.Button(panel, label="Simulate These Masses")
        self.Bind(wx.EVT_BUTTON, self.on_simulate, simbutton)

        sbs.Add(peakpop, 0, wx.EXPAND)
        sbs.Add(importbutton, 0, wx.EXPAND)
        sbs.Add(oligopop, 0, wx.EXPAND)
        sbs.Add(oligopop2, 0, wx.EXPAND)
        sbs.Add(addbutton, 0, wx.EXPAND)
        sbs.Add(clearbutt, 0, wx.EXPAND)
        self.masslistbox = MassListCtrlPanel(panel)

        sbs.Add(wx.StaticText(panel, label="Mass List"))
        sbs.Add(self.masslistbox)
        sbs.Add(simbutton, 0, wx.EXPAND)

        hbox.Add(sbs)

        sb2 = wx.StaticBox(panel, label='Oligomer Maker')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)

        clearbutt2 = wx.Button(panel, label="Clear Oligomer List")
        self.Bind(wx.EVT_BUTTON, self.on_clear_oligolist, clearbutt2)

        addbutton2 = wx.Button(panel, label="Add Oligomer Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_oligomer, addbutton2)

        importbutton2 = wx.Button(panel, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import_oligos, importbutton2)

        plotbutton = wx.Button(panel, label="View Autocorrelation Plot")
        self.Bind(wx.EVT_BUTTON, self.on_autocorr_window, plotbutton)
        buttonbox = wx.BoxSizer(wx.VERTICAL)
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        buttonbox.Add(importbutton2, 0, wx.EXPAND)
        buttonbox.Add(addbutton2, 0, wx.EXPAND)
        buttonbox.Add(clearbutt2, 0, wx.EXPAND)
        buttonbox.Add(plotbutton, 0, wx.EXPAND)
        hbox3.Add(buttonbox)
        textbox = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(panel,
                             label="  For i from Min # to Max #:\n      Mass(i)=Base Offset + Monomer Mass * i \n")
        font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        text.SetFont(font)
        textbox.Add(text)

        self.ctlmatcherror = wx.TextCtrl(panel, value=str(self.config.matchtolerance))
        textbox.Add(wx.StaticText(panel, label="Error Tolerance for Matching (Da)"), 0, wx.ALIGN_RIGHT)
        textbox.Add(self.ctlmatcherror, 0, wx.ALIGN_RIGHT)

        hbox3.Add(textbox)
        sbs2.Add(hbox3, 0, wx.EXPAND)

        self.oligomerlistbox = OligomerListCrtlPanel(panel)

        sbs2.Add(wx.StaticText(panel, label="Oligomer List"))
        sbs2.Add(self.oligomerlistbox)
        sbs2.Add(wx.StaticText(panel, label=""))

        sb4 = wx.StaticBox(panel, label="Match Peaks to Oligomers")
        sbs4 = wx.StaticBoxSizer(sb4, orient=wx.VERTICAL)
        match_iso_button = wx.Button(panel, label="Match to Isolated Oligomers")
        match_iso_button.SetToolTip(wx.ToolTip("Match peaks to isolated oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.on_match_isolated, match_iso_button)
        match_all_button = wx.Button(panel, label="Match to Mixed Oligomers")
        match_all_button.SetToolTip(
            wx.ToolTip("Match peaks to any possible combination of oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.on_match_all, match_all_button)
        self.matchlistbox = MatchListCrtlPanel(panel)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(match_iso_button, 1, wx.EXPAND)
        hbox2.Add(match_all_button, 1, wx.EXPAND)
        sbs4.Add(hbox2, 0, wx.EXPAND)
        sbs4.Add(self.matchlistbox)
        sbs2.AddStretchSpacer(prop=1)
        sbs2.Add(sbs4)
        hbox.Add(sbs2)

        sb5 = wx.StaticBox(panel, label='Common Masses')
        sbs5 = wx.StaticBoxSizer(sb5, orient=wx.VERTICAL)

        importbutton2 = wx.Button(panel, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_load_common_masses, importbutton2)

        savecommonbutton = wx.Button(panel, label="Save Common Masses")
        self.Bind(wx.EVT_BUTTON, self.on_save_common_masses, savecommonbutton)

        addbutton4 = wx.Button(panel, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_new_common_mass, addbutton4)

        sbs5.Add(importbutton2, 0, wx.EXPAND)
        sbs5.Add(savecommonbutton, 0, wx.EXPAND)
        sbs5.Add(addbutton4, 0, wx.EXPAND)
        self.commonmassespanel = CommonMassesPanel(panel)

        sbs5.Add(wx.StaticText(panel, label="Common Masses List"))
        sbs5.Add(self.commonmassespanel)

        hbox.Add(sbs5)

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

        self.masslistbox.list.populate(self.config.masslist)
        self.oligomerlistbox.list.populate(self.config.oligomerlist)

        defaultmatchlist = self.config.matchlist
        if len(defaultmatchlist) == 4:
            self.matchlistbox.list.populate(defaultmatchlist[0], defaultmatchlist[1], defaultmatchlist[2],
                                            defaultmatchlist[3])
        try:
            self.load_common_masses(self.config.masstablefile)
        except:
            print("Unable to load common masses")
        # self.load_common_masses(self.config.masstablefile)
        self.CenterOnParent()

    def on_common_to_oligo(self, e):
        outdata = self.commonmassespanel.list.get_selection()
        print(outdata)
        for o in outdata:
            print(o)
            self.oligomerlistbox.list.add_line(0, o[1], 0, 1, o[0])

    def on_add_new_common_mass(self, e):
        self.commonmassespanel.list.add_line()

    def on_save_common_masses(self, e):
        outdata = self.commonmassespanel.list.get_list()
        cmfilename = FileDialogs.save_file_dialog("Save Common Masses File", "*.csv*", self.config.masstablefile)
        if cmfilename is not None:
            np.savetxt(cmfilename, outdata, delimiter=",", fmt="%s")
        print("Saving Common Masses to:", cmfilename)

    def on_load_common_masses(self, e):
        cmfilename = FileDialogs.open_file_dialog("Open Common Masses File", file_types="*.csv*",
                                                  default=self.config.masstablefile)
        if cmfilename is not None:
            self.load_common_masses(cmfilename)

    def load_common_masses(self, file):
        print("Loading Common Masses From File:", file)
        self.commonmasses = np.genfromtxt(file, delimiter=',', usecols=(0, 1, 2), dtype=str, skip_header=1)
        self.commonmassespanel.list.populate(self.commonmasses[:, 0], self.commonmasses[:, 1], self.commonmasses[:, 2])

    def on_simulate(self, e):
        """
        Replaces the peaks in self.pks with the masses defined in the self.masslistbox.
        The intensity information comes from interpolating the mass distribution in self.massdat.
        It flags self.pks.changes as 1, which is picked up downstream to say that these new peaks should be replotted.
        :param e: Unused event
        :return: None
        """
        newmasslist = self.masslistbox.list.get_list()
        newpeaks = []
        f = interp1d(self.massdat[:, 0], self.massdat[:, 1], bounds_error=False, fill_value=0)
        for mass in newmasslist:
            intensity = f(mass)
            newpeaks.append([mass, intensity])
        newpeaks = np.array(newpeaks)
        self.pks.__init__()
        self.pks.add_peaks(newpeaks)
        self.pks.default_params(cmap=self.config.peakcmap)
        self.pks.changed = 1

    def on_autocorr_window(self, e):
        """
        Opens up an autocorrelation window for the purpose of figuring out mass differences.
        :param e: Unused event
        :return: None
        """
        if not ud.isempty(self.massdat):
            dlg = AutocorrWindow(self)
            dlg.initalize_dialog(self.config, self.massdat)
            dlg.ShowModal()

    def on_match_isolated(self, e):
        """
        Match the peaks in self.pks to possible combination of isolated rows of oligomers in the oligomerlist.
        Uses ud.make_isolated_match function.
        Populated the matchlistbox with the results of the matching.
        :param e: Unused event
        :return: None
        """
        oligos = self.oligomerlistbox.list.get_list()
        try:
            tolerance = float(self.ctlmatcherror.GetValue())
        except ValueError:
            tolerance = None
        oligomasslist, oligonames = ud.make_isolated_match(oligos)
        matchlist = ud.match(self.pks, oligomasslist, oligonames, tolerance=tolerance)
        self.matchlistbox.list.populate(matchlist[0], matchlist[1], matchlist[2], matchlist[3])

    def on_match_all(self, e):
        """
        Match the peaks in self.pks to all possible combination of oligomers in the oligomerlist.
        Uses ud.make_all_matches function.
        Populated the matchlistbox with the results of the matching.
        :param e: Unused event
        :return: None
        """
        oligos = self.oligomerlistbox.list.get_list()
        try:
            tolerance = float(self.ctlmatcherror.GetValue())
        except ValueError:
            tolerance = None
        oligomasslist, oligonames = ud.make_all_matches(oligos)
        if ud.isempty(oligomasslist):
            print("ERROR: Need to specify the Potential Oligomers")
            return
        matchlist = ud.match(self.pks, oligomasslist, oligonames, tolerance=tolerance)
        self.matchlistbox.list.populate(matchlist[0], matchlist[1], matchlist[2], matchlist[3])

    def on_close(self, e):
        """
        Close the dialog and apply the changes.

        Sets self.config.masslist to the defined mass list.
        Sets self.config.oligomerlist to the defined oligomers.
        Sets self.config.matchlist to the matched values
        :param e: Unused event
        :return: None
        """
        try:
            tolerance = float(self.ctlmatcherror.GetValue())
        except ValueError:
            tolerance = None
        self.config.matchtolerance = tolerance

        newmasslist = self.masslistbox.list.get_list()
        # print newmasslist
        if not ud.isempty(newmasslist):
            newmasslist = np.array(newmasslist)
            newmasslist = newmasslist[newmasslist > 0]
            # print newmasslist
            if len(newmasslist) > 0:
                self.config.masslist = newmasslist
            else:
                self.config.masslist = []
        oligos = self.oligomerlistbox.list.get_list()
        if not ud.isempty(oligos):
            oligos = np.array(oligos)
            oligoshort = oligos[:, :2]
            oligoshort = oligoshort.astype(np.float)
            oligos = oligos[np.any([oligoshort != 0], axis=2)[0], :]
            # oligos=oligos[::-1]
        self.config.oligomerlist = oligos

        matchlist = np.transpose(self.matchlistbox.list.get_list())
        if not ud.isempty(matchlist):
            self.config.matchlist = matchlist

        self.Destroy()
        try:
            self.EndModal(0)
        except Exception as e:
            pass

    def on_close_cancel(self, e):
        """
        Close the dialog but do not modify any of the values.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def pop_from_peaks(self, e):
        """
        Populate the mass list from the peak masses in self.pks
        :param e: Unused event
        :return: None
        """
        try:
            self.masslistbox.list.populate(self.pks.masses)
        except AttributeError:
            print("Pick peaks first")

    def pop_oligo_iso(self, e):
        """
        Populate the mass list with isolated oligomer rows.
        :param e: Unused event
        :return: None
        """
        oligos = self.oligomerlistbox.list.get_list()
        oligomasslist, oligonames = ud.make_isolated_match(oligos)
        oligomasslist = np.unique(oligomasslist)
        self.masslistbox.list.populate(oligomasslist)

    def pop_oligo_all(self, e):
        """
        Populates the mass list with all possible oligomers.
        :param e: Unused event
        :return: None
        """
        oligos = self.oligomerlistbox.list.get_list()
        oligomasslist, oligonames = ud.make_all_matches(oligos)
        oligomasslist = np.unique(oligomasslist)
        self.masslistbox.list.populate(oligomasslist)

    def on_clear_masslist(self, e):
        """
        Clears the mass list.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.clear()

    def on_clear_oligolist(self, e):
        """
        Clears the oligomer list.
        :param e: Unused Event
        :return: None
        """
        self.oligomerlistbox.list.clear()

    def on_add_mass(self, e):
        """
        Adds a blank line to the mass list.
        :param e: Unused Event
        :return: None
        """
        self.masslistbox.list.add_line()

    def on_add_oligomer(self, e):
        """
        Adds a blank line to the oligomer list.
        :param e: Unused event
        :return: None
        """
        self.oligomerlistbox.list.add_line()

    def on_import_masses(self, e):
        """
        Opens a dialog to import mass list files.
        :param e: Unused event
        :return: None
        """
        mfilename = FileDialogs.open_file_dialog("Open Mass File (mfile)", file_types="*.*")
        if mfilename is not None:
            importmass = np.loadtxt(mfilename)
            if importmass[0] > 0:
                self.masslistbox.list.populate(importmass)

    def on_import_oligos(self, e):
        """
        Open a file dialog to import oligomer files.
        :param e: Unused event
        :return: None
        """
        # TODO: This fails when the output ofile isn't square (some are unnamed and other are named)
        ofilename = FileDialogs.open_file_dialog("Open Oligomer File (ofile)", file_types="*.*")
        if ofilename is not None:
            importolig = np.genfromtxt(ofilename, dtype='str')
            if np.shape(importolig) == (5,) or np.shape(importolig) == (4,):
                importolig = [importolig]
            self.oligomerlistbox.list.populate(importolig)


# Main App Execution
if __name__ == "__main__":
    import unidec

    eng = unidec.UniDec()
    app = wx.App(False)
    frame = MassSelection(None)
    frame.init_dialog(eng.config, eng.pks, massdat=None)
    frame.ShowModal()
    app.MainLoop()

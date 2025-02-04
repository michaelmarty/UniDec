import string
import wx.lib.mixins.listctrl as listmix
import wx
import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from unidec.modules.AutocorrWindow import AutocorrWindow
from unidec.modules.isolated_packages import FileDialogs
from unidec.modules.isolated_packages import biopolymer_calculator, spreadsheet
from unidec.modules import matchtools
import unidec.tools as ud

'''
Module for window defining the oligomers, the expected masses, and for matching peaks to the defined oligomers.
'''

white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0, 0, 0)


class MassListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, panel, pos=wx.DefaultPosition, size=wx.DefaultSize, style=wx.LC_REPORT,
                 coltitle="Mass (Da)"):
        """
        Create the mass list ctrl with one column.
        :param parent: Passed to wx.ListCtrl
        :param panel: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        id_value = wx.NewIdRef()
        wx.ListCtrl.__init__(self, panel, id_value, pos, size, style)
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

    def on_right_click_masslist(self, event=None):
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


class OligomerListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, panel, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a 5 column list ctrl.
        :param parent: Passed to wx.ListCrtl
        :param panel: Passed to wx.ListCrtl
        :param pos: Passed to wx.ListCrtl
        :param size: Passed to wx.ListCrtl
        :param style: Passed to wx.ListCrtl
        :return: None
        """
        id_value = wx.NewIdRef()
        wx.ListCtrl.__init__(self, panel, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.parent = parent

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
        self.popupID4 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_oligo_delete, id=self.popupID2)
        self.Bind(wx.EVT_MENU, self.on_biopolymer, id=self.popupID3)
        self.Bind(wx.EVT_MENU, self.on_add_to_common, id=self.popupID4)
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
        indexes = range(0, count)
        return self.get_items(indexes)

    def get_items(self, indexes):
        list_out = []
        for i in indexes:
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
        menu.Append(self.popupID4, "Add to Common Masses")
        menu.AppendSeparator()
        menu.Append(self.popupID2, "Delete")
        self.PopupMenu(menu)
        menu.Destroy()

    def on_add_to_common(self, event):
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
        olist = self.get_items(selection)
        for o in olist:
            self.parent.add_to_common_masses(o[4], o[1], "User")
            print('Added', o[4], 'to common mass table')

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


class MatchListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, panel, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a four column uneditable list control.
        :param parent: Parent panel or window
        :param panel: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        id_value = wx.NewIdRef()
        wx.ListCtrl.__init__(self, panel, id_value, pos, size, style=wx.LC_REPORT)
        # wx.ListCtrl.__init__(self, pos=wx.DefaultPosition, size=size,
        # style=wx.LC_REPORT | wx.BORDER_SUNKEN)
        self.parent = parent
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
        self.mode = 0

        self.popupID1 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_view_alt, id=self.popupID1)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self)

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

    def color_item(self, i, color):
        # Adjust background color
        try:
            color = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)),
                              int(round(color[2] * 255)), alpha=255)
        except Exception:
            color = wx.Colour(255, 255, 255, alpha=255)
        self.SetItemBackgroundColour(i, col=color)

        # Adjust Text Color
        luminance = ud.get_luminance(color, type=2)
        # print(wx.Colour(colout), luminance)
        if luminance < ud.luminance_cutoff:
            self.SetItemTextColour(i, col=white_text)
        else:
            self.SetItemTextColour(i, col=black_text)

    def on_view_alt(self, e=None):
        item = self.GetFirstSelected()
        self.parent.on_view_alt(item)

    def on_right_click(self, event):
        """
        Create a right click menu.
        :param event: Unused event
        :return: None
        """
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            if self.mode == 0:
                menu.Append(self.popupID1, "View Alternates")
            else:
                menu.Append(self.popupID1, "Select Match")
            self.PopupMenu(menu)
            menu.Destroy()


class DummyPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        self.parent = parent


class CommonMasses(wx.ListCtrl,  # listmix.ListCtrlAutoWidthMixin,
                   listmix.TextEditMixin, listmix.ColumnSorterMixin):
    def __init__(self, parent, panel, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        """
        Create a four column uneditable list control.
        :param parent: Parent panel or window
        :param panel: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :return: None
        """
        id_value = wx.NewIdRef()
        self.parent = parent
        wx.ListCtrl.__init__(self, panel, id_value, pos, size, style)
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
        self.popupsitemode = wx.NewIdRef()
        self.popupID2 = wx.NewIdRef()
        self.popupID3 = wx.NewIdRef()
        self.popupID4 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_transfer, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_transfer_site, id=self.popupsitemode)
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
        menu.Append(self.popupsitemode, "Add to Site Mode")
        menu.AppendSeparator()
        menu.Append(self.popupID2, "Delete")
        menu.Append(self.popupID3, "Delete All")
        menu.Append(self.popupID4, "Repopulate")
        self.PopupMenu(menu)
        menu.Destroy()

    def on_transfer(self, e):
        self.parent.on_common_to_oligo(e)

    def on_transfer_site(self, e):
        self.parent.on_add_to_ss(e)

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

    def add_line(self, name=None, mass="0", type="User"):
        """
        Add a line to the list.
        :return: None
        """
        if name is None:
            name = string.ascii_uppercase[self.index]
        index = self.InsertItem(10000, str(name))
        self.SetItem(index, 1, str(mass))
        self.SetItem(index, 2, str(type))
        self.index += 1

    def repopulate(self, e):
        self.populate(self.parent.commonmasses[:, 0], self.parent.commonmasses[:, 1],
                      self.parent.commonmasses[:, 2])
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


class MassSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Initialize the dialog
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.matchtype = None
        self.oligoindexes = None
        self.basemass = 0
        self.sdf = None
        self.SetSize((1200, 700))
        self.SetTitle("Mass and Oligomer Tools")

        self.pks = None
        self.config = None
        self.massdat = None
        self.commonmasses = None

        self.masslistbox = None
        self.oligomerlistbox = None
        self.matchlistbox = None
        self.ctlmatcherror = None
        self.commonmassespanel = None
        self.diffmatrix = None
        self.matchlist = None
        self.oligomasslist = None
        self.oligonames = None
        self.peakindex = None
        self.peakmass = None
        self.mnames = None
        self.mmasses = None
        self.diffs = None
        self.tolerance = 1000
        self.notebook = None

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

        oligopop3 = wx.Button(panel, label="Populate from Site Mode")
        self.Bind(wx.EVT_BUTTON, self.pop_oligo_sites, oligopop3)

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
        sbs.Add(oligopop3, 0, wx.EXPAND)
        sbs.Add(addbutton, 0, wx.EXPAND)
        sbs.Add(clearbutt, 0, wx.EXPAND)
        self.masslistbox = MassListCtrl(self, panel, coltitle="Mass (Da)", size=(210, 380), style=wx.LC_REPORT)

        sbs.Add(wx.StaticText(panel, label="Mass List"))
        sbs.Add(self.masslistbox, 1, wx.EXPAND)
        sbs.Add(simbutton, 0, wx.EXPAND)

        hbox.Add(sbs, 0, wx.EXPAND)

        vmidsizer = wx.BoxSizer(wx.VERTICAL)

        self.notebook = wx.Notebook(panel)
        tab1 = wx.Panel(self.notebook)

        p1 = tab1

        sb2 = wx.StaticBox(p1, label='Oligomer Maker')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)

        clearbutt2 = wx.Button(p1, label="Clear Oligomer List")
        self.Bind(wx.EVT_BUTTON, self.on_clear_oligolist, clearbutt2)

        addbutton2 = wx.Button(p1, label="Add Oligomer Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_oligomer, addbutton2)

        importbutton2 = wx.Button(p1, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import_oligos, importbutton2)

        plotbutton = wx.Button(p1, label="View Autocorrelation Plot")
        self.Bind(wx.EVT_BUTTON, self.on_autocorr_window, plotbutton)
        buttonbox = wx.BoxSizer(wx.VERTICAL)
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        buttonbox.Add(importbutton2, 0, wx.EXPAND)
        buttonbox.Add(addbutton2, 0, wx.EXPAND)
        buttonbox.Add(clearbutt2, 0, wx.EXPAND)
        buttonbox.Add(plotbutton, 0, wx.EXPAND)
        hbox3.Add(buttonbox)
        textbox = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(p1,
                             label=" Oligomer Mode\n    Either isolated or mixed combinations\n    "
                                   "For i from Min # to Max #:\n       Mass(i)=Base Offset + Monomer Mass * i \n")
        font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        text.SetFont(font)
        textbox.Add(text)
        hbox3.Add(textbox)
        sbs2.Add(hbox3, 0, wx.EXPAND)

        self.oligomerlistbox = OligomerListCtrl(self, p1, size=(500, 200),
                                                style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)

        sbs2.Add(wx.StaticText(p1, label="Oligomer List"))
        sbs2.Add(self.oligomerlistbox, 1, wx.EXPAND)
        sbs2.Add(wx.StaticText(p1, label=""))

        tab1.SetSizerAndFit(sbs2)
        self.notebook.AddPage(tab1, "Oligomer Mode")

        tab2 = wx.Panel(self.notebook)
        sb2b = wx.StaticBox(tab2, label='Site Constructor')
        sbs2b = wx.StaticBoxSizer(sb2b, orient=wx.VERTICAL)

        buttonbox = wx.BoxSizer(wx.VERTICAL)
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)

        importbuttonsite = wx.Button(tab2, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import_sites, importbuttonsite)
        buttonbox.Add(importbuttonsite, 0, wx.EXPAND)

        addbutton3 = wx.Button(tab2, label="Add Row")
        self.Bind(wx.EVT_BUTTON, self.on_add_site, addbutton3)
        buttonbox.Add(addbutton3, 0, wx.EXPAND)

        addbutton4 = wx.Button(tab2, label="Add Column")
        self.Bind(wx.EVT_BUTTON, self.on_add_col, addbutton4)
        buttonbox.Add(addbutton4, 0, wx.EXPAND)

        clearbutt3 = wx.Button(tab2, label="Clear All")
        self.Bind(wx.EVT_BUTTON, self.on_clear_sites, clearbutt3)
        buttonbox.Add(clearbutt3, 0, wx.EXPAND)

        textbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlsitebasemass = wx.TextCtrl(tab2, value="")
        textbox2.Add(wx.StaticText(tab2, label="Base Mass (Da):"), 0, wx.ALIGN_CENTER_VERTICAL)
        textbox2.Add(self.ctlsitebasemass, 0)
        buttonbox.Add(textbox2, 0, wx.EXPAND)

        hbox3.Add(buttonbox)
        textbox = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(tab2,
                             label=" Site Mode\n    Only one bound per site\n    For any non-zero i at each site:\n "
                                   "      Mass(i, ...)=Base + Site1_i + ... \n")
        font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        text.SetFont(font)
        textbox.Add(text)
        hbox3.Add(textbox)

        sbs2b.Add(hbox3, 0, wx.EXPAND)

        self.ss = spreadsheet.SpreadsheetPanel(self, tab2, 4, 3)
        self.ss.set_col_labels(["Name", "Mass", "Site 1"])
        self.add_line_to_ss()

        sbs2b.Add(self.ss.ss, 1, wx.EXPAND)

        tab2.SetSizerAndFit(sbs2b)
        self.notebook.AddPage(tab2, "Site Mode")
        vmidsizer.Add(self.notebook, 1, wx.EXPAND)
        # self.notebook.SetSelection(1)

        p3 = panel
        sb4 = wx.StaticBox(p3, label="Match Peaks")
        sbs4 = wx.StaticBoxSizer(sb4, orient=wx.VERTICAL)

        match_iso_button = wx.Button(p3, label="Match Isolated Oligomers")
        match_iso_button.SetToolTip(wx.ToolTip("Match peaks to isolated oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.on_match_isolated, match_iso_button)
        match_all_button = wx.Button(p3, label="Match Mixed Oligomers")
        match_all_button.SetToolTip(
            wx.ToolTip("Match peaks to any possible combination of oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.on_match_all, match_all_button)
        match_sites_button = wx.Button(p3, label="Site Matching")
        match_sites_button.SetToolTip(
            wx.ToolTip(
                "Match peaks to any possible combination of rows from the Site Mode. Each site can only have one row."))
        self.Bind(wx.EVT_BUTTON, self.on_match_sites, match_sites_button)
        check_alt_button = wx.Button(p3, label="Check Alternates")
        check_alt_button.SetToolTip(
            wx.ToolTip("Check for alternative matches. Yellow indicates possible alternates within tolerance."))
        self.Bind(wx.EVT_BUTTON, self.on_check_for_alt_match, check_alt_button)
        self.matchlistbox = MatchListCtrl(self, p3, size=(500, 200), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(match_iso_button, 1, wx.EXPAND)
        hbox2.Add(match_all_button, 1, wx.EXPAND)
        hbox2.Add(match_sites_button, 1, wx.EXPAND)
        hbox2.Add(check_alt_button, 1, wx.EXPAND)
        sbs4.Add(hbox2, 0, wx.EXPAND)

        textbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmatcherror = wx.TextCtrl(p3, value=str(self.config.matchtolerance))
        textbox2.Add(wx.StaticText(p3, label="Error Tolerance for Matching (Da):"), 0, wx.ALIGN_CENTER_VERTICAL)
        textbox2.Add(self.ctlmatcherror, 0)
        sbs4.Add(textbox2, 0, wx.EXPAND)

        sbs4.Add(self.matchlistbox, 1, wx.EXPAND)
        vmidsizer.Add(sbs4, 1, wx.EXPAND)

        hbox.Add(vmidsizer, 1, wx.EXPAND)

        p4 = panel
        sb5 = wx.StaticBox(p4, label='Common Masses')
        sbs5 = wx.StaticBoxSizer(sb5, orient=wx.VERTICAL)

        importbutton2 = wx.Button(p4, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_load_common_masses, importbutton2)

        savecommonbutton = wx.Button(p4, label="Save Common Masses")
        self.Bind(wx.EVT_BUTTON, self.on_save_common_masses, savecommonbutton)

        addbutton4 = wx.Button(p4, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_new_common_mass, addbutton4)

        sbs5.Add(importbutton2, 0, wx.EXPAND)
        sbs5.Add(savecommonbutton, 0, wx.EXPAND)
        sbs5.Add(addbutton4, 0, wx.EXPAND)
        self.commonmassespanel = CommonMasses(self, p4, size=(400, 500), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)

        sbs5.Add(wx.StaticText(p4, label="Common Masses List"))
        sbs5.Add(self.commonmassespanel, 1, wx.EXPAND)

        hbox.Add(sbs5, 0, wx.EXPAND)

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

        self.masslistbox.populate(self.config.masslist)
        self.oligomerlistbox.populate(self.config.oligomerlist)

        defaultmatchlist = self.config.matchlist
        if len(defaultmatchlist) == 4:
            self.matchlistbox.populate(defaultmatchlist[0], defaultmatchlist[1], defaultmatchlist[2],
                                       defaultmatchlist[3])
            self.matchlist = defaultmatchlist
        try:
            self.load_common_masses(self.config.masstablefile)
        except:
            print("Unable to load common masses")
        # self.load_common_masses(self.config.masstablefile)
        self.CenterOnParent()
        self.Bind(wx.EVT_CLOSE, self.on_close_cancel)

    def get_site_df(self):
        try:
            self.basemass = float(self.ctlsitebasemass.GetValue())
        except Exception as e:
            self.basemass = 0
        self.ss.ss.remove_empty("Mass")
        self.sdf = self.ss.ss.get_df()
        typedict = {i: float for i in self.ss.ss.get_col_headers()}
        typedict["Name"] = str

        self.sdf = self.sdf.astype(typedict)
        print("Site Data Frame:", self.sdf)

    def get_site_masses_names(self):
        self.get_site_df()
        indexes, masses, probs, names = matchtools.get_sitematch_list(self.sdf, sites=self.get_site_headers(),
                                                                      masscolumn="Mass", namecolumn="Name", sort=None)
        oligomasslist = masses + self.basemass
        return oligomasslist, indexes, names

    def get_site_headers(self):
        col_labels = self.ss.ss.get_col_headers()
        sites = []
        for s in col_labels:
            if "Site" in s:
                sites.append(s)
        return sites

    def on_add_col(self, event=None, newcolindex=None):
        if newcolindex is None:
            newcolindex = self.ss.ss.GetNumberCols()
            self.ss.ss.InsertCols(newcolindex)
        col_labels = self.ss.ss.get_col_headers()
        maxsite = 1
        for i, s in enumerate(col_labels):
            if "Site" in col_labels[i - 1]:
                sitenumber = int(col_labels[i - 1].split()[1])
                if sitenumber > maxsite:
                    maxsite = sitenumber
        newsite = str(maxsite + 1)
        print(newsite, newcolindex)
        self.ss.ss.SetColLabelValue(newcolindex, "Site " + newsite)

    def on_add_to_ss(self, e=None):
        self.notebook.SetSelection(1)
        outdata = self.commonmassespanel.get_selection()
        sites = self.get_site_headers()
        sitedict = {sites[i]: str(1) for i in range(len(sites))}
        for o in outdata:
            print(o)
            row = self.ss.ss.add_line(Name=o[0], Mass=o[1], **sitedict)
            # self.ss.ss.set_row_values(row, **sitedict)

    def add_line_to_ss(self, name="Free", mass=0):
        sites = self.get_site_headers()
        sitedict = {sites[i]: str(1) for i in range(len(sites))}
        row = self.ss.ss.add_line(Name=name, Mass=str(mass), **sitedict)

    def on_common_to_oligo(self, e):
        self.notebook.SetSelection(0)
        outdata = self.commonmassespanel.get_selection()
        print(outdata)
        for o in outdata:
            print(o)
            self.oligomerlistbox.add_line(0, o[1], 0, 1, o[0])

    def on_add_new_common_mass(self, e):
        self.commonmassespanel.add_line()

    def add_to_common_masses(self, name, mass, type):
        self.commonmassespanel.add_line(name, mass, type)

    def on_save_common_masses(self, e):
        outdata = self.commonmassespanel.get_list()
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
        self.commonmassespanel.populate(self.commonmasses[:, 0], self.commonmasses[:, 1], self.commonmasses[:, 2])

    def on_simulate(self, e):
        """
        Replaces the peaks in self.pks with the masses defined in the self.masslistbox.
        The intensity information comes from interpolating the mass distribution in self.massdat.
        It flags self.pks.changes as 1, which is picked up downstream to say that these new peaks should be replotted.
        :param e: Unused event
        :return: None
        """
        newmasslist = self.masslistbox.get_list()
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
        self.notebook.SetSelection(0)
        self.on_match(type="isolated")

    def on_match_sites(self, e=None):
        """
        Match the peaks in self.pks to all possible combination of oligomers in the oligomerlist.
        Uses ud.make_all_matches function.
        Populated the matchlistbox with the results of the matching.
        :param e: Unused event
        :return: None
        """
        self.notebook.SetSelection(1)
        self.on_match(type="site")

    def on_match_all(self, e=None):
        """
        Match the peaks in self.pks to all possible combination of oligomers in the oligomerlist.
        Uses ud.make_all_matches function.
        Populated the matchlistbox with the results of the matching.
        :param e: Unused event
        :return: None
        """
        self.notebook.SetSelection(0)
        self.on_match(type="all")

    def on_match(self, type="all"):
        self.matchtype = type
        self.oligos = self.oligomerlistbox.get_list()
        self.oligos = np.array(self.oligos)
        self.read_tolerance()

        if type == "all":
            self.oligomasslist, self.oligonames = ud.make_all_matches(self.oligos)
        elif type == "isolated":
            self.oligomasslist, self.oligonames = ud.make_isolated_match(self.oligos)
        elif type == "site":
            self.oligomasslist, self.oligoindexes, self.oligonames = self.get_site_masses_names()
        else:
            print("Match type not recognized")
            return
        if ud.isempty(self.oligomasslist):
            print("ERROR: Need to specify the Potential Oligomers")
            return
        if type != "site":
            self.matchlist = ud.match(self.pks, self.oligomasslist, self.oligonames, self.oligos,
                                      tolerance=self.tolerance)
        else:
            self.matchlist = matchtools.site_match(self.pks, self.oligomasslist, self.oligoindexes, self.oligonames,
                                                   self.get_site_headers(), tolerance=self.tolerance)
        self.matchlistbox.populate(self.matchlist[0], self.matchlist[1], self.matchlist[2], self.matchlist[3])
        self.matchlistbox.mode = 0
        self.diffmatrix = None

    def on_check_for_alt_match(self, e=None):
        self.matchlistbox.populate(self.matchlist[0], self.matchlist[1], self.matchlist[2], self.matchlist[3])
        self.matchlistbox.mode = 0
        self.read_tolerance()
        peakmasses = [p.mass for p in self.pks.peaks]
        if self.oligomasslist is None:
            self.on_match()
        self.diffmatrix = np.subtract.outer(peakmasses, self.oligomasslist)
        matchcount = np.sum(np.abs(self.diffmatrix) < self.tolerance, axis=1)
        for i, c in enumerate(matchcount):
            if c == 0:
                color = [1, 0, 0]
            elif c == 1:
                color = [0, 1, 0]
            elif c > 1:
                color = [1, 1, 0]
            else:
                color = [0, 0, 1]
            self.matchlistbox.color_item(i, color)

    def on_view_alt(self, index=None):
        if self.matchlistbox.mode == 0:
            self.read_tolerance()
            if self.diffmatrix is None:
                self.on_check_for_alt_match()
            if index is None:
                index = 0

            self.peakindex = index
            self.peakmass = self.pks.peaks[index].mass
            b1 = np.abs(self.diffmatrix[index]) < self.tolerance
            l = np.sum(b1)
            peakmass = [self.peakmass for i in range(l)]
            self.diffs = self.diffmatrix[index][b1]
            self.mmasses = np.array(self.oligomasslist)[b1]
            if self.matchtype != "site":
                indexes = np.array(self.oligonames)[b1]
                self.mnames = np.array(
                    [ud.index_to_oname(i, self.oligos[:, 2].astype(int), self.oligos[:, 4]) for i in indexes])
            else:
                indexes = self.oligoindexes[b1]
                self.mnames = np.array(
                    [matchtools.index_to_sname(i, self.oligonames, self.get_site_headers()) for i in indexes])

            self.matchlistbox.populate(peakmass, self.mmasses, self.diffs, self.mnames)

            for i in range(l):
                absdiff = np.abs(self.diffs)
                min = np.amin(absdiff)
                max = np.amax(absdiff)  # self.tolerance
                if max != min:
                    try:
                        r = np.sqrt(np.abs(np.abs(self.diffs[i]) - min) / (max - min))
                    except Exception:
                        r = 1
                else:
                    r = 1
                color = [r, r, 1]
                self.matchlistbox.color_item(i, color)

            self.matchlistbox.mode = 1

        elif self.matchlistbox.mode == 1:
            selected_diff = self.diffs[index]
            selected_name = self.mnames[index]
            selected_match = self.mmasses[index]
            self.matchlist[1][self.peakindex] = selected_match
            self.matchlist[2][self.peakindex] = selected_diff
            self.matchlist[3][self.peakindex] = selected_name

            for p in self.pks.peaks:
                if p.mass == self.peakmass:
                    p.label = selected_name
                    p.match = selected_match
                    p.matcherror = selected_diff

            self.matchlistbox.populate(self.matchlist[0], self.matchlist[1], self.matchlist[2], self.matchlist[3])

            self.matchlistbox.mode = 0

    def read_tolerance(self):
        try:
            self.tolerance = float(self.ctlmatcherror.GetValue())
        except ValueError:
            self.tolerance = None

    def set_tolerance(self, tol=1000):
        self.ctlmatcherror.SetValue(str(tol))

    def on_close(self, e):
        """
        Close the dialog and apply the changes.

        Sets self.config.masslist to the defined mass list.
        Sets self.config.oligomerlist to the defined oligomers.
        Sets self.config.matchlist to the matched values
        :param e: Unused event
        :return: None
        """
        self.read_tolerance()
        self.config.matchtolerance = self.tolerance

        newmasslist = self.masslistbox.get_list()
        # print newmasslist
        if not ud.isempty(newmasslist):
            newmasslist = np.array(newmasslist)
            newmasslist = newmasslist[newmasslist > 0]
            # print newmasslist
            if len(newmasslist) > 0:
                self.config.masslist = newmasslist
            else:
                self.config.masslist = []
        oligos = self.oligomerlistbox.get_list()
        if not ud.isempty(oligos):
            oligos = np.array(oligos)
            oligoshort = oligos[:, :2]
            oligoshort = oligoshort.astype(float)
            oligos = oligos[np.any([oligoshort != 0], axis=2)[0], :]
            # oligos=oligos[::-1]
        self.config.oligomerlist = oligos

        matchlist = np.transpose(self.matchlistbox.get_list())
        if not ud.isempty(matchlist):
            self.config.matchlist = matchlist

        self.Destroy()
        try:
            self.EndModal(wx.ID_OK)
        except Exception as e:
            pass

    def on_close_cancel(self, e):
        """
        Close the dialog but do not modify any of the values.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(wx.ID_CANCEL)

    def pop_from_peaks(self, e):
        """
        Populate the mass list from the peak masses in self.pks
        :param e: Unused event
        :return: None
        """
        try:
            self.masslistbox.populate(self.pks.masses)
        except AttributeError:
            print("Pick peaks first")

    def pop_oligo_iso(self, e):
        """
        Populate the mass list with isolated oligomer rows.
        :param e: Unused event
        :return: None
        """
        self.notebook.SetSelection(0)
        oligos = self.oligomerlistbox.get_list()
        oligomasslist, oligonames = ud.make_isolated_match(oligos)
        oligomasslist = np.unique(oligomasslist)
        self.masslistbox.populate(oligomasslist)

    def pop_oligo_all(self, e):
        """
        Populates the mass list with all possible oligomers.
        :param e: Unused event
        :return: None
        """
        self.notebook.SetSelection(0)
        oligos = self.oligomerlistbox.get_list()
        oligomasslist, oligonames = ud.make_all_matches(oligos)
        oligomasslist = np.unique(oligomasslist)
        self.masslistbox.populate(oligomasslist)

    def pop_oligo_sites(self, e):
        """
        Populates the mass list with all possible oligomers.
        :param e: Unused event
        :return: None
        """
        self.notebook.SetSelection(1)
        self.oligomasslist, self.oligoindexes, self.oligonames = self.get_site_masses_names()
        self.masslistbox.populate(np.unique(self.oligomasslist))

    def on_clear_masslist(self, e):
        """
        Clears the mass list.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.clear()

    def on_clear_oligolist(self, e):
        """
        Clears the oligomer list.
        :param e: Unused Event
        :return: None
        """
        self.oligomerlistbox.clear()

    def on_clear_sites(self, e):
        """
        Clears the site list.
        :param e: Unused Event
        :return: None
        """
        self.ss.ss.clear_all()

    def on_add_mass(self, e):
        """
        Adds a blank line to the mass list.
        :param e: Unused Event
        :return: None
        """
        self.masslistbox.add_line()

    def on_add_oligomer(self, e):
        """
        Adds a blank line to the oligomer list.
        :param e: Unused event
        :return: None
        """
        self.oligomerlistbox.add_line()

    def on_add_site(self, e):
        """
        Adds a blank line to the site list.
        :param e: Unused event
        :return: None
        """
        self.ss.ss.add_line()

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
                self.masslistbox.populate(importmass)

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
            self.oligomerlistbox.populate(importolig)

    def on_import_sites(self, e=None):
        """
        Open a file dialog to import sites files.
        :param e: Unused event
        :return: None
        """
        path = FileDialogs.open_file_dialog("Open Sites File", file_types="*.*")
        print("Opening Path:", path)
        self.import_sites(path)

    def import_sites(self, path):
        if path is not None:
            df = spreadsheet.file_to_df(path)
            self.ss.ss.set_df(df)


# Main App Execution
if __name__ == "__main__":
    from unidec import engine

    eng = engine.UniDec()
    eng.open_file("C:\\Data\\NolanWashU\\20220429-S100B-Oligomer-1.txt")
    eng.unidec_imports()
    eng.pick_peaks()
    app = wx.App(False)
    frame = MassSelection(None)
    frame.init_dialog(eng.config, eng.pks, massdat=None)
    frame.set_tolerance(10)
    frame.on_match_all()
    frame.on_check_for_alt_match()
    frame.on_view_alt()
    frame.ShowModal()
    app.MainLoop()

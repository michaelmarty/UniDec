import sys
import string
import wx.lib.mixins.listctrl as listmix
import wx
import numpy as np
from scipy.interpolate import interp1d

from unidec_modules.AutocorrWindow import AutocorrWindow
from unidec_modules.isolated_packages import FileDialogs
import unidec_modules.unidectools as ud

'''
Module for window defining the oligomers, the expected masses, and for matching peaks to the defined oligomers.
'''

# TODO: Inherit a subclass with basic list control functions like delete


class MassListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Mass (Da)")
        self.SetColumnWidth(0, width=190)  # , wx.LIST_AUTOSIZE)

        self.popupID1 = wx.NewId()
        self.Bind(wx.EVT_MENU, self.on_masslist_delete, id=self.popupID1)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click_masslist, self)

    def populate(self, listctrldata):
        self.DeleteAllItems()
        for i in range(0, len(listctrldata)):
            index = self.InsertStringItem(sys.maxint, str(listctrldata[i]))
            self.SetItemData(index, i)

    def clear(self):
        self.DeleteAllItems()

    def add_line(self):
        self.InsertStringItem(sys.maxint, str(0))

    def get_list(self):
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            list_out.append(float(self.GetItemText(i)))
        return list_out

    def on_right_click_masslist(self, event):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_masslist_delete(self, event):
        item = self.GetFirstSelected()
        num = self.GetSelectedItemCount()
        selection = [item]
        for i in range(1, num):
            item = self.GetNextSelected(item)
            selection.append(item)
        for i in range(0, num):
            self.DeleteItem(selection[num - i - 1])


class MassListCtrlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = MassListCrtl(self, wx.NewId(), size=(210, 380), style=wx.LC_REPORT)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class OligomerListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
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

        self.popupID2 = wx.NewId()
        self.Bind(wx.EVT_MENU, self.on_oligo_delete, id=self.popupID2)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click_oligolist, self)

    def clear(self):
        self.DeleteAllItems()
        self.index = 0

    def add_line(self):
        index = self.InsertStringItem(sys.maxint, str(0))
        self.SetStringItem(index, 1, str(0))
        self.SetStringItem(index, 2, str(0))
        self.SetStringItem(index, 3, str(1))
        self.SetStringItem(index, 4, string.uppercase[self.index])
        self.index += 1

    def populate(self, data):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data)):
            index = self.InsertStringItem(sys.maxint, str(data[i][0]))
            try:
                self.SetStringItem(index, 1, str(data[i][1]))
                self.SetStringItem(index, 2, str(data[i][2]))
                self.SetStringItem(index, 3, str(data[i][3]))
            except (ValueError, IndexError):
                self.SetStringItem(index, 1, "")
                self.SetStringItem(index, 2, "")
                self.SetStringItem(index, 3, "")
            try:
                self.SetStringItem(index, 4, str(data[i][4]))
            except (ValueError, IndexError):
                self.SetStringItem(index, 4, "")
                # self.SetItemData(index, i)

    def get_list(self):
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       int(self.GetItem(i, col=2).GetText()), int(self.GetItem(i, col=3).GetText()),
                       self.GetItem(i, col=4).GetText()]
            list_out.append(sublist)
        return list_out

    def on_right_click_oligolist(self, event):
        menu = wx.Menu()
        menu.Append(self.popupID2, "Delete")
        self.PopupMenu(menu)
        menu.Destroy()

    def on_oligo_delete(self, event):
        item = self.GetFirstSelected()
        num = self.GetSelectedItemCount()
        selection = [item]
        for i in range(1, num):
            item = self.GetNextSelected(item)
            selection.append(item)
        for i in range(0, num):
            self.DeleteItem(selection[num - i - 1])


class OligomerListCrtlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = OligomerListCtrl(self, wx.NewId(), size=(500, 200), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class MatchListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
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
        self.DeleteAllItems()
        self.index = 0

    def populate(self, data1, data2, data3, data4):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data1)):
            index = self.InsertStringItem(sys.maxint, str(data1[i]))
            self.SetStringItem(index, 1, str(data2[i]))
            self.SetStringItem(index, 2, str(data3[i]))
            self.SetStringItem(index, 3, str(data4[i]))
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
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = MatchListCrtl(self, wx.NewId(), size=(500, 200), style=wx.LC_REPORT | wx.LC_SORT_ASCENDING)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


class MassSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((800, 700))
        self.SetTitle("Mass and Oligomer Tools")

        self.pks = None
        self.config = None
        self.massdat = None

        self.masslistbox = None
        self.oligomerlistbox = None
        self.matchlistbox = None

    def init_dialog(self, config, pks, massdat=None):
        # massbins = config.massbins
        self.massdat = massdat
        self.config = config
        self.pks = pks

        panel = wx.Panel(self)
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
        text = wx.StaticText(panel,
                             label="  For i from Min # to Max #:\n      Mass(i)=Base Offset + Monomer Mass * i ")
        font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        text.SetFont(font)
        hbox3.Add(text)
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
        match_all_button = wx.Button(panel, label="Matched to Mixed Oligomers")
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

        self.CenterOnParent()

    def on_simulate(self, e):
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
        if not ud.isempty(self.massdat):
            dlg = AutocorrWindow(self)
            dlg.initalize_dialog(self.config, self.massdat)
            dlg.ShowModal()

    def on_match_isolated(self, e):
        oligos = self.oligomerlistbox.list.get_list()
        oligomasslist, oligonames = ud.make_isolated_match(oligos)
        matchlist = ud.match(self.pks, oligomasslist, oligonames)
        self.matchlistbox.list.populate(matchlist[0], matchlist[1], matchlist[2], matchlist[3])

    def on_match_all(self, e):
        oligos = self.oligomerlistbox.list.get_list()
        oligomasslist, oligonames = ud.make_all_matches(oligos)
        matchlist = ud.match(self.pks, oligomasslist, oligonames)
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
        except Exception, e:
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
            print "Pick peaks first"

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
        ofilename = FileDialogs.open_file_dialog("Open Oligomer File (ofile)", file_types="*.*")
        if ofilename is not None:
            importolig = np.genfromtxt(ofilename, dtype='str')
            if np.shape(importolig) == (5,) or np.shape(importolig) == (4,):
                importolig = [importolig]
            self.oligomerlistbox.list.populate(importolig)

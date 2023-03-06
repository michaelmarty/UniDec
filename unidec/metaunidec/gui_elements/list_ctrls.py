import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np
import matplotlib.cm as cm

import unidec.tools as ud
from copy import deepcopy

luminance_cutoff = 135
white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0, 0, 0)


class YValueListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, config=None, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.config = config
        self.InsertColumn(0, "Index")
        self.InsertColumn(1, "Variable 1")
        self.InsertColumn(2, "Variable 2")
        self.InsertColumn(3, "Name")
        self.SetColumnWidth(0, width=75)
        self.SetColumnWidth(1, width=75)
        self.SetColumnWidth(2, width=75)
        self.SetColumnWidth(3, width=75)
        self.freeze_reorder = False
        self.parent = parent

    def populate(self, dataset, colors=None):
        self.DeleteAllItems()
        try:
            colormap = cm.get_cmap(self.config.spectracmap, dataset.len)
        except:
            print("Failed to get spectra cmap", self.config)
            colormap = cm.get_cmap('rainbow', dataset.len)
        peakcolors = colormap(np.arange(dataset.len))
        if colors is None:
            colors = peakcolors
        for i in range(0, dataset.len):
            s = dataset.spectra[i]
            index = self.InsertItem(dataset.len, str(s.index))
            try:
                self.SetItem(index, 1, str(s.var1))
            except:
                self.SetItem(index, 1, str(i))

            try:
                self.SetItem(index, 2, str(s.var2))
            except:
                self.SetItem(index, 2, str(0))
            try:
                self.SetItem(index, 3, s.name)
            except:
                self.SetItem(index, 3, "")
            self.SetItemData(index, i)
            if colors is not None:
                color = wx.Colour(int(round(colors[i][0] * 255)), int(round(colors[i][1] * 255)),
                                  int(round(colors[i][2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

                luminance = ud.get_luminance(color, type=2)
                # print(wx.Colour(colout), luminance)
                if luminance < luminance_cutoff:
                    self.SetItemTextColour(index, col=white_text)
                else:
                    self.SetItemTextColour(index, col=black_text)

                s.color = colors[i]
        self.data = dataset
        self.colors = colors
        self.rename_column(1, dataset.v1name)
        self.rename_column(2, dataset.v2name)
        self.freeze_reorder = False

    def recolor(self):
        dataset = self.data
        try:
            colormap = cm.get_cmap(self.config.spectracmap, dataset.len)
        except:
            print("Failed to get spectra cmap", self.config)
            colormap = cm.get_cmap('rainbow', dataset.len)
        peakcolors = colormap(np.arange(dataset.len))
        colors = peakcolors
        for i in range(0, dataset.len):
            s = dataset.spectra[i]
            index = int(s.index)
            # print(s.index, s.var1)
            if colors is not None:
                color = wx.Colour(int(round(colors[i][0] * 255)), int(round(colors[i][1] * 255)),
                                  int(round(colors[i][2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

                luminance = ud.get_luminance(color, type=2)
                # print(wx.Colour(colout), luminance)
                if luminance < luminance_cutoff:
                    self.SetItemTextColour(index, col=white_text)
                else:
                    self.SetItemTextColour(index, col=black_text)
                s.color = colors[i]
        self.colors = colors

    def clear_list(self):
        self.DeleteAllItems()
        self.freeze_reorder = False

    def add_line(self, var1="count", var2=0):
        if var1 == "count":
            var1 = self.GetItemCount()
        index = self.InsertItem(10000, str(self.GetItemCount()))
        self.SetItem(index, 1, str(var1))
        self.SetItem(index, 2, str(var2))
        self.SetItem(index, 3, str(""))

    def get_list(self):
        count = self.GetItemCount()
        # colormap = cm.get_cmap('rainbow', count)
        # peakcolors = colormap(np.arange(count))
        peakcolors = self.get_colors()
        list_output = []
        for i in range(0, count):
            it1 = self.GetItem(i, col=1).GetText()
            it2 = self.GetItem(i, col=2).GetText()
            try:
                it1 = float(it1)
            except:
                pass
            try:
                it2 = float(it2)
            except:
                pass

            sublist = [int(self.GetItem(i, col=0).GetText()), it1, it2, self.GetItem(i, col=3).GetText(),
                       peakcolors[i][0],
                       peakcolors[i][1], peakcolors[i][2]]
            list_output.append(sublist)

        indexes = np.array([i[0] for i in list_output])
        if np.any(indexes != np.arange(0, len(list_output))) and not self.freeze_reorder:
            list_output = self.reorder(list_output)
        return list_output

    def reorder(self, list):
        print("Reordering...")
        newlist = []
        indexes = np.array([i[0] for i in list])
        sind = np.sort(indexes)
        newspectra = []
        for i in sind:
            index = ud.nearestunsorted(indexes, i)
            s = self.data.spectra[index]
            s.index = i
            newspectra.append(s)
            newlist.append(list[index])
        self.data.spectra = newspectra
        self.repopulate()
        self.parent.pres.eng.data.export_hdf5()
        return newlist

    def repopulate(self):
        self.populate(self.data, self.colors)
        self.freeze_reorder = False
        pass

    def rename_column(self, num, text):
        col = self.GetColumn(num)
        col.SetText(text)
        self.SetColumn(num, col)

    def get_colors(self):
        colors = []
        count = self.GetItemCount()
        for i in range(0, count):
            c = self.GetItemBackgroundColour(i)
            colors.append(c)
        return colors


class ListCtrlPanel(wx.Panel):
    def __init__(self, parent, pres, size=(200, 400)):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        id_value = wx.ID_ANY
        self.selection = []
        self.pres = pres
        self.config = self.pres.eng.config
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.list = YValueListCtrl(self, id_value, config=self.config, size=size, style=wx.LC_REPORT | wx.BORDER_NONE)

        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self.list)

        self.popupID1 = wx.NewIdRef()
        self.popupID2 = wx.NewIdRef()
        self.popupID3 = wx.NewIdRef()
        self.popupID4 = wx.NewIdRef()
        self.popupID5 = wx.NewIdRef()
        self.popupID6 = wx.NewIdRef()
        self.popupID7 = wx.NewIdRef()
        self.popupID8 = wx.NewIdRef()
        self.popupID10 = wx.NewIdRef()
        self.popupID11 = wx.NewIdRef()

        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_popup_two, id=self.popupID2)
        self.Bind(wx.EVT_MENU, self.on_popup_three, id=self.popupID3)
        self.Bind(wx.EVT_MENU, self.on_popup_four, id=self.popupID4)
        self.Bind(wx.EVT_MENU, self.on_popup_five, id=self.popupID5)
        self.Bind(wx.EVT_MENU, self.on_popup_six, id=self.popupID6)
        self.Bind(wx.EVT_MENU, self.on_popup_seven, id=self.popupID7)
        self.Bind(wx.EVT_MENU, self.on_popup_eight, id=self.popupID8)
        self.Bind(wx.EVT_MENU, self.on_popup_ten, id=self.popupID10)
        self.Bind(wx.EVT_MENU, self.on_popup_eleven, id=self.popupID11)

    def on_right_click(self, event):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID4, "Ignore")
            menu.Append(self.popupID5, "Isolate")
            menu.Append(self.popupID6, "Repopulate")
            menu.AppendSeparator()
            menu2 = wx.Menu()

            menu2.Append(self.popupID10, "Autocorrelation")
            menu2.Append(self.popupID11, "FFT Window")
            menu.Append(wx.ID_ANY, "Analysis Tools", menu2)
            menu.AppendSeparator()
            menu.Append(self.popupID7, "Change Color")
            menu.Append(self.popupID2, "Make Top")
            if self.pres.chrommode:
                menu.Append(self.popupID8, "Make Selection")
            menu.Append(self.popupID3, "Fill Down Variable 2")
            menu.AppendSeparator()
            menu.Append(self.popupID1, "Delete")

            self.PopupMenu(menu)
            menu.Destroy()

    def on_popup_one(self, event):
        # Delete
        item = self.list.GetFirstSelected()
        num = self.list.GetSelectedItemCount()
        self.selection = []
        self.selection.append(item)
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(0, num):
            self.list.DeleteItem(self.selection[num - i - 1])
        self.pres.on_delete_spectrum(indexes=self.selection)

    def on_popup_two(self, event):
        item = self.list.GetFirstSelected()
        self.pres.make_top(item)

    def on_popup_eight(self, event):
        item = self.list.GetFirstSelected()
        self.pres.make_selection(item)

    def on_popup_three(self, event):
        item = self.list.GetFirstSelected()
        val = self.list.GetItem(item, col=2).GetText()
        count = self.list.GetItemCount()
        for i in range(0, count):
            self.list.SetItem(i, 2, val)

    def on_popup_four(self, event):
        item = self.list.GetFirstSelected()
        num = self.list.GetSelectedItemCount()
        self.selection = []
        self.selection.append(item)
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(0, num):
            self.list.DeleteItem(self.selection[num - i - 1])
        self.list.freeze_reorder = True
        self.pres.on_ignore(self.selection)

        if self.list.GetItemCount() < 1:
            print("Ignored everything; repopulating...")
            self.on_popup_six()

    def on_popup_five(self, event):
        item = self.list.GetFirstSelected()
        num = self.list.GetSelectedItemCount()
        tot = self.list.GetItemCount()
        self.selection = []
        self.selection.append(item)
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(tot - 1, -1, -1):
            if not np.any(np.array(self.selection) == i):
                self.list.DeleteItem(i)
        self.list.freeze_reorder = True
        self.pres.on_isolate(self.selection)

    def on_popup_six(self, event=None):
        self.list.get_list()
        self.list.repopulate()
        self.pres.on_repopulate()

    def on_popup_seven(self, event=None):
        """
        Spawns a dialog for the first selected item to select the color.
        Redraws the list control with the new colors and then triggers an EVT_DELETE_SELECTION_2.
        :param event: Unused Event
        :return: None
        """
        # Change Color
        item = self.list.GetFirstSelected()
        col = self.list.GetItemBackgroundColour(item)
        print("Color In:", col)
        col = wx.Colour(int(col[0]), int(col[1]), int(col[2]), alpha=int(col.alpha))
        col2 = wx.ColourData()
        col2.SetColour(col)
        colout = col2
        dlg = wx.ColourDialog(None, data=col2)
        if dlg.ShowModal() == wx.ID_OK:
            colout = dlg.GetColourData()
            colout = deepcopy(colout.GetColour())
            print("Color Out", colout)
        dlg.Destroy()
        self.list.SetItemBackgroundColour(item, col=colout)
        colout = colout.Get(includeAlpha=True)
        luminance = ud.get_luminance(wx.Colour(colout), type=2)
        colout = ([colout[0] / 255., colout[1] / 255., colout[2] / 255., colout[3] / 255.])

        if luminance < luminance_cutoff:
            self.list.SetItemTextColour(item, col=white_text)
        else:
            self.list.SetItemTextColour(item, col=black_text)

        self.pres.on_color_change(item, colout)

    def on_popup_ten(self, event=None):
        item = self.list.GetFirstSelected()
        self.pres.on_autocorr2(item)

    def on_popup_eleven(self, event=None):
        item = self.list.GetFirstSelected()
        self.pres.on_fft_window2(item)

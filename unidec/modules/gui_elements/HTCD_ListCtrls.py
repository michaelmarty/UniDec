import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np
import matplotlib as mpl

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
        self.InsertColumn(0, "")
        self.InsertColumn(1, "Label")
        self.InsertColumn(2, "Sum")
        # self.InsertColumn(3, "Name")
        self.SetColumnWidth(0, width=30)
        self.SetColumnWidth(1, width=175)
        self.SetColumnWidth(2, width=50)
        # self.SetColumnWidth(3, width=75)
        self.freeze_reorder = False
        self.parent = parent

    def populate(self, cc):
        self.DeleteAllItems()
        cc.len = len(cc.chromatograms)
        self.data = cc
        for i in range(0, cc.len):
            c = cc.chromatograms[i]
            if c.ignore:
                continue
            index = self.InsertItem(cc.len, str(c.index))
            try:
                self.SetItem(index, 1, str(c.label))
            except:
                self.SetItem(index, 1, "")

            try:
                self.SetItem(index, 2, str(c.sum))
            except:
                self.SetItem(index, 2, str(-1))
            '''
            try:
                self.SetItem(index, 3, s.name)
            except:
                self.SetItem(index, 3, "")'''
            self.SetItemData(index, i)
            color = c.color
            if color is not None:
                try:
                    color = mpl.colors.to_rgba(color)
                except Exception as e:
                    if type(color) is str or type(color) is np.str_:
                        color = str(color)
                        if "," in color:
                            color = np.fromstring(color[1:-1], sep=",", dtype=float)
                        elif " " in color:
                            color = np.fromstring(color[1:-1], sep=" ", dtype=float)
                        else:
                            print("Color not be recognized:", color)
                            color = [0, 0, 0, 1]
                    else:
                        print("Color not recognized:", color)
                        print(e)
                        print(type(color))
                        color = [0, 0, 0, 1]
                    c.color = color

                color = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)),
                                  int(round(color[2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

                luminance = ud.get_luminance(color, type=2)
                # print(wx.Colour(colout), luminance)
                if luminance < luminance_cutoff:
                    self.SetItemTextColour(index, col=white_text)
                else:
                    self.SetItemTextColour(index, col=black_text)

    def clear_list(self):
        self.DeleteAllItems()

    def repopulate(self, reset=True):
        if reset:
            for c in self.data.chromatograms:
                c.ignore = False
        self.populate(self.data)

    def get_data(self):
        return self.data


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
        self.popupID4 = wx.NewIdRef()
        self.popupID5 = wx.NewIdRef()
        self.popupID6 = wx.NewIdRef()
        self.popupID7 = wx.NewIdRef()
        self.popupID10 = wx.NewIdRef()
        self.popupID11 = wx.NewIdRef()

        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_popup_four, id=self.popupID4)
        self.Bind(wx.EVT_MENU, self.on_popup_five, id=self.popupID5)
        self.Bind(wx.EVT_MENU, self.on_popup_six, id=self.popupID6)
        self.Bind(wx.EVT_MENU, self.on_popup_seven, id=self.popupID7)
        self.Bind(wx.EVT_MENU, self.on_popup_ten, id=self.popupID10)
        self.Bind(wx.EVT_MENU, self.on_popup_eleven, id=self.popupID11)

    def on_right_click(self, event):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID4, "Ignore")
            menu.Append(self.popupID5, "Isolate")
            menu.Append(self.popupID6, "Repopulate")
            menu.AppendSeparator()

            menu.Append(self.popupID10, "Autocorrelation")

            menu.AppendSeparator()
            menu.Append(self.popupID7, "Change Color")

            menu.AppendSeparator()
            menu.Append(self.popupID1, "Delete")
            menu.Append(self.popupID11, "Delete All")

            self.PopupMenu(menu)
            menu.Destroy()

    def get_selected(self):
        item = self.list.GetFirstSelected()
        num = self.list.GetSelectedItemCount()
        self.selection = []
        self.selection.append(item)
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.selection.append(item)
        print("Selection:", self.selection)
        return self.selection

    def on_popup_one(self, event):
        # Delete
        self.selection = self.get_selected()
        for i in range(0, len(self.selection))[::-1]:
            self.list.data.chromatograms = np.delete(self.list.data.chromatograms, self.selection[i])
        self.list.repopulate()
        self.pres.on_ignore_repopulate()

    def on_popup_eleven(self, event):
        # Delete All
        self.list.data.chromatograms = []
        self.list.repopulate()
        self.pres.on_ignore_repopulate()

    def on_popup_four(self, event):
        self.selection = self.get_selected()
        for i in range(0, len(self.selection)):
            self.list.data.chromatograms[self.selection[i]].ignore = True

        self.list.repopulate(reset=False)
        self.pres.on_ignore_repopulate()

        if self.list.GetItemCount() < 1:
            print("Ignored everything; repopulating...")
            self.on_popup_six()

    def on_popup_five(self, event):
        self.selection = self.get_selected()
        for c in self.list.data.chromatograms:
            c.ignore = True
        for i in range(0, len(self.selection)):
            self.list.data.chromatograms[self.selection[i]].ignore = False

        self.list.repopulate(reset=False)
        self.pres.on_ignore_repopulate()

        if self.list.GetItemCount() < 1:
            print("Isolated nothing; repopulating...")
            self.on_popup_six()

    def on_popup_six(self, event=None):
        self.list.repopulate()
        self.pres.on_ignore_repopulate()

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

        self.list.data.chromatograms[item].color = colout

        self.pres.on_ignore_repopulate()

    def on_popup_ten(self, event=None):
        item = self.list.GetFirstSelected()
        self.pres.on_autocorr2(item)

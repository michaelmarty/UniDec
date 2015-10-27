import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np

from unidec_modules import unidectools as ud


class PeakListCtrlPanel(wx.Panel, listmix.ColumnSorterMixin):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)

        self.index = 0

        self.list_ctrl = wx.ListCtrl(self, pos=wx.DefaultPosition, size=(300, 1100),
                                     style=wx.LC_REPORT | wx.BORDER_SUNKEN
                                     # |wx.LC_SORT_ASCENDING
                                     )
        self.list_ctrl.InsertColumn(0, " ", width=25)
        self.list_ctrl.InsertColumn(1, "Mass (Da)", wx.LIST_FORMAT_RIGHT, width=70)
        self.list_ctrl.InsertColumn(2, "Intensity", width=65)
        self.list_ctrl.InsertColumn(3, "Area", width=45)
        self.list_ctrl.InsertColumn(4, "Name", width=80)

        listmix.ColumnSorterMixin.__init__(self, 4)
        self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.list_ctrl)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick, self.list_ctrl)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.list_ctrl, 0, wx.ALL | wx.EXPAND)
        self.SetSizer(sizer)

        self.EVT_DELETE_SELECTION_2 = wx.PyEventBinder(wx.NewEventType(), 1)
        self.EVT_CHARGE_STATE = wx.PyEventBinder(wx.NewEventType(), 1)
        self.EVT_DIFFERENCES = wx.PyEventBinder(wx.NewEventType(),1)

    def Clear(self):
        self.list_ctrl.DeleteAllItems()
        self.remove = []

    def AddData(self, pks, show="area", collab1="Mass (Da)"):
        self.list_ctrl.DeleteAllItems()
        self.pks = pks

        col = self.list_ctrl.GetColumn(3)
        col.SetText("Area")
        self.list_ctrl.SetColumn(3, col)
        self.list_ctrl.SetColumnWidth(3, 50)

        try:
            col = self.list_ctrl.GetColumn(1)
            col.SetText(collab1)
            self.list_ctrl.SetColumn(1, col)
            self.list_ctrl.SetColumnWidth(1, -2)
        except Exception, e:
            pass

        for i in xrange(0, self.pks.plen):
            p = pks.peaks[i]
            self.list_ctrl.InsertStringItem(i, p.textmarker)
            self.list_ctrl.SetStringItem(i, 1, str(p.mass))
            self.list_ctrl.SetStringItem(i, 2, "%.2f" % p.height)
            try:
                if show == "area":
                    self.list_ctrl.SetStringItem(i, 3, str(p.area))
                elif show == "integral":
                    self.list_ctrl.SetStringItem(i, 3, "%.2f" % p.integral)
                elif show == "diff":
                    self.list_ctrl.SetStringItem(i, 3, str(p.diff))
                else:
                    self.list_ctrl.SetStringItem(i, 3, "")
            except (ValueError, AttributeError, TypeError):
                self.list_ctrl.SetStringItem(i, 3, "")
            self.list_ctrl.SetStringItem(i, 4, str(p.label))
            self.list_ctrl.SetItemData(i, i)
            color = wx.Colour(round(p.color[0] * 255), round(p.color[1] * 255), round(p.color[2] * 255), alpha=255)
            self.list_ctrl.SetItemBackgroundColour(i, col=color)
        self.remove = []
        listmix.ColumnSorterMixin.__init__(self, 4)

    # ----------------------------------------------------------------------
    # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
    def GetListCtrl(self):
        return self.list_ctrl

    # ----------------------------------------------------------------------
    def OnColClick(self, event):
        # print "column clicked"
        event.Skip()

    def OnRightClick(self, event):
        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewId()
            self.popupID2 = wx.NewId()
            self.popupID3 = wx.NewId()
            self.popupID4 = wx.NewId()
            self.popupID5 = wx.NewId()
            self.popupID6 = wx.NewId()

            self.Bind(wx.EVT_MENU, self.OnPopupOne, id=self.popupID1)
            self.Bind(wx.EVT_MENU, self.OnPopupTwo, id=self.popupID2)
            self.Bind(wx.EVT_MENU, self.OnPopupThree, id=self.popupID3)
            self.Bind(wx.EVT_MENU, self.OnPopupFour, id=self.popupID4)
            self.Bind(wx.EVT_MENU, self.OnPopupFive, id=self.popupID5)
            self.Bind(wx.EVT_MENU, self.OnPopupSix, id=self.popupID6)
        menu = wx.Menu()
        menu.Append(self.popupID1, "Delete")
        menu.Append(self.popupID2, "Isolate")
        menu.Append(self.popupID3, "Repopulate")
        menu.Append(self.popupID4, "Label Charge States")
        menu.Append(self.popupID5, "Color Select")
        menu.Append(self.popupID6, "Display Differences")
        self.PopupMenu(menu)
        menu.Destroy()

    def OnItemSelected(self, event):
        self.currentItem = event.m_itemIndex
        print self.currentItem

    def OnPopupOne(self, event):
        # Delete
        item = self.list_ctrl.GetFirstSelected()
        num = self.list_ctrl.GetSelectedItemCount()
        self.selection = []
        self.selection.append(item)
        self.remove.append(float(self.list_ctrl.GetItem(item, col=1).GetText()))
        for i in range(1, num):
            item = self.list_ctrl.GetNextSelected(item)
            self.remove.append(float(self.list_ctrl.GetItem(item, col=1).GetText()))
            self.selection.append(item)
        for i in range(0, num):
            self.list_ctrl.DeleteItem(self.selection[num - i - 1])
        for p in self.pks.peaks:
            if p.mass in self.remove:
                p.ignore = 1
        newevent = wx.PyCommandEvent(self.EVT_DELETE_SELECTION_2._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

    def OnPopupTwo(self, event):
        # Isolate
        item = self.list_ctrl.GetFirstSelected()
        num = self.list_ctrl.GetSelectedItemCount()
        tot = self.list_ctrl.GetItemCount()
        self.selection = []
        self.selection.append(item)
        for i in range(1, num):
            item = self.list_ctrl.GetNextSelected(item)
            self.selection.append(item)
        self.selection = np.array(self.selection)
        for i in range(tot - 1, -1, -1):
            if not np.any(self.selection == i):
                self.remove.append(float(self.list_ctrl.GetItem(i, col=1).GetText()))
                self.list_ctrl.DeleteItem(i)
        for p in self.pks.peaks:
            if p.mass in self.remove:
                p.ignore = 1
        newevent = wx.PyCommandEvent(self.EVT_DELETE_SELECTION_2._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

    def OnPopupSix(self, event):
        # Show Differences
        item = self.list_ctrl.GetFirstSelected()
        num = self.list_ctrl.GetSelectedItemCount()
        self.selection2 = float(self.list_ctrl.GetItem(item, col=1).GetText())

        for p in self.pks.peaks:
            p.diff=p.mass-self.selection2
        self.list_ctrl.DeleteAllItems()

        self.AddData(self.pks,show="diff")

        col = self.list_ctrl.GetColumn(3)
        col.SetText("Diff.")
        self.list_ctrl.SetColumn(3, col)
        self.list_ctrl.SetColumnWidth(3, 65)


        newevent = wx.PyCommandEvent(self.EVT_DIFFERENCES._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

    def OnPopupThree(self, event):
        # Repopluate
        self.list_ctrl.DeleteAllItems()
        for p in self.pks.peaks:
            p.ignore = 0
        self.AddData(self.pks)
        self.remove = []
        listmix.ColumnSorterMixin.__init__(self, 4)

        newevent = wx.PyCommandEvent(self.EVT_DELETE_SELECTION_2._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

    def OnPopupFour(self, event):
        # Label Charge State
        item = self.list_ctrl.GetFirstSelected()
        num = self.list_ctrl.GetSelectedItemCount()
        self.selection2 = []
        self.selection2.append(float(self.list_ctrl.GetItem(item, col=1).GetText()))
        for i in range(1, num):
            item = self.list_ctrl.GetNextSelected(item)
            self.selection2.append(float(self.list_ctrl.GetItem(item, col=1).GetText()))
        newevent = wx.PyCommandEvent(self.EVT_CHARGE_STATE._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

    def OnPopupFive(self, event):
        # Change Color
        item = self.list_ctrl.GetFirstSelected()
        col = self.list_ctrl.GetItemBackgroundColour(item)
        print "Color In:", col
        col = wx.Colour(col[0], col[1], col[2], alpha=col.alpha)
        col2 = wx.ColourData()
        col2.SetColour(col)
        colout = col2
        dlg = wx.ColourDialog(None, data=col2)
        if dlg.ShowModal() == wx.ID_OK:
            colout = dlg.GetColourData()
            colout = colout.GetColour()
            print "Color Out", colout
        dlg.Destroy()
        self.list_ctrl.SetItemBackgroundColour(item, col=colout)
        peak = float(self.list_ctrl.GetItem(item, col=1).GetText())
        i = ud.nearest(self.pks.masses, peak)
        self.pks.peaks[i].color = ([colout[0] / 255., colout[1] / 255., colout[2] / 255])
        newevent = wx.PyCommandEvent(self.EVT_DELETE_SELECTION_2._getEvtType(), self.GetId())
        self.GetEventHandler().ProcessEvent(newevent)

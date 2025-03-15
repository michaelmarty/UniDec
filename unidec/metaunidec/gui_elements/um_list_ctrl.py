import os
import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np
import matplotlib as mpl


class XValueListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Marker")
        self.InsertColumn(1, "Peak Mass")
        self.InsertColumn(2, "Label")
        self.InsertColumn(3, "Color")
        self.SetColumnWidth(0, width=50)  # , wx.LIST_AUTOSIZE)
        self.SetColumnWidth(1, width=100)  # , wx.LIST_AUTOSIZE)
        self.SetColumnWidth(2, width=150)  # , wx.LIST_AUTOSIZE)
        self.SetColumnWidth(3, width=100)  # , wx.LIST_AUTOSIZE)
        self.parent = parent

    def populate(self, listctrldata, colors=None):
        self.DeleteAllItems()
        for i in range(0, len(listctrldata)):
            try:
                index = self.InsertItem(10000, listctrldata[i][0])
                self.SetItem(index, 1, listctrldata[i][1])
                try:
                    self.SetItem(index, 2, listctrldata[i][2])
                    self.SetItem(index, 3, listctrldata[i][3])
                except:
                    pass
                self.SetItemData(index, i)
            except (ValueError, TypeError):
                index = self.InsertItem(10000, listctrldata[i])

    def clear_list(self):
        self.DeleteAllItems()

    def add_line(self, val=0, marker='\u25CB'):
        try:
            index = self.InsertItem(10000, marker)
            self.SetItem(index, 1, str(int(val)))
            self.SetItem(index, 2, str(int(val)))
            self.SetItem(index, 3, "From File Above")
        except:
            index = self.InsertItem(10000, str('\u25CB'))
            self.SetItem(index, 1, str(0))
            self.SetItem(index, 2, str(0))
            self.SetItem(index, 3, "From File Above")

    def get_list(self):
        count = self.GetItemCount()
        list_output = []
        for i in range(0, count):
            sublist = [self.GetItem(i, col=0).GetText(), self.GetItem(i, col=1).GetText(),
                       self.GetItem(i, col=2).GetText(),  self.GetItem(i, col=3).GetText()]
            if sublist[0] != "":
                list_output.append(sublist)
        return list_output


class YValueListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "File Name")
        self.InsertColumn(1, "Color")
        self.InsertColumn(2, "Line Style")
        self.InsertColumn(3, "Label")
        self.SetColumnWidth(0, width=500)
        self.SetColumnWidth(1, width=50)
        self.SetColumnWidth(2, width=50)
        self.SetColumnWidth(3, width=100)

    def populate(self, listctrldata, colors=None):
        self.DeleteAllItems()
        for i in range(0, len(listctrldata)):
            index = self.InsertItem(i, str(listctrldata[i][0]))
            self.SetItem(index, 1, str(listctrldata[i][1]))
            self.SetItem(index, 2, str(listctrldata[i][2]))
            try:
                self.SetItem(index, 3, str(listctrldata[i][3]))
            except (ValueError, TypeError):
                self.SetItem(index, 3, "")
            self.SetItemData(index, i)
            if colors is not None:
                color = wx.Colour(int(round(colors[i][0] * 255)), int(round(colors[i][1] * 255)),
                                  int(round(colors[i][2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

    def clear_list(self):
        self.DeleteAllItems()

    def add_line(self, file_name="file.txt", var1=None, var2="-"):
        if var1 is None:
            var1 = "k"
        index = self.InsertItem(10000, str(file_name))
        self.SetItem(index, 1, str(var1))
        self.SetItem(index, 2, str(var2))
        self.SetItem(index, 3, str(""))

    def get_list(self):
        count = self.GetItemCount()
        # colormap = cm.get_cmap('rainbow', count)
        colormap = mpl.colormaps['rainbow'].resampled(count)
        peakcolors = colormap(np.arange(count))
        list_output = []
        for i in range(0, count):
            sublist = [str(self.GetItem(i, col=0).GetText()), str(self.GetItem(i, col=1).GetText()),
                       str(self.GetItem(i, col=2).GetText()), str(self.GetItem(i, col=3).GetText()), peakcolors[i][0],
                       peakcolors[i][1], peakcolors[i][2]]
            list_output.append(sublist)

        return list_output


class ListCtrlPanel(wx.Panel):
    def __init__(self, parent, pres=None, markers = None, list_type="X", size=(200, 300)):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        id_value = wx.NewIdRef()
        self.list_type = list_type
        self.parent = parent
        self.pres = pres
        self.markers = markers
        sizer = wx.BoxSizer(wx.VERTICAL)
        if list_type == "X":
            self.list = XValueListCtrl(self, id_value, size=size, style=wx.LC_REPORT | wx.BORDER_NONE)
        elif list_type == "Y":
            self.list = YValueListCtrl(self, id_value, size=size, style=wx.LC_REPORT | wx.BORDER_NONE)
        else:
            print("Error making ListCtrlPanel")
            exit()
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

        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_popup_two, id=self.popupID2)
        self.Bind(wx.EVT_MENU, self.on_popup_three, id=self.popupID3)
        self.Bind(wx.EVT_MENU, self.on_popup_four, id=self.popupID4)
        self.Bind(wx.EVT_MENU, self.on_popup_five, id=self.popupID5)
        self.Bind(wx.EVT_MENU, self.on_popup_six, id=self.popupID6)
        self.Bind(wx.EVT_MENU, self.on_popup_seven, id=self.popupID7)
        self.Bind(wx.EVT_MENU, self.on_popup_eight, id=self.popupID8)

    def on_right_click(self, event):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            if self.list_type == "Y":
                menu.Append(self.popupID3, "Fill Down Color")
                menu.Append(self.popupID4, "Fill Down Line Style")
                menu.Append(self.popupID5, "Fill Down Label")
                menu.AppendSeparator()
                menu.Append(self.popupID6, "Get Peak Indexes")
                menu.AppendSeparator()
                menu.Append(self.popupID8, "Open File in MetaUniDec")
                menu.AppendSeparator()
            if self.list_type == "X":
                menu.Append(self.popupID7, "Change Marker")
                menu.AppendSeparator()
            menu.Append(self.popupID1, "Delete")
            menu.Append(self.popupID2, "Delete All")
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

    def on_popup_two(self, event):
        self.list.DeleteAllItems()

    def on_popup_three(self, event):
        item = self.list.GetFirstSelected()
        val = self.list.GetItem(item, col=1).GetText()
        num = self.list.GetSelectedItemCount()
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.list.SetItem(item, 1, val)

    def on_popup_four(self, event):
        item = self.list.GetFirstSelected()
        val = self.list.GetItem(item, col=2).GetText()
        num = self.list.GetSelectedItemCount()
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.list.SetItem(item, 2, val)

    def on_popup_five(self, event):
        item = self.list.GetFirstSelected()
        val = self.list.GetItem(item, col=3).GetText()
        num = self.list.GetSelectedItemCount()
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.list.SetItem(item, 3, val)

    def on_popup_six(self, event):
        item = self.list.GetFirstSelected()
        self.pres.load_x_from_peaks(index=item)

    def on_popup_seven(self, event):
        item = self.list.GetFirstSelected()
        dlg = SelectMarker(self)
        dlg.initialize_interface(self.markers)
        self.list.SetItem(item, 0, dlg.textmarker)
        num = self.list.GetSelectedItemCount()
        for i in range(1, num):
            item = self.list.GetNextSelected(item)
            self.list.SetItem(item, 0, dlg.textmarker)

        self.pres.on_run()

    def on_popup_eight(self, event):
        item = self.list.GetFirstSelected()
        self.pres.open_file_in_meta(index=item)

class DCDropTarget(wx.FileDropTarget):
    """"""

    def __init__(self, window):
        """Constructor"""
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, filenames):
        """
        When files are dropped, either open a single file or run in batch.
        """
        if len(filenames) == 1:
            # Open a single file
            path = filenames[0]
            directory, fname = os.path.split(path)
            if os.path.splitext(fname)[1] == ".json":
                print("Loading .json file:", fname)
                self.window.load(path)
                return 0
        elif len(filenames) > 1:
            pass
        else:
            print("Error in file drop", filenames)
            return 1
        for f in filenames:
            self.window.ypanel.list.add_line(file_name=f)
        return 0


class SelectMarker(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create a dialog for setting some obscure additional parameters.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((285, 125))
        self.SetTitle("Pick the Peak Marker")

    def initialize_interface(self, mdkeys):
        """
        :return: None
        """
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        self.mdkeys = mdkeys

        sb = wx.StaticBox(pnl, label='Marker Type')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        hbox5 = wx.BoxSizer(wx.HORIZONTAL)

        for i, m in enumerate(mdkeys):
            button = wx.Button(pnl, i, m, size=(35, 35))
            hbox5.Add(button, 0)
            button.Bind(wx.EVT_BUTTON, self.on_close)

        sbs.Add(hbox5, 0)

        pnl.SetSizer(sbs)
        self.Bind(wx.EVT_CLOSE, self.on_close_cancel)
        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        self.SetSizer(vbox)
        self.ShowModal()

    def on_close(self, e):
        """
        Close the window.
        :param e:  Event
        :return: None
        """
        id = e.GetId()
        self.textmarker = self.mdkeys[id]
        self.Destroy()
        self.EndModal(wx.ID_OK)

    def on_close_cancel(self, e):
        """
        Close the window.
        :param e:  Event
        :return: None
        """
        self.Destroy()
        self.EndModal(wx.ID_CANCEL)
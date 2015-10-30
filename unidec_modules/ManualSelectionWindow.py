import sys
from copy import deepcopy
import wx.lib.mixins.listctrl as listmix
import wx
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from scipy.spatial.distance import euclidean

from unidec_modules import plot1d, plot2d
from unidec_modules.isolated_packages import FileDialogs
import unidec_modules.unidectools as ud

__author__ = 'Michael.Marty'


def closest(x, y, manlist):
    u = np.array([x, y])
    mindist = sys.maxint
    out = None
    pos = 0
    for i, l in enumerate(manlist):
        v = np.array([l[0], l[2]])
        dist = euclidean(u, v)
        if dist < mindist:
            mindist = dist
            out = l
            pos = i
    return out, pos


# NEED FIXING
def correctassignments(manlist, xdat, ydat):
    manlist2 = []
    for l in manlist:
        row = [l[0], l[0], l[2], l[2], l[4]]
        manlist2.append(row)
    manlist2 = np.array(manlist2)

    for i, x in enumerate(xdat):
        for j, y in enumerate(ydat):
            l, pos = closest(x, y, manlist)
            if abs(x - l[0]) < l[1] and abs(y - l[2]) < l[3]:
                if x < manlist2[pos, 0]:
                    manlist2[pos, 0] = x
                if x > manlist2[pos, 1]:
                    manlist2[pos, 1] = x
                if y < manlist2[pos, 2]:
                    manlist2[pos, 2] = y
                if y > manlist2[pos, 3]:
                    manlist2[pos, 3] = y
    manlist3 = []
    for l in manlist2:
        xwin = (l[1] - l[0]) / 2.
        ywin = (l[3] - l[2]) / 2.
        row = [l[0] + xwin, xwin, l[2] + ywin, ywin, l[4]]
        manlist3.append(row)
    return np.array(manlist3)


def range_overlap(a_min, a_max, b_min, b_max):
    return not ((a_min > b_max) or (b_min > a_max))


def checkoverlap(l1, l2):
    return range_overlap(l1[0], l1[1], l2[0], l2[1]) and range_overlap(l1[2], l1[3], l2[2], l2[3])


def detectoverlap(manlist):
    manlist2 = []
    for l in manlist:
        row = [l[0] - l[1], l[0] + l[1], l[2] - l[3], l[2] + l[3], l[4]]
        manlist2.append(row)
    manlist2 = np.array(manlist2)

    for l1 in manlist2:
        for l2 in manlist2:
            if np.all(np.array(l1) != np.array(l2)) and checkoverlap(l1, l2):
                print "Overlap of rectangles detected"
                print l1, l2
                return True
    return False


class ManualListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0, imflag=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        self.imflag = imflag
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        if self.imflag == 0:
            self.InsertColumn(0, "m/z value (Th)")
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Charge")
            self.SetColumnWidth(0, 100)
            self.SetColumnWidth(1, 100)
            self.SetColumnWidth(2, 100)
        else:
            self.InsertColumn(0, "m/z value (Th)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Arrival Time (ms)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(3, "+/-")
            self.InsertColumn(4, "Charge")
            self.SetColumnWidth(0, 90)
            self.SetColumnWidth(1, 50)
            self.SetColumnWidth(2, 105)
            self.SetColumnWidth(3, 50)
            self.SetColumnWidth(4, 75)

    def clear(self):
        self.DeleteAllItems()

    def add_line(self, line=None):
        if not line:
            line = [0, 0, 0, 0]
        index = self.InsertStringItem(sys.maxint, str(line[0]))
        self.SetStringItem(index, 1, str(line[1]))
        self.SetStringItem(index, 2, str(line[2]))
        if self.imflag == 1:
            self.SetStringItem(index, 3, str(line[3]))
            self.SetStringItem(index, 4, str(0))

    def populate(self, data, colors=None):
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertStringItem(sys.maxint, str(data[i][0]))
            self.SetStringItem(index, 1, str(data[i][1]))
            self.SetStringItem(index, 2, str(data[i][2]))
            if self.imflag == 1:
                self.SetStringItem(index, 3, str(data[i][3]))
                self.SetStringItem(index, 4, str(data[i][4]))
            if colors is not None:
                color = colors[i]
                color = wx.Colour(round(color[0] * 255), round(color[1] * 255), round(color[2] * 255), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

    def get_list(self):
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            if self.imflag == 0:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                           float(self.GetItem(i, col=2).GetText())]
            else:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                           float(self.GetItem(i, col=2).GetText()), float(self.GetItem(i, col=3).GetText()),
                           float(self.GetItem(i, col=4).GetText())]
            list_out.append(sublist)
        return list_out


class ManualListCtrlPanel(wx.Panel):
    def __init__(self, parent, imflag=0):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = ManualListCrtl(self, wx.NewId(), size=(300 + 100 * imflag, 150), style=wx.LC_REPORT, imflag=imflag)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self.list)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        self.selection = []
        self.popupID1 = wx.NewId()
        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)

    def on_right_click(self, e=None):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_popup_one(self, e=None):
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


class ManualSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetTitle("Manually Assign Masses")
        self.newtrunclist = []
        self.config = None
        self.defaulttrunclist = []
        self.data = []
        self.plot1 = None
        self.masslistbox = None

    def initiate_dialog(self, config, data):
        self.data = data
        self.SetSize((550 + config.imflag * 50, 600))
        self.config = config
        self.defaulttrunclist = deepcopy(self.config.manuallist)
        self.newtrunclist = deepcopy(self.defaulttrunclist)

        panel = wx.Panel(self)

        vbox = wx.BoxSizer(wx.VERTICAL)

        size = (7, 4.5)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        if self.config.imflag == 0:
            self.plot1 = plot1d.Plot1d(panel, figsize=size)
        else:
            self.plot1 = plot2d.Plot2d(panel, figsize=size)

        vbox2.Add(self.plot1, 0, wx.EXPAND)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        sb = wx.StaticBox(panel, label='Manually Set Masses')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        importbutton = wx.Button(panel, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import, importbutton)

        clearbutt = wx.Button(panel, label="Clear List")
        self.Bind(wx.EVT_BUTTON, self.on_clear, clearbutt)

        addbutton = wx.Button(panel, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add, addbutton)

        addbutton2 = wx.Button(panel, label="Add from Plot Zoom Range")
        self.Bind(wx.EVT_BUTTON, self.on_add_from_plot, addbutton2)

        plotbutton = wx.Button(panel, label="Plot Manual Assignments")
        self.Bind(wx.EVT_BUTTON, self.on_plot, plotbutton)

        sbs.Add(importbutton, 0, wx.EXPAND)
        sbs.Add(addbutton, 0, wx.EXPAND)
        sbs.Add(addbutton2, 0, wx.EXPAND)
        sbs.Add(clearbutt, 0, wx.EXPAND)
        sbs.Add(plotbutton, 0, wx.EXPAND)
        hbox.Add(sbs)

        sb2 = wx.StaticBox(panel, label='Manual List')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)
        self.masslistbox = ManualListCtrlPanel(panel, imflag=self.config.imflag)
        # sbs2.Add(wx.StaticText(panel, label="Manual List"))
        sbs2.Add(self.masslistbox)
        hbox.Add(sbs2)
        vbox2.Add(hbox, 1)
        panel.SetSizer(vbox2)

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

        self.masslistbox.list.populate(self.defaulttrunclist)
        if self.config.imflag == 0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], "", "m/z (Th)", "Normalized Intensity", "Data",
                                      self.config)
        else:
            self.plot1.contourplot(self.data, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="")
        self.CenterOnParent()

    def on_close(self, e):
        self.newtrunclist = self.masslistbox.list.get_list()
        if not ud.isempty(self.newtrunclist):
            self.newtrunclist = np.array(self.newtrunclist)
            self.newtrunclist = self.newtrunclist[np.any([self.newtrunclist != 0], axis=2)[0], :]

        if self.newtrunclist != "":
            if len(self.newtrunclist) > 0:
                self.config.manuallist = self.newtrunclist
            else:
                self.config.manuallist = []
        else:
            self.config.manuallist = []

        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        self.newtrunclist = self.defaulttrunclist
        self.Destroy()
        self.EndModal(1)

    def on_clear(self, e):
        self.masslistbox.list.clear()

    def on_add(self, e):
        self.masslistbox.list.add_line()

    def on_add_from_plot(self, e):
        x0, x1 = self.plot1.subplot1.get_xlim()
        xwin = (x1 - x0) / 2.
        if self.config.imflag == 0:
            line = [x0 + xwin, xwin, 0, 0]
        else:
            y0, y1 = self.plot1.subplot1.get_ylim()
            ywin = (y1 - y0) / 2.
            line = [x0 + xwin, xwin, y0 + ywin, ywin]
        self.masslistbox.list.add_line(line=line)
        pass

    def on_plot(self, e):
        if self.config.imflag == 0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], "", "m/z (Th)", "Normalized Intensity", "Data",
                                      self.config)
        else:
            self.plot1.contourplot(self.data, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="")
        self.newtrunclist = self.masslistbox.list.get_list()
        if not ud.isempty(self.newtrunclist):
            self.newtrunclist = np.array(self.newtrunclist)
            self.newtrunclist = self.newtrunclist[np.any([self.newtrunclist != 0], axis=2)[0], :]
            colormap = cm.get_cmap('rainbow', len(self.newtrunclist))
            xcolors = colormap(np.arange(len(self.newtrunclist)))
            if self.config.imflag == 1:
                try:
                    if detectoverlap(self.newtrunclist):
                        self.newtrunclist = correctassignments(self.newtrunclist, np.unique(self.data[:, 0]),
                                                               np.unique(self.data[:, 1]))
                except Exception, e:
                    print "Error with overlapping assignments. Try making sure regions don't intersect.", e
            for i, l in enumerate(self.newtrunclist):
                if self.config.imflag == 0:
                    y0 = np.amin(self.data[:, 1])
                    ywidth = np.amax(self.data[:, 1]) - y0
                else:
                    y0 = l[2] - l[3]
                    ywidth = l[3] * 2.
                self.plot1.subplot1.add_patch(
                    Rectangle((l[0] - l[1], y0), l[1] * 2., ywidth, alpha=0.5, facecolor=xcolors[i], edgecolor='black',
                              fill=True))
            self.plot1.repaint()
            self.masslistbox.list.populate(self.newtrunclist, colors=xcolors)
        pass

    def on_import(self, e):
        truncfilename = FileDialogs.open_file_dialog("Open File", file_types="*.*")
        if truncfilename is not None:
            importtrunc = np.loadtxt(truncfilename)
            if self.config.imflag == 0:
                if importtrunc.shape == (3,):
                    importtrunc = [importtrunc]
                if len(importtrunc[0]) == 3:
                    self.masslistbox.list.populate(importtrunc)
                    # print importtrunc
            else:
                if importtrunc.shape == (5,):
                    importtrunc = [importtrunc]
                if len(importtrunc[0]) == 5:
                    self.masslistbox.list.populate(importtrunc)
                    # print importtrunc

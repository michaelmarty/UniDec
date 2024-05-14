import sys
import wx.lib.mixins.listctrl as listmix
import wx
import numpy as np
from matplotlib.patches import Rectangle
# import matplotlib.cm as cm
import matplotlib as mpl
from scipy.spatial.distance import euclidean

from unidec.modules import PlottingWindow
from unidec.modules.isolated_packages import FileDialogs
import unidec.tools as ud

__author__ = 'Michael.Marty'

"""
Window for defining the manual assignments.
"""


def closest(x, y, manlist):
    """
    For manlist of manual assignments and an (x,y) point, finds the nearest manual assignment.
    Used in correcting overlapping rectanges.
    :param x: Float
    :param y: Float
    :param manlist: Array of manual assignment for IMMS (N x 5)
    :return: out, pos

    out is the line in manlist that is closest.
    pos is the position of the line in manlist.
    """
    u = np.array([x, y])
    mindist = sys.maxsize
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
    """
    Correct the manual assignments so that they are not overlapping in IM-MS.
    :param manlist: List of manual assignments in (N x 5) array.
    :param xdat: x-axis (m/z)
    :param ydat: y-axis (arrival time)
    :return: manlist3, a new array of values corrected to eliminate overlap.
    """
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
    """
    Checks whether two specific ranges of [a_min,a_max] and [b_min,b_max] overlap.
    :param a_min:
    :param a_max:
    :param b_min:
    :param b_max:
    :return: True if overlapping, False if not.
    """
    return not ((a_min > b_max) or (b_min > a_max))


def checkoverlap(l1, l2):
    """
    Check whether two lines of manual assignments are overlapping.
    :param l1: Line 1 (1 x 5)
    :param l2: Line 2 (1 x 5)
    :return: True if overlapping, False if not.
    """
    return range_overlap(l1[0], l1[1], l2[0], l2[1]) and range_overlap(l1[2], l1[3], l2[2], l2[3])


def detectoverlap(manlist):
    """
    For a list of IM-MS manual assignments, check if any of the rectangles are overlapping.
    :param manlist: Array of manual assignments (N x 5)
    :return: True if overlapping, False if not.
    """
    manlist2 = []
    for l in manlist:
        row = [l[0] - l[1], l[0] + l[1], l[2] - l[3], l[2] + l[3], l[4]]
        manlist2.append(row)
    manlist2 = np.array(manlist2)

    for l1 in manlist2:
        for l2 in manlist2:
            if np.all(np.array(l1) != np.array(l2)) and checkoverlap(l1, l2):
                print("Overlap of rectangles detected")
                print(l1, l2)
                return True
    return False


class ManualListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0, imflag=0):
        """
        Create a listctrl with either three or five columns depending on whether it is IM or MS mode.
        :param parent: Passed to wx.ListCtrl
        :param id_value: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :param imflag: Whether the listctrl should have 3 columns for MS (0) or 5 columns for IMMS (1)
        :return: None
        """
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
        elif self.imflag == 1:
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
        else:
            self.InsertColumn(0, "m/z value (Th)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Charge", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(3, "+/-")
            self.InsertColumn(4, "Charge")
            self.SetColumnWidth(0, 90)
            self.SetColumnWidth(1, 50)
            self.SetColumnWidth(2, 105)
            self.SetColumnWidth(3, 50)
            self.SetColumnWidth(4, 75)

    def clear(self):
        """
        Clear the list.
        :return: None
        """
        self.DeleteAllItems()

    def add_line(self, line=None):
        """
        Insert a new line into the listctrl.
        :param line: Optional list of values for the new line. Default is [0,0,0,0]
        :return: None
        """
        if not line:
            line = [0, 0, 0, 0]
        index = self.InsertItem(10000, str(line[0]))
        self.SetItem(index, 1, str(line[1]))
        self.SetItem(index, 2, str(line[2]))
        if self.imflag >= 1:
            self.SetItem(index, 3, str(line[3]))
            self.SetItem(index, 4, str(0))

    def populate(self, data, colors=None):
        """
        Add data from an array to the listctrl.
        :param data: Data array to be added
        :param colors: Optional list of background colors.
        :return: None
        """
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertItem(i, str(data[i][0]))
            self.SetItem(index, 1, str(data[i][1]))
            self.SetItem(index, 2, str(data[i][2]))
            if self.imflag >= 1:
                self.SetItem(index, 3, str(data[i][3]))
                self.SetItem(index, 4, str(data[i][4]))
            if colors is not None:
                color = colors[i]
                color = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)), int(round(color[2] * 255)),
                                  alpha=255)
                self.SetItemBackgroundColour(index, col=color)

    def get_list(self):
        """
        Return a nested list of the values in the listctrl.
        :return: Nested list of the outputs.
        """
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
        """
        Create the panel for the ManualListCtrl. Binds right click options.
        :param parent: Parent window or panel.
        :param imflag: Whether it should be IM-MS mode (1) or just MS (0).
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = ManualListCrtl(self, wx.NewIdRef(), size=(300 + 100 * np.amin([imflag, 1]), 150),
                                   style=wx.LC_REPORT,
                                   imflag=imflag)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self.list)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        self.selection = []
        self.popupID1 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)

    def on_right_click(self, e=None):
        """
        Open right click menu.
        :param e: Unused event
        :return: None
        """
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_popup_one(self, e=None):
        """
        Delete the selected items.
        :param e: Unused event
        :return: None
        """
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
        """
        Create the initial dialog.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetTitle("Manually Assign Masses")
        self.config = None
        self.data = []
        self.plot1 = None
        self.masslistbox = None
        self.tdflag = 0

    def initiate_dialog(self, config, data):
        """
        Initiate the dialog window, lay everything out, and plot the intial results
        :param config: UniDecConfig object
        :param data: Data to plot (either MS or IM-MS)
        :return: None
        """
        self.data = data

        if config.imflag == 1:
            self.tdflag = 1
        elif config.cdmsflag == 1:
            self.tdflag = 2
        else:
            self.tdflag = 0

        self.SetSize((700 + np.amin([self.tdflag, 1]) * 50, 650))
        self.config = config

        panel = wx.Panel(self)

        vbox = wx.BoxSizer(wx.VERTICAL)

        size = (7, 4)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        if self.tdflag == 0:
            self.plot1 = PlottingWindow.Plot1d(panel, figsize=size)
        else:
            self.plot1 = PlottingWindow.Plot2d(panel, figsize=size)

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
        self.masslistbox = ManualListCtrlPanel(panel, imflag=self.tdflag)
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

        self.masslistbox.list.populate(self.config.manuallist)
        self.on_plot(0)
        self.CenterOnParent()

    def on_close(self, e):
        """
        Get manuallist from listctrl, clean up, write it to self.config.manuallist, and then destroy the window.
        :param e: Unused event
        :return: None
        """
        manuallist = self.masslistbox.list.get_list()
        if not ud.isempty(manuallist):
            manuallist = np.array(manuallist)
            manuallist = manuallist[np.any([manuallist != 0], axis=2)[0], :]

        if len(manuallist) > 0:
            self.config.manuallist = manuallist
        else:
            self.config.manuallist = []

        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Destroy the dialog without updating the manual assignment list in UniDecConfig.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def on_clear(self, e):
        """
        Clear the listctrl.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.clear()

    def on_add(self, e):
        """
        Add a new blank line to the listctrl.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.add_line()

    def on_add_from_plot(self, e):
        """
        Grab the limits of the plot and make that the limits of the manual assignment.
        :param e: Unused event
        :return: None
        """
        x0, x1 = self.plot1.subplot1.get_xlim()
        xwin = (x1 - x0) / 2.
        if self.tdflag == 0:
            line = [x0 + xwin, xwin, 0, 0]
        else:
            y0, y1 = self.plot1.subplot1.get_ylim()
            ywin = (y1 - y0) / 2.
            line = [x0 + xwin, xwin, y0 + ywin, ywin]
        self.masslistbox.list.add_line(line=line)
        pass

    def on_plot(self, e):
        """
        Attempts to correct and plot the manual assignments.
        :param e: Unused event
        :return: None
        """
        # Make initial plot
        if self.tdflag == 0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], "", "m/z (Th)", "Normalized Intensity", "Data",
                                      self.config)
        else:
            self.plot1.contourplot(self.data, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="")

        # Get manual assignments from list
        manuallist = self.masslistbox.list.get_list()
        if not ud.isempty(manuallist):
            # Clean up list to remove empties
            manuallist = np.array(manuallist)
            manuallist = manuallist[np.any([manuallist != 0], axis=2)[0], :]
            # Make the color map
            # colormap = cm.get_cmap('rainbow', len(manuallist))
            colormap = mpl.colormaps['rainbow'].resampled(len(manuallist))
            xcolors = colormap(np.arange(len(manuallist)))
            if self.tdflag == 1:
                # Try to correct for overlapping regions in IM-MS.
                try:
                    if detectoverlap(manuallist):
                        manuallist = correctassignments(manuallist, np.unique(self.data[:, 0]),
                                                        np.unique(self.data[:, 1]))
                except Exception as e:
                    print("Error with overlapping assignments. Try making sure regions don't intersect.", e)
            # Plot the appropriate regions on the plot.
            for i, l in enumerate(manuallist):
                if self.tdflag == 0:
                    y0 = np.amin(self.data[:, 1])
                    ywidth = np.amax(self.data[:, 1]) - y0
                else:
                    y0 = l[2] - l[3]
                    ywidth = l[3] * 2.
                self.plot1.subplot1.add_patch(
                    Rectangle((l[0] - l[1], y0), l[1] * 2., ywidth, alpha=0.5, facecolor=xcolors[i], edgecolor='black',
                              fill=True))
            self.plot1.repaint()
            # Refill the list with the corrected values
            self.masslistbox.list.populate(manuallist, colors=xcolors)
        pass

    def on_import(self, e):
        """
        Open a file dialog and import a N x 3 (MS) or N x 5 (IMMS) array of manual assignment.
        Add the array to the listctrl.
        :param e: Unused event
        :return: None
        """
        import_file_name = FileDialogs.open_file_dialog("Open File", file_types="*.*")
        if import_file_name is not None:
            importtrunc = np.loadtxt(import_file_name)
            if self.tdflag == 0:
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


class SmashListCrtl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0, imflag=0):
        """
        Create a listctrl with either three or five columns depending on whether it is IM or MS mode.
        :param parent: Passed to wx.ListCtrl
        :param id_value: Passed to wx.ListCtrl
        :param pos: Passed to wx.ListCtrl
        :param size: Passed to wx.ListCtrl
        :param style: Passed to wx.ListCtrl
        :param imflag: Whether the listctrl should have 3 columns for MS (0) or 5 columns for IMMS (1)
        :return: None
        """
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        self.imflag = imflag
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        if self.imflag == 0:
            self.InsertColumn(0, "m/z value (Th)")
            self.InsertColumn(1, "+/-")
            self.SetColumnWidth(0, 100)
            self.SetColumnWidth(1, 100)
        elif self.imflag == 1:
            self.InsertColumn(0, "m/z value (Th)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Arrival Time (ms)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(3, "+/-")
            self.SetColumnWidth(0, 90)
            self.SetColumnWidth(1, 50)
            self.SetColumnWidth(2, 105)
            self.SetColumnWidth(3, 50)
        else:
            self.InsertColumn(0, "m/z value (Th)", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Charge", wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(3, "+/-")
            self.SetColumnWidth(0, 90)
            self.SetColumnWidth(1, 50)
            self.SetColumnWidth(2, 105)
            self.SetColumnWidth(3, 50)

    def clear(self):
        """
        Clear the list.
        :return: None
        """
        self.DeleteAllItems()

    def add_line(self, line=None):
        """
        Insert a new line into the listctrl.
        :param line: Optional list of values for the new line. Default is [0,0,0,0]
        :return: None
        """
        if not line:
            line = [0, 0, 0, 0]
        index = self.InsertItem(10000, str(line[0]))
        self.SetItem(index, 1, str(line[1]))
        if self.imflag >= 1:
            self.SetItem(index, 2, str(line[2]))
            self.SetItem(index, 3, str(line[3]))

    def populate(self, data, colors=None):
        """
        Add data from an array to the listctrl.
        :param data: Data array to be added
        :param colors: Optional list of background colors.
        :return: None
        """
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertItem(i, str(data[i][0]))
            self.SetItem(index, 1, str(data[i][1]))
            if self.imflag >= 1:
                self.SetItem(index, 2, str(data[i][2]))
                self.SetItem(index, 3, str(data[i][3]))
            if colors is not None:
                color = colors[i]
                color = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)), int(round(color[2] * 255)),
                                  alpha=255)
                self.SetItemBackgroundColour(index, col=color)

    def get_list(self):
        """
        Return a nested list of the values in the listctrl.
        :return: Nested list of the outputs.
        """
        count = self.GetItemCount()
        list_out = []
        for i in range(0, count):
            if self.imflag == 0:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText())]
            else:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                           float(self.GetItem(i, col=2).GetText()), float(self.GetItem(i, col=3).GetText())]
            list_out.append(sublist)
        return list_out


class SmashListCtrlPanel(wx.Panel):
    def __init__(self, parent, imflag=0):
        """
        Create the panel for the SmashListCtrl. Binds right click options.
        :param parent: Parent window or panel.
        :param imflag: Whether it should be IM-MS mode (1) or just MS (0).
        :return: None
        """
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = SmashListCrtl(self, wx.NewIdRef(), size=(300 + 100 * np.amin([imflag, 1]), 150), style=wx.LC_REPORT,
                                  imflag=imflag)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.on_right_click, self.list)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        self.selection = []
        self.popupID1 = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)

    def on_right_click(self, e=None):
        """
        Open right click menu.
        :param e: Unused event
        :return: None
        """
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def on_popup_one(self, e=None):
        """
        Delete the selected items.
        :param e: Unused event
        :return: None
        """
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


class SmashSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create the initial dialog.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetTitle("Select Noise Peaks to Smash")
        self.config = None
        self.data = []
        self.plot1 = None
        self.masslistbox = None
        self.tdflag = 0

    def initiate_dialog(self, config, data):
        """
        Initiate the dialog window, lay everything out, and plot the intial results
        :param config: UniDecConfig object
        :param data: Data to plot (either MS or IM-MS)
        :return: None
        """
        if config.imflag == 1:
            self.tdflag = 1
        elif config.cdmsflag == 1:
            self.tdflag = 2
        else:
            self.tdflag = 0

        self.data = data
        self.SetSize((700 + np.amin([self.tdflag, 1]) * 50, 650))
        self.config = config

        panel = wx.Panel(self)

        vbox = wx.BoxSizer(wx.VERTICAL)

        size = (7, 4)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        if self.tdflag == 0:
            self.plot1 = PlottingWindow.Plot1d(panel, figsize=size)
        else:
            self.plot1 = PlottingWindow.Plot2d(panel, figsize=size)

        vbox2.Add(self.plot1, 0, wx.EXPAND)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        sb = wx.StaticBox(panel, label='Define Smash Range')
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

        sb2 = wx.StaticBox(panel, label='Smash List')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)
        self.masslistbox = SmashListCtrlPanel(panel, imflag=self.tdflag)
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

        self.masslistbox.list.populate(self.config.smashlist)
        if self.data is not None:
            if len(self.data) > 1:
                self.on_plot(0)
        self.CenterOnParent()

    def on_close(self, e):
        """
        Get smashlist from listctrl, clean up, write it to self.config.smashlist, and then destroy the window.
        :param e: Unused event
        :return: None
        """
        smashlist = self.masslistbox.list.get_list()
        if not ud.isempty(smashlist):
            smashlist = np.array(smashlist)
            smashlist = smashlist[np.any([smashlist != 0], axis=2)[0], :]

        if type(smashlist) is not str:
            if len(smashlist) > 0:
                self.config.smashlist = smashlist
            else:
                self.config.smashlist = []
        else:
            self.config.smashlist = []

        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Destroy the dialog without updating the manual assignment list in UniDecConfig.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def on_clear(self, e):
        """
        Clear the listctrl.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.clear()

    def on_add(self, e):
        """
        Add a new blank line to the listctrl.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.add_line()

    def on_add_from_plot(self, e):
        """
        Grab the limits of the plot and make that the limits of the manual assignment.
        :param e: Unused event
        :return: None
        """
        x0, x1 = self.plot1.subplot1.get_xlim()
        xwin = (x1 - x0) / 2.
        if self.tdflag == 0:
            line = [x0 + xwin, xwin, 0, 0]
        else:
            y0, y1 = self.plot1.subplot1.get_ylim()
            ywin = (y1 - y0) / 2.
            line = [x0 + xwin, xwin, y0 + ywin, ywin]
        self.masslistbox.list.add_line(line=line)
        pass

    def on_plot(self, e):
        """
        Attempts to correct and plot the smash assignments.
        :param e: Unused event
        :return: None
        """
        # Make initial plot
        if self.tdflag == 0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], "", "m/z (Th)", "Normalized Intensity", "Data",
                                      self.config)
        elif self.tdflag == 1:
            self.plot1.contourplot(self.data, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)", title="")
        else:
            self.plot1.contourplot(self.data, self.config, xlab="m/z (Th)", ylab="Charge", title="")

        # Get manual assignments from list
        manuallist = self.masslistbox.list.get_list()
        if not ud.isempty(manuallist):
            # Clean up list to remove empties
            manuallist = np.array(manuallist)
            manuallist = manuallist[np.any([manuallist != 0], axis=2)[0], :]
            # Make the color map
            # colormap = cm.get_cmap('rainbow', len(manuallist))
            colormap = mpl.colormaps['rainbow'].resampled(len(manuallist))
            xcolors = colormap(np.arange(len(manuallist)))
            if self.tdflag == 1:
                # Try to correct for overlapping regions in IM-MS.
                try:
                    if detectoverlap(manuallist):
                        manuallist = correctassignments(manuallist, np.unique(self.data[:, 0]),
                                                        np.unique(self.data[:, 1]))
                except Exception as e:
                    print("Error with overlapping assignments. Try making sure regions don't intersect.", e)
            # Plot the appropriate regions on the plot.
            for i, l in enumerate(manuallist):
                if self.tdflag == 0:
                    y0 = np.amin(self.data[:, 1])
                    ywidth = np.amax(self.data[:, 1]) - y0
                else:
                    y0 = l[2] - l[3]
                    ywidth = l[3] * 2.
                self.plot1.subplot1.add_patch(Rectangle((l[0] - l[1], y0), l[1] * 2., ywidth, alpha=0.5,
                                                        facecolor=xcolors[i], edgecolor='black', fill=True))
            self.plot1.repaint(resetzoom=True)
            #self.plot1.zoomout()

            # Refill the list with the corrected values
            self.masslistbox.list.populate(manuallist, colors=xcolors)
        pass

    def on_import(self, e):
        """
        Open a file dialog and import a N x 3 (MS) or N x 5 (IMMS) array of manual assignment.
        Add the array to the listctrl.
        :param e: Unused event
        :return: None
        """
        import_file_name = FileDialogs.open_file_dialog("Open File", file_types="*.*")
        if import_file_name is not None:
            importtrunc = np.loadtxt(import_file_name)
            if self.tdflag == 0:
                if importtrunc.shape == (2,):
                    importtrunc = [importtrunc]
                if len(importtrunc[0]) == 2:
                    self.masslistbox.list.populate(importtrunc)
                    # print importtrunc
            else:
                if importtrunc.shape == (4,):
                    importtrunc = [importtrunc]
                if len(importtrunc[0]) == 4:
                    self.masslistbox.list.populate(importtrunc)
                    # print importtrunc


if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\unidec\\bin\\Example Data\\CDMS\\GroEL_CDMS_1.RAW"
    from unidec.modules import CDEng

    eng = CDEng.UniDecCD()
    eng.open_file(path)
    eng.process_data()
    eng.config.discreteplot = 1
    # eng.config.cdmsflag = 1
    # Launch SmashSelection
    app = wx.App(False)
    frame = SmashSelection(None)
    frame.initiate_dialog(eng.config, eng.data.data3)
    frame.ShowModal()
    app.MainLoop()

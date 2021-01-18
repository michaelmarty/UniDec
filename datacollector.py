import os
import sys
import json
import time
import wx
import wx.lib.mixins.listctrl as listmix
import numpy as np
import matplotlib.cm as cm
from matplotlib.pyplot import colormaps
from matplotlib import rcParams
from matplotlib.patches import Rectangle
from pubsub import pub

import multiprocessing
from unidec_modules import UniFit, Extract2D, unidecstructure, PlotAnimations, plot1d, plot2d, miscwindows, \
    MassDefects, nativez, IM_functions, DoubleDec
from unidec_modules.PlottingWindow import PlottingWindow
import unidec_modules.unidectools as ud
from unidec_modules.AutocorrWindow import AutocorrWindow
import h5py
from unidec_modules.hdf5_tools import replace_dataset, get_dataset

__author__ = 'michael.marty'

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42

luminance_cutoff = 135
white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0,0,0)

# TODO: Rewrite this code to clean it up...


class XValueListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "X Values")
        self.InsertColumn(1, "# Prot")
        self.InsertColumn(2, "# Lig")
        self.SetColumnWidth(0, width=100)  # , wx.LIST_AUTOSIZE)
        self.SetColumnWidth(1, width=50)  # , wx.LIST_AUTOSIZE)
        self.SetColumnWidth(2, width=50)  # , wx.LIST_AUTOSIZE)
        self.currentItem = 0

    def populate(self, listctrldata, colors=None):
        self.DeleteAllItems()
        listctrldata = np.array(listctrldata)
        for i in range(0, len(listctrldata)):
            try:
                index = self.InsertItem(i, str(listctrldata[i, 0]))
                self.SetItem(index, 1, str(listctrldata[i, 1]))
                self.SetItem(index, 2, str(listctrldata[i, 2]))
                self.SetItemData(index, i)
            except (ValueError, TypeError):
                index = self.InsertItem(i, str(listctrldata[i]))

            if colors is not None:
                # print listctrldata[i],colors[i]
                color = wx.Colour(int(round(colors[i][0] * 255)), int(round(colors[i][1] * 255)),
                                  int(round(colors[i][2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

                luminance = ud.get_luminance(color, type=2)
                if luminance < luminance_cutoff:
                    self.SetItemTextColour(index, col=white_text)
                else:
                    self.SetItemTextColour(index, col=black_text)

            self.SetItemData(index, i)
        self.currentItem = 0
        return self.get_maxes()

    def clear_list(self):
        self.DeleteAllItems()
        return self.get_maxes()

    def add_line(self, val=0):
        index = self.InsertItem(10000, str(val))
        self.SetItem(index, 1, str(1))
        self.SetItem(index, 2, str(self.GetItemCount() - 1))
        return self.get_maxes()

    def get_list(self):
        count = self.GetItemCount()
        list_output = []
        for i in range(0, count):
            sublist = [str(self.GetItem(i, col=0).GetText()), ud.string_to_int(self.GetItem(i, col=1).GetText(), 0),
                       ud.string_to_int(self.GetItem(i, col=2).GetText(), 0)]
            if sublist[0] != "":
                list_output.append(sublist)
        return list_output

    def get_maxes(self):
        try:
            array = np.array(self.get_list())
            maxp = np.amax([int(thing[1]) for thing in array])
            maxl = np.amax([int(thing[2]) for thing in array])
            return [maxp, maxl]
        except (ValueError, TypeError):
            return [0, 0]


class YValueListCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.TextEditMixin):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, id_value, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "File Name")
        self.InsertColumn(1, "Variable 1 (Ligand)")
        self.InsertColumn(2, "Variable 2 (Protein)")
        self.InsertColumn(3, "Charge State")
        self.SetColumnWidth(0, width=300)
        self.SetColumnWidth(1, width=100)
        self.SetColumnWidth(2, width=100)
        self.SetColumnWidth(3, width=100)
        self.parent = parent

    def populate(self, listctrldata, colors=None):
        self.DeleteAllItems()
        for i in range(0, len(listctrldata)):
            index = self.InsertItem(i, str(listctrldata[i][0]))
            self.SetItem(index, 1, str(listctrldata[i][1]))
            self.SetItem(index, 2, str(listctrldata[i][2]))
            try:
                self.SetItem(index, 3, str(listctrldata[i][3]))
            except (ValueError, TypeError):
                self.SetItem(index, 3, "All")
            self.SetItemData(index, i)
            if colors is not None:
                color = wx.Colour(int(round(colors[i][0] * 255)), int(round(colors[i][1] * 255)),
                                  int(round(colors[i][2] * 255)), alpha=255)
                self.SetItemBackgroundColour(index, col=color)

                luminance = ud.get_luminance(color, type=2)
                if luminance < luminance_cutoff:
                    self.SetItemTextColour(index, col=white_text)
                else:
                    self.SetItemTextColour(index, col=black_text)

    def clear_list(self):
        self.DeleteAllItems()

    def add_line(self, file_name="file.txt", var1="count", var2=0):
        if var1 == "count":
            var1 = self.GetItemCount()
        index = self.InsertItem(10000, str(file_name))
        self.SetItem(index, 1, str(var1))
        self.SetItem(index, 2, str(var2))
        self.SetItem(index, 3, str("All"))

    def get_list(self, fcmap='rainbow'):
        count = self.GetItemCount()
        colormap = cm.get_cmap(fcmap, count)
        peakcolors = colormap(np.arange(count))
        list_output = []
        for i in range(0, count):
            sublist = [str(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),
                       float(self.GetItem(i, col=2).GetText()), self.GetItem(i, col=3).GetText(), peakcolors[i][0],
                       peakcolors[i][1], peakcolors[i][2]]
            list_output.append(sublist)

        return list_output


class ListCtrlPanel(wx.Panel):
    def __init__(self, parent, list_type="X", size=(200, 400)):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        id_value = wx.NewIdRef()
        self.selection = []
        self.list_type = list_type
        self.parent=parent
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
        if list_type == "Y":
            self.popupID3 = wx.NewIdRef()
            self.popupID4 = wx.NewIdRef()

        self.Bind(wx.EVT_MENU, self.on_popup_one, id=self.popupID1)
        self.Bind(wx.EVT_MENU, self.on_popup_two, id=self.popupID2)
        if list_type == "Y":
            self.Bind(wx.EVT_MENU, self.on_popup_three, id=self.popupID3)
            self.Bind(wx.EVT_MENU, self.on_popup_four, id=self.popupID4)

    def on_right_click(self, event):
        if hasattr(self, "popupID1"):
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            menu.Append(self.popupID2, "Delete All")
            if self.list_type == "Y":
                menu.Append(self.popupID3, "Fill Down Variable 2")
                menu.Append(self.popupID4, "Fill Down Charge State")
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
        val = self.list.GetItem(item, col=2).GetText()
        count = self.list.GetItemCount()
        for i in range(0, count):
            self.list.SetItem(i, 2, val)

    def on_popup_four(self, event):
        item = self.list.GetFirstSelected()
        val = self.list.GetItem(item, col=3).GetText()
        count = self.list.GetItemCount()
        for i in range(0, count):
            self.list.SetItem(i, 3, val)


class NetworkFrame(PlottingWindow):
    def __init__(self, *args, **kwargs):
        PlottingWindow.__init__(self, *args, **kwargs)
        self.axes = self.figure.add_axes(self._axes)
        self.flag = True

    def clear_frame(self):
        self.figure.clear()
        self.axes = self.figure.add_axes(self._axes)
        self.repaint()


datachoices = {0: "Raw Data", 1: "Processed Data", 2: "Zero Charge Mass Spectrum", 3: "CCS (Experimental)", 4: "DoubleDec (Experimental)"}
extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position",
                  5: "Center of Mass 50%", 6: "Center of Mass 10%", 11: "Estimated Area"}
extractmethods = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 11}
extractlabels = {0: "Intensity", 1: "Intensity", 2: "Area", 3: "Mass", 4: "Mass", 5: "Mass", 6: "Mass", 11: "Area"}
modelchoices = {"Simple Single KD": "one", "Parallel KD's Chained": "parallel", "All KD's Free": "free",
                "Test for Best Model": "test", "Series KD's Chained": "series"}


# noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment,PyArgumentList
class DataCollector(wx.Frame):
    def __init__(self, parent, title, config=None, pks=None, *args, **kwargs):
        wx.Frame.__init__(self, parent, title=title)  # ,size=(200,-1))

        if "directory" in kwargs:
            self.directory = kwargs["directory"]
        else:
            self.directory = ""

        self.config = config
        self.pks = pks
        self.gridparams = None

        if self.config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            print("Using Empty Structure")
            self.config.publicationmode = 1
            if "viridis" in colormaps():
                self.config.cmap = u"viridis"
            else:
                self.config.cmap = u"jet"

        self.CreateStatusBar(2)
        self.SetStatusWidths([-1, 150])
        pub.subscribe(self.on_motion, 'newxy')

        self.filemenu = wx.Menu()

        self.menuSave = self.filemenu.Append(wx.ID_SAVE, "Save", "Save Parameters")
        self.menuLoad = self.filemenu.Append(wx.ID_ANY, "Load", "Load Parameters")
        self.filemenu.AppendSeparator()

        self.menuSaveFigPNG = self.filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                                   "Save all figures as PNG in central directory")
        self.menuSaveFigPDF = self.filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                                   "Save all figures as PDF in central directory")

        self.Bind(wx.EVT_MENU, self.on_save, self.menuSave)
        self.Bind(wx.EVT_MENU, self.on_load, self.menuLoad)
        self.Bind(wx.EVT_MENU, self.on_save_fig, self.menuSaveFigPNG)
        self.Bind(wx.EVT_MENU, self.on_save_figPDF, self.menuSaveFigPDF)
        self.toolsmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()
        self.menuOpen = self.experimentalmenu.Append(wx.ID_ANY, "Open HDF5 File", "Open HDF5 file from MetaUniDec")
        self.Bind(wx.EVT_MENU, self.on_hdf5_open, self.menuOpen)
        self.experimentalmenu.AppendSeparator()
        self.menuAnimation = self.experimentalmenu.Append(wx.ID_ANY, "Animate Spectra",
                                                          "Animation from spectra in list")
        self.Bind(wx.EVT_MENU, self.on_animate, self.menuAnimation)
        self.menuAnimation2 = self.experimentalmenu.Append(wx.ID_ANY, "Animate 2D Plots",
                                                           "Animate mass v. charge or m/z v. charge")
        self.Bind(wx.EVT_MENU, self.on_animate2, self.menuAnimation2)
        self.menu2dGrid = self.experimentalmenu.Append(wx.ID_ANY, "2D Grid Extraction",
                                                       "Extract a 2d Grid of intensities")
        self.Bind(wx.EVT_MENU, self.on_2dgrid, self.menu2dGrid)
        self.menudefect = self.experimentalmenu.Append(wx.ID_ANY, "Mass Defect Tools", "Mass Defect Analysis")
        self.Bind(wx.EVT_MENU, self.on_defect, self.menudefect)
        self.menulocalpath = self.toolsmenu.Append(wx.ID_ANY, "Convert to Local Path",
                                                   "Change file name to reflect local path for portability")
        self.Bind(wx.EVT_MENU, self.on_local_path, self.menulocalpath)
        self.menuabsolutepath = self.toolsmenu.Append(wx.ID_ANY, "Convert to Absolute Path",
                                                      "Change file name to reflect absolute path")
        self.Bind(wx.EVT_MENU, self.on_absolute_path, self.menuabsolutepath)
        self.menumsmsnorm = self.experimentalmenu.Append(wx.ID_ANY, "Normalize to MSMS",
                                                         "Normalizes mass deconvolutions to MS1 scan for variable 1 +/- variable 2")
        self.Bind(wx.EVT_MENU, self.on_msms_norm, self.menumsmsnorm)
        self.menudd = self.experimentalmenu.Append(wx.ID_ANY, "Batch DoubleDec",
                                                         "Run DoubleDec on All Files")
        self.Bind(wx.EVT_MENU, self.on_doubledec, self.menudd)
        self.menuautocorr = self.experimentalmenu.Append(wx.ID_ANY, "View Autocorrelation of Sum",
                                                         "Shows autocorelation plot of sum")
        self.Bind(wx.EVT_MENU, self.on_autocorr, self.menuautocorr)
        self.toolsmenu.AppendSeparator()
        self.menuylabel = self.toolsmenu.Append(wx.ID_ANY, "Specify Var. 1 Label", "Adds Var. 1 axis label to plot")
        self.Bind(wx.EVT_MENU, self.on_ylabel, self.menuylabel)

        self.menuplotx = self.toolsmenu.Append(wx.ID_ANY, "Plot X Ranges", "Plot X Ranges")
        self.Bind(wx.EVT_MENU, self.shade_plots, self.menuplotx)

        self.toolsmenu.AppendSeparator()

        ### CMap drop down menu
        self.scalemenu = wx.Menu()
        for i in range(0, len(self.config.cmaps)):
            idval = 500 + i
            self.scalemenu.Append(idval, self.config.cmaps[i], "", wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.menu_fcmaps, id=idval)
            if i%20==19:
                self.scalemenu.Break()
        self.toolsmenu.AppendSubMenu(self.scalemenu, 'Files Color Map')

        ### CMap drop down menu
        self.scalemenu2 = wx.Menu()
        for i in range(0, len(self.config.cmaps)):
            idval = 1500 + i
            self.scalemenu2.Append(idval, self.config.cmaps[i], "", wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.menu_xcmaps, id=idval)
            if i%20==19:
                self.scalemenu2.Break()
        self.toolsmenu.AppendSubMenu(self.scalemenu2, 'X Value Color Map')

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        self.menuBar.Append(self.experimentalmenu, "Experimental")
        self.SetMenuBar(self.menuBar)

        self.panel = wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.inputsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ypanelsizer = wx.BoxSizer(wx.VERTICAL)
        self.ypanel = ListCtrlPanel(self.panel, list_type="Y", size=(600, 500))
        self.ypanel.SetDropTarget(DCDropTarget(self))
        self.ypanelsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.addybutton = wx.Button(self.panel, label="Add Files")
        self.Bind(wx.EVT_BUTTON, self.on_add_y, self.addybutton)
        self.ypanelsizer2.Add(self.addybutton, 0, wx.EXPAND)
        self.ypanelsizer2.Add(wx.StaticText(self.panel, label="Directory:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.dirinput = wx.TextCtrl(self.panel, value="", size=(100, 20))
        self.ypanelsizer2.Add(self.dirinput, 1, wx.EXPAND)
        self.dirbutton = wx.Button(self.panel, label="...", size=(20, 20))
        self.Bind(wx.EVT_BUTTON, self.on_choose_dir, self.dirbutton)
        self.ypanelsizer2.Add(self.dirbutton, 0, wx.EXPAND)
        self.ypanelsizer.Add(self.ypanelsizer2, 0, wx.EXPAND)
        self.ypanelsizer.Add(self.ypanel, 0, wx.EXPAND)
        self.inputsizer.Add(self.ypanelsizer, 0, wx.EXPAND)

        self.xpanel = ListCtrlPanel(self.panel, list_type="X", size=(200, 500))
        self.xpanelsizer = wx.BoxSizer(wx.VERTICAL)
        self.addxbutton = wx.Button(self.panel, label="Add X Value")
        self.Bind(wx.EVT_BUTTON, self.on_add_x, self.addxbutton)
        self.xpanelsizer.Add(self.addxbutton, 0, wx.EXPAND)
        self.xpanelsizer.Add(self.xpanel, 0, wx.EXPAND)
        self.inputsizer.Add(self.xpanelsizer, 0, wx.EXPAND)

        self.plotwindow = wx.Notebook(self.panel)
        self.tab1 = wx.Panel(self.plotwindow)
        self.tab2 = wx.Panel(self.plotwindow)
        self.tab3 = wx.Panel(self.plotwindow)

        size = [5, 4]
        self.plot1 = plot1d.Plot1d(self.tab1, figsize=size)
        self.plot2d = plot2d.Plot2d(self.tab2, figsize=size)
        self.plot4 = plot1d.Plot1d(self.tab3, figsize=size)

        miscwindows.setup_tab_box(self.tab1, self.plot1)
        miscwindows.setup_tab_box(self.tab2, self.plot2d)
        miscwindows.setup_tab_box(self.tab3, self.plot4)
        self.plotwindow.AddPage(self.tab1, "1D Stack")
        self.plotwindow.AddPage(self.tab2, "2D Grid")
        self.plotwindow.AddPage(self.tab3, "Summation")

        self.plot2 = plot1d.Plot1d(self.panel)

        self.plotwindow2 = wx.Notebook(self.panel)
        self.tab12 = wx.Panel(self.plotwindow2)
        self.tab22 = wx.Panel(self.plotwindow2)

        self.plot3 = NetworkFrame(self.tab12)
        self.plot3h = plot1d.Plot1d(self.tab22)

        miscwindows.setup_tab_box(self.tab12, self.plot3)
        miscwindows.setup_tab_box(self.tab22, self.plot3h)
        self.plotwindow2.AddPage(self.tab12, "Network Model")
        self.plotwindow2.AddPage(self.tab22, "histogram")

        # self.plot3=NetworkFrame(self.panel)
        self.inputsizer.Add(self.plotwindow, 0, wx.EXPAND)
        self.sizer.Add(self.inputsizer, 1, wx.EXPAND)

        self.runsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer.Add(wx.StaticText(self.panel, label=" What to extract: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctldata = wx.ComboBox(self.panel, value="Raw Data", choices=list(datachoices.values()),
                                   style=wx.CB_READONLY)
        self.runsizer.Add(self.ctldata, 0, wx.EXPAND)

        self.runsizer.Add(wx.StaticText(self.panel, label=" Range:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlmin = wx.TextCtrl(self.panel, value="", size=(50, 20))
        self.runsizer.Add(self.ctlmin, 0, wx.ALIGN_CENTER_VERTICAL)
        self.runsizer.Add(wx.StaticText(self.panel, label=" to "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlmax = wx.TextCtrl(self.panel, value="", size=(50, 20))
        self.runsizer.Add(self.ctlmax, 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlnorm = wx.CheckBox(self.panel, label="Normalize Data")
        # self.Bind(wx.EVT_CHECKBOX, self.on_norm, self.ctlnorm)
        self.runsizer.Add(self.ctlnorm, 0, wx.EXPAND)
        self.runbutton = wx.Button(self.panel, label="Run Extraction")
        self.Bind(wx.EVT_BUTTON, self.on_run, self.runbutton)
        self.runsizer.Add(self.runbutton, 0, wx.EXPAND)
        self.ctlnorm2 = wx.CheckBox(self.panel, label="Normalize Extraction")
        self.runsizer.Add(self.ctlnorm2, 0, wx.EXPAND)
        self.runsizer.Add(wx.StaticText(self.panel, label=" How to extract: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlextract = wx.ComboBox(self.panel, value="Height", choices=list(extractchoices.values()),
                                      style=wx.CB_READONLY)
        self.runsizer.Add(self.ctlextract, 0, wx.EXPAND)
        self.runsizer.Add(wx.StaticText(self.panel, label=" Window:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlwindow = wx.TextCtrl(self.panel, value="", size=(50, 20))
        self.runsizer.Add(self.ctlwindow, 0, wx.EXPAND)

        self.runsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.runsizer2.Add(wx.StaticText(self.panel, label="Number of Proteins:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlprot = wx.TextCtrl(self.panel, value="", size=(50, 20))
        self.runsizer2.Add(self.ctlprot, 0, wx.ALIGN_CENTER_VERTICAL)

        self.runsizer2.Add(wx.StaticText(self.panel, label="Number of Ligands:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctllig = wx.TextCtrl(self.panel, value="", size=(50, 20))
        self.runsizer2.Add(self.ctllig, 0, wx.ALIGN_CENTER_VERTICAL)

        self.runsizer2.Add(wx.StaticText(self.panel, label="Number of Bootstraps:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlbootstrap = wx.TextCtrl(self.panel, value="0", size=(50, 20))
        self.runsizer2.Add(self.ctlbootstrap, 0, wx.ALIGN_CENTER_VERTICAL)

        self.kdbutton = wx.Button(self.panel, label="Fit KDs")
        self.Bind(wx.EVT_BUTTON, self.on_kd_fit, self.kdbutton)
        self.runsizer2.Add(self.kdbutton, 0, wx.EXPAND)

        choices = sorted(modelchoices.keys())
        self.ctlprotmodel = wx.ComboBox(self.panel, value=choices[0], choices=choices, style=wx.CB_READONLY)
        self.ctlligmodel = wx.ComboBox(self.panel, value=choices[0], choices=choices, style=wx.CB_READONLY)
        self.runsizer2.Add(wx.StaticText(self.panel, label=" Protein Model: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.runsizer2.Add(self.ctlprotmodel, 0, wx.EXPAND)
        self.runsizer2.Add(wx.StaticText(self.panel, label=" Ligand Model: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.runsizer2.Add(self.ctlligmodel, 0, wx.EXPAND)

        self.runsizer2.Add(wx.StaticText(self.panel, label="Binding Sites Per Protein:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlmaxsites = wx.TextCtrl(self.panel, value="None", size=(50, 20))
        self.runsizer2.Add(self.ctlmaxsites, 0, wx.ALIGN_CENTER_VERTICAL)

        # self.runsizer2.Add(wx.StaticText(self.panel,label="Binding Sites Per Protein:"),0,wx.ALIGN_CENTER_VERTICAL)
        self.ctloutliers = wx.CheckBox(self.panel, label="Remove Outliers")
        self.runsizer2.Add(self.ctloutliers, 0, wx.ALIGN_CENTER_VERTICAL)

        self.sizer.Add(self.runsizer, 0, wx.EXPAND)
        self.sizer.Add(self.runsizer2, 0, wx.EXPAND)

        self.plotsizer = wx.BoxSizer(wx.HORIZONTAL)

        # self.plot3=NetworkFrame(self.panel)
        self.plotsizer.Add(self.plot2, 0, wx.EXPAND)
        self.plotsizer.Add(self.plotwindow2, 1, wx.EXPAND)
        # self.plotsizer.Add(self.plot3,0,wx.EXPAND)
        self.sizer.Add(self.plotsizer, 0, wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)
        self.xvals = []
        self.yvals = []
        self.range = []
        self.extract = []
        self.normflag = True
        self.normflag2 = True
        self.protflag = "free"
        self.ligflag = "free"
        self.datachoice = 2
        self.numprot = 0
        self.numlig = 0
        self.bootstrap = 0
        self.window = ''
        self.maxsites = ''
        self.extractchoice = 0
        self.savename = "collection1.json"
        self.localpath = 0
        self.molig = None
        self.xcolors = []
        self.data = []
        self.grid = []
        self.var1 = []
        self.xlabel = "Mass"
        self.ylabel = ""
        self.fcmap = "rainbow"
        self.xcmap = "rainbow"
        self.offsets = []
        self.check_cmaps()
        self.hdf5_file = ""
        self.kernelfile=""
        self.filetype = 0
        self.update_set(0)
        self.Centre()
        self.Show(True)

        if "hdf_file" in kwargs:
            self.open_hdf5(kwargs["hdf_file"])

        try:
            self.load_x_from_peaks(0)
        except (ValueError, TypeError, AttributeError):
            print("Failed to load from peak list")

        if __name__ == "__main__":
            # self.directory="C:\\cprog\\Georg"
            # self.directory="C:\\Data\\AmtB_POPC"
            # self.load(os.path.join(self.directory,"AmtB_04_test.json"))
            # self.directory = "C:\\Data\\AmtB_DMPC"
            # self.load(os.path.join(self.directory, "AmtB_07.json"))
            if False:
                self.directory = "C:\\Data\\Others\\Miranda"
                self.load(os.path.join(self.directory, "collection1.json"))
                self.on_kd_fit(0)
            try:
                # testdir = "C:\Python\UniDec\unidec_src\UniDec\\x64\Release"
                # testfile = "JAW.hdf5"
                # self.open_hdf5(os.path.join(testdir, testfile))
                # self.directory="C:\\Data\\AmtB_POPC"
                # self.directory="C:\\cprog\\Shane_ND3"
                # self.directory="C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3"

                # self.load(os.path.join(self.directory,"AmtB_04_test.json"))

                # self.directory = "C:\\cprog\\Jon\\Jon Titration data\\Man9 141015"
                # self.load(os.path.join(self.directory, "collection1.json"))
                # self.on_run(0)
                # self.on_kd_fit(0)
                # self.on_animate2(0)
                pass
            except Exception as e:
                print(e)
                pass

    def check_cmaps(self):
        for i in range(len(self.config.cmaps)):
            if self.config.cmaps[i] == self.fcmap:
                self.scalemenu.Check(id=500 + i, check=True)

        for i in range(len(self.config.cmaps)):
            if self.config.cmaps[i] == self.xcmap:
                self.scalemenu2.Check(id=1500 + i, check=True)

    def menu_fcmaps(self, event=0):
        event_id = event.GetId() - 500
        print(event_id)
        self.fcmap = self.config.cmaps[event_id]
        print(self.fcmap)

    def menu_xcmaps(self, event=0):
        event_id = event.GetId() - 1500
        print(event_id)
        self.xcmap = self.config.cmaps[event_id]
        print(self.xcmap)

    def load_x_from_peaks(self, e):
        try:
            if not ud.isempty(self.pks.peaks):
                for p in self.pks.peaks:
                    maxes = self.xpanel.list.add_line(val=p.mass)
            self.ctlprot.SetValue(str(maxes[0]))
            self.ctllig.SetValue(str(maxes[1]))
        except Exception as ex:
            print("Unable to detect max # protein and ligands", ex)

    def on_hdf5_open(self, e):
        dlg = wx.FileDialog(self, "Open HDF5 File", self.directory, self.hdf5_file, "*.hdf5")
        if dlg.ShowModal() == wx.ID_OK:
            self.hdf5_file = dlg.GetPath()
            self.open_hdf5(self.hdf5_file)
        dlg.Destroy()

    def open_hdf5(self, path):
        self.topname = "ms_dataset"
        self.hdf5_file = path
        hdf = h5py.File(path)
        msdata = hdf.require_group(self.topname)
        keys = list(msdata.keys())
        self.indexes = []
        for k in keys:
            try:
                self.indexes.append(int(k))
            except:
                pass
        self.indexes = np.array(self.indexes)
        self.indexes = sorted(self.indexes)
        self.len = len(self.indexes)

        for f in self.indexes:
            msdata = hdf.get(self.topname + "/" + str(f))
            self.attrs = dict(list(msdata.attrs.items()))
            if "var1" in list(self.attrs.keys()):
                var1 = self.attrs["var1"]
            elif "collision_voltage" in list(self.attrs.keys()):
                var1 = self.attrs["collision_voltage"]
            else:
                var1 = f
            self.ypanel.list.add_line(file_name=f, var1=str(var1))

        pdataset = hdf.require_group("/peaks")
        peaks = get_dataset(pdataset, "peakdata")
        for p in peaks:
            maxes = self.xpanel.list.add_line(val=p[0])
        self.ctlprot.SetValue(str(maxes[0]))
        self.ctllig.SetValue(str(maxes[1]))

        hdf.close()
        self.update_get(0)
        self.directory = os.path.join(os.path.split(path)[0], "UniDec_Figures_and_Files")
        self.update_set(0)
        self.filetype = 1

    def on_save(self, e):
        self.update_get(e)
        try:
            exout = np.array(self.extract).tolist()
        except AttributeError:
            exout = []
        # print "Saved: ",self.gridparams
        outdict = {"x": self.xvals, "y": self.yvals, "dir": self.directory, "extractselection": self.extractchoice,
                   "window": self.window,
                   "selection": self.datachoice, "range": self.range, "normalize": self.normflag,
                   "normalize2": self.normflag2,
                   "numprot": self.numprot, "numlig": self.numlig, "bootstrap": self.bootstrap, "extract": exout,
                   "protflag": self.protflag,
                   "ligflag": self.ligflag, "maxsites": self.maxsites, "gridparams": self.gridparams,
                   "molig": self.molig, "filetype": self.filetype}

        dlg = wx.FileDialog(self, "Save Collection in JSON Format", self.directory, self.savename, "*.json")
        if dlg.ShowModal() == wx.ID_OK:
            self.savename = dlg.GetPath()
            with open(self.savename, "w+") as outfile:
                json.dump(outdict, outfile)
            print("Saved: ", self.savename)
        dlg.Destroy()

    def on_load(self, e):
        dlg = wx.FileDialog(self, "Load JSON Collection", self.directory, self.savename, "*.json")
        if dlg.ShowModal() == wx.ID_OK:
            self.savename = dlg.GetPath()
            self.load(self.savename)
        dlg.Destroy()

    def load(self, savename):
        with open(savename, "r") as outfile:
            indict = json.load(outfile)
        if "x" in indict:
            self.xvals = indict["x"]
        if "y" in indict:
            self.yvals = indict["y"]
        if "dir" in indict:
            self.directory = indict["dir"]
        if "extractselection" in indict:
            self.extractchoice = indict["extractselection"]
        if "selection" in indict:
            self.datachoice = indict["selection"]
        if "range" in indict:
            self.range = indict["range"]
        if "normalize" in indict:
            self.normflag = indict["normalize"]
        if "normalize2" in indict:
            self.normflag2 = indict["normalize2"]
        if "numprot" in indict:
            self.numprot = indict["numprot"]
        if "numlig" in indict:
            self.numlig = indict["numlig"]
        if "bootstrap" in indict:
            self.bootstrap = indict["bootstrap"]
        if "window" in indict:
            self.window = indict["window"]
        if "protflag" in indict:
            self.protflag = indict["protflag"]
        if "ligflag" in indict:
            self.ligflag = indict["ligflag"]
        if "maxsites" in indict:
            self.maxsites = indict["maxsites"]
        if "gridparams" in indict:
            self.gridparams = indict["gridparams"]
            # print "Loaded: ",self.gridparams
        if "molig" in indict:
            self.molig = indict["molig"]
        if "filetype" in indict:
            self.filetype = indict["filetype"]
        else:
            self.filetype = 0
        self.update_set(0)
        print("Loaded: ", savename)
        self.on_run(0)

    def on_add_x(self, e):
        maxes = self.xpanel.list.add_line()
        try:
            self.ctlprot.SetValue(str(maxes[0]))
            self.ctllig.SetValue(str(maxes[1]))
        except (ValueError, TypeError):
            print("Failed to autoupdate total # of protein and ligand, update manually in boxes.")

    def on_add_y(self, e):
        self.update_get(e)
        dlg = wx.FileDialog(self, "Load Files", self.directory, "", "*.*", wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetPaths()
            for f in filenames:
                self.ypanel.list.add_line(file_name=f)
            self.filetype = 0
        dlg.Destroy()
        self.localpath = 0

    def on_choose_dir(self, e):
        dlg = wx.DirDialog(None, "Choose Top Directory", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.directory = dlg.GetPath()
            self.dirinput.SetValue(self.directory)
            # print self.directory
        dlg.Destroy()

    def on_motion(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=1)
        except:
            pass

    def update_get(self, e):
        self.xvals = self.xpanel.list.get_list()
        try:
            maxes = self.xpanel.list.get_maxes()
            self.ctlprot.SetValue(str(maxes[0]))
            self.ctllig.SetValue(str(maxes[1]))
        except (ValueError, TypeError):
            print("Failed to autoupdate total # of protein and ligand, update manually in boxes.")
        self.yvals = self.ypanel.list.get_list(fcmap=self.fcmap)
        self.directory = self.dirinput.GetValue()
        self.normflag = self.ctlnorm.GetValue()
        self.normflag2 = self.ctlnorm2.GetValue()
        self.datachoice = self.ctldata.GetSelection()
        self.extractchoice = self.ctlextract.GetSelection()
        self.protflag = modelchoices[self.ctlprotmodel.GetStringSelection()]
        self.ligflag = modelchoices[self.ctlligmodel.GetStringSelection()]
        try:
            self.window = float(self.ctlwindow.GetValue())
        except ValueError:
            pass
        try:
            self.bootstrap = int(self.ctlbootstrap.GetValue())
        except ValueError:
            self.bootstrap = 1
            pass
        self.range = []
        try:
            self.range.append(float(self.ctlmin.GetValue()))
            self.range.append(float(self.ctlmax.GetValue()))
        except ValueError:
            pass
        try:
            self.numprot = int(self.ctlprot.GetValue())
            self.numlig = int(self.ctllig.GetValue())
        except ValueError:
            pass

        try:
            self.maxsites = int(self.ctlmaxsites.GetValue())
        except ValueError:
            self.maxsites = ''
            # print self.xvals
            # print self.yvals
            # print self.directory

    def update_set(self, e):
        self.ctlprotmodel.SetValue(
            next((label for label, flag in list(modelchoices.items()) if flag == self.protflag), "test"))
        self.ctlligmodel.SetValue(
            next((label for label, flag in list(modelchoices.items()) if flag == self.ligflag), "test"))
        self.dirinput.SetValue(self.directory)
        self.xpanel.list.populate(self.xvals)
        self.ypanel.list.populate(self.yvals)
        self.ctlnorm.SetValue(self.normflag)
        self.ctlnorm2.SetValue(self.normflag2)
        self.ctlprot.SetValue(str(self.numprot))
        self.ctllig.SetValue(str(self.numlig))
        self.ctlwindow.SetValue(str(self.window))
        self.ctlbootstrap.SetValue(str(self.bootstrap))
        self.ctldata.SetSelection(self.datachoice)
        self.ctlextract.SetSelection(self.extractchoice)
        self.ctlmaxsites.SetValue(str(self.maxsites))
        if not ud.isempty(self.range):
            self.ctlmin.SetValue(str(self.range[0]))
            self.ctlmax.SetValue(str(self.range[1]))


    def on_run(self, e, vals=None):
        tstart = time.perf_counter()
        self.update_get(e)
        # os.chdir(self.directory)
        self.data = []
        self.extract = []
        self.plot1.clear_plot()
        self.plot2.clear_plot()
        self.var1 = []
        self.grid = []
        if self.filetype == 1:
            hdf = h5py.File(self.hdf5_file)
        ycolors = []
        print("Directory:", self.directory)
        for k, l in enumerate(self.yvals):
            if self.filetype == 1:
                msdata = hdf.get(self.topname + "/" + str(l[0]))

            filename = l[0]
            header = os.path.splitext(filename)[0]
            if self.localpath == 1 or not os.path.isabs(filename):
                header = os.path.join(self.directory, header)
                filename = os.path.join(self.directory, filename)
            print("Loading:", filename)
            subheader = os.path.split(header)[1]
            if self.numlig == 0 and self.numprot > 1:
                print("Using Variable 2 rather than Variable 1:", l[2])
                self.var1.append(l[2])
            else:
                self.var1.append(l[1])
            ycolors.append(l[4:7])
            fcolor = np.array(l[4:7])
            if self.datachoice == 0:
                if self.filetype == 1:
                    data = get_dataset(msdata, "raw_data")
                else:
                    data = np.loadtxt(filename)
                self.xlabel = "m/z (Th)"
            elif self.datachoice == 1:
                filename = os.path.join(header + "_unidecfiles", subheader + "_input.dat")
                if self.filetype == 1:
                    data = get_dataset(msdata, "processed_data")
                else:
                    data = np.loadtxt(filename)
                self.xlabel = "m/z (Th)"
            elif self.datachoice == 2:
                self.xlabel = "Mass (Da)"
                zstate = l[3]
                filename = os.path.join(header + "_unidecfiles", subheader + "_mass.txt")
                if self.filetype == 1:
                    data = get_dataset(msdata, "mass_data")
                else:
                    data = np.loadtxt(filename)
                if not zstate == 'All':
                    try:
                        filename = os.path.join(header + "_unidecfiles", subheader + "_massgrid.bin")
                        massgrid = np.fromfile(filename, dtype=self.config.dtype)
                        configfile = os.path.join(header + "_unidecfiles", subheader + "_conf.dat")
                        f = open(configfile, 'r')
                        for line in f:
                            if line.startswith("startz"):
                                startz = int(line.split()[1])
                            if line.startswith("endz"):
                                endz = int(line.split()[1])
                        zlength = endz - startz + 1

                        massgrid = np.reshape(massgrid, (len(data), zlength))

                        if zstate[0] == "N":
                            srange = zstate[1:]
                            nzstart, nzend = srange.split(':')
                            nzstart = float(nzstart)
                            nzend = float(nzend)
                            if nzstart > nzend:
                                nzstart, nzend = nzend, nzstart
                            ztab = np.arange(startz, endz + 1)
                            mgrid, zgrid = np.meshgrid(data[:, 0], ztab, indexing='ij')

                            offsetgrid = nativez.GetOffset(mgrid, zgrid)
                            bool1 = offsetgrid >= nzstart
                            bool2 = offsetgrid <= nzend
                            bool3 = np.all([bool1, bool2], axis=0)
                            mout = np.sum(massgrid * bool3, axis=1)

                            data[:, 1] = mout
                        else:
                            zindex = int(zstate) - startz
                            data[:, 1] = massgrid[:, zindex]
                    except Exception as e:
                        print("FAILED TO EXTRACT CHARGE STATE\nUsing total for all charge states\n", e)
                        l[3] = "All"
                        try:
                            self.ypanel.list.populate(self.yvals, colors=ycolors)
                        except (ValueError, TypeError, AttributeError):
                            self.ypanel.list.populate(self.yvals)

            elif self.datachoice == 3:
                filename = os.path.join(header + "_unidecfiles", subheader + "_ccs.txt")
                if self.filetype == 1:
                    print("ERROR: HDF5 Files not supported with IMMS data")
                else:
                    data = np.loadtxt(filename)
                self.xlabel = "CCS (A^2)"

            elif self.datachoice == 4:
                self.xlabel = "Mass (Da)"
                filename = os.path.join(header + "_unidecfiles", subheader + "_massdd.txt")
                if self.filetype == 1:
                    data = get_dataset(msdata, "mass_data_dd")
                else:
                    data = np.loadtxt(filename)

            if not ud.isempty(self.range):
                bool1 = data[:, 0] >= self.range[0]
                bool2 = data[:, 0] <= self.range[1]
                bool3 = np.all([bool1, bool2], axis=0)
                data = data[bool3]
            if self.normflag:
                data[:, 1] = data[:, 1] / np.amax(data[:, 1])
            if vals is not None:
                data[:, 1] = data[:, 1] * vals[k]
            self.data.append(data)
            if not self.plot1.flag:
                self.plot1.plotrefreshtop(data[:, 0], data[:, 1], "Extracted Data", "", "", filename, None,
                                          nopaint=True,
                                          color=fcolor, test_kda=True)
                self.grid.append(data)
            else:
                self.plot1.plotadd(data[:, 0], data[:, 1], fcolor, filename)
                try:
                    self.grid.append(ud.mergedata(self.grid[0], data))
                except (ValueError, TypeError, AttributeError):
                    print("Error combining data")

            xext = []
            for i in range(0, len(self.xvals)):
                s = self.xvals[i][0]
                try:
                    window = float(self.window)
                except (ValueError, TypeError):
                    window = None

                if ";" not in s:
                    val = ud.data_extract(data, float(s), extractmethods[self.extractchoice], window=window)
                else:
                    xs = s.split(';')
                    val = 0
                    for xval in xs:
                        val += ud.data_extract(data, float(xval), extractmethods[self.extractchoice], window=window)
                xext.append(val)
            self.extract.append(xext)
        self.ypanel.list.populate(self.yvals, colors=ycolors)
        self.plot1.repaint()
        tend = time.perf_counter()
        print("Extraction time: %.2gs" % (tend - tstart))
        print("Plotting...")
        self.extract = np.array(self.extract)
        if len(self.xvals) != 0:
            if self.normflag2 == 1:
                sums = np.sum(self.extract, axis=1)
                self.extract = [self.extract[i] / sums[i] for i in range(0, len(self.yvals))]
                self.extract = np.array(self.extract)

            colormap = cm.get_cmap(self.xcmap, len(self.xvals))
            self.xcolors = colormap(np.arange(len(self.xvals)))
            self.makeplot2()
            # print np.mean(self.extract,axis=0)
            np.savetxt(os.path.join(self.directory, "extracts.txt"), self.extract)
            self.xpanel.list.populate(self.xvals, colors=self.xcolors)

        self.make_grid_plots(e)
        if self.filetype == 1:
            hdf.close()
        print("Extraction Complete")

    def shade_plots(self, e):
        self.update_get(e)
        for i in range(0, len(self.xvals)):
            color = self.xcolors[i]

            s = self.xvals[i][0]
            try:
                window = float(self.window)
            except (ValueError, TypeError):
                window = None

            if ";" not in s:
                val = float(s)
                self.make_shade_plot(val, window, color)
            else:
                xs = s.split(';')
                val = 0
                for xval in xs:
                    val += float(xval)
                    self.make_shade_plot(val,  window, color)
        self.plot1.repaint()

    def make_shade_plot(self, val, window, color):
        y0 = np.amin(self.data[0][:, 1])
        ywidth = np.amax(self.data[0][:, 1]) - y0
        if self.plot1.kdnorm == 1000:
            val = val/self.plot1.kdnorm
            window = window/self.plot1.kdnorm

        print(val, window)
        self.plot1.subplot1.add_patch(
            Rectangle((val - window, y0), window * 2., ywidth, alpha=0.5, facecolor=color, edgecolor='black',
                      fill=True))

    def make_grid_plots(self, e):
        self.grid = np.array(self.grid)
        if not ud.isempty(self.grid):
            try:
                self.plot4.plotrefreshtop(np.unique(self.grid[0, :, 0]), np.sum(self.grid[:, :, 1], axis=0), "Total",
                                          "", "Sum", "", self.config, test_kda=True)
                np.savetxt(os.path.join(self.directory, "sums.txt"),
                           np.transpose([np.unique(self.grid[0, :, 0]), np.sum(self.grid[:, :, 1], axis=0)]))
            except (ValueError, TypeError, AttributeError):
                print("Failed total plot")
                self.plot4.clear_plot()
            try:
                x, y = np.meshgrid(self.grid[0, :, 0], self.var1, indexing='ij')
                dat = np.transpose([np.ravel(x), np.ravel(y), np.ravel(np.transpose(self.grid[:, :, 1]))])
                self.plot2d.contourplot(dat, self.config, xlab=self.xlabel, ylab=self.ylabel, title="Extracted Data")
                self.on_export(e)
            except (ValueError, TypeError, AttributeError):
                print("Failed to make 2D plot")
                self.plot2d.clear_plot()

    def makeplot2(self):
        # This is a bit funny but the ylabel from this plot is actually the xlabel from the others or the intensity...
        ylabel = extractlabels[extractmethods[self.extractchoice]]
        if ylabel == "mass":
            ylabel = self.xlabel

        self.plot2.clear_plot()
        self.plot2._axes = [0.15, 0.1, 0.75, 0.8]
        for i in range(0, len(self.xvals)):
            color = self.xcolors[i]
            if not self.plot2.flag:

                self.plot2.plotrefreshtop(self.var1, self.extract[:, i], title="Extracted Data", xlabel=self.ylabel
                                          , ylabel=ylabel, color=color, test_kda=False)
                self.plot2.plotadddot(self.var1, self.extract[:, i], color, "o")
            else:
                self.plot2.plotadd(self.var1, self.extract[:, i], color)
                self.plot2.plotadddot(self.var1, self.extract[:, i], color, "o")
        if self.normflag2 == 1:
            self.plot2.subplot1.set_ylim([0, 1])
        self.plot2.repaint()

    def on_save_fig(self, e):
        name1 = os.path.join(self.directory, "Figure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print(name1)
        name2 = os.path.join(self.directory, "Figure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print(name2)
        name3 = os.path.join(self.directory, "Figure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
            print(name3)
        name1 = os.path.join(self.directory, "Figure1_2D.png")
        if self.plot2d.flag:
            self.plot2d.on_save_fig(e, name1)
            print(name1)
        name3 = os.path.join(self.directory, "Figure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name3)
            print(name3)

    def on_save_figPDF(self, e):
        name1 = os.path.join(self.directory, "Figure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print(name1)
        name2 = os.path.join(self.directory, "Figure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print(name2)
        name3 = os.path.join(self.directory, "Figure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
            print(name3)
        name1 = os.path.join(self.directory, "Figure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name1)
            print(name1)
        name1 = os.path.join(self.directory, "Figure1_2D.pdf")
        if self.plot2d.flag:
            self.plot2d.on_save_fig(e, name1)
            print(name1)

    def on_kd_fit(self, e):
        outlierflag = self.ctloutliers.GetValue()
        self.update_get(e)
        self.plot3.clear_frame()
        # self.makeplot2()
        nodelist = []
        for i in range(0, len(self.xvals)):
            nodelist.append([self.xvals[i][1], self.xvals[i][2]])
        nodelist = np.array(nodelist)
        print("Nodes: ", nodelist)
        self.yvals = np.array(self.yvals)
        self.plot2.clear_plot()
        # noinspection PyProtectedMember
        self.plot2.subplot1 = self.plot2.figure.add_axes(self.plot2._axes)
        try:
            maxsites = int(self.maxsites)
        except (ValueError, TypeError):
            maxsites = 0
        model = UniFit.KDmodel(np.transpose(self.extract), self.yvals[:, 2].astype(np.float64),
                               self.yvals[:, 1].astype(np.float64), nodelist, os.path.join(self.directory, "fits"),
                               numtotprot=self.numprot, numtotlig=self.numlig, removeoutliers=outlierflag,
                               plot1=self.plot2.subplot1,
                               plot2=self.plot3.axes, plot3=self.plot3h, bootnum=self.bootstrap, maxsites=maxsites,
                               prot=self.protflag, lig=self.ligflag, label=self.ylabel, cmap=self.xcmap)
        try:
            if self.bootstrap > 0:
                np.savetxt(os.path.join(self.directory, "fits_boots.txt"), model.randfit)
        except Exception as e:
            print(e)
        self.plot2.repaint()
        self.plot2.flag = True
        datalim = [np.amin(nodelist[:, 1]) - 0.5, np.amin(nodelist[:, 0] - 0.5), np.amax(nodelist[:, 1] + 0.5),
                   np.amax(nodelist[:, 0] + 0.5)]
        self.plot3.axes.set_xlim(datalim[0], datalim[2])
        self.plot3.axes.set_ylim(datalim[1], datalim[3])
        self.plot3.setup_zoom([self.plot3.axes], 'box', data_lims=datalim)
        self.plot3.repaint()

    def on_animate(self, e):
        self.yvals = np.array(self.yvals)

        if self.datachoice == 2:
            mode = "Mass"
        elif self.datachoice == 3:
            mode = "CCS"
        else:
            mode = "mz"

        PlotAnimations.AnimationWindow(self, self.grid, self.config, yvals=self.yvals[:, 1], pksmode=mode)

    def on_animate2(self, e):
        self.update_get(e)
        dlg = miscwindows.SingleInputDialog(self)
        dlg.initialize_interface(title="Set Compression", message="Number of x values to compress:", defaultvalue="10")
        dlg.ShowModal()
        try:
            compress = int(dlg.value)
            if compress > 1:
                print("Compressing Data by:", compress)
        except (ValueError, TypeError, AttributeError):
            print("Unrecognized compression value")
            compress = 0

        print("Loading 2D Data...")
        data2 = []
        for i, l in enumerate(self.yvals):
            filename = l[0]
            header = os.path.splitext(filename)[0]
            if self.localpath == 1:
                header = os.path.join(self.directory, header)
            subheader = os.path.split(header)[1]
            if self.datachoice < 2:
                filename = os.path.join(header + "_unidecfiles", subheader + "_grid.bin")
                file2 = os.path.join(header + "_unidecfiles", subheader + "_input.dat")
            elif self.datachoice == 2:
                filename = os.path.join(header + "_unidecfiles", subheader + "_massgrid.bin")
                file2 = os.path.join(header + "_unidecfiles", subheader + "_mass.txt")
            elif self.datachoice == 3:
                filename = os.path.join(header + "_unidecfiles", subheader + "_massccs.bin")
                file2 = os.path.join(header + "_unidecfiles", subheader + "_ccs.txt")
            else:
                file2 = filename
                print("Undefined choice for data type", self.datachoice, type(self.datachoice))
            data = np.loadtxt(file2)
            massgrid = np.fromfile(filename, dtype=self.config.dtype)
            configfile = os.path.join(header + "_unidecfiles", subheader + "_conf.dat")
            f = open(configfile, 'r')
            for line in f:
                if line.startswith("startz"):
                    startz = int(line.split()[1])
                if line.startswith("endz"):
                    endz = int(line.split()[1])
                if line.startswith("numz"):
                    numz = int(line.split()[1])
            f.close()
            try:
                zlength = endz - startz + 1
            except Exception as e1:
                try:
                    zlength = numz
                    endz = startz + numz - 1
                except Exception as e2:
                    print(e1)
                    print(e2)

            massgrid = np.reshape(massgrid, (len(data), zlength))
            ztab = np.arange(startz, endz + 1)
            mgrid, zgrid = np.meshgrid(data[:, 0], ztab, indexing='ij')

            if not ud.isempty(self.range):
                bool1 = data[:, 0] >= self.range[0]
                bool2 = data[:, 0] <= self.range[1]
                bool3 = np.all([bool1, bool2], axis=0)
                mgrid = mgrid[bool3]
                zgrid = zgrid[bool3]
                massgrid = massgrid[bool3]
            if self.normflag:
                massgrid /= np.amax(massgrid)

            if compress > 1:
                m, z, d = IM_functions.compress_2d(mgrid, zgrid, massgrid, compress)
                dat = np.transpose([np.ravel(m), np.ravel(z), np.ravel(d)])
            else:
                dat = np.transpose([np.ravel(mgrid), np.ravel(zgrid), np.ravel(massgrid)])

            data2.append(dat)
            print(i, end=' ')
        print("Loaded 2D Data")
        self.yvals = np.array(self.yvals)
        data2 = np.array(data2)
        print(data2.shape)

        if self.datachoice == 2:
            mode = "Mass"
        elif self.datachoice == 3:
            mode = "CCS"
        else:
            mode = "mz"

        PlotAnimations.AnimationWindow(self, data2, self.config, mode="2D", yvals=self.yvals[:, 1], pksmode=mode)

    def on_2dgrid(self, e):
        self.yvals = np.array(self.yvals)
        exwindow = Extract2D.Extract2DPlot(self, self.grid, self.config, yvals=self.yvals[:, 1], params=self.gridparams,
                                           directory=self.directory)
        self.gridparams = exwindow.params

    def on_defect(self, e):
        self.yvals = np.array(self.yvals)
        MassDefects.MassDefectWindow(self, self.grid, self.config, yvals=self.yvals[:, 1],
                                     directory=self.directory, value=self.molig)
        pass

    def on_autocorr(self, e):
        masssum = np.transpose([np.unique(self.grid[0, :, 0]), np.sum(self.grid[:, :, 1], axis=0)])
        print("Sum shape:", masssum.shape)
        autocorrwindow = AutocorrWindow(self)
        autocorrwindow.initalize_dialog(self.config, masssum)
        autocorrwindow.ShowModal()

    def on_local_path(self, e):
        self.update_get(0)
        for i, l in enumerate(self.yvals):
            filename = l[0]
            localpath = os.path.relpath(filename, self.directory)
            l[0] = localpath
        self.update_set(0)
        self.localpath = 1

    def on_absolute_path(self, e):
        self.update_get(0)
        for i, l in enumerate(self.yvals):
            filename = l[0]
            abspath = os.path.abspath(os.path.join(self.directory, filename))
            l[0] = abspath
        self.update_set(0)
        self.localpath = 0

    def on_export(self, e):
        if not ud.isempty(self.grid):
            try:
                x, y = np.meshgrid(self.grid[0, :, 0], self.var1, indexing='ij')
                dat = np.transpose([np.ravel(x), np.ravel(y), np.ravel(np.transpose(self.grid[:, :, 1]))])
                path = os.path.join(self.directory, "ExtractFull2D_" + self.xlabel[-3:-1] + ".txt")
                np.savetxt(path, dat)
                print("Saved to: ", path)
            except (ValueError, TypeError):
                print("Failed to export data")
        else:
            print("Grid is empty")

    def on_msms_norm(self, e):
        dlg = wx.FileDialog(self, "Choose MS1 data file in x y list format", '', "", "*.*", )
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            path = os.path.join(dirname, filename)

            dlg.Destroy()
            print("Openening: ", path)
            msdat = np.loadtxt(path)
            mids = np.array([y[1] for y in self.yvals]).astype(np.float)
            wins = np.array([y[2] for y in self.yvals]).astype(np.float)

            vals = []
            for i, m in enumerate(mids):
                w = float(wins[i])
                if w == 0:
                    index = ud.nearest(msdat[:, 0], m)
                    val = msdat[index, 1]
                    pass
                else:
                    limits = [m - w, m + w]
                    boo1 = msdat[:, 0] < limits[1]
                    boo2 = msdat[:, 0] > limits[0]
                    intdat = msdat[np.all([boo1, boo2], axis=0)]
                    val = np.trapz(intdat[:, 1], x=intdat[:, 0])
                    pass
                vals.append(val)
            vals = np.array(vals)
            vals /= np.amax(vals)
            print(vals)
            self.on_run(0, vals=vals)

    def on_doubledec(self, e):
        self.update_get(e)
        if self.filetype == 1:
            hdf = h5py.File(self.hdf5_file)
        self.paths = []
        self.headers = []
        for k, l in enumerate(self.yvals):
            if self.filetype == 1:
                msdata = hdf.get(self.topname + "/" + str(l[0]))
            filename = l[0]
            header = os.path.splitext(filename)[0]
            if self.localpath == 1 or not os.path.isabs(filename):
                header = os.path.join(self.directory, header)
                filename = os.path.join(self.directory, filename)
            subheader = os.path.split(header)[1]

            filename = os.path.join(header + "_unidecfiles", subheader + "_mass.txt")
            if self.filetype == 1:
                data = get_dataset(msdata, "mass_data")
            else:
                data = np.loadtxt(filename)

            self.paths.append(filename)
            self.headers.append(header)
        print(self.paths)

        dlg = wx.FileDialog(self, "Open Kernel Text File", self.directory, self.kernelfile, "*.txt")
        if dlg.ShowModal() == wx.ID_OK:
            self.kernelfile = dlg.GetPath()
        dlg.Destroy()
        print(self.kernelfile)
        DoubleDec.batch_dd(self.paths, self.kernelfile)



    def on_ylabel(self, e):
        dlg = miscwindows.SingleInputDialog(self)
        dlg.initialize_interface(title="Set Variable 1 Label", message="Variable 1 axis label:", defaultvalue="")
        dlg.ShowModal()
        self.ylabel = dlg.value
        print("New  var. 1 axis label:", self.ylabel)
        try:
            self.make_grid_plots(e)
        except Exception as ex:
            print("Could not plot grid:", ex)
        try:
            self.makeplot2()
        except Exception as ex:
            print("Could not plot extract:", ex)
        pass


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


# Main App Execution
if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = wx.App(False)
    frame = DataCollector(None, "Collect Data")
    app.MainLoop()

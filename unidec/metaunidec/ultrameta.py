import json
import time
from matplotlib.pyplot import colormaps
from matplotlib import rcParams
from matplotlib import colors as mplcol
import matplotlib as mpl

from pubsub import pub
import multiprocessing
from unidec.modules import unidecstructure, PlottingWindow
from unidec.modules import miscwindows
import unidec.tools as ud
import h5py
from unidec.modules.hdf5_tools import replace_dataset, get_dataset
from unidec.metaunidec import mudeng
from unidec import MetaUniDec as mudpres
from unidec.metaunidec.gui_elements.um_list_ctrl import *
import threading
import numpy as np
import wx

def read_attr(thing, string, config):
    if string in list(config.attrs.keys()):
        val = config.attrs.get(string)
        if isinstance(val, np.ndarray):
            return val[0]
        else:
            return val
    else:
        return thing


__author__ = 'michael.marty'

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42

extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position", 5: "Estimated Area"}
extractchoicesz = {0: "Maximum", 1: "Weighted Average"}
extractlabels = {0: "Intensity", 1: "Intensity", 2: "Area", 3: "Mass", 4: "Mass", 5: "Area"}

markdict = {'\u25CB': 'o', '\u25BD': 'v', '\u25B3': '^', '\u25B7': '>', '\u25A2': 's', '\u2662': 'd',
            '\u2606': '*'}
mdkeys = ['\u25CB', '\u25BD', '\u25B3', '\u25B7', '\u25A2', '\u2662', '\u2606']


# noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment,PyArgumentList
class DataCollector(wx.Frame):
    def __init__(self, parent, title, config=None, *args, **kwargs):
        wx.Frame.__init__(self, parent, title=title)  # ,size=(200,-1))
        if "directory" in kwargs:
            self.directory = kwargs["directory"]
        else:
            self.directory = ""

        self.config = config
        self.mudwin = None

        if self.config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.initialize_system_paths()
            print("Using Empty Structure")
            self.config.publicationmode = 1
            if "viridis" in colormaps():
                self.config.cmap = u"viridis"
            else:
                self.config.cmap = u"jet"

        self.CreateStatusBar(3)
        self.SetStatusWidths([-1, 150, 400])

        pub.subscribe(self.on_motion2, 'newxy')

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
        self.menulocalpath = self.toolsmenu.Append(wx.ID_ANY, "Convert to Local Path",
                                                   "Change file name to reflect local path for portability")
        self.Bind(wx.EVT_MENU, self.on_local_path, self.menulocalpath)
        self.menuabsolutepath = self.toolsmenu.Append(wx.ID_ANY, "Convert to Absolute Path",
                                                      "Change file name to reflect absolute path")
        self.Bind(wx.EVT_MENU, self.on_absolute_path, self.menuabsolutepath)
        self.toolsmenu.AppendSeparator()
        self.menuylabel = self.toolsmenu.Append(wx.ID_ANY, "Specify Var. 1 Label", "Adds Var. 1 axis label to plot")
        self.Bind(wx.EVT_MENU, self.on_ylabel, self.menuylabel)
        self.toolsmenu.AppendSeparator()
        self.menuexpfit = self.toolsmenu.Append(wx.ID_ANY, "Exponential Decay Fit",
                                                "Fit all plots to exponential decays")
        self.Bind(wx.EVT_MENU, self.on_exp_fit, self.menuexpfit)
        self.menulinfit = self.toolsmenu.Append(wx.ID_ANY, "Linear Fit",
                                                "Fit all plots to line")
        self.Bind(wx.EVT_MENU, self.on_lin_fit, self.menulinfit)
        self.menusigfit = self.toolsmenu.Append(wx.ID_ANY, "Logistic Fit",
                                                "Fit all plots to logistic equation")
        self.Bind(wx.EVT_MENU, self.on_sig_fit, self.menusigfit)
        self.toolsmenu.AppendSeparator()
        self.menuplots1 = self.toolsmenu.Append(wx.ID_ANY, "Plot Mass Distributions", "Plots Mass Distributions of All")
        self.Bind(wx.EVT_MENU, self.on_plot_all, self.menuplots1)
        self.menuplots2 = self.toolsmenu.Append(wx.ID_ANY, "Plot Mass Defects", "Plots Mass Defects of All")
        self.Bind(wx.EVT_MENU, self.on_plot_all_MD, self.menuplots2)

        self.toolsmenu.AppendSeparator()
        self.menulabframe = self.toolsmenu.Append(wx.ID_ANY, "Multiply X-axis by Average Charge State",
                                                  "Multiply the X-axis (Collision Voltage) by the average charge state to get the lab frame of energy")
        self.Bind(wx.EVT_MENU, self.on_lab_frame, self.menulabframe)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        self.SetMenuBar(self.menuBar)

        self.panel = wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.inputsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.displaysize = wx.GetDisplaySize()
        print(self.displaysize)
        self.ypanelsizer = wx.BoxSizer(wx.VERTICAL)
        self.ypanel = ListCtrlPanel(self.panel, pres=self, list_type="Y", size=(700, self.displaysize[1] - 300))
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

        self.xpanel = ListCtrlPanel(self.panel, pres=self, markers=mdkeys, list_type="X", size=(400, 200))
        self.xpanelsizer = wx.BoxSizer(wx.VERTICAL)
        self.addxbutton = wx.Button(self.panel, label="Add X Value")
        self.Bind(wx.EVT_BUTTON, self.on_add_x, self.addxbutton)
        self.xpanelsizer.Add(self.addxbutton, 0, wx.EXPAND)
        self.xpanelsizer.Add(self.xpanel, 0, wx.EXPAND)

        self.runsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer.Add(self.xpanelsizer, 0, wx.EXPAND)
        self.ctlnorm = wx.RadioBox(self.panel, label="Mass/Intensity Normalization",
                                   choices=["None", "Max", "Sum", "Peak Max", "Peak Sum"], majorDimension=1,
                                   style=wx.RA_SPECIFY_COLS)
        self.ctlnormz = wx.RadioBox(self.panel, label="Charge Normalization",
                                   choices=["None", "Max", "Sum", "Peak Max", "Peak Sum"], majorDimension=1,
                                   style=wx.RA_SPECIFY_COLS)
        self.ctlextractwindow = wx.TextCtrl(self.panel, value="", size=(100, 20))
        self.ctlextractthresh = wx.TextCtrl(self.panel, value="", size=(100, 20))
        self.ctlextract = wx.ComboBox(self.panel, value="Height", size=(150, 30), choices=list(extractchoices.values()),
                                      style=wx.CB_READONLY)
        self.ctlextractz = wx.ComboBox(self.panel, value="Weighted Average", size=(150, 30), choices=list(extractchoicesz.values()),
                                      style=wx.CB_READONLY)
        self.ctlextractwindow.SetValue("0")
        self.ctlextractthresh.SetValue("10")
        self.runbutton = wx.Button(self.panel, label="Run Extraction", size=(150, 200))
        self.runbutton.SetBackgroundColour((0, 200, 0))
        self.Bind(wx.EVT_BUTTON, self.on_run, self.runbutton)
        self.runsizer.Add(self.runbutton, 0, wx.EXPAND)
        self.runsizer.Add(self.ctlnorm, 0, wx.EXPAND)
        self.runsizer.Add(self.ctlnormz, 0, wx.EXPAND)
        self.runsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer2.Add(wx.StaticText(self.panel, label="Extraction Window: "), wx.FIXED_MINSIZE)
        self.runsizer2.Add(self.ctlextractwindow, 0, wx.FIXED_MINSIZE)

        self.runsizer2b = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer2b.Add(wx.StaticText(self.panel, label="Intensity Threshold: "), wx.FIXED_MINSIZE)
        self.runsizer2b.Add(self.ctlextractthresh, 0, wx.FIXED_MINSIZE)
        self.runsizer2b.Add(wx.StaticText(self.panel, label="%"), wx.FIXED_MINSIZE)

        self.runsizer3 = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer3.Add(wx.StaticText(self.panel, label="How to Extract Mass/Intensity: "), wx.FIXED_MINSIZE)
        self.runsizer3.Add(self.ctlextract, 0, wx.FIXED_MINSIZE)
        self.runsizer3b = wx.BoxSizer(wx.HORIZONTAL)
        self.runsizer3b.Add(wx.StaticText(self.panel, label="How to Extract Charge: "), wx.FIXED_MINSIZE)
        self.runsizer3b.Add(self.ctlextractz, 0, wx.FIXED_MINSIZE)
        self.runsizer4 = wx.BoxSizer(wx.VERTICAL)
        self.runsizer4.Add(self.runsizer2)
        self.runsizer4.Add(self.runsizer2b)
        self.runsizer4.Add(wx.StaticText(self.panel, label=" "), wx.FIXED_MINSIZE)
        self.runsizer4.Add(self.runsizer3)
        self.runsizer4.Add(self.runsizer3b)
        self.runsizer4.Add(wx.StaticText(self.panel, label=" "), wx.FIXED_MINSIZE)
        self.ctlzeros = wx.CheckBox(self.panel, label="Remove Zero Points")
        self.runsizer4.Add(self.ctlzeros)
        self.runsizer.Add(self.runsizer4)

        self.ypanelsizer.Add(self.runsizer, 0, wx.EXPAND)

        self.inputsizer.Add(self.ypanelsizer, 0, wx.EXPAND)
        self.sizer.Add(self.inputsizer, 1, wx.EXPAND)

        figsize = (7, 5)
        axes = [0.12, 0.12, 0.65, 0.8]
        self.plot1 = PlottingWindow.Plot1d(self.panel, figsize=figsize, axes=axes)
        self.plot2 = PlottingWindow.Plot1d(self.panel, figsize=figsize, axes=axes)
        self.plotsizer = wx.BoxSizer(wx.VERTICAL)
        self.plotsizer.Add(self.plot1, 0, wx.EXPAND)
        self.plotsizer.Add(self.plot2, 0, wx.EXPAND)
        self.inputsizer.Add(self.plotsizer, 0, wx.EXPAND)
        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)
        self.xvals = []
        self.yvals = []
        self.normflag2 = True
        self.datachoice = 2
        self.window = ''
        self.extractchoice = 0
        self.savename = "collection1.json"
        self.localpath = 0
        self.xcolors = []
        self.data = []
        self.grid = []
        self.var1 = []
        self.xlabel = "Mass"
        self.ylabel = ""
        self.v1name = "Variable 1"
        self.hdf5_file = ""
        self.topname = "ms_dataset"
        self.configname = "config"
        self.removezeros = False

        self.update_set(0)
        self.Centre()
        self.Show(True)

        if __name__ == "__main__":
            # testdir = "C:\Python\unidec\unidec_src\unidec\\x64\Release"
            if True:
                try:
                    testdir = "C:\\Python\\UniDec3\\TestSpectra"
                    testfile = "um_collection1.json"
                    self.load(os.path.join(testdir, testfile))
                    # self.on_plot_all()
                    pass
                except Exception as e:
                    print(e)
                    pass

    def open_file_in_meta(self, e=None, index=0):
        y = self.yvals[index]
        path = y[0]
        if self.localpath == 1:
            path = os.path.join(self.directory, path)
        if os.path.isfile(path):
            print(path)
            if self.mudwin is None:
                self.mudwin = mudpres.UniDecApp()
            try:
                self.mudwin.open_file(path)
            except:
                try:
                    self.mudwin = mudpres.UniDecApp()
                    self.mudwin.open_file(path)
                except Exception as e:
                    print("ERROR: MetaUniDec:", e)

    def load_x_from_peaks(self, e=None, index=0):
        self.update_get(e)
        try:
            y = self.yvals[index]
            path = y[0]
            if self.localpath == 1:
                path = os.path.join(self.directory, path)
            if os.path.isfile(path):
                print(path)
                self.hdf5_file = path

                hdf = h5py.File(path, 'r')
                pdataset = hdf.require_group("/peaks")
                ultrapeaks = get_dataset(pdataset, "ultrapeakdata")
                peaks = get_dataset(pdataset, "peakdata")[:, 0]
                hdf.close()

                if not ud.isempty(ultrapeaks):
                    peaks = ultrapeaks
                indexes = np.arange(0, len(peaks))
                self.xpanel.list.clear_list()
                if not ud.isempty(indexes):
                    for i in indexes:
                        marker = mdkeys[i % len(mdkeys)]
                        self.xpanel.list.add_line(val=peaks[i], marker=marker)
        except Exception as ex:
            print("Unable to detect peaks", ex)
        self.load_params_from_hdf5()

    def load_params_from_hdf5(self, e=None, index=0):
        y = self.yvals[index]
        path = y[0]
        if self.localpath == 1:
            path = os.path.join(self.directory, path)
        self.hdf5_file = path
        self.update_config()

        if not os.path.isfile(path):
            print("Error1: File Not Found", path)
            self.update_gui()
            return

        hdf = h5py.File(path, 'r')
        self.config.hdf_file = path
        h5_config = hdf.require_group("config")
        self.config.exnorm = h5_config.attrs["exnorm"]
        try:
            self.config.exnormz = h5_config.attrs["exnormz"]
        except:

            self.config.exnormz = 0

        self.config.exchoice = h5_config.attrs["exchoice"]

        try:
            self.config.exchoicez = h5_config.attrs["exchoicez"]
        except:
            self.config.exchoicez = 1

        if self.config.exchoice == 6:
            self.config.exchoice = 0
            self.config.exchoicez = 1

        if self.config.exchoice == 7:
            self.config.exchoice = 0
            self.config.exchoicez = 0

        self.config.exwindow = h5_config.attrs["exwindow"]

        try:
            self.config.exthresh = h5_config.attrs["exthresh"]
        except:
            self.config.exthresh = 10

        try:
            msdata = hdf.require_group("ms_dataset")
            self.v1name = msdata.attrs["v1name"]
        except:
            self.v1name = "Variable 1"
            pass
        hdf.close()

        self.update_gui()

    def on_save(self, e):
        self.update_get(e)
        # print "Saved: ",self.gridparams
        outdict = {"x": self.xvals, "y": self.yvals, "dir": self.directory, "zeros": self.removezeros}
        dlg = wx.FileDialog(self, "Save Collection in JSON Format", self.directory, self.savename, "*.json", wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.savename = dlg.GetPath()
            with open(self.savename, "w+") as outfile:
                json.dump(outdict, outfile)
            print("Saved: ", self.savename)
        dlg.Destroy()

    def on_load(self, e):
        dlg = wx.FileDialog(self, "Load JSON Collection", self.directory, self.savename, "*.json", wx.FD_OPEN)
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
            if not os.path.isdir(self.directory):
                self.directory = os.path.dirname(savename)
        if "zeros" in indict:
            self.removezeros = indict["zeros"]
        self.update_set(0)
        # self.load_x_from_peaks(0)
        self.load_params_from_hdf5()
        print("Loaded: ", savename)
        self.on_run(0)

    def on_add_x(self, e):
        self.update_get()
        index = len(self.xvals)
        print(index, self.xvals)
        marker = mdkeys[index % len(mdkeys)]
        self.xpanel.list.add_line(val=index, marker=marker)

    def on_add_y(self, e):
        self.update_get(e)
        dlg = wx.FileDialog(self, "Load Files", self.directory, "", "*.hdf5", wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            filenames = dlg.GetPaths()
            for f in filenames:
                self.ypanel.list.add_line(file_name=f)
            if len(self.yvals) == 0:
                self.load_x_from_peaks()
            self.localpath = 0
        dlg.Destroy()

    def update_get(self, e=None):
        self.xvals = self.xpanel.list.get_list()
        self.yvals = self.ypanel.list.get_list()
        self.directory = self.dirinput.GetValue()
        self.removezeros = self.ctlzeros.GetValue()

    def update_set(self, e):
        self.dirinput.SetValue(self.directory)
        self.xpanel.list.populate(self.xvals)
        self.ypanel.list.populate(self.yvals)
        self.ctlzeros.SetValue(self.removezeros)

    def on_exp_fit(self, e=None):
        self.on_run(fit="exp")

    def on_lin_fit(self, e=None):
        self.on_run(fit="lin")

    def on_sig_fit(self, e=None):
        self.on_run(fit="sig")

    def update_hdf5(self, path):
        try:
            if not os.path.isfile(path):
                print("Error2: File Not Found:", path)
                return
            hdf = h5py.File(path, 'a')
            h5_config = hdf.require_group("config")
            h5_config.attrs.modify("exnorm", self.config.exnorm)
            h5_config.attrs.modify("exnormz", self.config.exnormz)
            h5_config.attrs.modify("exchoice", self.config.exchoice)
            h5_config.attrs.modify("exchoicez", self.config.exchoicez)
            h5_config.attrs.modify("exwindow", self.config.exwindow)
            h5_config.attrs.modify("exthresh", self.config.exthresh)

            pdataset = hdf.require_group("/peaks")
            ultrapeakdata = np.array([int(x[1]) for x in self.xvals])
            replace_dataset(pdataset, "ultrapeakdata", data=ultrapeakdata)
            hdf.close()

            self.run_hdf5(path)
        except Exception as e:
            self.SetStatusText("ERROR with File: " + path, number=2)
            print("ERROR Python: File", path, e)

    def run_hdf5(self, path):
        if not os.path.isfile(path):
            print("Error3: File Not Found:", path)
            return
        print("Running MetaUniDec -ultraextract", path)
        out = mudeng.metaunidec_call(self.config, "-ultraextract", path=path)
        if out != 0:
            self.SetStatusText("ERROR with File: " + path, number=2)
            print("ERROR C: File", path, out)

    """
    Gets called when "Run Extraction" button gets clicked. For each file given in the list,
    runs the Unidec.exe -extract functions on the data. Then, plots the peaks from the hdf5
    extracts group. If the labels are the same on multiple files, then error bars are displayed.
    """

    def on_run(self, e=None, fit=None, labframe=False):
        self.update_config()
        tstart = time.perf_counter()
        self.plot1.clear_plot()
        self.plot2.clear_plot()
        self.update_get(e)
        print("Running\n")
        # print self.yvals
        labels = np.array(self.yvals)[:, 3]
        temp, indexes = np.unique(labels, return_index=True)
        uniquelabels = labels[np.sort(indexes)]

        output = []

        if ud.isempty(self.xvals):
            self.xvals = [['\u25CB', 0, "0", "From File"]]

        # run extraction on all files with parameters

        paths = []
        for y in self.yvals:
            path = y[0]
            if self.localpath == 1:
                path = os.path.join(self.directory, path)
            paths.append(path)

        if fit is None:
            if False:
                threads = []
                for p in paths:
                    t = threading.Thread(target=self.update_hdf5, args=(p,))
                    threads.append(t)
                    t.start()
                for t in threads:
                    t.join()
            else:
                for p in paths:
                    self.update_hdf5(p)
        print("Total Execution Time: %.2gs" % (time.perf_counter() - tstart))
        for x in self.xvals:
            print(x)
            try:
                index = int(x[1])
                marker = markdict[x[0]]
            except:
                print("Error with peak index:", x)
                index = 0
                marker = "o"
            try:
                xlabel = str(x[2])
            except:
                xlabel = str(index)

            try:
                if mplcol.is_color_like(x[3]):
                    xcolor = x[3]
                else:
                    xcolor = None
            except:
                xcolor = None

            bargraphfits = []
            bargraphlabels = []
            bargrapherrors = []
            bargraphcolors = []

            for u in uniquelabels:
                extracts = []
                zexts = []
                xvalsall = []
                for y in self.yvals:
                    label = y[3]
                    if label == u:
                        path = y[0]
                        if self.localpath == 1:
                            path = os.path.join(self.directory, path)
                        if not os.path.isfile(path):
                            path2 = os.path.join(self.directory, path)
                            if os.path.isfile(path2):
                                path = path2
                                self.localpath = 1
                                print("Switching to local path mode")
                            else:
                                print("Error4: File Not Found:", path)
                                return
                        color = y[1]
                        linestyle = y[2]
                        print(path, u, index, marker, linestyle)
                        self.hdf5_file = path

                        hdf = h5py.File(path, "r")

                        msdata1 = hdf.require_group(self.topname)
                        self.len = msdata1.attrs["num"]

                        xvals = []
                        zdat = []
                        for f in np.arange(0, self.len):
                            msdata = hdf.get(self.topname + "/" + str(f))
                            self.attrs = dict(list(msdata.attrs.items()))
                            if self.v1name in list(self.attrs.keys()):
                                var1 = self.attrs[self.v1name]
                                self.ylabel = self.v1name
                            elif "var1" in list(self.attrs.keys()):
                                var1 = self.attrs["var1"]
                            elif "Variable 1" in list(self.attrs.keys()):
                                var1 = self.attrs["Variable 1"]
                            elif "collision_voltage" in list(self.attrs.keys()):
                                var1 = self.attrs["collision_voltage"]
                                if self.ylabel == "":
                                    self.ylabel = "Collision Voltage"
                            elif "Collision Voltage" in list(self.attrs.keys()):
                                var1 = self.attrs["Collision Voltage"]
                                if self.ylabel == "":
                                    self.ylabel = "Collision Voltage"
                            else:
                                var1 = f

                            try:
                                var1 = float(var1)
                            except:
                                var1 = f

                            zdata = get_dataset(msdata, "charge_data")
                            try:
                                zdata[:, 1] /= np.amax(zdata[:, 1])
                                zval = ud.center_of_mass(zdata)[0]

                            except:
                                print("ERROR with Zdata")
                                print(zdata)
                                zval = 0
                                pass
                            zdat.append(zval)
                            xvals.append(var1)
                        xvalsall.append(xvals)
                        pdataset = hdf.require_group("/peaks")
                        # Get the peaks back in
                        ultrapeaks = get_dataset(pdataset, "ultrapeakdata")
                        peak = np.argwhere(ultrapeaks == index)[0][0]

                        try:
                            exz = get_dataset(pdataset, "ultrazextracts")[:, peak]
                        except Exception as e:
                            print(e)
                            exz = zdat
                        zexts.append(exz)

                        ex = get_dataset(pdataset, "ultraextracts")[:, peak]
                        extracts.append(ex)
                        hdf.close()

                if not ud.isempty(xvalsall) and not ud.isempty(extracts):
                    extracts = np.array(extracts)
                    zexts = np.array(zexts)
                    # xvals = np.array(xvals)
                    xvalsall = np.array(xvalsall)
                    l = len(xvalsall)
                    unix, counts = np.unique(np.hstack(xvalsall), return_counts=True)
                    b1 = counts == l
                    unix = unix[b1]

                    for i, x in enumerate(xvalsall):
                        x = np.array(x)
                        b2 = np.in1d(x, unix)
                        test = np.all(b2)
                        if not test:
                            print("Removing some data points to take only consensus x values")
                            print("X Values:", x)
                            print("Consensus:", unix)
                            extracts[i] = extracts[i][b2]
                            zexts[i] = zexts[i][b2]
                    xvals = unix

                    avg = np.mean(extracts, axis=0)
                    std = np.std(extracts, axis=0)
                    zdat = np.mean(zexts, axis=0)
                    zstd = np.std(zexts, axis=0)
                    lab = u + " " + xlabel
                    fits = []
                    zfits = []

                    if labframe:
                        xvals = xvals * zdat
                        self.ylabel = "Collision Energy (eV)"

                    if xcolor is not None:
                        color = xcolor
                        elab = None
                    else:
                        elab = lab
                        lab = None

                    if fit is None:
                        if self.removezeros:
                            boo1 = avg > 0
                            xvals = xvals[boo1]
                            avg = avg[boo1]
                            zdat = zdat[boo1]
                            std = std[boo1]
                            zstd = zstd[boo1]

                        if not self.plot1.flag:
                            self.plot1.plotrefreshtop(xvals, avg, "Mass Extracts", self.ylabel, self.ylabel2,
                                                      label=lab, marker=marker,
                                                      nopaint=True, color=color, test_kda=False, linestyle=linestyle)
                            pass
                        else:
                            self.plot1.plotadd(xvals, avg, color, linestyle=linestyle, newlabel=lab, marker=marker)
                            pass
                        if not self.plot2.flag:
                            self.plot2.plotrefreshtop(xvals, zdat, "Charge Extracts", self.ylabel, "Charge",
                                                      label=lab, marker=marker,
                                                      nopaint=True, color=color, test_kda=False, linestyle=linestyle)
                        else:
                            self.plot2.plotadd(xvals, zdat, color, linestyle=linestyle, newlabel=lab, marker=marker)
                    else:
                        print(zdat, zexts, avg)
                        if fit == "exp":
                            fits, fitdat = ud.exp_fit(xvals-xvals[0], avg)
                            zfits, zfitdat = ud.exp_fit(xvals-xvals[0], zdat)
                        elif fit == "lin":
                            fits, fitdat = ud.lin_fit(xvals, avg)
                            zfits, zfitdat = ud.lin_fit(xvals, zdat)
                        elif fit == "sig":
                            fits, fitdat = ud.sig_fit(xvals, avg)
                            zfits, zfitdat = ud.sig_fit(xvals, zdat)
                        else:
                            print("ERROR: Unsupported fit type")
                            break

                        print("Fits:", fits)
                        print("Charge Fits:", zfits)

                        bargraphfits.append(fits)
                        bargraphlabels.append(u)
                        bargraphcolors.append(color)

                        errors = [[], [], [], []]
                        sums = [0, 0, 0, 0]
                        for x in extracts:
                            if fit == "exp":
                                tmpfits, tmpfitdat = ud.exp_fit(xvals, x)
                                errors[0].append(tmpfits[0])
                                errors[1].append(tmpfits[1])
                            elif fit == "lin":
                                tmpfits, tmpfitdat = ud.lin_fit(xvals, x)
                                errors[0].append(tmpfits[0])
                                errors[1].append(tmpfits[1])
                            elif fit == "sig":
                                tmpfits, tmpfitdat = ud.sig_fit(xvals, x)
                                errors[0].append(tmpfits[0])
                                errors[1].append(tmpfits[1])
                                errors[2].append(tmpfits[2])
                                errors[3].append(tmpfits[3])
                        if fit == "sig":
                            tmp = []
                            for x in range(0, 4):
                                tmp.append(np.std(errors[x]))
                            bargrapherrors.append(tmp)
                        elif fit == "exp" or fit == "lin":
                            tmp = []
                            for x in range(0, 2):
                                tmp.append(np.std(errors[x]))
                            bargrapherrors.append(tmp)

                        try:
                            print("Standard Deviations of Fist:", bargrapherrors)
                            print("Average Fits:", bargraphfits)
                        except:
                            pass

                        if not self.plot1.flag:
                            self.plot1.plotrefreshtop(xvals, fitdat, "Mass Extracts", self.ylabel, self.ylabel2,
                                                      nopaint=True, color=color, test_kda=False, label=lab,
                                                      marker=None, linestyle=linestyle)
                            pass
                        else:
                            self.plot1.plotadd(xvals, fitdat, color, linestyle=linestyle, newlabel=lab,
                                               marker=None)
                            pass
                        if not self.plot2.flag:
                            self.plot2.plotrefreshtop(xvals, zfitdat, "Charge Extracts", self.ylabel, "Charge",
                                                      label=lab, marker=None,
                                                      nopaint=True, color=color, test_kda=False, linestyle=linestyle)
                        else:
                            self.plot2.plotadd(xvals, zfitdat, color, linestyle=linestyle, newlabel=lab,
                                               marker=None)

                        if fit == "sig":
                            self.plot1.addtext("", fits[0], (fits[3] + fits[2]) / 0.95, ymin=fits[3], color=color)
                            self.plot2.addtext("", zfits[0], (zfits[3] + zfits[2]) / 0.95, ymin=zfits[3], color=color)

                    if self.removezeros:
                        boo1 = avg > 0
                        xvals = xvals[boo1]
                        avg = avg[boo1]
                        zdat = zdat[boo1]
                        std = std[boo1]
                        zstd = zstd[boo1]



                    self.plot1.errorbars(xvals, avg, yerr=std, color=color, linestyle=" ", marker=marker, newlabel=elab)

                    self.plot2.errorbars(xvals, zdat, yerr=zstd, color=color, linestyle=" ",
                                         marker=marker, newlabel=elab)
                    out = [[lab], xvals, avg, std, zdat, zstd, fits, zfits]
                    output.append(out)

            if fit is not None:
                self.on_bar_graphs(bargraphfits, bargraphlabels, bargrapherrors, fit=fit, colors=bargraphcolors)
        print("Plotting Done")

        #output = np.array(output)
        self.data = output
        try:
            hdf = h5py.File(os.path.join(self.directory, "Extracts.hdf5"), 'a')
            for n, l in enumerate(output):
                name = l[0][0]
                try:
                    name = str(name)
                except:
                    name = str(n)
                dataset = hdf.require_group("/" + name)
                data = np.array([l[i] for i in range(1, 6)])
                fits = l[6]
                zfits = l[7]
                replace_dataset(dataset, "ultraextracts", data)
                if not ud.isempty(fits):
                    replace_dataset(dataset, "fits", fits)
                if not ud.isempty(zfits):
                    replace_dataset(dataset, "zfits", zfits)
            hdf.close()
            # outfile = open(os.path.join(self.directory, "Extracts.pkl"), "wb")
            # pickle.dump(output, outfile)
            # outfile.close()
        except Exception as ex:
            print("Failed to Export Output:", ex)
        print("Exports Done")
        self.plot1.add_legend(anchor=(1.35, 1))
        self.plot2.add_legend(anchor=(1.35, 1))
        self.plot2.repaint()
        self.plot1.repaint()

    def on_save_fig(self, e):
        """Finished"""
        self.update_get(e)
        name1 = os.path.join(self.directory, "Figure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print(name1)
        name2 = os.path.join(self.directory, "Figure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print(name2)

    def on_save_figPDF(self, e):
        """Finished"""
        self.update_get(e)
        name1 = os.path.join(self.directory, "Figure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print(name1)
        name2 = os.path.join(self.directory, "Figure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print(name2)

    def on_local_path(self, e):
        """Finished"""
        self.update_get(0)
        paths = [l[0] for l in self.yvals]
        for i, l in enumerate(self.yvals):
            filename = l[0]
            toppath = ud.commonprefix(paths)
            localpath = os.path.relpath(filename, toppath)  # self.directory)
            l[0] = localpath
        self.update_set(0)
        self.localpath = 1

    def on_absolute_path(self, e):
        """Finished"""
        self.update_get(0)
        for i, l in enumerate(self.yvals):
            filename = l[0]
            abspath = os.path.abspath(os.path.join(self.directory, filename))
            l[0] = abspath
        self.update_set(0)
        self.localpath = 0

    def on_ylabel(self, e):
        """Finished"""
        dlg = miscwindows.SingleInputDialog(self)
        dlg.initialize_interface(title="Set Variable 1 Label", message="Variable 1 axis label:", defaultvalue="")
        dlg.ShowModal()
        self.ylabel = dlg.value
        print("New  var. 1 axis label:", self.ylabel)
        try:
            self.on_run()
        except Exception as ex:
            print("Could not plot extract:", ex)
        pass

    def on_choose_dir(self, e):
        """Finished"""
        dlg = wx.DirDialog(None, "Choose Top Directory", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.directory = dlg.GetPath()
            self.dirinput.SetValue(self.directory)
            # print self.directory
        dlg.Destroy()

    def on_motion2(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=1)
        except:
            pass

    def on_bar_graphs(self, fits=None, labels=None, errors=None, fit=None, colors=None):
        if fit == "exp":
            bargraph = BarGraphWindow(self, title="Exponential Decay Fit", directory=self.directory)
            bargraph.on_exp_plot(fits, labels, errors, colors=colors)
        elif fit == "lin":
            bargraph = BarGraphWindow(self, title="Linear Fit", directory=self.directory, colors=colors)
            bargraph.on_lin_plot(fits, labels, errors, colors)
        elif fit == "sig":
            bargraph = BarGraphWindow(self, title="Logistic Fit", directory=self.directory, colors=colors)
            bargraph.on_sig_plot(fits, labels, errors, colors)
        else:
            print("Error: bad fit given")

    def update_config(self):
        self.config.exnorm = self.ctlnorm.GetSelection()
        self.config.exnormz = self.ctlnormz.GetSelection()
        self.config.exchoice = self.ctlextract.GetSelection()
        self.config.exchoicez = self.ctlextractz.GetSelection()
        self.ylabel2 = extractlabels[self.config.exchoice]
        self.removezeros = self.ctlzeros.GetValue()
        try:
            self.config.exwindow = float(self.ctlextractwindow.GetValue())
        except ValueError:
            self.config.exwindow = 0
        try:
            self.config.exthresh = float(self.ctlextractthresh.GetValue())
        except ValueError:
            self.config.exthresh = 0

    def update_gui(self):
        self.ctlnorm.SetSelection(self.config.exnorm)
        self.ctlnormz.SetSelection(self.config.exnormz)
        if self.config.exchoice > len(extractchoices):
            print("Error: Extract Choice out of range")
            if self.config.exchoice == 5:
                self.config.exchoice = 3
                self.config.exthresh = 50
        self.ctlextract.SetSelection(self.config.exchoice)
        self.ctlextractz.SetSelection(self.config.exchoicez)
        self.ctlextractwindow.SetValue(str(self.config.exwindow))
        self.ctlextractthresh.SetValue(str(self.config.exthresh))
        self.ctlzeros.SetValue(self.removezeros)

    def on_plot_all(self, e=None, type="dist"):
        tstart = time.perf_counter()
        self.update_get(e)
        print("Running Plot All")
        # print self.yvals
        self.yvals = np.array(self.yvals)
        labels = self.yvals[:, 3]
        uniquelabels = np.unique(labels)
        linestyles = self.yvals[:, 2]
        uniquelinestyles = np.unique(linestyles)
        colors = self.yvals[:, 1]
        uniquecolors = np.unique(colors)
        output = []

        xdim = len(uniquelinestyles)
        ydim = len(uniquecolors)

        plotwindow = BarGraphWindow(self, title="Plots", directory=self.directory)
        if type == "massdefect":
            plotwindow.setup_window_generic(xdim, ydim, type="2D")
        else:
            plotwindow.setup_window_generic(xdim, ydim)

        # run extraction on all files with parameters
        for y in self.yvals:
            path = y[0]
            if self.localpath == 1:
                path = os.path.join(self.directory, path)
            self.hdf5_file = path
            self.config.hdf_file = path
            out = mudeng.metaunidec_call(self.config, "-grids")

        for u in uniquelabels:
            extracts = []
            gridextracts = []
            xvals = []
            for y in self.yvals:
                label = y[3]
                if label == u:

                    path = y[0]
                    if self.localpath == 1:
                        path = os.path.join(self.directory, path)
                    color = y[1]
                    linestyle = y[2]
                    print(path, u, linestyle)
                    self.hdf5_file = path
                    hdf = h5py.File(path, "r")
                    msdata1 = hdf.require_group(self.topname)
                    config = hdf.get("config")
                    molig = read_attr(0, "molig", config)

                    self.len = msdata1.attrs["num"]

                    # Get the x values
                    xvals = []
                    for f in np.arange(0, self.len):
                        msdata = hdf.get(self.topname + "/" + str(f))
                        self.attrs = dict(list(msdata.attrs.items()))
                        if "var1" in list(self.attrs.keys()):
                            var1 = self.attrs["var1"]
                        elif "collision_voltage" in list(self.attrs.keys()):
                            var1 = self.attrs["collision_voltage"]
                            if self.ylabel == "":
                                self.ylabel = "Collision Voltage"
                        elif "Collision Voltage" in list(self.attrs.keys()):
                            var1 = self.attrs["Collision Voltage"]
                            if self.ylabel == "":
                                self.ylabel = "Collision Voltage"
                        else:
                            var1 = f
                        xvals.append(var1)

                    # Get the peaks back in
                    msdataset = hdf.require_group("/ms_dataset")
                    massaxis = get_dataset(msdataset, "mass_axis")
                    masssum = get_dataset(msdataset, "mass_sum")
                    massgrid = get_dataset(msdataset, "mass_grid")
                    extracts.append(np.transpose([massaxis, masssum]))
                    gridextracts.append(massgrid)
                    hdf.close()

            if not ud.isempty(xvals) and not ud.isempty(extracts):

                extracts = np.array(extracts)
                print(np.shape(extracts))
                xvals = np.array(xvals)
                # avg = np.mean(extracts, axis=0)
                # std = np.std(extracts, axis=0)

                edat = extracts[0]
                ypos = np.where(uniquecolors == color)[0][0]
                xpos = np.where(uniquelinestyles == linestyle)[0][0]

                if type == "massdefect":
                    data1d, data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(edat, molig)
                    if molig != 0:
                        # noinspection PyUnresolvedReferences
                        plotwindow.plots[xpos][ypos].contourplot(data2d, self.config, xlab="Mass", ylab="Mass Defect",
                                                                 normflag=1, title=u, test_kda=True, repaint=False)
                        plotwindow.plots[xpos][ypos].setup_zoom([plotwindow.plots[xpos][ypos].subplot1], 'box')
                        plotwindow.plots[xpos][ypos].subplot1.set_title(u)
                        plotwindow.plots[xpos][ypos].repaint()
                    else:
                        print("Could not create mass defect plots, oligomer mass is 0.")
                        print("To define oligomer mass, reopen the HDF5 files and set the Mass Difference.")
                else:
                    plotwindow.plots[xpos][ypos].plotrefreshtop(
                        edat[:, 0], edat[:, 1], u,
                        self.ylabel,
                        "Mass", None, None,
                        nopaint=False, color=color, test_kda=True,
                        linestyle=linestyle)

        plotwindow.Show()

    def on_plot_all_MD(self, e=None):
        self.on_plot_all(type="massdefect")

    def on_lab_frame(self, e=None):
        self.on_run(labframe=True)


class BarGraphWindow(wx.Frame):
    def __init__(self, parent, title, directory, *args, **kwargs):
        wx.Frame.__init__(self, parent, size=(700, 700), title=title)  # ,size=(200,-1))
        self.directory = directory
        self.plotname = None
        self.plotmenu = wx.Menu()
        self.parent = parent

        self.menuPNG = self.plotmenu.Append(wx.ID_ANY, "Save Figures as PNG", "Save Figures as PNG")
        self.menuPDF = self.plotmenu.Append(wx.ID_ANY, "Save Figures as PDF", "Save Figures as PDF")
        self.Bind(wx.EVT_MENU, self.on_save_figure_png, self.menuPNG)
        self.Bind(wx.EVT_MENU, self.on_save_figure_pdf, self.menuPDF)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(self.menuBar)

    def on_exp_plot(self, fits=None, labels=None, errors=None, colors=None):
        self.setup_window(2)
        self.colors = colors
        barlabels = []
        midpoint = []
        slope = []
        num = 0
        if self.colors is None:
            # colormap = cm.get_cmap('rainbow', len(fits))
            colormap = mpl.colormaps['rainbow'].resampled(len(fits))
            cols = colormap(np.arange(len(fits)))
        else:
            cols = self.colors
            print(self.colors)
        midpointerr = []
        slopeerr = []
        plotcols = []
        for x, f in enumerate(fits):
            midpoint.append(0)
            midpoint.append(f[0])
            midpointerr.append(0)
            midpointerr.append(errors[x][0])
            slope.append(0)
            slope.append(f[1])
            slopeerr.append(0)
            slopeerr.append((errors[x][1]))
            barlabels.append("")
            barlabels.append(labels[x])
            num += 2
            plotcols.append("black")
            plotcols.append(cols[x])
        midpoint.append(0)
        slope.append(0)
        midpointerr.append(0)
        slopeerr.append(0)
        barlabels.append("")
        plotcols.append("black")
        num += 1
        xvals = list(range(0, num))
        self.p1.barplottoperrors(xarr=xvals, yarr=midpoint, yerr=midpointerr,
                                 peakval=barlabels, colortab=plotcols, title="Rate")
        self.p2.barplottoperrors(xarr=xvals, yarr=slope, yerr=slopeerr,
                                 peakval=barlabels, colortab=plotcols, title="Amplitude")
        self.p1.repaint()
        self.p2.repaint()
        self.plotname = "Exp_Fit"
        self.Show()

    def on_lin_plot(self, fits=None, labels=None, errors=None, colors=None):
        self.setup_window(2)
        barlabels = []
        midpoint = []
        slope = []
        num = 0
        if colors is None:
            # colormap = cm.get_cmap('rainbow', len(fits))
            colormap = mpl.colormaps['rainbow'].resampled(len(fits))
            cols = colormap(np.arange(len(fits)))
        else:
            cols = colors

        midpointerr = []
        slopeerr = []
        plotcols = []
        for x, f in enumerate(fits):
            midpoint.append(0)
            midpoint.append(f[0])
            midpointerr.append(0)
            midpointerr.append(errors[x][0])
            slope.append(0)
            slope.append(f[1])
            slopeerr.append(0)
            slopeerr.append((errors[x][1]))
            barlabels.append("")
            barlabels.append(labels[x])
            num += 2
            plotcols.append("black")
            plotcols.append(cols[x])
        midpoint.append(0)
        slope.append(0)
        midpointerr.append(0)
        slopeerr.append(0)
        barlabels.append("")
        plotcols.append("black")
        num += 1
        xvals = list(range(0, num))
        self.p1.barplottoperrors(xarr=xvals, yarr=midpoint, yerr=midpointerr,
                                 peakval=barlabels, colortab=plotcols, title="Slope")
        self.p2.barplottoperrors(xarr=xvals, yarr=slope, yerr=slopeerr,
                                 peakval=barlabels, colortab=plotcols, title="Intercept")
        self.p1.repaint()
        self.p2.repaint()
        self.plotname = "Lin_Fit"
        self.Show()

    def on_sig_plot(self, fits=None, labels=None, errors=None, colors=None):
        self.setup_window(4)
        barlabels = []
        baseline = []
        slope = []
        amplitude = []
        midpoint = []
        num = 0
        if colors is None:
            # colormap = cm.get_cmap('rainbow', len(fits))
            colormap = mpl.colormaps['rainbow'].resampled(len(fits))
            cols = colormap(np.arange(len(fits)))
        else:
            cols = colors
        midpointerr = []
        slopeerr = []
        amplitudeerr = []
        baselineerr = []
        plotcols = []
        for x, f in enumerate(fits):
            midpoint.append(0)
            midpoint.append(f[0])
            midpointerr.append(0)
            midpointerr.append(errors[x][0])
            slope.append(0)
            slope.append(f[1])
            slopeerr.append(0)
            slopeerr.append((errors[x][1]))
            amplitude.append(0)
            amplitude.append(f[2])
            amplitudeerr.append(0)
            amplitudeerr.append(errors[x][2])
            baseline.append(0)
            baseline.append(f[3])
            baselineerr.append(0)
            baselineerr.append(errors[x][3])
            barlabels.append("")
            barlabels.append(labels[x])
            num += 2
            plotcols.append("black")
            plotcols.append(cols[x])
        midpoint.append(0)
        slope.append(0)
        amplitude.append(0)
        baseline.append(0)
        midpointerr.append(0)
        slopeerr.append(0)
        amplitudeerr.append(0)
        baselineerr.append(0)
        barlabels.append("")
        plotcols.append("black")
        num += 1
        xvals = list(range(0, num))
        self.p1.barplottoperrors(xarr=xvals, yarr=midpoint, yerr=midpointerr,
                                 peakval=barlabels, colortab=plotcols, title="Midpoint")
        self.p2.barplottoperrors(xarr=xvals, yarr=slope, yerr=slopeerr,
                                 peakval=barlabels, colortab=plotcols, title="Slope Parameter")
        self.p3.barplottoperrors(xarr=xvals, yarr=amplitude, yerr=amplitudeerr,
                                 peakval=barlabels, colortab=plotcols, title="Amplitude")
        self.p4.barplottoperrors(xarr=xvals, yarr=baseline, yerr=baselineerr,
                                 peakval=barlabels, colortab=plotcols, title="Baseline")
        self.p1.repaint()
        self.p2.repaint()
        self.p3.repaint()
        self.p4.repaint()
        self.plotname = "Sig_Fit"
        self.Show()

    def setup_window(self, numplots=2):
        # self.panel = wx.Panel(self)
        # noinspection PyUnresolvedReferences
        self.panel = wx.lib.scrolledpanel.ScrolledPanel(self)
        self.panel.SetupScrolling()
        plotsizer = wx.GridBagSizer()
        figsize = (6, 5)
        self.plots = []
        self.p1 = PlottingWindow.Plot1d(self.panel, figsize=figsize)
        self.p2 = PlottingWindow.Plot1d(self.panel, figsize=figsize)
        self.plots.append(self.p1)
        self.plots.append(self.p2)
        plotsizer.Add(self.p1, (0, 0), span=(1, 1), flag=wx.EXPAND)
        plotsizer.Add(self.p2, (0, 1), span=(1, 1), flag=wx.EXPAND)
        if numplots == 4:
            self.p3 = PlottingWindow.Plot1d(self.panel, figsize=figsize)
            self.p4 = PlottingWindow.Plot1d(self.panel, figsize=figsize)
            self.plots.append(self.p3)
            self.plots.append(self.p4)
            plotsizer.Add(self.p3, (1, 0), span=(1, 1), flag=wx.EXPAND)
            plotsizer.Add(self.p4, (1, 1), span=(1, 1), flag=wx.EXPAND)
        self.panel.SetSizer(plotsizer)
        plotsizer.Fit(self)

    def setup_window_generic(self, xdim=1, ydim=1, type="1D"):
        # self.panel = wx.Panel(self)
        self.panel = wx.Panel(self)  # wx.lib.scrolledpanel.ScrolledPanel(self)
        # self.panel.SetupScrolling()
        plotsizer = wx.GridBagSizer()
        figsize = (18 / ydim, 10 / xdim)
        self.plots = []
        for x in range(0, xdim):
            ptemp = []
            for y in range(0, ydim):
                if type == "2D":
                    p = PlottingWindow.Plot2d(self.panel, figsize=figsize)
                else:
                    p = PlottingWindow.Plot1d(self.panel, figsize=figsize)
                plotsizer.Add(p, (x, y), span=(1, 1), flag=wx.EXPAND)
                ptemp.append(p)
            self.plots.append(ptemp)
        self.panel.SetSizer(plotsizer)
        self.plotname = "Plot"
        plotsizer.Fit(self)

    def save_all_figures(self, extension, e=0, header=None, **kwargs):
        """
        Save All of the Figures. Will name as header+extension2+_FigureX.+exetension
        :param extension: Figure type (pdf, eps, png). Anything accepted by matplotlib
        :param extension2: Additional text to include in the figure header.
        :param e: Dummy wx Event
        :param header: Option to add different header. Default of none yields self.outfname as the path header
        :param kwargs: Any keywards to pass to the matplotlib savefig command such as Transparent or DPI
        :return: figureflags, files (the figures that were successfully saved and the files that they were saved to)
        """
        figureflags = []
        files = []
        if header is None:
            header = self.plotname
        else:
            header += self.plotname

        for i, plot in enumerate(np.ravel(self.plots)):
            name1 = header + "_" + str(i) + "." + extension
            path = os.path.join(self.directory, name1)
            if plot.flag:
                plot.on_save_fig(e, path, **kwargs)
                figureflags.append(i + 1)
                files.append([i + 1, name1])
        return figureflags, files

    def on_save_figure_png(self, e, **kwargs):
        """
        Save all figures as PNG
        :param e: Dummy wx event
        :param kwargs: keywards to pass to matplotlib savefig
        :return: None
        """
        self.save_all_figures("png", **kwargs)
        pass

    def on_save_figure_pdf(self, e):
        self.save_all_figures("pdf")


# Main App Execution
if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = wx.App(False)
    frame = DataCollector(None, "Ultra Meta Data Collector")
    app.MainLoop()

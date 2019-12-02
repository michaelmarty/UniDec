import os

import matplotlib.pyplot as plt
import numpy as np
import wx
from pubsub import pub

import GUniDec
import unidec_modules.unidectools as ud
from unidec_modules import plot1d
from unidec_modules.mzMLimporter import mzMLimporter

try:
    from unidec_modules.data_reader import DataImporter
except:
    print("Could not open Data Reader")

try:
    from unidec_modules.waters_importer import Importer
except:
    print("Could not open Waters Importer")


class ChromEngine:
    """
    Imports mzML data files.
    """

    def __init__(self):
        self.filename = None
        self.dirname = None
        self.path = None
        self.cdata = None
        self.tic = None
        self.ticdat = None

    def load_mzml(self, path, *args, **kwargs):
        if os.path.splitext(path)[1] == ".mzML":
            self.cdata = mzMLimporter(path)
            try:
                self.tic = self.cdata.get_tic()
                self.ticdat = np.transpose([self.tic.time, self.tic.i])
            except:
                print("Error getting TIC in mzML; trying to make it...")
                t = self.cdata.times
                tic = [np.sum(d[:,1]) for d in self.cdata.data]
                self.ticdat = np.transpose([t, tic])
        elif os.path.splitext(path)[1] == ".raw" and os.path.isdir(path):
            self.cdata = Importer.WatersDataImporter(path)
            self.tic = self.cdata.get_tic()
            self.ticdat = np.array(self.tic)
        else:
            self.cdata = DataImporter(path)
            self.tic = self.cdata.get_tic()
            self.ticdat = np.array(self.tic)
        print(self.ticdat.shape)
        if False:
            plt.plot(self.ticdat[:, 0], self.ticdat[:, 1])
            plt.show()

    def get_scans(self, scan_range=None):
        self.mzdata = self.cdata.get_data(scan_range)
        return self.mzdata


class ChromWindow(wx.Frame):
    def __init__(self, parent, title, path=None, *args, **kwargs):
        wx.Frame.__init__(self, parent, title=title)  # ,size=(200,-1))
        self.CreateStatusBar(4)
        self.SetStatusWidths([-1, 150, 150, 150])
        pub.subscribe(self.on_motion, 'newxy')
        pub.subscribe(self.on_selection, 'scans_selected')

        self.filemenu = wx.Menu()
        self.menuOpen = self.filemenu.Append(wx.ID_SAVE, "Open File", "Open File")
        self.Bind(wx.EVT_MENU, self.on_open, self.menuOpen)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.SetMenuBar(self.menuBar)

        self.panel = wx.Panel(self)
        self.panel.SetDropTarget(ChromDropTarget(self))
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.inputsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plotc = plot1d.Plot1d(self.panel)
        self.plotm = plot1d.Plot1d(self.panel)

        self.inputsizer.Add(self.plotc)
        self.inputsizer.Add(self.plotm)

        self.sizer.Add(self.inputsizer, 1, wx.EXPAND)

        self.ctlsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.runbutton = wx.Button(self.panel, label="Run UniDec")
        self.Bind(wx.EVT_BUTTON, self.on_unidec_run, self.runbutton)
        self.ctlsizer.Add(self.runbutton)

        self.exportbutton = wx.Button(self.panel, label="Export m/z data as:")
        self.Bind(wx.EVT_BUTTON, self.export_mz_file, self.exportbutton)
        self.ctlsizer.Add(self.exportbutton)

        self.ctlfname = wx.TextCtrl(self.panel, value="", size=(100, 20))
        self.ctlsizer.Add(self.ctlfname, 1, wx.EXPAND)

        self.autobutton = wx.Button(self.panel, label="Autoexport in steps of:")
        self.Bind(wx.EVT_BUTTON, self.on_autoexport, self.autobutton)
        self.ctlsizer.Add(self.autobutton)

        self.ctlmin = wx.TextCtrl(self.panel, value="2", size=(100, 20))
        self.ctlsizer.Add(self.ctlmin, 0)
        self.ctlsizer.Add(wx.StaticText(self.panel, label="minutes. Start label at:"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlstart = wx.TextCtrl(self.panel, value="40", size=(100, 20))
        self.ctlsizer.Add(self.ctlstart, 0)
        self.ctlsizer.Add(wx.StaticText(self.panel, label=" and increment label by:"), 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlincr = wx.TextCtrl(self.panel, value="10", size=(100, 20))
        self.ctlsizer.Add(self.ctlincr, 0)

        self.sizer.Add(self.ctlsizer, 0, wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)

        self.Centre()
        self.Show(True)

        self.eng = ChromEngine()
        self.outputfile = None
        self.udwin = None
        self.scans = None
        self.autoparams = None

        if path is not None:
            self.open_file(path)

    def get_from_gui(self):
        self.outputfile = self.ctlfname.GetValue()

        minutes = float(self.ctlmin.GetValue())
        increment = float(self.ctlincr.GetValue())
        start = float(self.ctlstart.GetValue())

        self.autoparams = [minutes, start, increment]
        pass

    def load_to_gui(self):
        self.ctlfname.SetValue(self.eng.filename)
        pass

    def on_open(self, e=None):
        """
        Open dialog for file opening
        :param e: unused space for event
        :return: None
        """
        dlg = wx.FileDialog(self, "Choose an MS data file", '', "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            self.open_file(dlg.GetPath())
        dlg.Destroy()

    def open_file(self, path):
        self.eng.dirname, self.eng.filename = os.path.split(path)
        self.SetStatusText("Opening", number=0)
        print("Openening: ", self.eng.filename)
        self.eng.path = path
        self.eng.load_mzml(self.eng.path)
        self.SetStatusText(self.eng.filename, number=0)
        self.load_to_gui()
        self.plotc.plotrefreshtop(self.eng.ticdat[:, 0], self.eng.ticdat[:, 1], title="Total Ion Chromatogram",
                                  xlabel="Time (min)", ylabel="Intensity", zoom="fixed_span")

    def plot_mz(self):
        self.plotm.plotrefreshtop(self.eng.mzdata[:, 0], self.eng.mzdata[:, 1])

    def on_selection(self, min, max):
        minscan = ud.nearest(self.eng.ticdat[:, 0], min)
        if self.eng.ticdat[minscan, 0] < min:
            minscan += 1
        maxscan = ud.nearest(self.eng.ticdat[:, 0], max)
        if self.eng.ticdat[maxscan, 0] > max:
            maxscan -= 1
        if maxscan <= minscan:
            maxscan = minscan + 1
        self.scans = [minscan, maxscan, min, max]
        self.get_scan_data()

    def get_scan_data(self):
        self.SetStatusText("Time: %.2f to %.2f" % (self.scans[2], self.scans[3]), number=1)
        self.SetStatusText("Scans: " + str(self.scans[0]) + " to " + str(self.scans[1]), number=2)
        self.eng.get_scans(scan_range=self.scans[:2])
        try:
            self.plot_mz()
        except (TypeError, IndexError):
            print("Error Plotting Scans")

    def on_unidec_run(self, e=None):
        savefile = self.export_mz_file()
        if self.udwin is None:
            self.udwin = GUniDec.UniDecApp()
        self.udwin.on_open_file(savefile, self.eng.dirname)
        pass

    def export_mz_file(self, e=None, append=None, directory=None):
        extension = ".txt"
        if append is not None:
            extension = "_" + append + extension
        if directory is None:
            directory = self.eng.dirname

        self.get_from_gui()
        if self.outputfile != "":
            savefile = self.outputfile
        else:
            savefile = self.eng.filename
        savefile = os.path.splitext(savefile)[0] + extension
        savepath = os.path.join(directory, savefile)
        print(self.outputfile, savepath)
        np.savetxt(savepath, self.eng.mzdata)
        return savefile

    def autoexport(self, time_window, start_count, increment):
        times = np.arange(0, np.amax(self.eng.ticdat[:, 0]), time_window)
        labels = np.arange(start_count, start_count + len(times) * increment, increment)
        directory = os.path.join(self.eng.dirname, os.path.splitext(self.eng.filename)[0] + "_extracts")
        if not os.path.isdir(directory):
            os.mkdir(directory)
        for i, t in enumerate(times):
            self.on_selection(t, t + time_window)
            self.export_mz_file(0, append="{0:g}".format(labels[i]), directory=directory)

    def on_autoexport(self, e=None):
        self.get_from_gui()
        self.autoexport(self.autoparams[0], self.autoparams[1], self.autoparams[2])

    def on_motion(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=3)
        except:
            pass


class ChromDropTarget(wx.FileDropTarget):
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
            path = filenames[0]
            self.window.open_file(path)
        elif len(filenames) > 1:
            for f in filenames:
                self.window.open_file(f)
                self.window.on_autoexport(0)
        else:
            print("Error in drag and drop.")
        return 0


if __name__ == "__main__":
    # dir = "C:\\cprog\\UniDecDemo"
    # file = "20150304_myo_70r_5u_full-ALL.mzML"
    # dir = "C:\\Users\\michael.marty\\Hugh\Data\\290415\\mzML"
    # file = "20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_2_35000.mzML"
    # file = "20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_1_35000.mzML"
    # file="20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_24_35000.mzML"
    dir = "C:\Data\EmptyNanodiscs\MixedBelts"
    file = "160302_MTM_ND_MixPO_Ramp40_190_2_Test.mzML"

    path = os.path.join(dir, file)
    path = None
    # data = ChromWindow(os.path.join(dir, file))
    app = wx.App(False)
    frame = ChromWindow(None, "Chromatogram Viewer", path)
    app.MainLoop()
    # print data

    # TODO: Control click to add multiple windows

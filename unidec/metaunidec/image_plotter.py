import time
import os
import numpy as np
from pubsub import pub
import wx
from scipy.spatial.distance import cosine
from unidec.modules import PlottingWindow, unidecstructure
import unidec.tools as ud
from unidec.metaunidec import mudstruct

__author__ = 'Michael.Marty'

'''
Window for viewing imaging results
'''


class ImagingWindow(wx.Frame):
    """
    Dialog window for Imaging.
    """

    def __init__(self, *args, **kwargs):
        """
        Creates a dialog.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Frame.__init__(self, *args, **kwargs)
        self.SetSize((1600, 1000))
        self.SetTitle("Imaging Plots")
        self.config = None
        self.plot1 = None
        self.plot2 = None
        self.plot3 = None
        self.plot4 = None
        self.plot5 = None
        self.plot6 = None
        self.exchoice = 2
        self.mass_extracted = None
        self.mz_extracted = None

    def init(self, data, config=None):
        """
        """
        self.config = config
        if self.config is None:
            self.config = unidecstructure.UniDecConfig()
        self.init_data(data)

        # Make the menu
        filemenu = wx.Menu()
        menu_open = filemenu.Append(wx.ID_ANY, "Open HDF5 File",
                                    "Open an HDF5 to view")
        self.Bind(wx.EVT_MENU, self.on_open_hdf5, menu_open)
        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")

        # Add an analysis menu
        analysismenu = wx.Menu()
        menu_analyze = analysismenu.Append(wx.ID_ANY, "Cosine Similarity Mass vs. m/z",
                                           "Print the Cosine Similarity of the Mass vs. m/z images")
        self.Bind(wx.EVT_MENU, self.on_cosine, menu_analyze)
        menu_bar.Append(analysismenu, "&Analysis")

        self.SetMenuBar(menu_bar)

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.plot1 = PlottingWindow.Plot1d(panel, smash=2)
        self.Bind(self.plot1.EVT_MZLIMITS, self.extract, self.plot1)

        self.plot2 = PlottingWindow.Plot2d(panel, integrate=1)
        pub.subscribe(self.sum_region, 'integrate')

        self.plot3 = PlottingWindow.Plot1d(panel, smash=2)
        self.Bind(self.plot3.EVT_MZLIMITS, self.extract2, self.plot3)

        self.plot4 = PlottingWindow.Plot1d(panel, smash=2)
        self.Bind(self.plot4.EVT_MZLIMITS, self.extract3, self.plot4)

        self.plot5 = PlottingWindow.Plot2d(panel, integrate=1)
        pub.subscribe(self.sum_region, 'integrate')

        self.plot6 = PlottingWindow.Plot1d(panel, smash=2)
        self.Bind(self.plot6.EVT_MZLIMITS, self.extract4, self.plot6)

        hbox.Add(self.plot1)
        hbox.Add(self.plot2)
        hbox.Add(self.plot3)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(self.plot4)
        hbox2.Add(self.plot5)
        hbox2.Add(self.plot6)

        vbox.Add(hbox, 1, wx.EXPAND)
        vbox.Add(hbox2, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer.Add(wx.StaticText(panel, label="Extraction Type: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctletype = wx.ComboBox(panel, value="Area", choices=["Midpoint Height", "Local Max", "Area"],
                                    style=wx.CB_READONLY | wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctletype, 0, wx.ALIGN_CENTER_VERTICAL)

        zbutton = wx.Button(panel, label="Plot Charge States")
        controlsizer.Add(zbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.plot_mz_values, zbutton)

        # controlsizer.Add(wx.StaticText(panel, label="   Window (Da):"), 0, wx.ALIGN_CENTER_VERTICAL)
        # self.ctlwindow = wx.TextCtrl(panel)
        # controlsizer.Add(self.ctlwindow)

        vbox.Add(controlsizer, 0, wx.EXPAND)

        panel.SetSizer(vbox)
        vbox.Fit(self)
        self.CenterOnParent()
        try:
            self.load()
            self.extract(no_face=True)
        except:
            pass
        self.Show()

    def init_data(self, data):
        self.data = data
        self.x = np.array(self.data.var1)
        self.y = np.array(self.data.var2)
        self.mass_extracted = None
        self.mz_extracted = None

    def update(self, e=None):
        self.exchoice = self.ctletype.GetSelection()

    def load(self, e=None):
        self.load_plot1()
        self.load_plot4()

    def load_plot4(self, e=None):
        self.plot4.plotrefreshtop(self.data.mzdat[:, 0], self.data.mzdat[:, 1], xlabel="m/z (Th)",
                                  ylabel="Intensity")

    def load_plot1(self, e=None):
        self.plot1.plotrefreshtop(self.data.massdat[:, 0], self.data.massdat[:, 1], xlabel="Mass (Da)",
                                  ylabel="Intensity")

    def extract(self, e=None, erange=None, dtype="mass", plot=None, no_face=False):
        self.update()
        if erange is None:
            erange = self.plot1.subplot1.get_xlim()
            self.load_plot1()
            if not no_face:
                self.plot1.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.data.massdat[:, 1]),
                                    facecolor="y")
        if dtype == "mass":
            grid = self.data.massgrid[:, :, 1]
            dat = self.data.massdat[:, 0]
            self.massrange = erange
        elif dtype == "mz":
            grid = self.data.mzgrid[:, :, 1]
            dat = self.data.mzdat[:, 0]
            self.mzrange = erange
        else:
            return

        print("Extracting in Range:", erange)
        midpoint = np.mean(erange)
        window = midpoint - erange[0]
        self.exdata = ud.extract_from_data_matrix(dat, grid, midpoint=midpoint, window=window,
                                                  extract_method=self.exchoice)
        self.image_plot(plot=plot)

        if dtype == "mass":
            self.mass_extracted = self.exdata
        elif dtype == "mz":
            self.mz_extracted = self.exdata

    def extract2(self, e=None):
        erange = self.plot3.subplot1.get_xlim()
        self.load_plot1()
        self.plot1.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.data.massdat[:, 1]), facecolor="y")
        self.plot3.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.region_massdat[:, 1]), facecolor="y")
        self.extract(erange=erange)

    def sum_region(self, e=None):
        xrange = self.plot2.subplot1.get_xlim()
        yrange = self.plot2.subplot1.get_ylim()
        print("Extracting in Shape Range:", xrange, yrange)

        b1 = self.x >= xrange[0]
        b2 = self.x <= xrange[1]
        b3 = self.y >= yrange[0]
        b4 = self.y <= yrange[1]

        ball = b1 * b2 * b3 * b4

        regiondat = np.sum(self.data.massgrid[ball, :, 1], axis=0)
        regiondat = np.transpose([self.data.massdat[:, 0], regiondat])
        self.plot3.plotrefreshtop(regiondat[:, 0], regiondat[:, 1], xlabel="Mass (Da)", ylabel="Intensity")
        self.region_massdat = regiondat

        regiondat2 = np.sum(self.data.mzgrid[ball, :, 1], axis=0)
        regiondat2 = np.transpose([self.data.mzdat[:, 0], regiondat2])
        self.plot6.plotrefreshtop(regiondat2[:, 0], regiondat2[:, 1], xlabel="m/z (Th)", ylabel="Intensity")
        self.region_mzdat = regiondat2

    def image_plot(self, e=None, data=None, plot=None):
        if data is None:
            data = self.exdata
        if plot is None:
            plot = self.plot2

        dat = np.transpose([self.x, self.y, data])
        plot.contourplot(dat=dat, normflag=1, config=self.config,
                         xlab="x", ylab="y", discrete=1, order=None)

    def extract3(self, e=None):
        erange = self.plot4.subplot1.get_xlim()
        self.load_plot4()
        self.plot4.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.data.mzdat[:, 1]), facecolor="y")
        self.extract(erange=erange, dtype="mz", plot=self.plot5)

    def extract4(self, e=None):
        erange = self.plot6.subplot1.get_xlim()
        self.load_plot4()
        self.plot4.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.data.mzdat[:, 1]), facecolor="y")
        self.plot6.add_rect(erange[0], 0, erange[1] - erange[0], np.amax(self.region_mzdat[:, 1]), facecolor="y")
        self.extract(erange=erange, dtype="mz", plot=self.plot5)

    def plot_mz_values(self, e=None):
        print("Plotting mz")
        erange = self.massrange
        midpoint = np.mean(erange)
        localmaxpos = ud.data_extract(self.data.massdat, midpoint, window=(erange[1] - erange[0]) / 2.,
                                      extract_method=4)
        print(localmaxpos, midpoint, erange)
        self.ztab = np.arange(self.config.startz, self.config.endz + 1)
        mztab = (localmaxpos + self.ztab * self.config.adductmass) / self.ztab
        b1 = mztab < np.amax(self.data.mzdat[:, 0])
        b2 = mztab > np.amin(self.data.mzdat[:, 0])
        b3 = b1 * b2
        ztab = self.ztab[b3]
        mztab = mztab[b3]
        self.plot1.addtext(str(localmaxpos), localmaxpos, 0.95 * np.amax(self.data.massdat[:, 1]), vlines=True)
        for i, z in enumerate(ztab):
            self.plot4.addtext(str(z), mztab[i], 0.95 * np.amax(self.data.mzdat[:, 1]), vlines=True)

    def on_open_hdf5(self, e=None):
        dlg = wx.FileDialog(self, "Choose a data file in HDF5 format", '', "", "*.hdf5*")
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print("Openening: ", path)
            if os.path.splitext(path)[1] != ".hdf5":
                print("Need HDF5 file")
                return
        dlg.Destroy()
        data = mudstruct.MetaDataSet(None)
        data.import_hdf5(path, speedy=True)
        data.import_grids()
        self.init_data(data)
        self.load()
        self.extract()

    def on_cosine(self, event=None):
        print("Calculating Cosine Similarity")
        d1 = self.mass_extracted
        d2 = self.mz_extracted

        if d1 is None or d2 is None:
            print("Need to select mass and m/z data first")
            return

        if np.shape(d1) != np.shape(d2):
            print("Mass and m/z data must be the same size")
            print(np.shape(d1), np.shape(d2))
            return

        cos = 1 - cosine(d1, d2)

        outtext = "Cosine Similarity: " + str(np.round(cos * 100, 2)) + "%"

        print(outtext)
        self.plot5.add_title(outtext)


if __name__ == "__main__":
    path = "C:\\Data\\Imaging\\UoB_RepairedRK_image_200x200um_proc\\RepairedRK_image_200x200um_proc - Copy.hdf5"
    os.chdir(os.path.dirname(path))
    s = time.perf_counter()
    data = mudstruct.MetaDataSet(None)
    data.import_hdf5(path, speedy=True)
    data.import_grids()
    print("Import Time:", time.perf_counter() - s)

    app = wx.App(False)
    frame = ImagingWindow(None)
    frame.init(data)
    app.MainLoop()

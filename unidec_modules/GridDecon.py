import os
from copy import deepcopy
import numpy as np
import wx
import time

import scipy.ndimage.filters as filt
import scipy.spatial.distance as distance
from scipy.sparse import csr_matrix
from scipy.signal import savgol_filter

import unidecstructure
import plot2d
import plot1d
import unidectools as ud

__author__ = 'Michael.Marty'

def make_pmat(mzgrid, fwhm):
    vec = np.ravel(mzgrid)
    g1, g2 = np.meshgrid(vec, vec)
    sig = fwhm / 2.35482
    pmat = np.exp(-(g1 - g2) ** 2. / (2. * sig * sig))
    thresh = 0.001
    pmat = np.clip(pmat, thresh, 1) - thresh
    pmat /= np.amax(pmat)
    return pmat


def conv(b, pmat):
    shape = b.shape
    dotproduct = np.dot(np.ravel(b), pmat)
    return np.reshape(dotproduct, shape)


def blur(b, msig, zsig):
    return filt.gaussian_filter(b, [msig, zsig])


def grid_unidec(mzgrid, igrid, numit=1000, fwhm=5, mode=0, msig=0, zsig=1):
    pmat = make_pmat(mzgrid, fwhm)
    b = deepcopy(igrid)
    b2 = b
    diff = 1
    i = 0
    while diff > 1e-9 and i < numit:
        # Blur
        b = blur(b,msig,zsig)
        # Deconvolve
        if mode == 0:
            b = ud.safedivide(b * igrid, conv(b, pmat))
        elif mode == 1:
            b = ud.safedivide(b * igrid * igrid, conv(b * igrid, pmat))
        elif mode == 2:
            b *= conv(ud.safedivide(igrid, conv(b, pmat)), pmat)
        else:
            print "Mode is invalid:", mode
            exit()
        # Convergence Test
        b /= np.amax(b)
        diff = np.sum((b2 - b) ** 2.)
        i += 1
        b2 = b
    print "Interations:", i
    return b
    pass


class GridDeconWindow(wx.Frame):
    def __init__(self, parent, data, config=None, directory=None):
        """

        masses = m0 + m1 * range(min m1, max m1 +1)

        :return: None
        """

        wx.Frame.__init__(self, parent, title="2D Grid Extraction", size=(-1, -1))
        # Make Menu
        self.filemenu = wx.Menu()
        self.menuSaveFigPNG = self.filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                                   "Save all figures as PNG in central directory")
        self.menuSaveFigPDF = self.filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                                   "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_save_fig, self.menuSaveFigPNG)
        self.Bind(wx.EVT_MENU, self.on_save_figPDF, self.menuSaveFigPDF)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.SetMenuBar(self.menuBar)

        # Initialize Parameters
        self.data = data

        self.config = config
        if self.config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.cmap="viridis"
            if directory is None:
                self.directory = os.getcwd()
        else:
            if directory is None:
                self.directory = self.config.dirname

        if self.config.griddecon is None:
            self.params = [98868, 760.076, 0, 120, 1, 25, 0, 20, 3, 1]
        else:
            self.params = self.config.griddecon

        if config is not None:
            self.params[0] = self.config.oligomerlist[0][0]
            self.params[1] = self.config.oligomerlist[0][1]
            self.params[2] = self.config.oligomerlist[0][2]
            self.params[3] = self.config.oligomerlist[0][3]
            self.params[4] = self.config.startz
            self.params[5] = self.config.endz
            self.params[7] = self.config.mzsig
            self.params[8] = self.config.msig
            self.params[9] = self.config.zzsig
            self.storediscrete = deepcopy(self.config.discreteplot)
        else:
            self.storediscrete = 1

        # Setup GUI
        panel = wx.Panel(self)
        topsizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot2d.Plot2d(panel)
        sizer.Add(self.plot1, 1, wx.EXPAND)
        sizer.Add(self.plot2, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.VERTICAL)
        controlsizer1 = wx.BoxSizer(wx.VERTICAL)

        self.ctlm0 = wx.TextCtrl(panel, value=str(self.params[0]))
        self.ctlm1 = wx.TextCtrl(panel, value=str(self.params[1]))
        self.ctlm1min = wx.TextCtrl(panel, value=str(self.params[2]))
        self.ctlm1max = wx.TextCtrl(panel, value=str(self.params[3]))
        self.ctlm2min = wx.TextCtrl(panel, value=str(self.params[4]))
        self.ctlm2max = wx.TextCtrl(panel, value=str(self.params[5]))
        self.ctlwindow = wx.TextCtrl(panel, value=str(self.params[6]))
        self.ctlfwhm = wx.TextCtrl(panel, value=str(self.params[7]))
        self.ctlmsig = wx.TextCtrl(panel, value=str(self.params[8]))
        self.ctlzzsig = wx.TextCtrl(panel, value=str(self.params[9]))

        controlsizer1.Add(wx.StaticText(panel, label="Base Mass"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm0, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Monomer Mass"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm1, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass Window"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlwindow, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Monomer Min #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm1min, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Monomer Max #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm1max, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Minimum Charge State"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm2min, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Maximum Charge State"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm2max, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Peak Full Width at Half Max"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlfwhm, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass Filter Width"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlmsig, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Charge Filter Width"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlzzsig, 0, wx.EXPAND)

        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        extractbutton = wx.Button(panel, label="Extract")
        deconbutton = wx.Button(panel, label="Deconvolute")

        controlsizer2.Add(extractbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.extract, extractbutton)
        controlsizer2.Add(deconbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.deconvolve, deconbutton)

        controlsizer.Add(controlsizer1, 0, wx.EXPAND)
        controlsizer.Add(controlsizer2, 0, wx.EXPAND)

        topsizer.Add(sizer, 0, wx.EXPAND)
        topsizer.Add(controlsizer, 1, wx.EXPAND)

        self.Bind(wx.EVT_CLOSE, self.on_close, self)

        panel.SetSizer(topsizer)
        topsizer.Fit(self)

        # Run initial extraction
        self.extract(0)
        self.deconvolve(0)

        self.Centre()
        self.Show(True)

    def getfromgui(self):
        """
        Extract parameters from GUI to self.params
        :return:
        """
        try:
            self.params[0] = float(self.ctlm0.GetValue())
            self.params[1] = float(self.ctlm1.GetValue())
            self.params[2] = float(self.ctlm1min.GetValue())
            self.params[3] = float(self.ctlm1max.GetValue())
            self.params[4] = float(self.ctlm2min.GetValue())
            self.params[5] = float(self.ctlm2max.GetValue())
            self.params[7] = float(self.ctlfwhm.GetValue())
            self.params[8] = float(self.ctlmsig.GetValue())
            self.params[9] = float(self.ctlzzsig.GetValue())
            try:
                self.params[6] = float(self.ctlwindow.GetValue())
            except ValueError:
                self.params[6] = 0
        except ValueError:
            print "Failed to get from gui"

    def extract(self, e):
        self.getfromgui()
        self.makegrids()
        self.makeplot()

    def makegrids(self):
        """
        Make grid of mass values for potential combinations of m1 and m2 + m0.
        :return: None
        """
        self.mrange = np.arange(self.params[2], self.params[3] + 1)
        self.zrange = np.arange(self.params[4], self.params[5] + 1)
        self.mgrid, self.zgrid = np.meshgrid(self.mrange, self.zrange, indexing='ij')
        self.mzgrid = (self.params[0] + self.mgrid * self.params[1] + self.zgrid) / self.zgrid
        self.igrid = ud.data_extract_grid(self.data, self.mzgrid, window=self.params[6])

    def deconvolve(self, e):
        tstart = time.clock()
        self.igrid = grid_unidec(self.mzgrid, self.igrid, fwhm=self.params[7], msig=self.params[8], zsig=self.params[9])
        tend = time.clock()
        print "Deconvolution time: %.2gs" % (tend - tstart)
        self.makeplot(title="Deconvolved Data")
        pass

    def makeplot(self, title="Extracted Data"):
        """
        Make the 2D and 1D plots for element self.pos in data_list.
        :return: None
        """
        dat = np.transpose([np.ravel(self.mgrid), np.ravel(self.zgrid), np.ravel(self.igrid)])
        self.config.discreteplot = 1
        # self.config.cmap="jet"
        try:
            self.plot2.contourplot(dat, self.config, xlab="Oligomer Number", ylab="Charge", title=title, normflag=True,
                                   normrange=[0, 1])
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.plotrefreshtop(np.unique(self.mgrid), np.sum(self.igrid, axis=1), title, "Oligomer Number",
                                      "Total Intensity", "", self.config, test_kda=False, nopaint=False)
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

    def on_close(self, e=None):
        """
        Close the window. Return self.config.discreteplot to its original value.
        :param e:
        :return:
        """
        print "Closing"
        self.config.discreteplot = self.storediscrete
        self.Destroy()

    def on_save_fig(self, e):
        """
        Save figures as a PNG in self.directory.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Extract2DFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Extract2DFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2

    # noinspection PyPep8Naming
    def on_save_figPDF(self, e):
        """
        Save figures as a PDF in self.directory.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Extract2DFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Extract2DFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2


if __name__ == "__main__":
    dir_name = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\"
    file_name = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit100.txt"

    path = os.path.join(dir_name, file_name)

    data = np.loadtxt(path)

    app = wx.App(False)
    frame = GridDeconWindow(None, data)
    app.MainLoop()

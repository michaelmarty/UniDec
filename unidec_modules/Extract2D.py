import os
from copy import deepcopy
import numpy as np
import wx

import unidecstructure
import plot2d
import plot1d
import unidectools as ud

__author__ = 'Michael.Marty'


# def __init__(self,datalist,config=None,*args,**kwargs):
#    wx.Window.__init__(self, *args, **kwargs)
class Extract2DPlot(wx.Frame):
    def __init__(self, parent, datalist, config=None, yvals=None, directory=None, header=None, params=None, *args,
                 **kwargs):
        wx.Frame.__init__(self, parent, title="2D Grid Extraction", size=(-1, -1))

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

        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config

        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        if header is None:
            self.header = "Extract"
        else:
            self.header = header

        if params is None:
            self.params = [98868, 760.076, 22044, 0, 90, 0, 2, 0]
        else:
            self.params = params

        self.datalist = datalist
        self.dlen = len(datalist)
        self.pos = -1
        self.yvals = np.array(yvals).astype(np.float)
        if ud.isempty(yvals):
            self.yvals = np.arange(0, len(datalist))
        self.storediscrete = deepcopy(self.config.discreteplot)

        self.panel = wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = plot1d.Plot1d(self.panel)
        self.plot2 = plot2d.Plot2d(self.panel)
        self.sizer.Add(self.plot1, 1, wx.EXPAND)
        self.sizer.Add(self.plot2, 1, wx.EXPAND)

        self.controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.controlsizer1 = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlm0 = wx.TextCtrl(self.panel, value=str(self.params[0]))
        self.ctlm1 = wx.TextCtrl(self.panel, value=str(self.params[1]))
        self.ctlm2 = wx.TextCtrl(self.panel, value=str(self.params[2]))
        self.ctlm1min = wx.TextCtrl(self.panel, value=str(self.params[3]))
        self.ctlm1max = wx.TextCtrl(self.panel, value=str(self.params[4]))
        self.ctlm2min = wx.TextCtrl(self.panel, value=str(self.params[5]))
        self.ctlm2max = wx.TextCtrl(self.panel, value=str(self.params[6]))
        self.ctlwindow = wx.TextCtrl(self.panel, value=str(self.params[7]))

        self.controlsizer.Add(wx.StaticText(self.panel, label="Mass 0"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlm0, 0, wx.EXPAND)
        self.controlsizer.Add(wx.StaticText(self.panel, label="Mass 1"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlm1, 0, wx.EXPAND)
        self.controlsizer.Add(wx.StaticText(self.panel, label="Mass 2"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlm2, 0, wx.EXPAND)
        self.controlsizer.Add(wx.StaticText(self.panel, label="Mass Window"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlwindow, 0, wx.EXPAND)
        if self.dlen > 1:
            self.ctlnorm = wx.CheckBox(self.panel, label="Normalize")
            self.controlsizer.Add(self.ctlnorm, 0, wx.EXPAND)
        self.controlsizer1.Add(wx.StaticText(self.panel, label="Mass 1 Min #"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer1.Add(self.ctlm1min, 0, wx.EXPAND)
        self.controlsizer1.Add(wx.StaticText(self.panel, label="Mass 1 Max #"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer1.Add(self.ctlm1max, 0, wx.EXPAND)
        self.controlsizer1.Add(wx.StaticText(self.panel, label="Mass 2 Min #"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer1.Add(self.ctlm2min, 0, wx.EXPAND)
        self.controlsizer1.Add(wx.StaticText(self.panel, label="Mass 2 Max #"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer1.Add(self.ctlm2max, 0, wx.EXPAND)

        self.controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.backbutton = wx.Button(self.panel, label="Back")
        self.nextbutton = wx.Button(self.panel, label="Next")
        if self.dlen > 1:
            self.totalbutton = wx.Button(self.panel, label="Total")
        else:
            self.totalbutton = wx.Button(self.panel, label="Replot")
        self.wapbutton = wx.Button(self.panel, label="WAP")
        if self.dlen > 1:
            self.controlsizer2.Add(self.backbutton, 0, wx.EXPAND)
            self.controlsizer2.Add(self.nextbutton, 0, wx.EXPAND)
        self.controlsizer2.Add(self.totalbutton, 0, wx.EXPAND)
        if self.dlen > 1:
            self.controlsizer2.Add(self.wapbutton, 0, wx.EXPAND)

        self.Bind(wx.EVT_BUTTON, self.on_back, self.backbutton)
        self.Bind(wx.EVT_BUTTON, self.on_next, self.nextbutton)
        self.Bind(wx.EVT_BUTTON, self.on_total, self.totalbutton)
        self.Bind(wx.EVT_BUTTON, self.on_wap, self.wapbutton)

        self.sizer.Add(self.controlsizer, 0, wx.EXPAND)
        self.sizer.Add(self.controlsizer1, 0, wx.EXPAND)
        self.sizer.Add(self.controlsizer2, 0, wx.EXPAND)

        self.Bind(wx.EVT_CLOSE, self.on_close, self)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)
        try:
            self.on_total(0)
        except Exception, e:
            self.on_next(0)
            print e
        self.Centre()
        self.Show(True)

    def modparams(self):
        self.params[0] = self.m0
        self.params[1] = self.m1
        self.params[2] = self.m2
        self.params[3] = self.m1minmax[0]
        self.params[4] = self.m1minmax[1]
        self.params[5] = self.m2minmax[0]
        self.params[6] = self.m2minmax[1]
        self.params[7] = self.window

    def getfromgui(self):
        try:
            self.m0 = float(self.ctlm0.GetValue())
            self.m1 = float(self.ctlm1.GetValue())
            self.m2 = float(self.ctlm2.GetValue())
            self.m1minmax = [int(self.ctlm1min.GetValue()), int(self.ctlm1max.GetValue())]
            self.m2minmax = [int(self.ctlm2min.GetValue()), int(self.ctlm2max.GetValue())]
            try:
                self.window = float(self.ctlwindow.GetValue())
            except ValueError:
                self.window = 0
            try:
                self.normflag = int(self.ctlnorm.GetValue())
            except ValueError:
                self.normflag = 1
            self.modparams()
        except ValueError:
            print "Failed to get from gui"

    def makegrid(self):
        self.grids = []
        self.m1range = np.arange(self.m1minmax[0], self.m1minmax[1] + 1)
        self.m2range = np.arange(self.m2minmax[0], self.m2minmax[1] + 1)
        self.m1grid, self.m2grid = np.meshgrid(self.m1range, self.m2range, indexing='ij')
        self.massgrid = self.m0 + self.m1grid * self.m1 + self.m2grid * self.m2

    def extractall(self):
        self.igrid = np.zeros((len(self.datalist), len(self.m1range), len(self.m2range)))
        for i, data in enumerate(self.datalist):
            for j in xrange(0, len(self.m1range)):
                for k in xrange(0, len(self.m2range)):
                    mass = self.massgrid[j, k]

                    if not self.window > 0:
                        pos = ud.nearest(data[:, 0], mass)
                        if pos != 0 and pos != len(data) - 1:
                            intensity = data[pos, 1]
                        else:
                            intensity = 0
                    else:
                        mtest = np.abs(data[:, 0] - mass)
                        btest = mtest < self.window
                        pint = data[btest]
                        if not ud.isempty(pint):
                            intensity = np.amax(pint[:, 1])
                        else:
                            intensity = 0
                    self.igrid[i, j, k] = intensity
        self.igrid /= np.amax(self.igrid)

    def makeplot(self):
        i = self.pos
        dat = np.transpose([np.ravel(self.m1grid), np.ravel(self.m2grid), np.ravel(self.igrid[i])])
        self.config.discreteplot = 1
        # self.config.cmap="jet"

        title = str(self.yvals[i])
        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title=title, normflag=self.normflag,
                                   normrange=[0, 1])
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.sum(self.igrid[i], axis=1), title, "mass 1",
                                      "Total Intensity", "", self.config, test_kda=False, nopaint=False)
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

    def makeplottotal(self):
        grid = np.sum(self.igrid, axis=0)
        dat = np.transpose([np.ravel(self.m1grid), np.ravel(self.m2grid), np.ravel(grid)])
        self.config.discreteplot = 1
        # self.config.cmap="jet"
        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title="Total Extraction", normflag=0,
                                   normrange=[np.amin(grid), np.amax(grid)])
            outfile = os.path.join(self.directory, self.header + "_grid_2D_Extract.txt")
            np.savetxt(outfile, dat)
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.sum(grid, axis=1), "Total Projection", "mass 1",
                                      "Total Intensity", "", self.config, test_kda=False, nopaint=False)
            outputdata = np.transpose([np.unique(self.m1grid), np.sum(grid, axis=1)])
            outfile = os.path.join(self.directory, self.header + "_total_2D_Extract.txt")
            np.savetxt(outfile, outputdata)
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

    def makeplotwap(self):
        grid = np.zeros((len(self.m1range), len(self.m2range)))
        for j in xrange(0, len(self.m1range)):
            for k in xrange(0, len(self.m2range)):
                grid[j, k] = np.average(self.yvals, weights=self.igrid[:, j, k])
        dat = np.transpose([np.ravel(self.m1grid), np.ravel(self.m2grid), np.ravel(grid)])
        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title="Weighted Average of Position",
                                   normflag=0, normrange=[np.amin(grid), np.amax(grid)])
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.average(grid, axis=1), "Average Projection", "mass 1",
                                      "Average Position", "", self.config, test_kda=False, nopaint=False)
            outputdata = np.transpose([np.unique(self.m1grid), np.average(grid, axis=1)])
            outfile = os.path.join(self.directory, self.header + "_WAP_2D_Extract.txt")
            np.savetxt(outfile, outputdata)
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

    def on_close(self, e=None):
        print "Closing"
        self.config.discreteplot = self.storediscrete
        self.Destroy()

    def on_back(self, e=None):
        self.getfromgui()
        self.pos += -1
        self.pos %= len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_next(self, e=None):
        self.getfromgui()
        self.pos += 1
        self.pos %= len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_total(self, e=None):
        self.getfromgui()
        self.makegrid()
        self.extractall()
        self.makeplottotal()

    def on_wap(self, e):
        self.getfromgui()
        self.makegrid()
        self.extractall()
        self.makeplotwap()

    def on_save_fig(self, e):

        name1 = os.path.join(self.directory, "Extract2DFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print name1
        name2 = os.path.join(self.directory, "Extract2DFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print name2

    # noinspection PyPep8Naming
    def on_save_figPDF(self, e):

        name1 = os.path.join(self.directory, "Extract2DFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print name1
        name2 = os.path.join(self.directory, "Extract2DFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print name2


# Main App Execution
if __name__ == "__main__":
    dir_name = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit100_unidecfiles"
    file_name = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit100_mass.txt"

    path = os.path.join(dir_name, file_name)

    data = np.loadtxt(path)

    dir_name = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_unidecfiles"
    file_name = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_mass.txt"
    path = os.path.join(dir_name, file_name)
    data2 = np.loadtxt(path)

    dir_name = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit180_unidecfiles"
    file_name = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit180_mass.txt"
    path = os.path.join(dir_name, file_name)
    data3 = np.loadtxt(path)

    datalist = [data, data2, data3]

    app = wx.App(False)
    frame = Extract2DPlot(None, datalist, yvals=np.array([100., 120., 180.]))
    app.MainLoop()

import os
from copy import deepcopy
import numpy as np
import wx

import unidecstructure
import plot2d
import plot1d
import unidectools as ud

__author__ = 'Michael.Marty'

"""
Window for extracting intensity values from data to make 2D plots and 1D plots.
"""


class Extract2DPlot(wx.Frame):
    def __init__(self, parent, data_list, config=None, yvals=None, directory=None, header=None, params=None):
        """
        Create wx.Frame and initialzie components
        :param parent: Parent window or panel passed to wx.Frame
        :param data_list: Input data for extraction in a list of arrays (N x 2)
        :param config: UniDecConfig object. If None, will use defaults.
        :param yvals: Position values for corresponding data_list elements.
        For plots, these become the titles. For the Weighted-Average-of-Position (WAP) these are the position values.
        :param directory: Directory for saving files. Default is current working directory.
        :param header: Header for files that are written. Default is "Extract"
        :param params: List of 8 values that define the parameters for extraction.

        0=mass 0
        1=mass 1
        2=mass 2
        3=minimum oligomeric state of mass 1
        4=maximum oligomeric state of mass 1
        5=minimum oligomeric state of mass 2
        6=maximum oligomeric state of mass 2
        7=Error window for finding intensity value

        masses = m0 + m1 * range(min m1, max m1 +1) + m2 * range(min m2, max m2 +1)

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
        if config is None:
            # Default UniDecConfig object
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

        self.datalist = data_list
        self.dlen = len(data_list)
        self.pos = -1
        self.yvals = np.array(yvals).astype(np.float)
        if ud.isempty(yvals):
            self.yvals = np.arange(0, len(data_list))
        self.storediscrete = deepcopy(self.config.discreteplot)

        # Setup GUI
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot2d.Plot2d(panel)
        sizer.Add(self.plot1, 1, wx.EXPAND)
        sizer.Add(self.plot2, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer1 = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlm0 = wx.TextCtrl(panel, value=str(self.params[0]))
        self.ctlm1 = wx.TextCtrl(panel, value=str(self.params[1]))
        self.ctlm2 = wx.TextCtrl(panel, value=str(self.params[2]))
        self.ctlm1min = wx.TextCtrl(panel, value=str(self.params[3]))
        self.ctlm1max = wx.TextCtrl(panel, value=str(self.params[4]))
        self.ctlm2min = wx.TextCtrl(panel, value=str(self.params[5]))
        self.ctlm2max = wx.TextCtrl(panel, value=str(self.params[6]))
        self.ctlwindow = wx.TextCtrl(panel, value=str(self.params[7]))

        controlsizer.Add(wx.StaticText(panel, label="Mass 0"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm0, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass 1"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm1, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass 2"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm2, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass Window"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwindow, 0, wx.EXPAND)
        if self.dlen > 1:
            self.ctlnorm = wx.CheckBox(panel, label="Normalize")
            controlsizer.Add(self.ctlnorm, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass 1 Min #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm1min, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass 1 Max #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm1max, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass 2 Min #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm2min, 0, wx.EXPAND)
        controlsizer1.Add(wx.StaticText(panel, label="Mass 2 Max #"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlm2max, 0, wx.EXPAND)

        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        backbutton = wx.Button(panel, label="Back")
        nextbutton = wx.Button(panel, label="Next")
        if self.dlen > 1:
            totalbutton = wx.Button(panel, label="Total")
        else:
            totalbutton = wx.Button(panel, label="Replot")
        wapbutton = wx.Button(panel, label="WAP")
        if self.dlen > 1:
            controlsizer2.Add(backbutton, 0, wx.EXPAND)
            controlsizer2.Add(nextbutton, 0, wx.EXPAND)
        controlsizer2.Add(totalbutton, 0, wx.EXPAND)
        if self.dlen > 1:
            controlsizer2.Add(wapbutton, 0, wx.EXPAND)

        self.Bind(wx.EVT_BUTTON, self.on_back, backbutton)
        self.Bind(wx.EVT_BUTTON, self.on_next, nextbutton)
        self.Bind(wx.EVT_BUTTON, self.on_total, totalbutton)
        self.Bind(wx.EVT_BUTTON, self.on_wap, wapbutton)

        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer1, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)

        self.Bind(wx.EVT_CLOSE, self.on_close, self)

        panel.SetSizer(sizer)
        sizer.Fit(self)
        # Run initial extraction
        try:
            self.on_total(0)
        except Exception, e:
            self.on_next(0)
            print e
        self.Centre()
        self.Show(True)
        self.normflag = 1

    def getfromgui(self):
        """
        Extract parameters from GUI to self.params
        :return:
        """
        try:
            self.params[0] = float(self.ctlm0.GetValue())
            self.params[1] = float(self.ctlm1.GetValue())
            self.params[2] = float(self.ctlm2.GetValue())
            self.params[3] = int(self.ctlm1min.GetValue())
            self.params[4] = int(self.ctlm1max.GetValue())
            self.params[5] = int(self.ctlm2min.GetValue())
            self.params[6] = int(self.ctlm2max.GetValue())
            try:
                self.params[7] = float(self.ctlwindow.GetValue())
            except ValueError:
                self.params[7] = 0
            try:
                self.normflag = int(self.ctlnorm.GetValue())
            except (ValueError, AttributeError):
                self.normflag = 1
        except ValueError:
            print "Failed to get from gui"

    def makegrid(self):
        """
        Make grid of mass values for potential combinations of m1 and m2 + m0.
        :return: None
        """
        self.m1range = np.arange(self.params[3], self.params[4] + 1)
        self.m2range = np.arange(self.params[5], self.params[6] + 1)
        self.m1grid, self.m2grid = np.meshgrid(self.m1range, self.m2range, indexing='ij')
        self.massgrid = self.params[0] + self.m1grid * self.params[1] + self.m2grid * self.params[2]

    def extractall(self):
        """
        Extract intensity values from self.datalist for mass values in self.massgrid.
        :return: None
        """
        # TODO: Remove this functionality to unidectools to merge with datacollector.
        # TODO: Optimize function. It currently preforms the extraction on all every time but doesn't need to...
        self.igrid = np.zeros((len(self.datalist), len(self.m1range), len(self.m2range)))
        for i, data in enumerate(self.datalist):
            for j in xrange(0, len(self.m1range)):
                for k in xrange(0, len(self.m2range)):
                    mass = self.massgrid[j, k]
                    if not self.params[7] > 0:
                        pos = ud.nearest(data[:, 0], mass)
                        if pos != 0 and pos != len(data) - 1:
                            intensity = data[pos, 1]
                        else:
                            intensity = 0
                    else:
                        mtest = np.abs(data[:, 0] - mass)
                        btest = mtest < self.params[7]
                        pint = data[btest]
                        if not ud.isempty(pint):
                            intensity = np.amax(pint[:, 1])
                        else:
                            intensity = 0
                    self.igrid[i, j, k] = intensity
        self.igrid /= np.amax(self.igrid)

    def makeplot(self):
        """
        Make the 2D and 1D plots for element self.pos in data_list.
        :return: None
        """
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
        """
        Make the 2D and 1D plots for the sum of all arrays in data_list. Write outputs to _grid_2D_extract.txt and
        _total_2D_extract.txt.
        :return:
        """
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
        """
        Calculates the weighted average of position (WAP) for each element in the intensity grid.
        Writes 1D output to _WAP_2D_Extract.txt.
        :return: None
        """
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
        """
        Close the window. Return self.config.discreteplot to its original value.
        :param e:
        :return:
        """
        print "Closing"
        self.config.discreteplot = self.storediscrete
        self.Destroy()

    def on_back(self, e=None):
        """
        PLot the extraction from the previous array in data_list.
        :param e: Unused event
        :return:
        """
        self.getfromgui()
        self.pos += -1
        self.pos %= len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_next(self, e=None):
        """
        Plot the extraction from the next array in data_list.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()
        self.pos += 1
        self.pos %= len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_total(self, e=None):
        """
        Extract all and plot total.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()
        self.makegrid()
        self.extractall()
        self.makeplottotal()

    def on_wap(self, e):
        """
        Extract all and plot weighted average of position.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()
        self.makegrid()
        self.extractall()
        self.makeplotwap()

    def on_save_fig(self, e):
        """
        Save figures as a PNG in self.directory.
        :param e: Unused event
        :return: None
        """
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
        """
        Save figures as a PDF in self.directory.
        :param e: Unused event
        :return: None
        """
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

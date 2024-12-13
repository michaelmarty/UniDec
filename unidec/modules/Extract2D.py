import os
from copy import deepcopy
import numpy as np
import wx
# import matplotlib.cm as cm
import matplotlib as mpl

from unidec.modules import unidecstructure
from modules.plotting import PlottingWindow
from unidec.modules import miscwindows
from unidec import tools as ud
from unidec.modules.MassFitter import MassFitter

__author__ = 'Michael.Marty'

"""
Window for extracting intensity values from data to make 2D plots and 1D plots.
"""


def get_kd(mutant, wt, n=0, l=1, temp=25, units="kJ"):
    if "kcal" in units:
        R = 1.98720425864083 / 1000.
    else:
        R = 8.31446261815324 / 1000.
    T = 273.15 + temp
    a = wt[n]
    al = wt[n + l]
    b = mutant[n]
    bl = mutant[n + l]
    kd = a * bl / (al * b)
    deltag = -R * T * np.log(kd)
    return kd, deltag, deltag / l


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
        8=Extract Method

        masses = m0 + m1 * range(min m1, max m1 +1) + m2 * range(min m2, max m2 +1)

        :return: None
        """
        wx.Frame.__init__(self, parent, title="2D Grid Extraction", size=(-1, -1))
        # Make Menu
        self.filemenu = wx.Menu()
        self.menuRunAll = self.filemenu.Append(wx.ID_ANY, "Run Next On All",
                                                   "Run Next on All")
        self.menuSaveFigPNG = self.filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                                   "Save all figures as PNG in central directory")
        self.menuSaveFigPDF = self.filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                                   "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_all_next, self.menuRunAll)
        self.Bind(wx.EVT_MENU, self.on_save_fig, self.menuSaveFigPNG)
        self.Bind(wx.EVT_MENU, self.on_save_figPDF, self.menuSaveFigPDF)

        self.plotmenu = wx.Menu()
        self.menufit = self.plotmenu.Append(wx.ID_ANY, "Fit Gaussians",
                                            "Fit total distribution to a series of Gaussians")
        self.menufit2 = self.plotmenu.Append(wx.ID_ANY, "Fit Poisson",
                                             "Fit total distribution to a Poisson Distribution")
        self.menufit3 = self.plotmenu.Append(wx.ID_ANY, "Fit Binomial",
                                             "Fit total distribution to a Binomial Distribution")
        self.menufit4 = self.plotmenu.Append(wx.ID_ANY, "Fit Multiple Poissons",
                                             "Fit total distribution to multiple Poisson distributions")
        self.Bind(wx.EVT_MENU, self.on_fit, self.menufit)
        self.Bind(wx.EVT_MENU, self.on_fit2, self.menufit2)
        self.Bind(wx.EVT_MENU, self.on_fit3, self.menufit3)
        self.Bind(wx.EVT_MENU, self.on_fit4, self.menufit4)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(self.menuBar)
        # Initialize Parameters
        if config is None:
            # Default UniDecConfig object
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config
            if directory is None:
                directory = self.config.udir

        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        if header is None:
            self.header = "Extract"
        else:
            self.header = header

        if params is None:
            self.params = [98868, 760.076, 22044, 0, 90, 0, 2, 0, 1]
            self.params = [0, 4493, 678, 1, 20, 0, 30, 10, 1]
            self.params = [98601, 1400, 352, 0, 7, 0, 1, 100, 1]
        else:
            self.params = params

        self.extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position",
                               5: "Estimated Area", 6: "Area 10% Threshold"}

        self.extractdict = {0:0, 1:1, 2:2, 3:3, 4:4, 5:11, 6:12}

        self.datalist = data_list
        self.dlen = len(data_list)
        self.pos = -1
        self.yvals = np.array(yvals).astype(float)
        if ud.isempty(yvals):
            self.yvals = np.arange(0, len(data_list))
        self.storediscrete = deepcopy(self.config.discreteplot)

        # Setup GUI
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = PlottingWindow.Plot1d(panel)
        self.plot2 = PlottingWindow.Plot2d(panel)
        self.plot3 = PlottingWindow.Plot1d(panel)
        self.plot4 = PlottingWindow.Plot1d(panel)
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2.Add(self.plot1, 1, wx.EXPAND)
        sizer2.Add(self.plot3, 1, wx.EXPAND)
        sizer.Add(sizer2, 1, wx.EXPAND)
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2.Add(self.plot2, 1, wx.EXPAND)
        sizer2.Add(self.plot4, 1, wx.EXPAND)
        sizer.Add(sizer2, 1, wx.EXPAND)

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
        self.ctlextract = wx.ComboBox(panel, value=self.extractchoices[self.params[8]], choices=list(self.extractchoices.values()), style=wx.CB_READONLY)

        controlsizer.Add(wx.StaticText(panel, label="Mass 0"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm0, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass 1"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm1, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass 2"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm2, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label="Mass Window"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwindow, 0, wx.EXPAND)
        controlsizer.Add(wx.StaticText(panel, label=" How to extract: "), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlextract, 0, wx.EXPAND)

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
        # wapbutton = wx.Button(panel, label="WAP")
        if self.dlen > 1:
            controlsizer2.Add(backbutton, 0, wx.EXPAND)
            controlsizer2.Add(nextbutton, 0, wx.EXPAND)
        controlsizer2.Add(totalbutton, 0, wx.EXPAND)
        # if self.dlen > 1:
        # controlsizer2.Add(wapbutton, 0, wx.EXPAND)

        self.Bind(wx.EVT_BUTTON, self.on_back, backbutton)
        self.Bind(wx.EVT_BUTTON, self.on_next, nextbutton)
        self.Bind(wx.EVT_BUTTON, self.on_total, totalbutton)
        # self.Bind(wx.EVT_BUTTON, self.on_wap, wapbutton)

        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer1, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)

        self.Bind(wx.EVT_CLOSE, self.on_close, self)

        panel.SetSizer(sizer)
        sizer.Fit(self)
        # Run initial extraction
        try:
            self.on_total(0)
        except Exception as e:
            try:
                self.on_next(0)
            except:
                pass
            print(e)
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
            self.params[8] = int(self.ctlextract.GetSelection())
            try:
                self.params[7] = float(self.ctlwindow.GetValue())
            except ValueError:
                self.params[7] = 0
            try:
                self.normflag = int(self.ctlnorm.GetValue())
            except (ValueError, AttributeError):
                self.normflag = 1
        except ValueError:
            print("Failed to get from gui")

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
        # TODO: Optimize function. It currently preforms the extraction on all every time but doesn't need to...
        self.igrid = np.zeros((len(self.datalist), len(self.m1range), len(self.m2range)))
        method = self.extractdict[self.params[8]]
        for i, data in enumerate(self.datalist):
            self.igrid[i] = ud.data_extract_grid(data, self.massgrid, extract_method=method, window=self.params[7])
        try:
            self.igrid /= np.amax(self.igrid)
        except Exception as e:
            print(e)
        self.makedataplot()

    def makedataplot(self):
        i = self.pos
        data = self.datalist[i]
        data[:, 1] = data[:, 1] / np.amax(data[:, 1])
        self.plot3.plotrefreshtop(data[:, 0], data[:, 1], "Data", "Mass", "Intensity", "", self.config, test_kda=True)
        window = self.params[7]
        # max = np.amax(data[:, 1])
        for n, m in enumerate(self.massgrid.ravel()):
            low = (m - window) / self.plot3.kdnorm
            maximum = np.ravel(self.igrid[i])[n]
            w = window * 2 / self.plot3.kdnorm
            self.plot3.add_rect(low, 0, w, maximum, alpha=0.35, nopaint=True)
        self.plot3.repaint()
        self.makebarchart()

    def makebarchart(self):
        i = self.pos
        mutant = self.igrid[i][:,0]
        wt = self.igrid[i][:,1]
        lval = len(mutant)-1
        units="kJ/mol"
        x = []
        y = []
        col = []
        lab = []
        for n in range(lval):
            value = get_kd(mutant, wt, n, 1, units=units)
            x.append(2*n)
            y.append(mutant[n] / np.amax(mutant))
            lab.append(str(n))
            col.append("b")
            x.append(2*n+1)
            y.append(wt[n] / np.amax(wt))
            lab.append(" ")
            col.append("g")
            print(n, " =", value, units)
        self.plot4.barplottop(x, y, lab, col, "Number of Bound Lipids", "Intensity", repaint=True)
        pass

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
        if type(self.directory) is not str and type(self.directory) is not os.PathLike:
            d = self.directory[i]
        else:
            d = self.directory

        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title=title, normflag=self.normflag,
                                   normrange=[0, 1])
            outfile = os.path.join(d, self.header + "_grid_2D_Extract"+ str(self.params[8])+".txt")
            np.savetxt(outfile, dat)
            print("Saved to:", outfile)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2a", e)

        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.sum(self.igrid[i], axis=1), title, "mass 1",
                                      "Total Intensity", "", self.config, test_kda=False, nopaint=False)
            self.data1d = np.transpose([np.unique(self.m1grid), np.sum(self.igrid[i], axis=1)])
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1a", e)

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
        if type(self.directory) is not str and type(self.directory) is not os.PathLike:
            newdir = self.directory[self.pos]
        else:
            newdir = self.directory


        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title="Total Extraction", normflag=0,
                                   normrange=[np.amin(grid), np.amax(grid)])
            outfile = os.path.join(newdir, self.header + "_grid_2D_Extract.txt")
            np.savetxt(outfile, dat)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)

        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.sum(grid, axis=1) / np.amax(np.sum(grid, axis=1)),
                                      "Total Projection", "mass 1",
                                      "Total Intensity", "", self.config, test_kda=False, nopaint=False)
            outputdata = np.transpose([np.unique(self.m1grid), np.sum(grid, axis=1)])
            self.data1d = outputdata

            for d in range(len(self.m2range)):
                if np.amax(grid[:, d]) != 0:
                    data = deepcopy(grid[:, d] / np.amax(grid[:, d]))
                else:
                    data = deepcopy(grid[:, d])
                # color = cm.get_cmap("viridis")(float(d / len(self.m2range)))
                color = mpl.colormaps["viridis"].resampled(float(d / len(self.m2range)))
                self.plot1.plotadd(np.unique(self.m1grid), data, colval=color)
            self.plot1.repaint()

            outfile = os.path.join(newdir, self.header + "_total_2D_Extract.txt")
            print(outfile)
            np.savetxt(outfile, outputdata)
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)

    def makeplotwap(self):
        """
        Calculates the weighted average of position (WAP) for each element in the intensity grid.
        Writes 1D output to _WAP_2D_Extract.txt.
        :return: None
        """
        grid = np.zeros((len(self.m1range), len(self.m2range)))
        for j in range(0, len(self.m1range)):
            for k in range(0, len(self.m2range)):
                grid[j, k] = np.average(self.yvals, weights=self.igrid[:, j, k])
        dat = np.transpose([np.ravel(self.m1grid), np.ravel(self.m2grid), np.ravel(grid)])
        try:
            self.plot2.contourplot(dat, self.config, xlab="mass 1", ylab="mass 2", title="Weighted Average of Position",
                                   normflag=0, normrange=[np.amin(grid), np.amax(grid)])
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)
        try:
            self.plot1.plotrefreshtop(np.unique(self.m1grid), np.average(grid, axis=1), "Average Projection", "mass 1",
                                      "Average Position", "", self.config, test_kda=False, nopaint=False)
            outputdata = np.transpose([np.unique(self.m1grid), np.average(grid, axis=1)])
            self.data1d = outputdata

            if type(self.directory) is not str and type(self.directory) is not os.PathLike:
                d = self.directory[self.pos]
            else:
                d = self.directory

            outfile = os.path.join(d, self.header + "_WAP_2D_Extract.txt")
            np.savetxt(outfile, outputdata)
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)

    def on_close(self, e=None):
        """
        Close the window. Return self.config.discreteplot to its original value.
        :param e:
        :return:
        """
        print("Closing")
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

    def on_all_next(self, e=None):
        self.getfromgui()
        for self.pos in range(len(self.datalist)):
            print(self.pos, self.params)
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
        if type(self.directory) is not str and type(self.directory) is not os.PathLike:
            d = self.directory[self.pos]
        else:
            d = self.directory
        name1 = os.path.join(d, "Extract2DFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(d, "Extract2DFigure2.png")
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
        if type(self.directory) is not str and type(self.directory) is not os.PathLike:
            d = self.directory[self.pos]
        else:
            d = self.directory
        name1 = os.path.join(d, "Extract2DFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(d, "Extract2DFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2

    def on_fit(self, e):
        peaks = ud.peakdetect(self.data1d, window=2)
        print("Peaks:", peaks[:, 0])
        peaks = np.concatenate((peaks, [[0, np.amin(self.data1d[:, 1])]]))
        fitdat, fits = MassFitter(self.data1d, peaks, 3, "microguess").perform_fit()
        print("Fits:", fits[:, 0])

        self.plot1.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)
        self.plot1.repaint()

    def on_fit2(self, e):
        fits, fitdat = ud.poisson_fit(self.data1d[:, 0], self.data1d[:, 1])
        print("Fits:", fits)
        self.plot1.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)
        self.plot1.repaint()

    def on_fit3(self, e):
        fits, fitdat = ud.binomial_fit(self.data1d[:, 0], self.data1d[:, 1])
        print("Fits:", fits)
        self.plot1.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)
        self.plot1.repaint()

    def on_fit4(self, e):
        dialog = miscwindows.SingleInputDialog(self)
        dialog.initialize_interface(title="Possible Oligomers", message="Potential Oligomers (comma separated): ")
        dialog.ShowModal()
        self.run_multip(dialog.value)

    def run_multip(self, string):
        try:
            array = string.split(",")
            array = np.array(array).astype(int)
            fits, fitdat, i1, i2 = ud.complex_poisson(self.data1d, array, background=True)
            print("Fits:", fits)
            self.plot1.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)
            self.plot1.subplot1.bar(array, i2 / np.amax(i2) * np.amax(self.data1d[:, 1]))
            self.plot1.subplot1.set_ylim(0, np.amax(self.data1d[:, 1]))
            self.plot1.repaint()
        except Exception as e:
            print(e)


# Main App Execution
if __name__ == "__main__":
    # data3 = np.loadtxt(
    #    "C:\UniDecPastedSpectra\PastedSpectrum_2017_Dec_11_11_30_45_unidecfiles\PastedSpectrum_2017_Dec_11_11_30_45_mass.txt")
    data3 = np.loadtxt(
        "Z:\\mtmarty\\Data\\Others\\Miranda\\sample_data\\MC_20170904_1to1_aB_C137S_11000_520_HCD300_CAL_unidecfiles\\CorrectedMassData.txt")
    #data3 = np.loadtxt("/TestSpectra/60_unidecfiles/60_mass.txt")
    #file = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL Project\\Final\\WT_R224A Titrations\WT_R224A_25C\\rep3\\20220303_HSJ_MAE_AqpZWT_R224A_160CL_TS_25C_3_unidecfiles\\20220303_HSJ_MAE_AqpZWT_R224A_160CL_TS_25C_3_mass.txt"
    #file = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\WT_R224A Titrations\WT_R224A_25C\\rep3\\20220308_HSJ_MAE_AqpZWT_R224A_160CL_TS_25_3_unidecfiles\\20220308_HSJ_MAE_AqpZWT_R224A_160CL_TS_25_3_mass.txt"
    # file = "Z:\Group Share\Hiruni Jayasekera\\HSJ_MAE AqpZ CL Project\\Final\\R224A_R224-75A Titrations All\\R224A_R224-75A_30C\\rep3\\20220317_HSJ_AqpZR224A_R224-75A_100CL_TS_30_3_unidecfiles\\20220317_HSJ_AqpZR224A_R224-75A_100CL_TS_30_3_mass.txt"
    #data3 = np.loadtxt(file)
    datalist = [data3]

    app = wx.App(False)
    frame = Extract2DPlot(None, datalist)
    # frame.run_multip("1,2,3,4,5,6")
    app.MainLoop()

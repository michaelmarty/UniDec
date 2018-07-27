import os
from copy import deepcopy
import numpy as np
import wx

from unidec_modules import unidecstructure, plot1d, plot2d
import unidec_modules.unidectools as ud
import matplotlib.cm as cm

from iFAMS import GuiTestFun


class iFAMS_Window(wx.Frame):
    def __init__(self, parent, data, config=None, directory=""):
        wx.Frame.__init__(self, parent, title="iFAMS")  # ,size=(-1,-1))

        self.parent = parent
        self.data = data
        self.xdata = data[:, 0]
        self.ydata = data[:, 1]

        self.directory=directory

        # Set up the config file
        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config

        # Make the menu
        filemenu = wx.Menu()
        menu_save_fig_png = filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                            "Save all figures as PNG in central directory")
        menu_save_fig_pdf = filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                            "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig_png)
        self.Bind(wx.EVT_MENU, self.on_save_fig_pdf, menu_save_fig_pdf)

        self.plotmenu = wx.Menu()

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the Plots
        panel = wx.Panel(self)

        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot1d.Plot1d(panel)
        self.plot3 = plot1d.Plot1d(panel)
        self.plot4 = plot1d.Plot1d(panel)

        sizer = wx.BoxSizer(wx.VERTICAL)
        plotsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer1.Add(self.plot1, 2, wx.EXPAND)
        plotsizer1.Add(self.plot4, 0, wx.EXPAND)

        plotsizer2.Add(self.plot2, 2, wx.EXPAND)
        plotsizer2.Add(self.plot3, 0, wx.EXPAND)

        sizer.Add(plotsizer1, 1, wx.EXPAND)
        sizer.Add(plotsizer2, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Set up the controls
        self.ctlm0 = wx.TextCtrl(panel, value=str(0))
        self.ctlwindow = wx.TextCtrl(panel, value=str(0))
        controlsizer.Add(wx.StaticText(panel, label="Dummy"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm0, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label="Blank"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwindow, 0, wx.ALIGN_CENTER_VERTICAL)

        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        fftbutton = wx.Button(panel, label="FFT")
        controlsizer2.Add(fftbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_fft, fftbutton)

        calcbutton = wx.Button(panel, label="Calc. Z and subunit")
        controlsizer2.Add(calcbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_calc, calcbutton)

        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()

        self.makeplot()

    def on_close(self, e):
        self.Destroy()
        # self.MakeModal(False)

    def getfromgui(self):
        """
        Update parameters from GUI.
        :return: None
        """
        try:
            self.m0 = float(self.ctlm0.GetValue())
            try:
                self.nbins = float(self.ctlwindow.GetValue())
            except ValueError:
                self.nbins = 0
        except ValueError:
            print("Failed to get from gui")

    def makeplot(self):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.getfromgui()

        self.plot1.plotrefreshtop(self.xdata, self.ydata, "Data", "m/z (Th)",
                                  "Intensity", "", color="b", config=self.config)

    def on_fft(self, e=None):
        self.getfromgui()

        self.yfull, self.expandedspan, self.xnew, self.ftspacing, self.paddedxnew, self.ynew = GuiTestFun.plot_function(
            self.xdata, self.ydata)

        self.maxfreq = GuiTestFun.maxfreq(self, self.yfull, self.expandedspan)
        self.ftx, self.ABFT, self.FT = GuiTestFun.Fourier(self, self.maxfreq, self.yfull)
        self.refmaxtab = GuiTestFun.findmax(self.expandedspan, self.ABFT, self.ftx, 0.001, 5, 10)
        self.plot2.plotrefreshtop(self.ftx, self.ABFT, "FFT", "Frequency", "Amplitude", color='k', config=self.config)
        self.plot2.plotadddot(np.array(self.refmaxtab)[:, 0], np.array(self.refmaxtab)[:, 1],"r","o")
        self.plot2.repaint()

    def on_calc(self, e=None):

        self.newcalcX = []
        self.newcalcY = []
        self.refmaxtabCalc = np.array(self.refmaxtab)
        self.numchar = GuiTestFun.spacing(self, self.refmaxtab)
        self.omega = GuiTestFun.omega(self.refmaxtab, self.numchar)
        self.chargestates, self.chargestatesr = GuiTestFun.charge(self.refmaxtab, self.numchar, self.omega)
        self.chargestateints = [int(self.chargestatesr[i]) for i in range(0, len(self.chargestatesr))]
        for i in range(0, len(self.chargestatesr)):
            self.newcalcX.append(self.refmaxtabCalc[i, 0])
            self.newcalcY.append(self.refmaxtabCalc[i, 1])
        print(self.newcalcX)
        print("Need to add plotting")

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Figure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "Figure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Figure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "Figure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2


# Main App Execution
if __name__ == "__main__":
    dir = "C:\\Python\\UniDec\\TestSpectra\\60_unidecfiles"
    file = "60_input.dat"

    path = os.path.join(dir, file)

    data = np.loadtxt(path)
    # data = ud.datachop(data, 0,20000)

    app = wx.App(False)
    frame = iFAMS_Window(None, data)
    app.MainLoop()

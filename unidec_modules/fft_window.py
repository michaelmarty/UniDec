import os
from copy import deepcopy
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows
from unidectools import win_fft_grid

__author__ = 'Michael.Marty'





class FFTWindow(wx.Frame):
    def __init__(self, parent, data, config=None, directory=None):

        wx.Frame.__init__(self, parent, title="FFT Analysis")

        # Setup initial values
        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.discreteplot = 0
        else:
            self.config = config

        self.config.publicationmode = 0

        self.window_fwhm = 500
        self.binsize = 0.5
        self.wbin = 200
        self.diffrange = [740, 770]

        self.rawdata = data
        self.rawdata[:, 1] /= np.amax(self.rawdata[:, 1])
        self.xlims = [min(data[:, 0]), max(data[:, 0])]

        # Make the menu
        filemenu = wx.Menu()
        menu_save_fig_png = filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                            "Save all figures as PNG in central directory")
        menu_save_fig_pdf = filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                            "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig_png)
        self.Bind(wx.EVT_MENU, self.on_save_fig_pdf, menu_save_fig_pdf)

        self.plotmenu = wx.Menu()
        self.menuaddline = self.plotmenu.Append(wx.ID_ANY, "Add Horizontal Line",
                                                "Add Horizontal Line at Specific Y Value")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the GUI
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot2d.Plot2d(panel)
        sizer.Add(self.plot1, 1, wx.EXPAND)
        sizer.Add(self.plot2, 1, wx.EXPAND)
        self.plot1._axes = [0.1, 0.1, 0.64, 0.8]

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlfwhm = wx.TextCtrl(panel, value=str(self.window_fwhm))
        self.ctlwbin = wx.TextCtrl(panel, value=str(self.wbin))
        self.ctlbinsize = wx.TextCtrl(panel, value=str(self.binsize))
        controlsizer.Add(wx.StaticText(panel, label=" Window FWHM:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlfwhm, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Window every:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwbin, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Linearization Bin Size:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlbinsize, 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlmin = wx.TextCtrl(panel, value=str(self.diffrange[0]))
        self.ctlmax = wx.TextCtrl(panel, value=str(self.diffrange[1]))
        controlsizer2.Add(wx.StaticText(panel, label=" Min Difference:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlmin, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(wx.StaticText(panel, label=" Max Difference:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlmax, 0, wx.ALIGN_CENTER_VERTICAL)



        label = "Replot"
        replotbutton = wx.Button(panel, label=label)
        controlsizer2.Add(replotbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.makeplot, replotbutton)

        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.makeplot()
        self.Centre()
        # self.MakeModal(True)
        self.Show(True)

    def on_close(self, e):
        """
        Close the window. Will try to set self.config.molig as self.m0.
        :param e: Unused event
        :return: None
        """
        try:
            self.config.molig = self.window_fwhm
        except AttributeError:
            pass
        self.Destroy()
        # self.MakeModal(False)

    def getfromgui(self):
        """
        Update parameters from GUI.
        :return: None
        """
        try:
            self.window_fwhm = float(self.ctlfwhm.GetValue())
            try:
                self.binsize = float(self.ctlbinsize.GetValue())
            except ValueError:
                self.binsize = 1
                self.ctlbinsize.SetValue(str(self.binsize))

            try:
                self.wbin = float(self.ctlwbin.GetValue())
            except ValueError:
                self.wbin = 500
                self.ctlwbin.SetValue(str(self.wbin))

            try:
                self.diffrange[0] = float(self.ctlmin.GetValue())
                self.diffrange[1] = float(self.ctlmax.GetValue())
            except ValueError:
                pass
        except ValueError:
            print "Failed to get from gui"

    def makeplot(self, e=None):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.getfromgui()
        data2d = win_fft_grid(self.rawdata, self.binsize, self.wbin, self.window_fwhm, self.diffrange)
        self.config.cmap = 'jet'
        try:
            self.plot2.contourplot(data2d, self.config, xlab="m/z (Th)", ylab="Mass Difference", title="", normflag=1)
            self.plot2.subplot1.set_xlim(self.xlims)
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.plotrefreshtop(self.rawdata[:, 0], self.rawdata[:, 1], "", "m/z (Th)",
                                      "Intensity", "", self.config)
            self.plot1.subplot1.set_xlim(self.xlims)
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "FFTFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "FFTFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "FFTFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "FFTFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2

    def on_add_line(self, e):
        """
        Add a horizontal line to the plot to visualize a predicted mass defect.
        Opens a dialog to input the value. Can be called more than once to add multiple lines.
        :param e: Unused event
        :return: None
        """
        dialog = miscwindows.SingleInputDialog(self)
        dialog.initialize_interface(title="Add Line", message="Add Line at Value: ")
        dialog.ShowModal()

        try:
            vval = float(dialog.value)
            xlim = self.plot2.subplot1.get_xlim()
            self.plot2.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot2.tickcolor)
            self.plot2.repaint()
            print xlim, vval
        except Exception, e:
            print "Failed: ", dialog.value, e
            pass


# Main App Execution
if __name__ == "__main__":
    datfile = "C:\\NDData\\PG25\\CG_07\\150820_CG_07_ramp90.txt"

    data2 = np.loadtxt(datfile)

    app = wx.App(False)
    frame = FFTWindow(None, data2)
    app.MainLoop()

import os
from copy import deepcopy
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows, fitting
from unidec_modules.unidectools import double_fft_diff, nearest, peakdetect, mergedata, datachop
import matplotlib.cm as cm

__author__ = 'Michael.Marty'


class FFTWindow(wx.Frame):
    def __init__(self, parent, datalist, xvals=None, config=None, directory=None):

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

        self.binsize = 0.5
        self.diffrange = [740, 770]

        if xvals is None:
            xvals = np.arange(0, len(datalist))
        self.xvals = xvals

        self.datalist = np.array(datalist)
        mins = []
        maxes = []
        for d in self.datalist:
            d[:, 1] /= np.amax(d[:, 1])
            mins.append(np.amin(d[:, 0]))
            maxes.append(np.amax(d[:, 0]))
        self.xlims = [min(mins), max(maxes)]
        self.colors = self.config.get_colors(len(self.datalist))

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

        #self.menupeaks = self.plotmenu.Append(wx.ID_ANY, "Get Peaks",
        #                                      "Get Peaks from Spectrum")
        #self.Bind(wx.EVT_MENU, self.on_get_peaks, self.menupeaks)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the GUI
        displaysize = wx.GetDisplaySize()
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        plotpanel1 = wx.BoxSizer(wx.HORIZONTAL)
        plotpanel2 = wx.BoxSizer(wx.HORIZONTAL)
        if displaysize[0] < 1700:
            figsize = (4, 3)
        else:
            figsize = (5, 4)
        self.plot1 = plot1d.Plot1d(panel, figsize=figsize)
        self.plot2 = plot2d.Plot2d(panel, figsize=figsize)
        self.plot3 = plot1d.Plot1d(panel, figsize=figsize)
        self.plot4 = plot1d.Plot1d(panel, figsize=figsize)
        plotpanel1.Add(self.plot1, 1, flag=wx.EXPAND)
        plotpanel2.Add(self.plot2, 1, flag=wx.EXPAND)
        plotpanel1.Add(self.plot3, 1, flag=wx.EXPAND)
        plotpanel2.Add(self.plot4, 1, flag=wx.EXPAND)
        self.plot1._axes = [0.1, 0.1, 0.64, 0.8]

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlminmz = wx.TextCtrl(panel, value=str(self.xlims[0]))
        self.ctlmaxmz = wx.TextCtrl(panel, value=str(self.xlims[1]))
        self.ctlbinsize = wx.TextCtrl(panel, value=str(self.binsize))
        controlsizer.Add(wx.StaticText(panel, label=" Min m/z:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlminmz, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Max m/z:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlmaxmz, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Linearization Bin Size:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlbinsize, 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlmin = wx.TextCtrl(panel, value=str(self.diffrange[0]))
        self.ctlmax = wx.TextCtrl(panel, value=str(self.diffrange[1]))
        controlsizer2.Add(wx.StaticText(panel, label=" Min Difference:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlmin, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(wx.StaticText(panel, label=" Max Difference:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlmax, 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlsep = wx.TextCtrl(panel, value=str(self.config.separation))
        controlsizer2.Add(wx.StaticText(panel, label=" Plot Separation:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlsep, 0, wx.ALIGN_CENTER_VERTICAL)

        label = "Replot"
        replotbutton = wx.Button(panel, label=label)
        controlsizer2.Add(replotbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.makeplot, replotbutton)

        sizer.Add(plotpanel1, 0, wx.EXPAND)
        sizer.Add(plotpanel2, 0, wx.EXPAND)
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
        self.Destroy()
        # self.MakeModal(False)

    def getfromgui(self):
        """
        Update parameters from GUI.
        :return: None
        """
        try:

            try:
                self.binsize = float(self.ctlbinsize.GetValue())
            except ValueError:
                self.binsize = 1
                self.ctlbinsize.SetValue(str(self.binsize))

            try:
                self.maxmz = float(self.ctlmaxmz.GetValue())
                self.minmz = float(self.ctlminmz.GetValue())
            except ValueError:
                self.minmz = self.xlims[0]
                self.maxmz = self.xlims[1]
                self.ctlmaxmz.SetValue(str(self.maxmz))
                self.ctlminmz.SetValue(str(self.minmz))

            try:
                self.diffrange[0] = float(self.ctlmin.GetValue())
                self.diffrange[1] = float(self.ctlmax.GetValue())
            except ValueError:
                pass

            try:
                self.config.separation = float(self.ctlsep.GetValue())
            except ValueError:
                pass
        except ValueError:
            print("Failed to get from gui")

    def fft_process(self, e=None):
        maxes = []
        self.fftdat = []
        for i, d in enumerate(self.datalist):
            maxpos, fft2 = double_fft_diff(datachop(d,self.minmz,self.maxmz), diffrange=self.diffrange, binsize=self.binsize)
            maxes.append(maxpos)
            fft2[:, 1] -= np.amin(fft2[:, 1])
            fft2[:, 1] /= np.amax(fft2[:, 1])
            if i > 0:
                fft2 = mergedata(self.fftdat[0], fft2)
            self.fftdat.append(fft2)
        self.peakvals = np.array(maxes)
        self.fftdat = np.array(self.fftdat)
        print(self.peakvals)

    def makeplot(self, e=None):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.getfromgui()

        try:
            for i, d in enumerate(self.datalist):
                d2=datachop(d,self.minmz,self.maxmz)
                if i == 0:
                    self.plot1.plotrefreshtop(d2[:, 0], d2[:, 1], xlabel="m/z (Th)",
                                              ylabel="Intensity", config=self.config, color=self.colors[i])
                else:
                    self.plot1.plotadd(d2[:, 0], d2[:, 1] - i * self.config.separation, colval=self.colors[i])
            #self.plot1.subplot1.set_xlim(self.xlims)
            self.plot1.repaint()
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)

        self.fft_process()
        try:
            for i, d in enumerate(self.fftdat):
                if i == 0:
                    self.plot3.plotrefreshtop(d[:, 0], d[:, 1], xlabel="Mass Difference (Da)",
                                              ylabel="Intensity", config=self.config, color=self.colors[i])
                else:
                    self.plot3.plotadd(d[:, 0], d[:, 1] - i * self.config.separation, colval=self.colors[i])
            self.plot3.repaint()
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot3", e)

        try:
            self.plot4.plotrefreshtop(self.xvals, self.peakvals, marker="o", config=self.config,
                                      xlabel="Variable 1", ylabel="Mass Difference (Da)")
        except Exception as e:
            self.plot4.clear_plot()
            print("Failed Plot4", e)

        try:
            x = self.xvals
            y = self.fftdat[0, :, 0]
            z = self.fftdat[:, :, 1]
            #print(x.shape, y.shape, z.shape)
            self.plot2.contourplot(xvals=x, yvals=y, zgrid=z,config=self.config, xlab="Variable 1", ylab="Mass Difference (Da)")
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2",e)

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
        name3 = os.path.join(self.directory, "FFTFigure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
        name4 = os.path.join(self.directory, "FFTFigure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name4)

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
        name3 = os.path.join(self.directory, "FFTFigure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
        name4 = os.path.join(self.directory, "FFTFigure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name4)

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
            print(xlim, vval)
        except Exception as e:
            print("Failed: ", dialog.value, e)
            pass



# Main App Execution
if __name__ == "__main__":
    datfile = "C:\\Data\\New\POPC_D1T0-2m_ISTRAP\\20170207_P1D_POPC_ND_D1T0-2m_ISTRAP_RAMP_0_275_25_1_200.0.txt"

    data2 = np.loadtxt(datfile)
    data3 = deepcopy(data2)
    data3[:, 0] = data3[:, 0] - 100

    datalist = [data2, data3]

    app = wx.App(False)
    frame = FFTWindow(None, datalist)

    app.MainLoop()

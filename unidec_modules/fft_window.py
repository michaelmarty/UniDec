import os
from copy import deepcopy
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows, fitting
from unidec_modules.unidectools import win_fft_grid, nearest, peakdetect
import matplotlib.cm as cm

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
        self.diffrange = [700, 800]
        self.norm=True

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

        self.menupeaks = self.plotmenu.Append(wx.ID_ANY, "Get Peaks",
                                              "Get Peaks from Spectrum")
        self.Bind(wx.EVT_MENU, self.on_get_peaks, self.menupeaks)

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

        self.ctlfwhm = wx.TextCtrl(panel, value=str(self.window_fwhm))
        self.ctlwbin = wx.TextCtrl(panel, value=str(self.wbin))
        self.ctlbinsize = wx.TextCtrl(panel, value=str(self.binsize))
        self.ctlnorm = wx.CheckBox(panel, label="Normalize FFT to Data")
        self.ctlnorm.SetValue(self.norm)
        controlsizer.Add(wx.StaticText(panel, label=" Window FWHM:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlfwhm, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Window every:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwbin, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label=" Linearization Bin Size:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlbinsize, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlnorm, 0, wx.ALIGN_CENTER_VERTICAL)


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
        comparebutton = wx.Button(panel, label="Compare")
        controlsizer2.Add(comparebutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_compare_regions, comparebutton)
        self.compareclicks = 0
        self.comparetext = wx.StaticText(panel, label="Click to activate compare mode")
        controlsizer2.Add(self.comparetext, 0, wx.ALIGN_CENTER_VERTICAL)

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

            self.norm=self.ctlnorm.GetValue()
        except ValueError:
            print("Failed to get from gui")

    def makeplot(self, e=None):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.getfromgui()
        data2d = win_fft_grid(self.rawdata, self.binsize, self.wbin, self.window_fwhm, self.diffrange, norm=self.norm)
        self.config.cmap = u'jet'
        try:
            self.plot2.contourplot(data2d, self.config, xlab="m/z (Th)", ylab="Mass Difference", title="", normflag=1)
            self.plot2.subplot1.set_xlim(self.xlims)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)
        try:
            self.plot1.plotrefreshtop(self.rawdata[:, 0], self.rawdata[:, 1], "", "m/z (Th)",
                                      "Intensity", "", self.config)
            self.plot1.subplot1.set_xlim(self.xlims)
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)
        try:
            indexdiff = ((self.diffrange[1] - self.diffrange[0]) / self.binsize) - 1
            rowsums = []
            rowdiff = []
            for row in range(0, int(indexdiff)):
                sum = 0
                for col in range(0, int(len(data2d) / indexdiff)):
                    sum += data2d[int(col * indexdiff + row), 2]
                rowsums.append(sum)
                rowdiff.append(self.diffrange[0] + self.binsize + (self.binsize * row))
            maxsum = np.amax(np.asarray(rowsums))
            tmp = [x / maxsum for x in rowsums]
            self.diffdat = np.transpose([rowdiff, tmp])
            self.plot3.plotrefreshtop(rowdiff, tmp, xlabel="Mass Difference", ylabel="Intensity")
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot3", e)

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

    def on_get_peaks(self, e=None, data=None):
        if data is None:
            data=self.diffdat
        print("Data range", np.amin(data[:,0]), "to", np.amax(data[:,0]))
        peaks = peakdetect(data, window=1000.)[:, 0]
        for p in peaks:
            index = nearest(data[:, 0], p)
            idiff = int(3 / self.binsize)
            print("Peak value:", p)
            fit, fitdat = fitting.isolated_peak_fit(data[index - idiff:index + idiff, 0],
                                                    data[index - idiff:index + idiff, 1], psfun=0)
            print("Peak Fit:", fit[1, 0], "+/-", fit[1, 1])

    def on_compare_regions(self, e=None):
        # First click
        if self.compareclicks == 0:
            self.plot2.zoom.comparemode = True
            self.compareclicks += 1
            self.comparetext.SetLabel("Drag box(es) and click again to compare")
        # Turn off if no boxes were drawn
        elif self.compareclicks == 1 and len(self.plot2.zoom.comparexvals) == 0:
            self.compareclicks = 0
            self.plot2.zoom.comparexvals = []
            self.plot2.zoom.compareyvals = []
            self.comparetext.SetLabel("Compare Mode Deactivated")
            self.plot2.zoom.comparemode = False
        # Admonish user for not using it correctly
        elif len(self.plot2.zoom.comparexvals) % 2 != 0:
            self.comparetext.SetLabel("Please select 2 or more regions to compare")
            self.compareclicks = 1
            self.plot2.zoom.comparexvals = []
            self.plot2.zoom.compareyvals = []
        # 2 boxes drawn, 2 clicks
        elif len(self.plot2.zoom.comparexvals) % 2 == 0:
            ret = self.createcompareplot()
            if ret == 1:
                self.comparetext.SetLabel("Plot Created")
            self.compareclicks = 0
            self.plot2.zoom.compareyvals = []
            self.plot2.zoom.comparexvals = []
            self.plot2.zoom.comparemode = False
        # Case i dont know
        else:
            self.comparetext.SetLabel("Unknown case found")

    def createcompareplot(self):
        xvals = self.plot2.zoom.comparexvals
        yvals = self.plot2.zoom.compareyvals
        data2d = win_fft_grid(self.rawdata, self.binsize, self.wbin, self.window_fwhm, self.diffrange, norm=self.norm)
        for x in range(0, len(xvals) // 2):
            if xvals[x * 2] > xvals[x * 2 + 1]: xvals[x * 2], xvals[x * 2 + 1] = xvals[x * 2 + 1], xvals[x * 2]
            if yvals[x * 2] > yvals[x * 2 + 1]: yvals[x * 2], yvals[x * 2 + 1] = yvals[x * 2 + 1], yvals[x * 2]
            # assure that x and y values are not equal
            if xvals[x * 2] == xvals[x * 2 + 1] or yvals[x * 2] == yvals[x * 2 + 1]:
                self.comparetext.SetLabel("Line or point drawn, please try again.")
                return 0
        # Round to nearest value, find index for that value
        nearestxvalind = []
        for val in xvals:
            roundedval = self.wbin * round(val / self.wbin)
            nearestxvalind.append(nearest(data2d[:, 0], roundedval))
        nearestyvalind = []
        indexdiff = ((self.diffrange[1] - self.diffrange[0]) / self.binsize) - 1
        for val in yvals:
            roundedval = self.binsize * round(val / self.binsize)
            nearestyvalind.append(nearest(data2d[0:int(indexdiff), 1], roundedval))
        # Sum up the rows in each box
        boxsums = []
        rowdiffs = []
        maxsum = 0
        for x in range(0, len(xvals) // 2):
            rowsums = []
            tmpdiff = []
            for row in range(nearestyvalind[x * 2], nearestyvalind[x * 2 + 1]):
                sum = 0
                for col in range(nearestxvalind[x * 2], nearestxvalind[x * 2 + 1], int(indexdiff)):
                    sum += data2d[int(col + row), 2]
                rowsums.append(sum)
                tmpdiff.append(self.diffrange[0] + self.binsize + (self.binsize * row))
                if sum > maxsum:
                    maxsum = sum
            boxsums.append(rowsums)
            rowdiffs.append(tmpdiff)
        colormap = cm.get_cmap('rainbow', len(xvals) / 2)
        cols = colormap(np.arange(len(xvals) / 2))
        tmp = [y / maxsum for y in boxsums[0]]
        self.on_get_peaks(data=np.transpose([rowdiffs[0],tmp]))
        self.plot4.plotrefreshtop(rowdiffs[0], tmp, color=cols[0],
                                  xlabel="Mass Difference", ylabel="Intensity", linestyle="solid")
        for x in range(1, len(xvals) // 2):
            tmp = [y / maxsum for y in boxsums[x]]
            self.on_get_peaks(data=np.transpose([rowdiffs[x], tmp]))
            self.plot4.plotadd(rowdiffs[x], tmp, cols[x], None)
        self.plot4.repaint()
        return 1


# Main App Execution
if __name__ == "__main__":
    datfile = "C:\\Data\\New\POPC_D1T0-2m_ISTRAP\\20170207_P1D_POPC_ND_D1T0-2m_ISTRAP_RAMP_0_275_25_1_200.0.txt"
    datfile = "C:\\Python\\UniDec3\\TestSpectra\\60.dat"

    data2 = np.loadtxt(datfile)

    app = wx.App(False)
    frame = FFTWindow(None, data2)

    frame.on_get_peaks()
    app.MainLoop()

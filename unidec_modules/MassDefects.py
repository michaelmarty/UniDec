import os
from copy import deepcopy
import numpy as np
import wx
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows
import unidec_modules.unidectools as ud
from MassFitter import MassFitter

__author__ = 'Michael.Marty'


class MassDefectWindow(wx.Frame):
    def __init__(self, parent, data_list, config=None, yvals=None, pks=None, value=None, directory=None):
        """
        Creates a window for visualizing high-mass mass defects.
        :param parent: Passed to wx.Frame
        :param data_list: List of mass distribution data.
        :param config: UniDecConfig object
        :param yvals: List of titles for each mass distribution in data_list
        :param pks: Peaks object
        :param value: Kendrick reference mass (default is 760.076, the mass of POPC)
        :param directory: Directory to save files to (default is os.getcwd())
        :return: None
        """
        wx.Frame.__init__(self, parent, title="Mass Defect")  # ,size=(-1,-1))

        # Setup initial values
        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        self.parent = parent
        self.m0 = deepcopy(value)
        if self.m0 is None:
            self.m0 = 760.076

        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.discreteplot = 0
            self.config.cmap = "jet"
        else:
            self.config = config

        self.config.publicationmode = 0
        self.pos = -1
        self.yvals = yvals
        self.ylab = "Normalized Mass Defect"
        self.nbins = 50
        self.transformmode = 1
        self.centermode = 1
        self.xtype = 1
        self.factor = 1
        self.xlab = ""

        self.datalist = [data_list[0]]
        for i in range(1, len(data_list)):
            self.datalist.append(ud.mergedata(data_list[0], data_list[i]))
        self.datalist = np.array(self.datalist)
        self.datasum = np.transpose([self.datalist[0, :, 0], np.sum(self.datalist[:, :, 1], axis=0)])
        print "Data list shape:", self.datalist.shape

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
        self.menufit = self.plotmenu.Append(wx.ID_ANY, "Fit Peaks",
                                            "Fit total mass defect peaks")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)
        self.Bind(wx.EVT_MENU, self.on_fit, self.menufit)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the GUI
        panel = wx.Panel(self)

        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot2d.Plot2d(panel)
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

        self.ctlm0 = wx.TextCtrl(panel, value=str(self.m0))
        self.ctlwindow = wx.TextCtrl(panel, value=str(self.nbins))
        controlsizer.Add(wx.StaticText(panel, label="Kendrick Mass"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlm0, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(wx.StaticText(panel, label="Number of Defect Bins"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwindow, 0, wx.ALIGN_CENTER_VERTICAL)

        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        if len(data_list) > 1:
            label = "Back"
        else:
            label = "Replot"
        backbutton = wx.Button(panel, label=label)
        controlsizer2.Add(backbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_back, backbutton)
        if len(data_list) > 1:
            nextbutton = wx.Button(panel, label="Next")
            controlsizer2.Add(nextbutton, 0, wx.EXPAND)
            self.Bind(wx.EVT_BUTTON, self.on_next, nextbutton)
            totalbutton = wx.Button(panel, label="Total")
            controlsizer2.Add(totalbutton, 0, wx.EXPAND)
            self.Bind(wx.EVT_BUTTON, self.makeplottotal, totalbutton)
        else:
            if not ud.isempty(pks):
                self.pks = pks
                peaksbutton = wx.Button(panel, label="Plot Peaks")
                controlsizer2.Add(peaksbutton, 0, wx.EXPAND)
                self.Bind(wx.EVT_BUTTON, self.on_peaks, peaksbutton)

        self.radiobox = wx.RadioBox(panel, choices=["Integrate", "Interpolate"], label="Type of Transform")
        controlsizer.Add(self.radiobox, 0, wx.EXPAND)
        self.radiobox.SetSelection(self.transformmode)

        self.radiobox2 = wx.RadioBox(panel, choices=["-0.5:0.5", "0:1"], label="Range")
        controlsizer.Add(self.radiobox2, 0, wx.EXPAND)
        self.radiobox2.SetSelection(self.centermode)

        self.radiobox3 = wx.RadioBox(panel, choices=["Mass Number", "Mass"], label="X-Axis")
        controlsizer.Add(self.radiobox3, 0, wx.EXPAND)
        self.radiobox3.SetSelection(self.xtype)

        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        try:
            self.makeplottotal(0)
        except Exception, e:
            self.on_next(0)
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
            self.config.molig = self.m0
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
            self.m0 = float(self.ctlm0.GetValue())
            try:
                self.nbins = float(self.ctlwindow.GetValue())
            except ValueError:
                self.nbins = 0
        except ValueError:
            print "Failed to get from gui"
        self.transformmode = self.radiobox.GetSelection()
        self.centermode = self.radiobox2.GetSelection()
        self.xtype = self.radiobox3.GetSelection()
        if self.xtype == 0:
            self.factor = 1
            self.xlab = "Nominal Mass Number"
        else:
            self.factor = self.m0
            self.xlab = "Mass"

    def makeplot(self):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.getfromgui()
        self.data1d, data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(self.datalist[self.pos], self.m0,
                                                                          centermode=self.centermode,
                                                                          nbins=self.nbins,
                                                                          transformmode=self.transformmode,
                                                                          xaxistype=self.xtype)
        if self.yvals is not None:
            title = str(self.yvals[self.pos])
            spacer = "_"
        else:
            title = ""
            spacer = ""
        try:
            save_path2d = os.path.join(self.directory, title + spacer + "2D_Mass_Defects.txt")
            np.savetxt(save_path2d, data2d)
            save_path1d = os.path.join(self.directory, title + spacer + "1D_Mass_Defects.txt")
            np.savetxt(save_path1d, self.data1d)
            print 'Saved: ', save_path2d, save_path1d
        except Exception, e:
            print "Failed save", e
        try:
            self.plot2.contourplot(data2d, self.config, xlab=self.xlab, ylab=self.ylab, title=title, normflag=1)
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], title, "Mass Defect",
                                      "Total Intensity", "", self.config)
        except Exception, e:
            self.plot3.clear_plot()
            print "Failed Plot 3", e
        try:
            if self.m0 == 0:
                return
            self.plot1.colorplotMD(self.datalist[self.pos, :, 0], self.datalist[self.pos, :, 1],
                                   self.datalist[self.pos, :, 0] / float(self.m0) % 1.0,
                                   title="Zero-Charge Mass Spectrum",
                                   xlabel="Mass", ylabel="Intensity")
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0] % 1.0,
                                   title="Total Projection Color", xlabel="Mass Defect", cmap="hsv",
                                   ylabel="Total Intensity", config=self.config)
        except Exception, e:
            self.plot4.clear_plot()
            print "Failed Plot 4", e

    def makeplottotal(self, e=None):
        """
        Runs the kendrick analysis on all of the mass distributions in self.datalist.
        Assumes they are all the same dimensions, which is why we run mergedata when it is loaded.
        Sums the results and plots the sums in 1D and 2D.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()
        # Run on all
        igrids = []
        m1grid, m2grid = None, None
        for i, dat in enumerate(self.datalist):
            self.data1d, data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(dat, self.m0,
                                                                              centermode=self.centermode,
                                                                              nbins=self.nbins,
                                                                              transformmode=self.transformmode,
                                                                              xaxistype=self.xtype)
            igrids.append(igrid)
        # Sum and reshape
        igrids = np.array(igrids)
        igrids /= np.amax(igrids)
        sumgrid = np.sum(igrids, axis=0)
        data2d = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(sumgrid) / np.amax(sumgrid)])
        self.data1d = np.transpose([np.unique(m2grid), np.sum(sumgrid, axis=0)])
        # Save Results
        save_path2d = os.path.join(self.directory, "Total_2D_Mass_Defects.txt")
        np.savetxt(save_path2d, data2d)
        save_path1d = os.path.join(self.directory, "Total_1D_Mass_Defects.txt")
        np.savetxt(save_path1d, self.data1d)
        print 'Saved: ', save_path2d, save_path1d
        # Plots
        try:
            self.plot2.contourplot(data2d, self.config, xlab=self.xlab, ylab=self.ylab, title="Total", normflag=1)
        except Exception, e:
            self.plot2.clear_plot()
            print "Failed Plot2", e
        try:
            self.plot1.colorplotMD(self.datasum[:, 0], self.datasum[:, 1], self.datasum[:, 0] / float(self.m0) % 1.0,
                                   title="Zero-Charge Mass Spectrum",
                                   xlabel="Mass", ylabel="Intensity")
        except Exception, e:
            self.plot1.clear_plot()
            print "Failed Plot1", e
        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], "Total Projection Black", "Mass Defect",
                                      "Total Intensity", "", self.config)
        except Exception, e:
            self.plot3.clear_plot()
            print "Failed Plot 3", e

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0] % 1.0,
                                   title="Total Projection Color", xlabel="Mass Defect", cmap="hsv",
                                   ylabel="Total Intensity", config=self.config)
        except Exception, e:
            self.plot4.clear_plot()
            print "Failed Plot 4", e

    def on_back(self, e):
        """
        Plot the mass defect of the previous thing in the mass defect list.
        :param e: Unused event
        :return: None
        """
        self.pos += -1
        self.pos %= len(self.datalist)
        self.makeplot()

    def on_next(self, e):
        """
        Plot the mass defect of the next thing in the mass defect list.
        :param e: Unused event
        :return: None
        """
        self.pos += 1
        self.pos %= len(self.datalist)
        self.makeplot()

    def on_peaks(self, e):
        """
        For each peak in self.pks, get the mass defects and plot them.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()
        flag = self.radiobox2.GetSelection()
        self.pks.get_mass_defects(self.m0, mode=flag)

        # Plot the mass defect peaks
        self.plot3.clear_plot()
        xvals = []
        yvals = []
        for p in self.pks.peaks:
            x = p.kendricknum * self.factor
            y = p.kendrickdefect
            if not self.plot3.flag:
                self.plot3.plotrefreshtop(x, y, "Mass Peaks", self.xlab, self.ylab, "", self.config, color=p.color,
                                          marker=p.marker)
            else:
                self.plot3.plotadddot(x, y, p.color, p.marker)
            xvals.append(x)
            yvals.append(y)
        datalims = [np.amin(xvals), np.amin(yvals), np.amax(xvals), np.amax(yvals)]
        self.plot3.setup_zoom([self.plot3.subplot1], "box", data_lims=datalims)
        self.plot3.subplot1.set_ylim(self.plot2.subplot1.get_ylim())
        self.plot3.subplot1.set_xlim(self.plot2.subplot1.get_xlim())
        self.plot3.repaint()

        # Save to txt file output
        save_path = os.path.join(self.directory, "Peaks_Mass_Defects.txt")
        np.savetxt(save_path, np.transpose([xvals, yvals]))

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "MassDefectFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "MassDefectFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "MassDefectFigure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "MassDefectFigure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "MassDefectFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "MassDefectFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "MassDefectFigure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "MassDefectFigure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

    def on_add_line(self, e):
        """
        Add a horizontal line to the plot to visualize a predicted mass defect.
        Opens a dialog to input the value. Can be called more than once to add multiple lines.
        :param e: Unused event
        :return: None
        """
        dialog = miscwindows.SingleInputDialog(self)
        dialog.initialize_interface(title="Add Line", message="Add Line at Defect Value: ")
        dialog.ShowModal()

        try:
            vval = float(dialog.value)
            xlim = self.plot2.subplot1.get_xlim()
            self.plot2.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot2.tickcolor)
            self.plot2.repaint()

            xlim2 = self.plot3.subplot1.get_xlim()
            if xlim2[1] > 1:
                self.plot3.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot3.tickcolor)
                self.plot3.repaint()
        except Exception, e:
            print "Failed: ", dialog.value, e
            pass

    def on_fit(self, e):
        peaks = ud.peakdetect(self.data1d, window=3)
        print "Peaks:", peaks[:, 0]
        peaks = np.concatenate((peaks, [[0, np.amin(self.data1d[:, 1])]]))
        fitdat, fits = MassFitter(self.data1d, peaks, 3, "microguess").perform_fit()
        print "Fits:", fits[:, 0]

        self.plot3.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)


# Main App Execution
if __name__ == "__main__":
    dir = "C:\\Python\\UniDec\\TestSpectra\\60_unidecfiles"
    file = "60_mass.txt"

    path = os.path.join(dir, file)

    data = np.loadtxt(path)

    # dir = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_unidecfiles"
    # file = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_mass.txt"

    # path = os.path.join(dir, file)

    # data2 = np.loadtxt(path)
    datalist = [data]  # , data2]

    app = wx.App(False)
    frame = MassDefectWindow(None, datalist)
    app.MainLoop()

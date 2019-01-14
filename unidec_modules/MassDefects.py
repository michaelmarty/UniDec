import os
from copy import deepcopy
import numpy as np
import wx
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows, MassDefectExtractor
import unidec_modules.unidectools as ud
from unidec_modules.MassFitter import MassFitter
import matplotlib.cm as cm

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
            self.config.cmap = u"jet"
            self.config.peakcmap = u"rainbow"
            self.config.separation = 0.025
        else:
            self.config = config

        self.config.publicationmode = 0
        self.pos = -1
        self.yvals = yvals

        self.ylab = "Normalized Mass Defect"
        if self.config.defectparams is not None:
            p = self.config.defectparams
            self.nbins = p[0]
            self.transformmode = p[1]
            self.centermode = p[2]
            self.xtype = p[3]
        else:
            self.nbins = 50
            self.transformmode = 1
            self.centermode = 1
            self.xtype = 1
        self.factor = 1
        self.xlab = ""
        self.outfname = os.path.splitext(self.config.filename)[0]
        if self.outfname is not "":
            self.outfname += "_"

        try:
            self.datalist = [data_list[0]]
        except:
            self.datalist = data_list[0]

        for i in range(1, len(data_list)):
            self.datalist.append(ud.mergedata(data_list[0], data_list[i]))
        self.datalist = np.array(self.datalist)

        if self.yvals is not None:
            try:
                self.yvals = np.array(self.yvals, dtype="float")
            except:
                self.yvals = np.arange(0, len(self.datalist))

        self.datasum = np.transpose([self.datalist[0, :, 0], np.sum(self.datalist[:, :, 1], axis=0)])
        print("Data list shape:", self.datalist.shape)

        # Make the menu
        filemenu = wx.Menu()
        if self.datalist.shape[0] > 1:
            extractwindow = filemenu.Append(wx.ID_ANY, "Extract Mass Defect Values",
                                            "Open Window to Extract Mass Defect Values")
            self.Bind(wx.EVT_MENU, self.on_extract_window, extractwindow)
            filemenu.AppendSeparator()

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
        self.menupeaks = self.plotmenu.Append(wx.ID_ANY, "Label Peaks",
                                              "Label peaks")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)
        self.Bind(wx.EVT_MENU, self.on_fit, self.menufit)
        self.Bind(wx.EVT_MENU, self.on_label_peaks, self.menupeaks)

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

        if self.datalist.shape[0] > 1:
            self.flag2 = True
            self.plot5 = plot1d.Plot1d(panel)
            self.plot6 = plot2d.Plot2d(panel)
        else:
            self.flag2 = False
            self.plot5 = None
            self.plot6 = None

        sizer = wx.BoxSizer(wx.VERTICAL)
        plotsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer1.Add(self.plot1, 2, wx.EXPAND)
        plotsizer1.Add(self.plot4, 0, wx.EXPAND)

        plotsizer2.Add(self.plot2, 2, wx.EXPAND)
        plotsizer2.Add(self.plot3, 0, wx.EXPAND)

        if self.flag2:
            plotsizer1.Add(self.plot5, 0, wx.EXPAND)
            plotsizer2.Add(self.plot6, 0, wx.EXPAND)
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
        except Exception as e:
            self.on_next(0)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()

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
            print("Failed to get from gui")
        self.transformmode = self.radiobox.GetSelection()
        self.centermode = self.radiobox2.GetSelection()
        self.xtype = self.radiobox3.GetSelection()
        if self.xtype == 0:
            self.factor = 1
            self.xlab = "Nominal Mass Number"
        else:
            self.factor = self.m0
            self.xlab = "Mass"
        self.config.defectparams = [self.nbins, self.transformmode, self.centermode, self.xtype]

    def make_list_plots(self):
        print(self.igrids.shape)
        self.colormap = cm.get_cmap(ud.smartdecode(self.config.peakcmap), len(self.datalist))
        if self.colormap is None:
            self.colormap = cm.get_cmap(u"rainbow", len(self.datalist))
        self.peakcolors = self.colormap(np.arange(len(self.datalist)))

        dat3 = []
        for i, dat in enumerate(self.igrids):
            dat2 = np.sum(dat, axis=0)
            dat2 = dat2 - np.amin(dat2)
            dat2 /= np.amax(dat2)
            dat3.append(dat2)
            try:
                if i == 0:
                    self.plot5.plotrefreshtop(self.data1d[:, 0], dat2 - self.config.separation * i, "All Data",
                                              "Mass Defect",
                                              "Total Intensity", "", color=self.peakcolors[i], config=self.config)
                else:
                    self.plot5.plotadd(self.data1d[:, 0], dat2 - self.config.separation * i, colval=self.peakcolors[i])
            except Exception as e:
                self.plot5.clear_plot()
                print("Failed Plot 5", e)
        self.plot5.repaint()

        dat3 = np.array(dat3)
        self.dat3 = dat3
        if self.yvals is None:
            yvals = np.arange(0, len(self.datalist))
            self.yvals = yvals
        else:
            yvals = self.yvals

        m1grid, m2grid = np.meshgrid(self.data1d[:, 0], yvals, indexing='ij')
        data2 = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(dat3.transpose())])
        try:
            save_path2d = os.path.join(self.directory, self.outfname + "Mass_Defect_Grid.txt")
            np.savetxt(save_path2d, dat3)
            print('Saved: ', save_path2d)
        except Exception as e:
            print("Failed Data Export 6", e)

        try:
            self.plot6.contourplot(data2, self.config, xlab="Mass Defect", ylab="Individual Spectra", title="",
                                   normflag=1)
        except Exception as e:
            self.plot6.clear_plot()
            print("Failed Plot2", e)

        pass

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
            title = self.outfname + str(self.yvals[self.pos])
            spacer = "_"
        else:
            title = self.outfname
            spacer = ""
        try:
            save_path2d = os.path.join(self.directory, title + spacer + "2D_Mass_Defects.txt")
            np.savetxt(save_path2d, data2d)
            save_path1d = os.path.join(self.directory, title + spacer + "1D_Mass_Defects.txt")
            np.savetxt(save_path1d, self.data1d)
            print('Saved: ', save_path2d, save_path1d)
        except Exception as e:
            print("Failed save", e)
        try:
            self.plot2.contourplot(data2d, self.config, xlab=self.xlab, ylab=self.ylab, title=title, normflag=1)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)
        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], title, "Mass Defect",
                                      "Total Intensity", "", self.config)
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot 3", e)
        try:
            if self.m0 == 0:
                return
            self.plot1.colorplotMD(self.datalist[self.pos, :, 0], self.datalist[self.pos, :, 1],
                                   self.datalist[self.pos, :, 0] / float(self.m0) % 1.0,
                                   title="Zero-Charge Mass Spectrum",
                                   xlabel="Mass", ylabel="Intensity")
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0] % 1.0,
                                   title="Total Projection Color", xlabel="Mass Defect", cmap="hsv",
                                   ylabel="Total Intensity", config=self.config)
        except Exception as e:
            self.plot4.clear_plot()
            print("Failed Plot 4", e)

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
        self.igrids = igrids
        data2d = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(sumgrid) / np.amax(sumgrid)])
        self.data1d = np.transpose([np.unique(m2grid), np.sum(sumgrid, axis=0)])
        # Save Results
        save_path2d = os.path.join(self.directory, self.outfname + "Total_2D_Mass_Defects.txt")
        np.savetxt(save_path2d, data2d)
        save_path1d = os.path.join(self.directory, self.outfname + "Total_1D_Mass_Defects.txt")
        np.savetxt(save_path1d, self.data1d)
        print('Saved: ', save_path2d, save_path1d)
        # Plots
        try:
            self.plot2.contourplot(data2d, self.config, xlab=self.xlab, ylab=self.ylab, title="Total", normflag=1)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)
        try:
            self.plot1.colorplotMD(self.datasum[:, 0], self.datasum[:, 1], self.datasum[:, 0] / float(self.m0) % 1.0,
                                   title="Zero-Charge Mass Spectrum",
                                   xlabel="Mass", ylabel="Intensity")
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)
        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], "Total Projection Black", "Mass Defect",
                                      "Total Intensity", "", self.config)
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot 3", e)

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0] % 1.0,
                                   title="Total Projection Color", xlabel="Mass Defect", cmap="hsv",
                                   ylabel="Total Intensity", config=self.config)
        except Exception as e:
            self.plot4.clear_plot()
            print("Failed Plot 4", e)

        if self.flag2:
            self.make_list_plots()

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

    def on_extract_window(self, e=None):
        print(self.yvals)
        self.config.kendrickmass = self.m0
        frame = MassDefectExtractor.MassDefectExtractorWindow(self, self.dat3, self.data1d[:, 0], self.yvals,
                                                              config=self.config)

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

        try:
            name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure5.png")
            if self.plot5.flag:
                self.plot5.on_save_fig(e, name1)
                # print name1
            name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure6.png")
            if self.plot6.flag:
                self.plot6.on_save_fig(e, name2)
                # print name2
        except:
            pass

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

        try:
            name1 = os.path.join(self.directory, self.outfname + "MassDefectFigure5.pdf")
            if self.plot5.flag:
                self.plot5.on_save_fig(e, name1)
                # print name1
            name2 = os.path.join(self.directory, self.outfname + "MassDefectFigure6.pdf")
            if self.plot6.flag:
                self.plot6.on_save_fig(e, name2)
                # print name2
        except:
            pass

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

            ylim4 = self.plot3.subplot1.get_ylim()
            self.plot3.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot3.tickcolor)
            self.plot3.repaint()

            ylim4 = self.plot4.subplot1.get_ylim()
            self.plot4.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot4.tickcolor)
            self.plot4.repaint()
        except Exception as e:
            print("Failed: ", dialog.value, e)
            pass

        if len(self.datalist) > 1:
            try:
                xval = float(dialog.value)
                ylim = self.plot6.subplot1.get_ylim()
                self.plot6.subplot1.plot((xval, xval), (ylim[0], ylim[1]), color=self.plot2.tickcolor)
                self.plot6.repaint()

                ylim4 = self.plot5.subplot1.get_ylim()
                self.plot5.subplot1.plot((xval, xval), (ylim4[0], ylim4[1]), color=self.plot5.tickcolor)
                self.plot5.repaint()
            except Exception as e:
                print("Failed: ", dialog.value, e)
                pass

    def on_fit(self, e):
        peaks = ud.peakdetect(self.data1d, window=3)
        print("Peaks:", peaks[:, 0])
        peaks = np.concatenate((peaks, [[0, np.amin(self.data1d[:, 1])]]))
        fitdat, fits = MassFitter(self.data1d, peaks, 3, "microguess").perform_fit()
        print("Fits:", fits[:, 0])

        self.plot3.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)

        for f in fits[:, 0]:
            if np.amin(self.data1d[:, 0]) <= f <= np.amax(self.data1d[:, 0]):
                y = np.amax(self.data1d[ud.nearest(self.data1d[:, 0], f), 1])
                self.plot3.addtext(str(np.round(f, 3)), f + 0.075, y * 0.95, vlines=False)
                self.plot3.addtext("", f, y, vlines=True)

    def on_label_peaks(self, e=None):
        peaks = ud.peakdetect(self.data1d, window=3)
        print("Peaks:", peaks[:, 0])

        for p in peaks:
            y = p[1]
            self.plot3.addtext(str(np.round(p[0], 3)), p[0] + 0.075, y * 0.95, vlines=False)
            self.plot3.addtext("", p[0], y, vlines=True)


# Main App Execution
if __name__ == "__main__":
    dir = "C:\\Python\\UniDec\\TestSpectra\\60_unidecfiles"
    file = "60_mass.txt"
    # file = "60_input.dat"

    path = os.path.join(dir, file)

    data = np.loadtxt(path)

    dir = "C:\Python\\UniDec\TestSpectra\\180_unidecfiles"
    file = "180_mass.txt"

    path = os.path.join(dir, file)

    data2 = np.loadtxt(path)
    # datalist = [data, data2]
    datalist = [data]

    app = wx.App(False)
    frame = MassDefectWindow(None, datalist)
    frame.on_label_peaks(0)
    app.MainLoop()

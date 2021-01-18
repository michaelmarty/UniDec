import os
from copy import deepcopy
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows, MassDefectExtractor
import unidec_modules.unidectools as ud
from unidec_modules.MassFitter import MassFitter
from unidec_modules.isolated_packages.MD_compare import MassDefectCompareWindow
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
        print(self.directory)

        self.parent = parent
        self.m0 = deepcopy(value)
        if self.m0 is None:
            self.m0 = 760.076

        self.defaultvalues = [44088, 3492, 0, 20, 1]

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
            self.xtype = 0
        self.factor = 1
        self.xlab = ""
        self.outfname = self.config.outfname
        if self.outfname != "":
            self.outfname += "_"
        print(self.outfname)
        self.total = False
        self.notchanged = False
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

        comparewindow = filemenu.Append(wx.ID_ANY, "Compare Mass Defect Files",
                                        "Open Window to Compare Mass Defect 2D files")
        self.Bind(wx.EVT_MENU, self.on_compare_window, comparewindow)
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

        self.menuaddmultiline = self.plotmenu.Append(wx.ID_ANY, "Add Multiple Lines",
                                                     "Add Multiple Horizontal Lines at Specific Y Values")

        self.menufit = self.plotmenu.Append(wx.ID_ANY, "Fit Peaks",
                                            "Fit total mass defect peaks")
        self.menupeaks = self.plotmenu.Append(wx.ID_ANY, "Label Peaks",
                                              "Label peaks")
        self.menulinreg = self.plotmenu.Append(wx.ID_ANY, "Linear Regression",
                                               "Linear Regression")
        self.menucom = self.plotmenu.Append(wx.ID_ANY, "Center of Mass", "Center of Mass")

        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)
        self.Bind(wx.EVT_MENU, self.on_add_multilines, self.menuaddmultiline)
        self.Bind(wx.EVT_MENU, self.on_fit, self.menufit)
        self.Bind(wx.EVT_MENU, self.on_label_peaks, self.menupeaks)
        self.Bind(wx.EVT_MENU, self.on_linear_regression, self.menulinreg)
        self.Bind(wx.EVT_MENU, self.on_com, self.menucom)

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
        controlsizer.Add(wx.StaticText(panel, label="Reference Mass"), 0, wx.ALIGN_CENTER_VERTICAL)
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

        self.radiobox3 = wx.RadioBox(panel, choices=["Normalized", "Mass (Da)"], label="Mass Defect Units")
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
        oldparams = deepcopy(self.config.defectparams)
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
            self.xlab = "Mass"
            self.ylab = "Normalized Mass Defect"
        else:
            self.factor = self.m0
            self.xlab = "Mass"
            self.ylab = "Mass Defect (Da)"
        self.config.defectparams = [self.nbins, self.transformmode, self.centermode, self.xtype]
        self.notchanged=np.all(oldparams == self.config.defectparams)

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
                                              self.ylab,
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
            self.plot6.contourplot(data2, self.config, xlab=self.ylab, ylab="Individual Spectra", title="",
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
        self.data1d, self.data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(self.datalist[self.pos], self.m0,
                                                                               centermode=self.centermode,
                                                                               nbins=self.nbins,
                                                                               transformmode=self.transformmode,
                                                                               xaxistype=self.xtype)
        if self.xtype == 0:
            factor = 1.0
        else:
            factor = self.m0

        if self.yvals is not None:
            title = str(self.yvals[self.pos])
            spacer = "_"
        else:
            title = str(self.pos)#self.outfname
            spacer = ""
        try:
            save_path2d = os.path.join(self.directory, title + spacer + "2D_Mass_Defects.txt")
            np.savetxt(save_path2d, self.data2d)

            if save_path2d != self.config.defectcomparefiles[0]:
                self.config.defectcomparefiles[1] = self.config.defectcomparefiles[0]
                self.config.defectcomparefiles[0]= save_path2d

            save_path1d = os.path.join(self.directory, title + spacer + "1D_Mass_Defects.txt")
            np.savetxt(save_path1d, self.data1d)
            print('Saved: ', save_path2d, save_path1d)
        except Exception as e:
            print("Failed save", e)

        try:
            self.plot2.contourplot(self.data2d, self.config, xlab=self.xlab, ylab=self.ylab, title=title, normflag=1,
                                   test_kda=True)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)

        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], title, self.ylab,
                                      "Total Intensity", "", self.config)
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot 3", e)
        try:
            if self.m0 == 0:
                return
            self.plot1.colorplotMD(self.datalist[self.pos, :, 0], self.datalist[self.pos, :, 1],
                                   self.datalist[self.pos, :, 0] / float(self.m0) % 1.0 * factor, max=factor,
                                   title="Zero-Charge Mass Spectrum",
                                   xlabel="Mass", ylabel="Intensity", test_kda=True)
            self.plotdat = self.datalist[self.pos]
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0],
                                   title="Total Projection Color", xlabel=self.ylab, cmap="hsv", max=factor,
                                   ylabel="Total Intensity", config=self.config)
        except Exception as e:
            self.plot4.clear_plot()
            print("Failed Plot 4", e)

        self.total=False

    def makeplottotal(self, e=None):
        """
        Runs the kendrick analysis on all of the mass distributions in self.datalist.
        Assumes they are all the same dimensions, which is why we run mergedata when it is loaded.
        Sums the results and plots the sums in 1D and 2D.
        :param e: Unused event
        :return: None
        """
        self.getfromgui()

        if not self.notchanged or not self.total:
            # Run on all
            igrids = []
            m1grid, m2grid = None, None
            for i, dat in enumerate(self.datalist):
                self.data1d, self.data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(dat, self.m0,
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
            self.data2d = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(sumgrid) / np.amax(sumgrid)])
            self.data1d = np.transpose([np.unique(m2grid), np.sum(sumgrid, axis=0)])
            # Save Results
            save_path2d = os.path.join(self.directory, self.outfname + "Total_2D_Mass_Defects.txt")
            np.savetxt(save_path2d, self.data2d)

            if save_path2d != self.config.defectcomparefiles[0]:
                self.config.defectcomparefiles[1] = self.config.defectcomparefiles[0]
                self.config.defectcomparefiles[0]= save_path2d

            save_path1d = os.path.join(self.directory, self.outfname + "Total_1D_Mass_Defects.txt")
            np.savetxt(save_path1d, self.data1d)
            print('Saved: ', save_path2d, save_path1d)
            self.total=True
        else:
            pass

        # Plots

        if self.xtype == 0:
            factor = 1.0
        else:
            factor = self.m0
        try:
            self.plot2.contourplot(self.data2d, self.config, xlab=self.xlab, ylab=self.ylab, title="Total", normflag=1,
                                   test_kda=True)
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot2", e)
        try:
            self.plot1.colorplotMD(self.datasum[:, 0], self.datasum[:, 1],
                                   self.datasum[:, 0] / float(self.m0) % 1.0 * factor,
                                   title="Zero-Charge Mass Spectrum", max=factor,
                                   xlabel="Mass", ylabel="Intensity", test_kda=True)
            self.plotdat = self.datasum
        except Exception as e:
            self.plot1.clear_plot()
            print("Failed Plot1", e)
        try:
            self.plot3.plotrefreshtop(self.data1d[:, 0], self.data1d[:, 1], "Total Projection Black", self.ylab,
                                      "Total Intensity", "", self.config)
        except Exception as e:
            self.plot3.clear_plot()
            print("Failed Plot 3", e)

        try:
            self.plot4.colorplotMD(self.data1d[:, 0], self.data1d[:, 1], self.data1d[:, 0],
                                   title="Total Projection Color", xlabel=self.ylab, cmap="hsv", max=factor,
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
                                                              config=self.config, xtype=self.xtype)
    def on_compare_window(self, e=None):
        frame= MassDefectCompareWindow(self, None, config=self.config)

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
        dialog.initialize_interface(title="Add Line", message="Add Line at Defect Value or Mass:")
        dialog.ShowModal()

        try:
            vval = float(dialog.value)

            if (vval > 1 and self.xtype == 0) or (vval > self.m0 and self.xtype == 0):
                vval = ud.simple_mass_defect(vval, refmass=self.m0, centermode=self.centermode, normtype=self.xtype)

            xlim = self.plot2.subplot1.get_xlim()
            self.plot2.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot2.tickcolor)
            self.plot2.repaint()

            ylim4 = self.plot3.subplot1.get_ylim()
            self.plot3.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot3.tickcolor)
            self.plot3.repaint()

            ylim4 = self.plot4.subplot1.get_ylim()
            self.plot4.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot4.tickcolor)
            self.plot4.repaint()

            if len(self.datalist) > 1:
                try:
                    ylim = self.plot6.subplot1.get_ylim()
                    self.plot6.subplot1.plot((vval, vval), (ylim[0], ylim[1]), color=self.plot2.tickcolor)
                    self.plot6.repaint()

                    ylim4 = self.plot5.subplot1.get_ylim()
                    self.plot5.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot5.tickcolor)
                    self.plot5.repaint()
                except Exception as e:
                    print("Failed 1: ", dialog.value, e)
                    pass

        except Exception as e:
            print("Failed 2: ", dialog.value, e)
            pass

    def on_add_multilines(self, e):
        """
        Add a horizontal line to the plot to visualize a predicted mass defect.
        Opens a dialog to input the value. Can be called more than once to add multiple lines.
        :param e: Unused event
        :return: None
        """
        dialog = miscwindows.MultiInputDialog(self)
        messages = ["Base Mass (Da):", "Oligomer Mass:", "Min #:", "Max #:", "# Step:"]
        defaults = self.defaultvalues
        dialog.initialize_interface(title="Add Lines For Oligomers", messages=messages, defaultvalues=defaults)
        dialog.ShowModal()

        if self.xtype == 0:
            factor = 1.0
        else:
            factor = self.m0

        try:
            values = dialog.values
            self.defaultvalues = values
            basemass = float(values[0])
            omass = float(values[1])
            minnum = int(values[2])
            maxnum = int(values[3])
            step = int(values[4])

            nums = np.arange(minnum, maxnum + 1, step=step)
            masses = basemass + omass * nums

        except:
            print("No values returned")
            return

        vvals = []
        shifts = []
        for i, m in enumerate(masses):
            if True:
                vval = float(m)

                label = str(nums[i])

                if (vval > 1 and self.xtype == 0) or (vval > self.m0 and self.xtype == 0):
                    vval = ud.simple_mass_defect(vval, refmass=self.m0, centermode=self.centermode, normtype=self.xtype)

                shift = 0.98
                boo1=(np.abs(np.array(vvals) - vval)) < factor * 0.05
                numclose = np.sum(boo1)
                if numclose > 0:
                    sarray = np.array(shifts)
                    pastshift = np.sum(sarray[boo1])
                    shift = shift - 0.04 - pastshift
                    shifts.append(0.04 + pastshift)
                else:
                    shifts.append(0)

                vvals.append(vval)

                xlim = self.plot2.subplot1.get_xlim()
                # self.plot2.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot2.tickcolor)
                self.plot2.addtext(label, xlim[1] * shift * self.plot2.kdnorm, vval, color=self.plot2.tickcolor,
                                   verticalalignment="center", vlines=False, hlines=True, nopaint=True)

                ylim4 = self.plot3.subplot1.get_ylim()
                # self.plot3.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot3.tickcolor)
                self.plot3.addtext(label, vval, ylim4[1] * shift, color=self.plot3.tickcolor, nopaint=True)

                ylim4 = self.plot4.subplot1.get_ylim()
                # self.plot4.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot4.tickcolor)
                self.plot4.addtext(label, vval, ylim4[1] * shift, color=self.plot4.tickcolor, nopaint=True)

                if len(self.datalist) > 1:
                    try:
                        ylim = self.plot6.subplot1.get_ylim()
                        # self.plot6.subplot1.plot((vval, vval), (ylim[0], ylim[1]), color=self.plot2.tickcolor)
                        self.plot6.addtext(label, vval, ylim[1] * shift, color=self.plot6.tickcolor, vlines=True,
                                           nopaint=True)

                        ylim4 = self.plot5.subplot1.get_ylim()
                        # self.plot5.subplot1.plot((vval, vval), (ylim4[0], ylim4[1]), color=self.plot5.tickcolor)
                        self.plot5.addtext(label, vval, ylim4[1] * shift, color=self.plot5.tickcolor, nopaint=True)

                    except Exception as e:
                        print("Failed Multiline 1: ", m, e)
                        pass

            #except Exception as e:
            #    print("Failed Multiline 2: ", m, e)
            #    pass

        self.plot2.repaint()
        self.plot3.repaint()
        self.plot4.repaint()
        if len(self.datalist) > 1:
            self.plot6.repaint()
            self.plot5.repaint()

    def on_fit(self, e):
        minval = np.amin(self.data1d[:, 1])
        self.data1d[:, 1] -= minval
        peaks = ud.peakdetect(self.data1d, window=3)
        print("Peaks:", peaks[:, 0])
        peaks = np.concatenate((peaks, [[0, np.amin(self.data1d[:, 1])]]))
        fitdat, fits = MassFitter(self.data1d, peaks, 3, "microguess").perform_fit()
        print("Fits:", fits)

        self.plot3.plotadd(self.data1d[:, 0], fitdat, "green", nopaint=False)

        for f in fits[:, 0]:
            if np.amin(self.data1d[:, 0]) <= f <= np.amax(self.data1d[:, 0]):
                y = np.amax(self.data1d[ud.nearest(self.data1d[:, 0], f), 1])
                self.plot3.addtext(str(np.round(f, 3)), f + 0.075, y * 0.95, vlines=False)
                self.plot3.addtext("", f, y, vlines=True)

    def on_label_peaks(self, e=None):
        self.peaks = ud.peakdetect(self.data1d, window=3)
        print("Peaks:", self.peaks[:, 0])

        for i, p in enumerate(self.peaks):
            y = p[1]
            label = str(np.round(p[0], 3))
            self.plot3.addtext(label, p[0] + 0.075, y * 0.95, vlines=False)
            self.plot3.addtext("", p[0], y, vlines=True)

    def on_linear_regression(self, e=None):
        print("Starting Linear Regression")
        dims = self.data2d.shape

        x = np.unique(self.data2d[:, 0])
        y = np.unique(self.data2d[:, 1])
        xl = len(x)
        yl = len(y)
        z = self.data2d[:, 2].reshape((xl, yl))

        xlin = np.arange(0, xl)
        ylin = np.array([np.argmax(d) / yl * self.m0 for d in z]) + xlin * self.m0
        zlin = np.array([np.max(d) for d in z])

        boo1 = zlin > self.config.peakthresh * np.amax(zlin)

        xlin = xlin[boo1]
        ylin = ylin[boo1]
        zlin = zlin[boo1]

        fit = np.polyfit(xlin, ylin, 1, w=zlin)
        slope = fit[0]
        intercept = fit[1]
        print(fit)

        # mnum=np.floor(self.plotdat[:,0]/self.m0)
        # masses = self.plotdat[:,0]
        # fit2 = np.polyfit(mnum, masses, 1, w=self.plotdat[:,1])
        # print(fit2)

        self.plot1.colorplotMD(xlin, ylin, zlin, max=np.amax(z), cmap="binary",
                               title="Linear Regression Plot", clabel="Intensity",
                               xlabel="Mass Number", ylabel="Mass", test_kda=False)

        # self.plot1.plotadd(xlin, xlin*fit[0]+fit[1], colval="r")
        sval = "Slope: " + str(np.round(slope, 3)) + "\nIntercept: " + str(np.round(intercept, 3)) \
               + "\nInt./RefMass: " + str(np.round(intercept / self.m0, 3))
        self.plot1.addtext(sval, xl * 0.7, np.amax(ylin) * 0.2, vlines=False)

    def on_com(self, e=None):
        print("Starting COM calculations")
        dims = self.data2d.shape

        x = np.unique(self.data2d[:, 0])
        y = np.unique(self.data2d[:, 1])
        xl = len(x)
        yl = len(y)
        z = self.data2d[:, 2].reshape((xl, yl))
        zlin = np.array([np.max(d) for d in z])

        self.on_label_peaks(e=0)

        coms = []
        for d in z.transpose():
            data = np.transpose([x, d])
            b1 = d > np.amax(d) * 0.1
            com, std = ud.center_of_mass(data[b1])
            coms.append(com)

        for p in self.peaks:
            index = ud.nearest(y, p[0])
            c = coms[index]
            carray = [c]
            try:
                carray.append(coms[index - 1])
            except:
                pass
            try:
                carray.append(coms[index + 1])
            except:
                pass

            com = np.average(carray)
            print(p[0], com)

        self.plot1.plotrefreshtop(y, coms, xlabel="Mass Defect", ylabel="Center of Mass (Da)")

        # print(zlin)
        # self.plot1.colorplotMD(y, coms, zlin, max=0, cmap="binary",
        #                       title="Linear Regression Plot", clabel="Intensity",
        #                       xlabel="Mass Number", ylabel="Mass", test_kda=False)


# Main App Execution
if __name__ == "__main__":
    dir = "C:\\Python\\UniDec3\\TestSpectra\\60_unidecfiles"
    file = "60_mass.txt"
    # file = "60_input.dat"

    path = os.path.join(dir, file)

    data = np.loadtxt(path)

    dir = "C:\Python\\UniDec3\TestSpectra\\180_unidecfiles"
    file = "180_mass.txt"

    path = os.path.join(dir, file)

    data2 = np.loadtxt(path)
    datalist = [data, data2]
    # datalist = [data]

    app = wx.App(False)
    frame = MassDefectWindow(None, datalist)
    # frame.on_com(0)
    app.MainLoop()

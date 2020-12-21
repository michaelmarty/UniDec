import matplotlib.pyplot as plt
import os
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows, MassDefectExtractor, ColorPlot
import unidec_modules.unidectools as ud


class MassDefectCompareWindow(wx.Frame):
    def __init__(self, parent, datalist=None, filelist=None, config=None, ):
        """
        """
        wx.Frame.__init__(self, parent, title="Mass Defect Comparison")  # ,size=(-1,-1))
        file1 = None
        file2 = None
        if datalist is None:
            if filelist is None:
                if config is None:
                    print("ERROR: No data supplied to Mass Defect Compare Window")
                    return
                else:
                    self.config = config
                    file1 = self.config.defectcomparefiles[0]
                    file2 = self.config.defectcomparefiles[1]
            else:
                file1 = filelist[0]
                file2 = filelist[1]

            self.get_data(file1, file2)
        self.file1 = file1
        self.file2 = file2

        self.parent = parent

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
        self.total = False
        self.notchanged = False



        '''
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
        self.SetMenuBar(menu_bar)'''

        # Setup the GUI
        panel = wx.Panel(self)

        self.plot1 = plot2d.Plot2d(panel)
        self.plot2 = plot2d.Plot2d(panel)
        self.plot3 = ColorPlot.ColorPlot2D(panel)
        # self.plot4 = plot1d.Plot1d(panel)

        sizer = wx.BoxSizer(wx.VERTICAL)
        # plotsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        plotsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        # plotsizer1.Add(self.plot4, 0, wx.EXPAND)
        plotsizer2.Add(self.plot1, 2, wx.EXPAND)
        plotsizer2.Add(self.plot2, 2, wx.EXPAND)
        plotsizer2.Add(self.plot3, 2, wx.EXPAND)

        # sizer.Add(plotsizer1, 1, wx.EXPAND)
        sizer.Add(plotsizer2, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlf1 = wx.TextCtrl(panel, value=str(self.file1), size=(1500, 25))
        self.f1button = wx.Button(panel, label="File 1: ")
        self.Bind(wx.EVT_BUTTON, self.on_f1button, self.f1button)
        controlsizer.Add(self.f1button, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlf1, 0, wx.ALIGN_CENTER_VERTICAL)

        controlsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlf2 = wx.TextCtrl(panel, value=str(self.file2), size=(1500, 25))
        self.f2button = wx.Button(panel, label="File 2: ")
        self.Bind(wx.EVT_BUTTON, self.on_f2button, self.f2button)
        controlsizer1.Add(self.f2button, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer1.Add(self.ctlf2, 0, wx.ALIGN_CENTER_VERTICAL)

        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        label = "Replot"
        ReplotButton = wx.Button(panel, label=label)
        controlsizer2.Add(ReplotButton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_replot, ReplotButton)

        '''
        self.radiobox = wx.RadioBox(panel, choices=["Integrate", "Interpolate"], label="Type of Transform")
        controlsizer.Add(self.radiobox, 0, wx.EXPAND)
        self.radiobox.SetSelection(self.transformmode)

        self.radiobox2 = wx.RadioBox(panel, choices=["-0.5:0.5", "0:1"], label="Range")
        controlsizer.Add(self.radiobox2, 0, wx.EXPAND)
        self.radiobox2.SetSelection(self.centermode)

        self.radiobox3 = wx.RadioBox(panel, choices=["Normalized", "Mass (Da)"], label="Mass Defect Units")
        controlsizer.Add(self.radiobox3, 0, wx.EXPAND)
        self.radiobox3.SetSelection(self.xtype)
        '''
        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer1, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)


        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.makeplot(0)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()

    def get_data(self, file1, file2):
        try:
            data1 = np.loadtxt(file1)
        except Exception as e:
            print("Error with File 1:", file1, e)
            data1 = None
        try:
            data2 = np.loadtxt(file2)
        except Exception as e:
            print("Error with File 2:", file2, e)
            data2 = None
        datalist = [data1, data2]
        self.datalist = datalist

    def on_replot(self, e=0):
        self.file1 = self.ctlf1.GetValue()
        self.file2 = self.ctlf2.GetValue()
        self.get_data(self.file1, self.file2)
        self.makeplot(e)

    def on_f1button(self, e=0):
        newfile = self.file_dialog()
        if newfile is not None:
            self.ctlf1.SetValue(newfile)
        self.on_replot()

    def on_f2button(self, e=0):
        newfile = self.file_dialog()
        if newfile is not None:
            self.ctlf2.SetValue(newfile)
        self.on_replot()

    def file_dialog(self, e=0):
        filename=None
        dlg = wx.FileDialog(self, "Choose a data file", '', "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            print("Opening: ", filename)
        dlg.Destroy()
        return filename

    def makeplot(self, e):
        self.plot1.contourplot(self.datalist[0], self.config, xlab=self.xlab, ylab=self.ylab, title="", normflag=1,
                               test_kda=True)
        if self.datalist[1] is not None:
            self.plot2.contourplot(self.datalist[1], self.config, xlab=self.xlab, ylab=self.ylab, title="", normflag=1,
                                   test_kda=True)
        if self.datalist[1] is not None:
            self.plot3.make_compare_color_plot(self.datalist[0], self.datalist[1], xlab=self.xlab, ylab=self.ylab,)

    def on_close(self, e):
        """
        Close the window. Will try to set self.config.molig as self.m0.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        # self.MakeModal(False)

# Main App Execution
if __name__ == "__main__":
    file1 = "Z:\Group Share\Julia\M2 protein\wildtype\\201116\\20201116_jat_am2_wt_dmpc_d1t0_frac1c_3to1_pc_40uMdrug_201116124724_unidecfiles\\20201116_jat_am2_wt_dmpc_d1t0_frac1c_3to1_pc_40uMdrug_201116124724_Total_2D_Mass_Defects.txt"
    file2 = "Z:\Group Share\Julia\M2 protein\wildtype\\201116\\20201116_jat_am2_wt_dmpc_frac1c_3to1_pc_unidecfiles\\20201116_jat_am2_wt_dmpc_frac1c_3to1_pc_Total_2D_Mass_Defects.txt"

    filelist=[file1,file2]

    app = wx.App(False)
    frame = MassDefectCompareWindow(None, filelist=filelist)
    # frame.on_com(0)
    app.MainLoop()

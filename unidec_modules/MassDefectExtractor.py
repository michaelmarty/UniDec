import os
import numpy as np
import wx
from unidec_modules import unidecstructure, plot1d, plot2d
import unidec_modules.unidectools as ud
import unidec_modules.masstools as masstools
import matplotlib.cm as cm
from unidec_modules.isolated_packages import FileDialogs
from matplotlib.ticker import FixedLocator

__author__ = 'Michael.Marty'

extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position",
                  5: "Center of Mass 50%", 6: "Center of Mass 10%"}
extractlabels = {0: "Intensity", 1: "Intensity", 2: "Area", 3: "Mass", 4: "Mass", 5: "Mass", 6: "Mass"}
markers = ['o', 'v', '^', '>', 's', 'd', '*']


class MassDefectExtractorWindow(wx.Frame):
    def __init__(self, parent, datalist, xarray, yarray, config=None):
        wx.Frame.__init__(self, parent, title="Mass Defect")  # ,size=(-1,-1))

        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.discreteplot = 0
            self.config.cmap = u"jet"
            self.config.peakcmap = u"rainbow"
            self.config.separation = 0.025
        else:
            self.config = config

        # Setup initial values
        try:
            self.directory = parent.directory
        except:
            self.directory = os.getcwd()

        self.outfname = os.path.splitext(self.config.filename)[0]
        if self.outfname is not "":
            self.outfname += "_"

        self.window = 0.05
        self.data = datalist
        self.xdat = xarray
        self.ydat = yarray
        defaultexchoice = "Local Max"

        # Make the menu
        filemenu = wx.Menu()
        menu_oligomer_tools = filemenu.Append(wx.ID_ANY, "Generate from Masses", "Oligomer and Mass Tools")
        filemenu.AppendSeparator()
        menu_save_fig_png = filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                            "Save all figures as PNG in central directory")
        menu_save_fig_pdf = filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                            "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_oligomer_tools, menu_oligomer_tools)
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig_png)
        self.Bind(wx.EVT_MENU, self.on_save_fig_pdf, menu_save_fig_pdf)

        self.plotmenu = wx.Menu()
        self.menuaddline = self.plotmenu.Append(wx.ID_ANY, "Add Horizontal Line",
                                                "Add Horizontal Line at Specific Y Value")
        # self.menufit = self.plotmenu.Append(wx.ID_ANY, "Fit Peaks",
        #                                    "Fit total mass defect peaks")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)
        # self.Bind(wx.EVT_MENU, self.on_fit, self.menufit)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the GUI
        panel = wx.Panel(self)

        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot2d.Plot2d(panel)
        self.plot5 = plot1d.Plot1d(panel)
        self.plot6 = plot2d.Plot2d(panel)

        sizer = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        sb = wx.StaticBox(panel, label='Set the Mass Defect Values to Extract')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        importbutton = wx.Button(panel, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.on_import_masses, importbutton)

        clearbutt = wx.Button(panel, label="Clear List")
        self.Bind(wx.EVT_BUTTON, self.on_clear_masslist, clearbutt)

        addbutton = wx.Button(panel, label="Add Species")
        self.Bind(wx.EVT_BUTTON, self.on_add_mass, addbutton)

        sbs.Add(importbutton, 0, wx.EXPAND)
        sbs.Add(addbutton, 0, wx.EXPAND)
        self.masslistbox = masstools.MassListCtrlPanel(panel, columntitle="Mass Defect Value", size=(210, 320))
        sbs.Add(wx.StaticText(panel, label="Mass Defect List"))
        sbs.Add(self.masslistbox)
        sbs.Add(clearbutt, 0, wx.EXPAND)
        # hbox.Add(sbs, 0, wx.EXPAND)

        plotsizer1 = wx.BoxSizer(wx.VERTICAL)
        plotsizer2 = wx.BoxSizer(wx.VERTICAL)
        plotsizer1.Add(self.plot1, 2, wx.EXPAND)
        plotsizer1.Add(self.plot2, 2, wx.EXPAND)

        plotsizer2.Add(self.plot5, 0, wx.EXPAND)
        plotsizer2.Add(self.plot6, 0, wx.EXPAND)

        hbox.Add(sbs, 0)
        hbox.Add(plotsizer1, 1, wx.EXPAND)
        hbox.Add(plotsizer2, 1, wx.EXPAND)

        sizer.Add(hbox, 1, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)

        totalbutton = wx.Button(panel, label="Extract")
        controlsizer.Add(totalbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_extract, totalbutton)

        controlsizer.Add(wx.StaticText(panel, label=" How to extract: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctlextract = wx.ComboBox(panel, value=defaultexchoice, choices=list(extractchoices.values()),
                                      style=wx.CB_READONLY | wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlextract, 0, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)

        self.ctlwindow = wx.TextCtrl(panel, value=str(self.window))
        controlsizer.Add(wx.StaticText(panel, label="Window:"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctlwindow, 0, wx.ALIGN_CENTER_VERTICAL)

        self.ctlnorm = wx.RadioBox(panel, label="Extract Normalization",
                                   choices=["None", "Max", "Sum", "Peak Max", "Peak Sum"])
        self.ctlnorm.SetSelection(2)
        # , majorDimension=3,style=wx.RA_SPECIFY_COLS)
        self.ctlnorm.SetSelection(2)
        controlsizer.Add(self.ctlnorm, 0, wx.ALIGN_CENTER_VERTICAL)

        sizer.Add(controlsizer, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()
        self.make_list_plots()

        # defaultmdlist = [0.0, 0.25, 0.55]
        # self.masslistbox.list.populate(defaultmdlist)
        # self.on_extract()

    def on_extract(self, e=None):
        print("Extracting...")
        self.extractchoice = self.ctlextract.GetSelection()
        self.norm = self.ctlnorm.GetSelection()
        self.mdlist = self.masslistbox.list.get_list()
        self.window = float(self.ctlwindow.GetValue())
        print(self.extractchoice, self.window, self.norm, self.mdlist)

        grid = []
        for i, x in enumerate(self.mdlist):
            ext = []
            for j, dat in enumerate(self.data):
                val = ud.data_extract(np.transpose([self.xdat, dat]), x, self.extractchoice, self.window,
                                      zero_edge=False)
                ext.append(val)
            grid.append(ext)
        self.grid = np.array(grid)
        ud.normalize_extracts(self.grid, self.norm)
        try:
            save_path2d = os.path.join(self.directory, self.outfname + "Mass_Defect_Extracts.txt")
            np.savetxt(save_path2d, self.grid)
            np.savetxt(os.path.join(self.directory, self.outfname + "Mass_Defect_Extracts_xvals.txt"), self.ydat)
            print('Saved: ', save_path2d)
        except Exception as e:
            print("Failed Data Export Extracts", e)
        self.make_ext_plots()

    def make_ext_plots(self):
        self.colormap = cm.get_cmap(ud.smartdecode(self.config.peakcmap), len(self.mdlist))
        if self.colormap is None:
            self.colormap = cm.get_cmap(u"rainbow", len(self.mdlist))
        self.peakcolors = self.colormap(np.arange(len(self.mdlist)))

        for i, x in enumerate(self.mdlist):
            print(self.grid[i])
            try:
                if i == 0:
                    self.plot1.plotrefreshtop(self.ydat, self.grid[i], title="Extracted Data",
                                              xlabel="Variable 1",
                                              ylabel=extractlabels[self.extractchoice],
                                              marker=markers[i % len(markers)],
                                              label=str(x), color=self.peakcolors[i], config=self.config)
                else:
                    self.plot1.plotadd(self.ydat, self.grid[i], newlabel=str(x), marker=markers[i % len(markers)],
                                       colval=self.peakcolors[i])
            except Exception as e:
                self.plot1.clear_plot()
                print("Failed Plot Ext1", e)
        self.plot1.add_legend()

        try:
            self.make_2d_plot()
        except Exception as e:
            print(e)

    def make_2d_plot(self):
        iarray = np.arange(0, len(self.mdlist))
        m1grid, m2grid = np.meshgrid(self.ydat, iarray, indexing='ij')
        data2 = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(self.grid.transpose())])
        try:
            self.plot2.contourplot(data2, self.config, xlab="Variable 1", ylab="Extracts", title="",
                                   normflag=1, discrete=True)
            self.plot2.subplot1.yaxis.set_major_locator(FixedLocator(iarray))
            self.plot2.subplot1.set_yticklabels(self.mdlist)  # , rotation=90)
            self.plot2.repaint()
        except Exception as e:
            self.plot2.clear_plot()
            print("Failed Plot 2", e)
        pass

    def make_list_plots(self):
        self.colormap = cm.get_cmap(ud.smartdecode(self.config.peakcmap), len(self.ydat))
        if self.colormap is None:
            self.colormap = cm.get_cmap(u"rainbow", len(self.ydat))
        self.peakcolors = self.colormap(np.arange(len(self.ydat)))

        for i, dat in enumerate(self.data):
            dat2 = dat
            try:
                if i == 0:
                    self.plot5.plotrefreshtop(self.xdat, dat2 - self.config.separation * i, "All Data",
                                              "Mass Defect",
                                              "Total Intensity", "", color=self.peakcolors[i], config=self.config)
                else:
                    self.plot5.plotadd(self.xdat, dat2 - self.config.separation * i, colval=self.peakcolors[i])
            except Exception as e:
                self.plot5.clear_plot()
                print("Failed Plot Ext5", e)
        self.plot5.repaint()

        m1grid, m2grid = np.meshgrid(self.xdat, self.ydat, indexing='ij')
        data2 = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(self.data.transpose())])
        try:
            self.plot6.contourplot(data2, self.config, xlab="Mass Defect", ylab="Individual Spectra", title="",
                                   normflag=1)
        except Exception as e:
            self.plot6.clear_plot()
            print("Failed Plot Ext6", e)
        pass

    def on_add_line(self, e):
        """
        Add a horizontal line to the plot to visualize a predicted mass defect.
        Opens a dialog to input the value. Can be called more than once to add multiple lines.
        :param e: Unused event
        :return: None
        """

        for x in self.mdlist:
            #print(x)
            try:
                ylim = self.plot6.subplot1.get_ylim()
                self.plot6.subplot1.plot((x, x), (ylim[0], ylim[1]), color=self.plot2.tickcolor)
                self.plot6.repaint()

                ylim4 = self.plot5.subplot1.get_ylim()
                self.plot5.subplot1.plot((x, x), (ylim4[0], ylim4[1]), color=self.plot5.tickcolor)
                self.plot5.repaint()
            except Exception as e:
                print("Failed: ", x, e)
                pass

    def on_oligomer_tools(self, e):
        dlg = masstools.MassSelection(self)
        dlg.init_dialog(self.config, None, massdat=None)
        result = dlg.ShowModal()
        if self.config.masslist != [] and result == 0 and self.config.kendrickmass is not None:
            print(self.config.masslist)
            kmass = self.config.masslist / self.config.kendrickmass
            if np.amin(self.xdat) < 0:
                centermode = 0
            else:
                centermode = 1

            if centermode == 1:
                nominalkmass = np.floor(kmass)
            else:
                nominalkmass = np.round(kmass)
            kmdefect = kmass - nominalkmass
            print(kmdefect)
            self.masslistbox.list.populate(kmdefect)

    def on_import_masses(self, e):
        """
        Opens a dialog to import mass list files.
        :param e: Unused event
        :return: None
        """
        mfilename = FileDialogs.open_file_dialog("Open Text File with List of Mass Defects", file_types="*.*")
        if mfilename is not None:
            importmass = np.loadtxt(mfilename)
            self.masslistbox.list.populate(importmass)

    def on_clear_masslist(self, e):
        """
        Clears the mass list.
        :param e: Unused event
        :return: None
        """
        self.masslistbox.list.clear()

    def on_add_mass(self, e):
        """
        Adds a blank line to the mass list.
        :param e: Unused Event
        :return: None
        """
        self.masslistbox.list.add_line()

    def on_close(self, e):
        """
        Close the window.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        # self.MakeModal(False)

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "MassDefectExtract1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1

        name1 = os.path.join(self.directory, self.outfname + "MassDefectExtract2.png")
        if self.plot5.flag:
            self.plot5.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectExtract3.png")
        if self.plot6.flag:
            self.plot6.on_save_fig(e, name2)
            # print name2
        name2 = os.path.join(self.directory, self.outfname + "MassDefectExtract4.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "MassDefectExtract1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1

        name1 = os.path.join(self.directory, self.outfname + "MassDefectExtract2.pdf")
        if self.plot5.flag:
            self.plot5.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, self.outfname + "MassDefectExtract3.pdf")
        if self.plot6.flag:
            self.plot6.on_save_fig(e, name2)
            # print name2
        name2 = os.path.join(self.directory, self.outfname + "MassDefectExtract4.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2


# Main App Execution
if __name__ == "__main__":
    fname = "C:\Python\\UniDec3\\unidec_modules\Mass_Defect_Grid.txt"
    data = np.loadtxt(fname)
    # print(data)
    xarray = np.arange(0, len(data[0]))
    xarray = xarray / np.amax(xarray)
    # print(xarray)
    yarray = np.arange(0, len(data))

    app = wx.App(False)
    frame = MassDefectExtractorWindow(None, data, xarray, yarray)
    frame.mdlist = [0.0, 0.1, 0.2, 0.3]
    frame.masslistbox.list.populate(frame.mdlist)
    frame.on_extract()
    frame.make_2d_plot()
    app.MainLoop()

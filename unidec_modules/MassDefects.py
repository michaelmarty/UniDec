import os
from copy import deepcopy

import numpy as np
import wx
from scipy.interpolate import interp1d

from unidec_modules import unidecstructure, plot1d, plot2d, miscwindows
import unidec_modules.unidectools as ud

__author__ = 'Michael.Marty'


class MassDefectWindow(wx.Frame):
    def __init__(self, parent, datalist, config=None, yvals=None, pks=None, value=None, dir=None, show=True, *args,
                 **kwargs):
        wx.Frame.__init__(self, parent, title="Mass Defect")  # ,size=(-1,-1))
        if dir is None:
            self.directory = os.getcwd()
        else:
            self.directory = dir

        self.parent = parent
        defaultvalue = deepcopy(value)
        if defaultvalue is None:
            defaultvalue = 760.076

        self.filemenu = wx.Menu()

        self.menuSaveFigPNG = self.filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                                   "Save all figures as PNG in central directory")
        self.menuSaveFigPDF = self.filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                                   "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_save_fig, self.menuSaveFigPNG)
        self.Bind(wx.EVT_MENU, self.on_save_figPDF, self.menuSaveFigPDF)

        self.plotmenu = wx.Menu()
        self.menuaddline = self.plotmenu.Append(wx.ID_ANY, "Add Horizontal Line",
                                                "Add Horizontal Line at Specific Y Value")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)

        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(self.menuBar)

        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
            self.config.discreteplot = 0
            self.config.cmap = "jet"
        else:
            self.config = config

        self.config.publicationmode = 0

        self.datalist = datalist
        self.pos = -1
        self.yvals = yvals

        self.panel = wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = plot1d.Plot1d(self.panel)
        self.plot2 = plot2d.Plot2d(self.panel)
        self.sizer.Add(self.plot1, 1, wx.EXPAND)
        self.sizer.Add(self.plot2, 1, wx.EXPAND)

        self.controlsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.ctlm0 = wx.TextCtrl(self.panel, value=str(defaultvalue))
        self.ctlwindow = wx.TextCtrl(self.panel, value="50")
        self.controlsizer.Add(wx.StaticText(self.panel, label="Kendrick Mass"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlm0, 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(wx.StaticText(self.panel, label="Number of Defect Bins"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.controlsizer.Add(self.ctlwindow, 0, wx.ALIGN_CENTER_VERTICAL)

        self.controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        if len(datalist) > 1:
            label = "Back"
        else:
            label = "Replot"
        self.backbutton = wx.Button(self.panel, label=label)
        self.controlsizer2.Add(self.backbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_back, self.backbutton)
        if len(datalist) > 1:
            self.nextbutton = wx.Button(self.panel, label="Next")
            self.controlsizer2.Add(self.nextbutton, 0, wx.EXPAND)
            self.Bind(wx.EVT_BUTTON, self.on_next, self.nextbutton)
            self.totalbutton = wx.Button(self.panel, label="Total")
            self.controlsizer2.Add(self.totalbutton, 0, wx.EXPAND)
            self.Bind(wx.EVT_BUTTON, self.on_total, self.totalbutton)
        else:
            if not ud.isempty(pks):
                self.pks = pks
                self.peaksbutton = wx.Button(self.panel, label="Plot Peaks")
                self.controlsizer2.Add(self.peaksbutton, 0, wx.EXPAND)
                self.Bind(wx.EVT_BUTTON, self.on_peaks, self.peaksbutton)

        self.radiobox = wx.RadioBox(self.panel, choices=["Integrate", "Interpolate"], label="Type of Transform")
        self.controlsizer.Add(self.radiobox, 0, wx.EXPAND)
        self.radiobox.SetSelection(1)

        self.radiobox2 = wx.RadioBox(self.panel, choices=["-0.5:0.5", "0:1"], label="Range")
        self.controlsizer.Add(self.radiobox2, 0, wx.EXPAND)
        self.radiobox2.SetSelection(1)

        self.radiobox3 = wx.RadioBox(self.panel, choices=["Mass Number", "Mass"], label="X-Axis")
        self.controlsizer.Add(self.radiobox3, 0, wx.EXPAND)
        self.radiobox3.SetSelection(1)

        self.sizer.Add(self.controlsizer, 0, wx.EXPAND)
        self.sizer.Add(self.controlsizer2, 0, wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.ylab = "Normalized Mass Defect"
        try:
            self.on_total(0)
        except:
            self.on_next(0)
        self.Centre()
        # self.MakeModal(True)
        self.Show(show)

    def OnClose(self, e):
        try:
            self.parent.molig = self.m0
        except:
            pass
        self.Destroy()
        # self.MakeModal(False)

    def getfromgui(self):
        try:
            self.m0 = float(self.ctlm0.GetValue())
            try:
                self.window = float(self.ctlwindow.GetValue())
            except:
                self.window = 0
        except:
            print "Failed to get from gui"
        self.ktype = self.radiobox3.GetSelection()
        if self.ktype == 0:
            self.factor = 1
            self.xlab = "Nominal Mass Number"
        else:
            self.factor = self.m0
            self.xlab = "Mass"

    def makegrid(self):
        self.xaxis = self.datalist[0][:, 0]
        self.kmass = self.xaxis * 1. / float(self.m0)
        flag = self.radiobox2.GetSelection()
        if flag == 1:
            self.nominalkmass = np.floor(self.kmass)
        else:
            self.nominalkmass = np.round(self.kmass)
        self.kmdefectexact = self.kmass - self.nominalkmass
        self.defects = np.linspace(np.amin(self.kmdefectexact), np.amax(self.kmdefectexact), self.window, endpoint=True)
        self.nominal = np.unique(self.nominalkmass)
        self.m1grid, self.m2grid = np.meshgrid(self.nominal, self.defects, indexing='ij')

    def extractall(self):
        self.igrid = np.zeros((len(self.datalist), len(self.nominal), len(self.defects)))
        flag = self.radiobox.GetSelection()
        for i, data in enumerate(self.datalist):
            if flag == 1:
                # Interpolation
                f = interp1d(data[:, 0], data[:, 1], bounds_error=False, fill_value=0)
                for j in xrange(0, len(self.nominal)):
                    for k in xrange(0, len(self.defects)):
                        nommass = self.nominal[j]
                        defect = self.defects[k]
                        mass = (nommass + defect) * self.m0
                        intensity = f(mass)
                        self.igrid[i, j, k] = intensity
            else:
                # Integration
                for j in xrange(0, len(self.kmass)):
                    nommass = self.nominalkmass[j]
                    defect = self.kmdefectexact[j]
                    pos = ud.nearest(self.defects, defect)
                    pos2 = ud.nearest(self.nominal, nommass)
                    try:
                        intensity = data[j, 1]
                    except:
                        intensity = 0
                    self.igrid[i, pos2, pos] += intensity

        self.igrid = self.igrid / np.amax(self.igrid)

    def makeplot(self):
        i = self.pos
        maxval = np.amax(self.igrid[i])
        dat = np.transpose(
            [np.ravel(self.m1grid) * self.factor, np.ravel(self.m2grid), np.ravel(self.igrid[i]) / maxval])

        if self.yvals is not None:
            title = str(self.yvals[i])
        else:
            title = ""
        try:
            self.plot2.contourplot(dat, self.config, xlab=self.xlab, ylab=self.ylab, title=title, normflag=1)
        except:
            self.plot2.clear_plot()
            print "Failed Plot2"
        try:
            self.plot1.plotrefreshtop(np.unique(self.m2grid), np.sum(self.igrid[i], axis=0), title, "Mass Defect",
                                      "Total Intensity", "", self.config)
        except:
            self.plot1.clear_plot()
            print "Failed Plot1"

    def makeplottotal(self):
        grid = np.sum(self.igrid, axis=0)
        maxval = np.amax(grid)
        dat = np.transpose([np.ravel(self.m1grid) * self.factor, np.ravel(self.m2grid), np.ravel(grid) / maxval])
        try:
            self.plot2.contourplot(dat, self.config, xlab=self.xlab, ylab=self.ylab, title="Total", normflag=1)
            path = os.path.join(self.directory, "Total_2D_Mass_Defects.txt")
            np.savetxt(path, dat)
            print 'Saved: ', path
        except:
            self.plot2.clear_plot()
            print "Failed Plot2"
        try:
            self.plot1.plotrefreshtop(np.unique(self.m2grid), np.sum(grid, axis=0), "Total Projection", "Mass Defect",
                                      "Total Intensity", "", self.config)

            outputdata = np.transpose([np.unique(self.m2grid), np.sum(grid, axis=0)])
            outfile = os.path.join(self.directory, "Total_1D_Mass_Defects.txt")
            np.savetxt(outfile, outputdata)
        except:
            self.plot1.clear_plot()
            print "Failed Plot1"

    def on_back(self, e):
        self.getfromgui()
        self.pos += -1
        self.pos = self.pos % len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_next(self, e):
        self.getfromgui()
        self.pos += 1
        self.pos = self.pos % len(self.datalist)
        self.makegrid()
        self.extractall()
        self.makeplot()

    def on_total(self, e):
        self.getfromgui()
        self.makegrid()
        self.extractall()
        self.makeplottotal()

    def on_peaks(self, e):

        self.getfromgui()
        flag = self.radiobox2.GetSelection()
        self.pks.get_mass_defects(self.m0, mode=flag)
        self.plot1.clear_plot()
        xvals = []
        yvals = []
        for p in self.pks.peaks:
            x = p.kendricknum * self.factor
            y = p.kendrickdefect
            if not self.plot1.flag:
                self.plot1.plotrefreshtop(x, y, "Mass Peaks", self.xlab, self.ylab, "", self.config, color=p.color,
                                          marker=p.marker)
            else:
                self.plot1.plotadddot(x, y, p.color, p.marker)
            xvals.append(x)
            yvals.append(y)
        # self.plot1.subplot1.set_xlim((np.amin(xvals)-1,np.amax(xvals)+1))
        # self.plot1.subplot1.set_ylim((-0.5,0.5))
        datalims = [np.amin(xvals), np.amin(yvals), np.amax(xvals), np.amax(yvals)]
        self.plot1.setup_zoom([self.plot1.subplot1], "box", data_lims=datalims)
        self.plot1.subplot1.set_ylim(self.plot2.subplot1.get_ylim())
        self.plot1.subplot1.set_xlim(self.plot2.subplot1.get_xlim())
        self.plot1.repaint()

        pdat = np.transpose([xvals, yvals])
        path = os.path.join(self.directory, "Total_2D_Mass_Defects.txt")
        np.savetxt(path, pdat)
        highmassdefects = pdat[pdat[:, 0] > 80000, 1]
        print np.average(highmassdefects), np.std(highmassdefects)

    def on_save_fig(self, e):

        name1 = os.path.join(self.directory, "MassDefectFigure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print name1
        name2 = os.path.join(self.directory, "MassDefectFigure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print name2

    def on_save_figPDF(self, e):

        name1 = os.path.join(self.directory, "MassDefectFigure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            print name1
        name2 = os.path.join(self.directory, "MassDefectFigure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            print name2

    def on_add_line(self, e):
        dialog = miscwindows.SingleInputDialog(self)
        dialog.InitUI(title="Add Line", message="Add Line at Defect Value: ")
        dialog.ShowModal()

        try:
            vval = float(dialog.value)
            xlim = self.plot2.subplot1.get_xlim()
            self.plot2.subplot1.plot((xlim[0], xlim[1]), (vval, vval), color=self.plot2.ticcol)
            self.plot2.repaint()
        except:
            print "Failed: ", dialog.value
            pass


# Main App Execution
if __name__ == "__main__":
    dir = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit100_unidecfiles"
    file = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit100_mass.txt"

    path = os.path.join(dir, file)

    data = np.loadtxt(path)

    dir = "C:\\MassLynx\\Mike.PRO\Data\\150521\\mzML\\Aqpz_05_Ramp3\\MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_unidecfiles"
    file = "MTM_150521_AqpZ_05_POPC_Ramp_1-5pbar_20mit120_mass.txt"

    path = os.path.join(dir, file)

    data2 = np.loadtxt(path)
    datalist = [data, data2]

    app = wx.App(False)
    frame = MassDefectWindow(None, datalist)
    app.MainLoop()

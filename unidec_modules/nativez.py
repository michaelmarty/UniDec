import sys
import time
import wx
import numpy as np
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from wx.lib.agw import ultimatelistctrl as ULC

import unidectools as ud
import plot1d
import plot2d
import MassFitter

__author__ = 'Michael.Marty'


def f0(x, array):
    # s2=s/2.35482
    out = 0
    for i in xrange(0, len(array) / 3):
        m = array[i * 3]
        s2 = array[i * 3 + 1]
        out += np.exp(-(x - m) * (x - m) / (2 * s2 * s2))
    return out


def simchargefit(x2):
    fit = [ud.predict_charge(i) for i in x2]
    return np.array(fit)


def localmax(Y, index, window):
    start = np.amax([0, index - window])
    end = np.amin([len(Y), index + window])
    return np.amax(Y[start:end])


class zoffset:
    def __init__(self, *args):
        self.offset = 0
        self.intensity = 0
        self.index = 0
        self.color = [255, 255, 255, 255]
        self.marker = "."
        self.id = 0
        self.width = 0
        self.nstate = 0
        self.extractwidth = 1

    def Make(self, offset, intensity, index, color, marker):
        self.offset = offset
        self.intensity = intensity
        self.index = index
        self.color = color
        self.marker = marker


class NativeZ(wx.Dialog):
    def __init__(self, *args, **kwargs):
        super(NativeZ, self).__init__(*args, **kwargs)
        self.SetSize((1400, 1000))
        self.SetTitle("Native Charge Tools")

    def InitUI(self, xvals, yvals, zdat, config, pks):
        self.config = config
        self.pks = pks
        self.xlen = len(xvals)
        self.ylen = len(yvals)
        self.xvals = np.array(xvals)
        self.yvals = np.array(yvals)
        self.newgrid = np.reshape(zdat, (self.xlen, self.ylen))

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Set Parameters to Plot Native Z')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        size = (5.5, 3.75)
        self.plot1 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot2 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot3 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot4 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot5 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot6 = plot1d.Plot1d(self.pnl, figsize=size)
        self.plot7 = plot2d.Plot2d(self.pnl, figsize=size)
        self.hbox0.Add(self.plot1, 0)
        self.hbox0.Add(self.plot3, 0)
        self.hbox0.Add(self.plot4, 0)

        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.addbutton = wx.Button(self.pnl, label="Add Line")
        self.fitbutton = wx.Button(self.pnl, label="Fit")
        self.resetbutton = wx.Button(self.pnl, label="Reset to Default")
        self.extractbutton = wx.Button(self.pnl, label="Extract")
        self.massoffset = wx.TextCtrl(self.pnl, value="0", size=(50, -1))
        self.ctlfilt = wx.RadioBox(self.pnl, label="Extract Shape", choices=["Box", "Gaussian"])
        self.savefigbutt = wx.Button(self.pnl, label="Save Figures")
        self.replotbutton = wx.Button(self.pnl, label="Replot")
        self.hbox1.Add(self.addbutton, 0)
        self.hbox1.Add(self.replotbutton, 0)
        self.hbox1.Add(self.resetbutton, 0)
        self.hbox1.Add(wx.StaticText(self.pnl, label="     "), 0, wx.ALIGN_CENTER_VERTICAL)

        self.hbox1.Add(self.fitbutton, 0)
        self.hbox1.Add(self.extractbutton, 0)
        self.hbox1.Add(wx.StaticText(self.pnl, label="     Monomer Mass: "), 0)  # , wx.ALIGN_CENTER_VERTICAL)
        self.hbox1.Add(self.massoffset, 0)
        self.hbox1.Add(self.ctlfilt, 0)

        self.hbox1.Add(self.savefigbutt, 0)

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.plot2, 0)
        self.hbox2.Add(self.plot5, 0)
        self.hbox2.Add(self.plot6, 0)

        self.list = ColorList(self.pnl)

        self.sbs.Add(self.hbox0, 0, wx.EXPAND)
        self.sbs.Add(self.hbox2, 0, wx.EXPAND)
        self.sbs.Add(self.hbox1, 0, wx.EXPAND)

        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.list, 1, wx.EXPAND)
        self.hbox3.Add(self.plot7, 0)

        self.sbs.Add(self.hbox3, 0, wx.EXPAND)

        self.pnl.SetSizer(self.sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        self.list.ultimateList.Bind(wx.EVT_BUTTON, self.on_delete)
        self.list.ultimateList.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.update)
        # self.list.ultimateList.Bind(ULC.EVT_LIST_END_LABEL_EDIT,self.on_edit)
        self.fitbutton.Bind(wx.EVT_BUTTON, self.fit)
        self.addbutton.Bind(wx.EVT_BUTTON, self.onadd)
        self.resetbutton.Bind(wx.EVT_BUTTON, self.OnReset)
        self.extractbutton.Bind(wx.EVT_BUTTON, self.Extract)
        self.savefigbutt.Bind(wx.EVT_BUTTON, self.SaveFig)
        self.replotbutton.Bind(wx.EVT_BUTTON, self.OnReplot)

        self.Center()

        # Set Range Here
        tstart = time.clock()
        self.MakeFArray(-50, 15)
        tend = time.clock()
        print "F Array Time: %.2gs" % (tend - tstart)
        if self.config.zoffs == []:
            self.GetMaxima()
            print self.maxes
        else:
            self.zoffs = self.config.zoffs
            self.PlotZoffs()
        self.PopulateList(0)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

    def OnReplot(self, e):
        self.update(e)
        self.PlotZoffs()
        self.UpdateList()

    def PlotZoffs(self):
        self.plot1.plotrefreshtop(self.ftot[:, 0], self.ftot[:, 1], title="Native Z Offset", xlabel="Offset",
                                  ylabel="Intensity", zoom="span")
        for i in xrange(0, len(self.zoffs)):
            self.plot1.plotadddot(self.zoffs[i].offset, self.zoffs[i].intensity,
                                  np.array(self.zoffs[i].color[:3]) / 255, self.zoffs[i].marker)
        self.plot1.repaint()
        self.make_plot_7(0)

    def make_plot_7(self, e):
        self.plot7.contourplot(xvals=self.xvals, yvals=self.yvals, zgrid=self.newgrid, config=self.config,
                               title="Mass vs. Charge", test_kda=True)

        try:
            for zoff in self.zoffs:
                zwidth = zoff.extractwidth
                self.eshape = self.ctlfilt.GetSelection()
                self.plot7.plot_native_z(zoff.offset, np.array(zoff.color) / 255, self.xvals, width=zwidth, alpha=0.5,
                                         shape=self.eshape)
        except:
            "Failed to plot"
            pass

    def on_edit(self, event):
        index = self.list.ultimateList.GetFirstSelected()
        self.list.ultimateList.SetStringItem(index, 1, str(0))

    def on_delete(self, event):
        self.list.ReturnData()
        btn = event.GetId()
        ids = [int(p.id) for p in self.list.zoffouts]
        # i=ids.index(int(btn))
        i = [i for i, j in enumerate(ids) if j == btn][0]
        index = self.list.zoffouts[i].index
        # print "click",index,ids,btn
        self.list.ultimateList.DeleteItem(index)
        # print "Deleted Test!!!"
        self.update(0)

    def fit(self, e):
        self.update(e)
        guess = []
        for z in self.zoffs:
            guess.append([z.offset, z.intensity * 1 / (0.5 * np.sqrt(2 * np.pi))])
        guess = np.array(guess)
        mf = MassFitter.MassFitter(self.ftot, guess, 0, "smallguess")
        fitspec, fitres = mf.perform_fit("nonorm", "smallguess")
        print "Output", fitres

        for i in xrange(0, len(self.zoffs)):
            fit = fitres[i]
            self.zoffs[i].offset = fit[0]
            self.zoffs[i].width = fit[1]
            self.zoffs[i].intensity = fit[2] * 1 / (fit[1] * np.sqrt(2 * np.pi))
        self.UpdateList()
        self.PlotZoffs()
        for z in self.zoffs:
            yvals = ud.ndis_std(z.offset, self.ftot[:, 0], z.width, a=z.intensity, norm_area=False)
            self.plot1.plotadd(self.ftot[:, 0], yvals, np.array(z.color) / 255, "Fit")
        self.plot1.repaint()

    def update(self, e):
        self.zoffs = self.list.ReturnData()
        for z in self.zoffs:
            if z.offset >= np.amin(self.ftot[:, 0]) and z.offset <= np.amax(self.ftot[:, 0]):
                z.intensity = self.f(z.offset)
            else:
                z.intensity = 0
            self.list.ultimateList.SetStringItem(z.index, 1, str(z.intensity))

    def OnReset(self, e):
        self.GetMaxima()
        self.PlotZoffs()
        self.list.SuperDelete()
        self.PopulateList(e)

    def UpdateList(self):
        for z in self.zoffs:
            index = z.index
            self.list.ultimateList.SetStringItem(index, 1, str(z.intensity))
            self.list.ultimateList.SetStringItem(index, 0, str(z.offset))
            self.list.ultimateList.SetStringItem(index, 2, str(z.width))

    def PopulateList(self, e):
        for i in self.zoffs:
            self.list.AddLine(i)

    def onadd(self, e):
        self.list.AddLineEmpty()

    def MakeFArray(self, min, max):
        zwidth = 1
        frange = np.arange(min, max, 0.5)
        ftot = []

        mgrid, zgrid = np.meshgrid(self.xvals, self.yvals, indexing='ij')
        self.offsetgrid = ud.get_z_offset(mgrid, zgrid)
        for f in frange:
            bool1 = self.offsetgrid >= f - zwidth
            bool2 = self.offsetgrid < f + zwidth
            bool3 = np.all([bool1, bool2], axis=0)
            ftot.append([f, np.sum(self.newgrid[bool3])])
        '''
        zvals=[simchargefit(self.xvals)+F for F in frange]
        for f in xrange(0,len(zvals)):
            intensity=0;
            for x in xrange(0,len(zvals[f])):
                zind=ud.nearest(self.yvals,zvals[f][x])
                start=np.amax([0,zind-zwidth])
                end=np.amin([zind+zwidth+1,self.ylen])
                for i in xrange(start,end):
                    intensity+=self.newgrid[x,i]
            ftot.append([frange[f],intensity])
        print ftot
        '''
        self.ftot = np.array(ftot)
        self.f = interp1d(self.ftot[:, 0], self.ftot[:, 1])

    def fastextract(self, f, width):
        bool1 = self.offsetgrid >= f - width
        bool2 = self.offsetgrid < f + width
        bool3 = np.all([bool1, bool2], axis=0)
        out = self.newgrid * bool3
        if self.eshape == 1:
            wgrid = np.exp(-((self.offsetgrid - f) ** 2.) / (2. * width * width))
            out = out * wgrid
        return out

    def Extract(self, e):
        self.zoffs = self.list.ReturnData()
        self.eshape = self.ctlfilt.GetSelection()
        self.make_plot_7(e)

        tstart = time.clock()
        frange = [z.offset for z in self.zoffs]
        colors = [z.color for z in self.zoffs]
        widths = [z.extractwidth for z in self.zoffs]
        zvals = [simchargefit(self.xvals) + F for F in frange]

        print frange
        offsetnum = [z.nstate for z in self.zoffs]
        self.offsetvalue = float(self.massoffset.GetValue())

        extracts = []
        extract = []
        f = 0

        for i, f in enumerate(frange):
            zwidth = widths[i]
            if self.eshape == 1:
                zbuff = zwidth * 3
            else:
                zbuff = zwidth
            intensity = self.fastextract(f, zwidth)
            extract = np.transpose([self.xvals + offsetnum[i] * self.offsetvalue, np.sum(intensity, axis=1)])
            extracts.append(extract)
        extracts = np.array(extracts)
        f = 0
        print extracts.shape

        self.plot2.plotrefreshtop(extracts[0, :, 0], extracts[0, :, 1], title="Extracted Intensities",
                                  xlabel="Mass (Da)",
                                  ylabel="Intensity", color=np.array(colors[f]) / 255, nticks=5, test_kda=True)
        for f in xrange(1, len(frange)):
            self.plot2.plotadd(extracts[f, :, 0], extracts[f, :, 1], np.array(colors[f]) / 255, "Offset=" + str(f))
        '''
        for x in xrange(0,len(zvals[f])):
            intensity=0;
            zind=ud.nearest(self.yvals,zvals[f][x])
            zwidth=widths[f]
            if self.eshape==1:
                zbuff=zwidth*3
            else:
                zbuff=zwidth
            start=max(0,zind-zbuff)
            end=min(zind+zbuff+1,self.ylen)
            zrange=xrange(start,end)
            weights=np.ones(len(zrange))

            for i,z in enumerate(zrange):
                if self.eshape==1:
                    weights[i]=np.exp(-((z-zind)**2.)/(2.*zwidth*zwidth))
                intensity+=weights[i]*self.newgrid[x,z]
            extract.append([self.xvals[x]+offsetnum[f]*self.offsetvalue,intensity])
        extract=np.array(extract)
        extracts.append(extract)
        self.plot2.plotrefreshtopbox(extract[:,0],extract[:,1], "Extracted Intensities", "Mass (Da)","Intensity",color=np.array(colors[f])/255,bins=5)
        for f in xrange(1,len(zvals)):
            extract=[]
            zwidth=widths[f]
            for x in xrange(0,len(zvals[f])):
                intensity=0;
                zind=ud.nearest(self.yvals,zvals[f][x])
                if self.eshape==1:
                    zbuff=zwidth*3
                else:
                    zbuff=zwidth
                start=max(0,zind-zbuff)
                end=min(zind+zbuff+1,self.ylen)
                zrange=xrange(start,end)
                weights=np.ones(len(zrange))
                for i,z in enumerate(zrange):
                    if self.eshape==1:
                        weights[i]=np.exp(-((z-zind)**2.)/(2.*zwidth*zwidth))
                    intensity+=weights[i]*self.newgrid[x,z]
                extract.append([self.xvals[x]+offsetnum[f]*self.offsetvalue,intensity])
            extract=np.array(extract)
            extracts.append(extract)
            self.plot2.plotadd(extract[:,0],extract[:,1],np.array(colors[f])/255,"Offset="+str(f))
        '''
        self.extracts = np.array(extracts)
        tend = time.clock()
        print "Extraction Time: %.2gs" % (tend - tstart)
        self.plot2.repaint()
        for f in xrange(0, len(zvals)):
            fileout = self.config.outfname + "_extracted_intensities" + str(offsetnum[f]) + ".txt"
            np.savetxt(fileout, self.extracts[f])
        self.PeakExtract()

    def GetMaxima(self):
        max = []
        window = 2
        threshold = 0.1
        defaultcolors = [[255, 0, 0, 255], [0, 0, 255, 255], [0, 255, 0, 255], [255, 0, 255, 255]]
        defaultmarkers = ['o', 'v', '^', '>', 's', 'd', '*']
        self.zoffs = []
        for i in xrange(len(self.ftot) - window - 1, window - 1, -1):
            if self.ftot[i, 1] > np.amax(
                    [np.amax(self.ftot[i - window:i, 1]), np.amax(self.ftot[i + 1:i + window + 1, 1])]):
                if self.ftot[i, 1] > np.amax(self.ftot[:, 1]) * threshold:
                    max.append([i, self.ftot[i, 0], self.ftot[i, 1]])
                    zoff = zoffset()
                    zoff.Make(self.ftot[i, 0], self.ftot[i, 1], i, defaultcolors[len(self.zoffs) % len(defaultcolors)],
                              defaultmarkers[len(self.zoffs)])
                    self.zoffs.append(zoff)
        self.maxes = np.array(max)
        self.PlotZoffs()

    def PeakExtract(self):
        self.peakextracts = np.zeros((len(self.extracts), self.pks.plen))
        self.peakextractsarea = np.zeros((len(self.extracts), self.pks.plen))
        try:
            X = self.pks.masses
        except:
            print "No Peaks to Extract"
            X = None
        if X is not None:
            for i in xrange(0, len(self.extracts)):
                for j in xrange(0, self.pks.plen):
                    p = self.pks.peaks[j]
                    extract = self.extracts[i]
                    xvals = extract[:, 0]
                    yvals = extract[:, 1]
                    index = ud.nearest(xvals, p.mass)
                    self.peakextracts[i, j] = localmax(yvals, index, int(self.config.peakwindow / self.config.massbins))
                    if not ud.isempty(p.integralrange):
                        boo1 = xvals < p.integralrange[1]
                        boo2 = xvals > p.integralrange[0]
                        intdat = extract[np.all([boo1, boo2], axis=0)]
                        integral = np.trapz(intdat[:, 1], x=intdat[:, 0])
                        self.peakextractsarea[i, j] = integral
            try:
                np.savetxt(self.config.outfname + "_extracted_heights.txt", self.peakextracts)
                np.savetxt(self.config.outfname + "_extracted_areas.txt", self.peakextractsarea)
            except:
                print "Error saving files"

            if self.offsetvalue > 0:
                X = np.round(X / self.offsetvalue)
                label = "Subunit Number"
            else:
                label = "Mass (Da)"
            sum1 = np.sum(self.peakextracts, axis=0)
            self.plot3.plotrefreshtop(X, sum1, title="Extracted Heights", xlabel=label, ylabel="Intensity",
                                      integerticks=True, test_kda=True)
            for i in xrange(0, len(self.extracts)):
                self.plot3.plotadd(X, self.peakextracts[i], np.array(self.zoffs[i].color) / 255., str(i))
            self.plot3.repaint()
            sum2 = np.sum(self.peakextractsarea, axis=0)

            X2 = np.arange(0, len(sum1))
            cmap = self.config.peakcmap
            colormap = cm.get_cmap(cmap, len(X))
            peakcolors = colormap(np.arange(len(X)))
            self.plot5.barplottop(X2, sum1 / np.amax(sum1), [int(i) for i in X], peakcolors, label,
                                  "Normalized Intensity", "Extracted Total Peak Heights")
            self.plot5.repaint()
            if np.amax(sum2) != 0:
                self.plot4.plotrefreshtop(X, sum2, title="Extracted Areas", xlabel=label, ylabel="Area",
                                          integerticks=True, test_kda=True)
                for i in xrange(0, len(self.extracts)):
                    self.plot4.plotadd(X, self.peakextractsarea[i], np.array(self.zoffs[i].color) / 255., str(i))
                self.plot4.repaint()
                self.plot6.barplottop(X2, sum2 / np.amax(sum2), [int(i) for i in X], peakcolors, label,
                                      "Normalized Intensity", "Extracted Total Peak Areas")
                self.plot6.repaint()
            else:
                print "No integration provided"

            try:
                np.savetxt(self.config.outfname + "_total_extracted_heights.txt",
                           np.transpose([self.pks.masses, X, sum1 / np.amax(sum1)]))
                if np.amax(sum2) != 0:
                    np.savetxt(self.config.outfname + "_total_extracted_areas.txt",
                               np.transpose([self.pks.masses, X, sum2 / np.amax(sum2)]))
            except:
                print "Error saving total files"

    def SaveFig(self, e):
        extraheader = "_NativeZ"
        name1 = self.config.outfname + extraheader + "_Figure1.pdf"
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
        name2 = self.config.outfname + extraheader + "_Figure2.pdf"
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
        name3 = self.config.outfname + extraheader + "_Figure3.pdf"
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
        name4 = self.config.outfname + extraheader + "_Figure4.pdf"
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name4)
        name5 = self.config.outfname + extraheader + "_Figure5.pdf"
        if self.plot5.flag:
            self.plot5.on_save_fig(e, name5)
        name6 = self.config.outfname + extraheader + "_Figure6.pdf"
        if self.plot6.flag:
            self.plot6.on_save_fig(e, name6)
        name7 = self.config.outfname + extraheader + "_Figure7.pdf"
        if self.plot7.flag:
            self.plot7.on_save_fig(e, name7)

    def OnClose(self, e):
        self.update(e)
        self.zoffouts = self.zoffs
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.zoffouts = []
        self.Destroy()
        self.EndModal(1)


class ColorList(wx.Panel):
    def __init__(self, parent):
        """Constructor"""
        wx.Panel.__init__(self, parent)

        # font = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        # boldfont = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        # boldfont.SetWeight(wx.BOLD)
        # boldfont.SetPointSize(12)

        self.ultimateList = ULC.UltimateListCtrl(self, size=(500, 310), agwStyle=wx.LC_REPORT
                                                                                 | wx.LC_VRULES | wx.LC_EDIT_LABELS
                                                                                 | wx.LC_HRULES | ULC.ULC_USER_ROW_HEIGHT
                                                                                 | ULC.ULC_HAS_VARIABLE_ROW_HEIGHT
                                                                                 | ULC.ULC_EDIT_LABELS)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._image = []
        info._format = 0
        info._kind = 1
        info._text = "Native Z Offset"
        self.ultimateList.InsertColumnInfo(0, info)

        info = ULC.UltimateListItem()
        info._format = wx.LIST_FORMAT_RIGHT
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT | ULC.ULC_MASK_CHECK
        info._image = []
        info._text = "Intensity"
        # info._font = boldfont
        self.ultimateList.InsertColumnInfo(1, info)

        info = ULC.UltimateListItem()
        info._format = wx.LIST_FORMAT_RIGHT
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT | ULC.ULC_MASK_CHECK
        info._image = []
        info._text = "Width"
        # info._font = boldfont
        self.ultimateList.InsertColumnInfo(2, info)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "Color"
        # info._font = font
        info._image = []
        self.ultimateList.InsertColumnInfo(3, info)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "ID"
        # info._font = font
        info._image = []
        self.ultimateList.InsertColumnInfo(4, info)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = ""
        # info._font = font
        info._image = []
        self.ultimateList.InsertColumnInfo(5, info)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "N"
        # info._font = font
        info._image = []
        self.ultimateList.InsertColumnInfo(6, info)

        info = ULC.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "Z Width"
        # info._font = font
        info._image = []
        self.ultimateList.InsertColumnInfo(7, info)

        self.ultimateList.SetColumnWidth(0, 100)
        self.ultimateList.SetColumnWidth(1, 100)
        self.ultimateList.SetColumnWidth(2, 100)
        self.ultimateList.SetColumnWidth(3, 100)
        self.ultimateList.SetColumnWidth(4, 20)
        self.ultimateList.SetColumnWidth(5, 100)
        self.ultimateList.SetColumnWidth(6, 50)
        self.ultimateList.SetColumnWidth(7, 50)
        self.ultimateList.SetUserLineHeight(25)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ultimateList, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.buttontot = 0

    def AddLineEmpty(self):
        index = self.ultimateList.InsertStringItem(sys.maxint, str(0))
        self.ultimateList.SetStringItem(index, 1, str(0))
        self.ultimateList.SetStringItem(index, 2, str(0))
        size = (50, 10)
        self.colorbox = wx.ColourPickerCtrl(self.ultimateList, size=size)
        self.ultimateList.SetItemWindow(index, col=3, wnd=self.colorbox, expand=True)
        self.ultimateList.SetStringItem(index, 4, str(self.buttontot))
        deletebutton = wx.Button(self.ultimateList, label="Delete", id=self.buttontot)
        self.ultimateList.SetItemWindow(index, col=5, wnd=deletebutton, expand=True)
        self.buttontot += 1
        textinput = wx.TextCtrl(self.ultimateList, value="0")
        self.ultimateList.SetItemWindow(index, col=6, wnd=textinput, expand=True)
        textinput2 = wx.TextCtrl(self.ultimateList, value="1")
        self.ultimateList.SetItemWindow(index, col=7, wnd=textinput2, expand=True)

    def AddLine(self, zoff):
        index = self.ultimateList.InsertStringItem(sys.maxint, str(zoff.offset))
        self.ultimateList.SetStringItem(index, 1, str(zoff.intensity))
        self.ultimateList.SetStringItem(index, 2, str(zoff.width))
        size = (50, 10)
        colorarray = list(zoff.color)
        if len(colorarray) < 4:
            colorarray.append(255)
        self.colorbox = wx.ColourPickerCtrl(self.ultimateList, size=size,
                                            col=wx.Colour(colorarray[0], colorarray[1], colorarray[2],
                                                          alpha=colorarray[3]))
        self.ultimateList.SetItemWindow(index, col=3, wnd=self.colorbox, expand=True)
        self.ultimateList.SetStringItem(index, 4, str(self.buttontot))
        deletebutton = wx.Button(self.ultimateList, label="Delete", id=self.buttontot)
        self.ultimateList.SetItemWindow(index, col=5, wnd=deletebutton, expand=True)
        self.buttontot += 1
        textinput = wx.TextCtrl(self.ultimateList, value=str(zoff.nstate))
        self.ultimateList.SetItemWindow(index, col=6, wnd=textinput, expand=True)
        textinput2 = wx.TextCtrl(self.ultimateList, value="1")
        self.ultimateList.SetItemWindow(index, col=7, wnd=textinput2, expand=True)

    def ReturnData(self):
        count = self.ultimateList.GetItemCount()
        self.zoffouts = []
        # print "Count",count
        defaultmarkers = ['o', 'v', '^', '>', 's', 'd', '*']
        for index in xrange(0, count):
            zout = zoffset()
            zout.index = index
            # print index
            colorwindow = self.ultimateList.GetItemWindow(index, col=3)
            zout.color = colorwindow.GetColour()
            zout.offset = float(self.ultimateList.GetItem(index, col=0).GetText())
            zout.width = float(self.ultimateList.GetItem(index, col=2).GetText())
            # zout.intensity=float(self.ultimateList.GetItem(index,col=1).GetText())
            zout.id = int(self.ultimateList.GetItem(index, col=4).GetText())
            zout.marker = defaultmarkers[len(self.zoffouts)]
            inputbox = self.ultimateList.GetItemWindow(index, col=6)
            zout.nstate = int(inputbox.GetValue())
            inputbox = self.ultimateList.GetItemWindow(index, col=7)
            zout.extractwidth = int(inputbox.GetValue())
            self.zoffouts.append(zout)
        return self.zoffouts

    def SuperDelete(self):
        count = self.ultimateList.GetItemCount()
        topcount = count
        num = 0
        while count > 0 and num <= topcount:
            try:
                self.ultimateList.DeleteItemWindow(0, col=3)
                self.ultimateList.DeleteItem(0)
                count = self.ultimateList.GetItemCount()
                num += 1
            except:
                num += 1
                pass
        self.ultimateList.DeleteAllItems()

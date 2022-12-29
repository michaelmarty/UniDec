import numpy as np
import wx
import unidec.tools as ud
from unidec.modules import PlottingWindow, peakstructure, unidecstructure


class SubDivFrame(wx.Frame):
    def __init__(self, parent, data, pks, config=None):
        wx.Frame.__init__(self, parent, title="Subtract and Divide Tools", size=(1100, 600))  # ,size=(-1,-1))
        self.parent = parent
        self.div = 760
        self.sub = 44088
        self.data = data
        self.pks = pks
        self.avgs = None
        self.nums = None
        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config
            self.sub = config.molig

        try:
            self.masses = []
            for p in self.pks.peaks:
                if p.ignore == 0:
                    self.masses.append(p.mass)
            defaultsub = np.amin(self.masses)
        except:
            self.masses = None
            defaultsub = 44088
        self.sub = defaultsub

        self.panel = wx.Panel(self)
        self.plot1 = PlottingWindow.Plot1d(self.panel)
        self.plot2 = PlottingWindow.Plot1d(self.panel)
        self.ctlsub = wx.TextCtrl(self.panel, value=str(self.sub))
        self.ctldiv = wx.TextCtrl(self.panel, value=str(self.div))
        self.ctlresult = wx.TextCtrl(self.panel, value="")
        self.ctltype = wx.RadioBox(self.panel, label="X-axis", choices=["Mass", "Number"])
        self.replotbutton = wx.Button(self.panel, -1, "Replot")
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        s1 = wx.BoxSizer(wx.HORIZONTAL)
        s1.Add(self.plot1, 1, wx.GROW)
        s1.Add(self.plot2, 1, wx.GROW)
        self.sizer.Add(s1, 1)

        s2 = wx.BoxSizer(wx.HORIZONTAL)
        s2.Add(wx.StaticText(self.panel, label="  Subtract (Da): "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlsub, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(wx.StaticText(self.panel, label="  Divide (Da): "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctldiv, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctltype, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        s2.Add(self.replotbutton, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        s2.Add(wx.StaticText(self.panel, label="  Result (Da): "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlresult, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        self.sizer.Add(s2, 0)

        self.Bind(wx.EVT_BUTTON, self.replot, self.replotbutton)
        self.Bind(wx.EVT_RADIOBOX, self.change_type, self.ctltype)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()
        self.Show(True)
        self.replot()

    def make_plot(self, data=None):
        if data is not None:
            self.data = data

        self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1],
                                  "Zero-charge Mass Spectrum", "Mass (Da)",
                                  "Intensity", "Mass Distribution", self.config, test_kda=True,
                                  nopaint=True)
        if self.pks.plen > 0:
            for p in self.pks.peaks:
                if p.ignore == 0:
                    self.plot1.plotadddot(p.mass, p.height, p.color, p.marker)
        self.plot1.repaint()

        if self.masses is not None:
            type = self.ctltype.GetSelection()
            if type == 0:
                x = self.masses
                xlabel = "Mass (Da)"
                test_kda = True

            elif type == 1:
                x = self.nums
                xlabel = "Number of Monomers"
                test_kda = False

            self.plot2.plotrefreshtop(x, self.avgs, "Average Masses", xlabel, "Average Monomer Mass",
                                      test_kda=test_kda)

            if self.pks.plen > 0:
                for p in self.pks.peaks:
                    if p.ignore == 0:
                        if type == 0:
                            x = p.mass
                        elif type == 1:
                            x = p.sdnum
                        self.plot2.plotadddot(x, p.sdval, p.color, p.marker)
            self.plot2.repaint()

    def subdiv(self):
        try:
            sub = float(self.sub)
            div = float(self.div)
        except:
            print("Error with Subtract and Divide Inputs:", self.sub, self.div)
            return 0
        '''
        try:
            sd_results = []
            for s in self.data.spectra:
                pks = peakstructure.Peaks()
                peaks = ud.peakdetect(s.massdat, self.config)
                pks.add_peaks(peaks, massbins=self.config.massbins)
                sd_result, avgs, ints = ud.subtract_and_divide(pks, sub, div, outputall=True)
                sd_results.append(sd_result)
            avg = np.mean(sd_results)
            self.ctlresult.SetValue(str(avg))

        except:'''

        sd_result, self.avgs, ints, self.nums, self.masses = ud.subtract_and_divide(self.pks, sub, div, outputall=True)
        self.ctlresult.SetValue(str(sd_result))

    def replot(self, e=None):
        sub = self.ctlsub.GetValue()
        try:
            sub = float(sub)
        except:
            sub = None
        self.sub = sub

        div = self.ctldiv.GetValue()
        try:
            div = float(div)
        except:
            div = None
        self.div = div
        self.subdiv()
        self.make_plot()

    def change_type(self, e=None):
        self.replot()


if __name__ == "__main__":
    path = "/unidec/bin\\Example Data\\POPC_Nanodiscs_unidecfiles\\POPC_Nanodiscs_mass.txt"
    data = np.loadtxt(path)

    pks = peakstructure.Peaks()
    peaks = ud.peakdetect(data, threshold=0.1)
    pks.add_peaks(peaks)

    app = wx.App()
    panel = SubDivFrame(None, data, pks)
    panel.replot()
    app.MainLoop()

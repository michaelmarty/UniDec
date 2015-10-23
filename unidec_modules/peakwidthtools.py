import wx
from scipy.optimize import curve_fit
import numpy as np

from unidec_modules import plot1d
from unidec_modules.unidectools import make_peak_shape, simp_string_to_value

__author__ = 'Michael.Marty'


class PeakTools1d(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((700, 700))
        self.SetTitle("Peak Fitting Tools")

    def InitUI(self, config, data):
        self.config = config
        self.data = data
        self.topmax = np.amax(self.data[:, 0])
        self.topmin = np.amin(self.data[:, 0])
        self.mzsig = config.mzsig
        self.length = len(self.data)
        self.psfun = config.psfun

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Peak Shape Tool')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.plot1 = plot1d.Plot1d(self.pnl)

        self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")
        # self.plot1.plotadd(self.x,self.peakshape, "blue", "Peak Shape Guess")
        self.sbs.Add(self.plot1, 1, wx.EXPAND)

        self.hbox11 = wx.BoxSizer(wx.HORIZONTAL)

        self.hbox10 = wx.BoxSizer(wx.VERTICAL)
        self.centerbutton = wx.Button(self.pnl, label='Reset Range')
        self.centerbutton.Bind(wx.EVT_BUTTON, self.OnReset)
        self.hbox10.Add(self.centerbutton, 0, wx.ALIGN_LEFT)
        self.hbox11.Add(self.hbox10, 0, wx.ALIGN_LEFT)

        self.hbox9 = wx.BoxSizer(wx.VERTICAL)
        self.ctlpsfun = wx.RadioBox(self.pnl, label="Peak Shape Fit Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        self.ctlpsfun.SetSelection(self.psfun)
        self.fitbutton = wx.Button(self.pnl, label='Fit Peak Shape')
        self.plotbutton = wx.Button(self.pnl, label='Plot Guess')
        self.hbox9.Add(self.ctlpsfun, 0, wx.ALIGN_RIGHT)
        self.hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmzsig = wx.TextCtrl(self.pnl, value=str(self.mzsig))
        self.hbox8.Add(wx.StaticText(self.pnl, label='Guess for Peak FWHM: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox8.Add(self.ctlmzsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox9.Add(self.hbox8, 0, wx.ALIGN_RIGHT)
        self.hbox9.Add(self.plotbutton, 0, wx.ALIGN_RIGHT)
        self.hbox9.Add(self.fitbutton, 0, wx.ALIGN_RIGHT)
        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.errorbox = wx.TextCtrl(self.pnl, value="", style=wx.TE_READONLY)
        self.hbox6.Add(wx.StaticText(self.pnl, label='Error in Fit: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox6.Add(self.errorbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox9.Add(self.hbox6, 0, wx.ALIGN_RIGHT)

        self.hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        self.resbox = wx.TextCtrl(self.pnl, value=str(""))
        self.hbox7.Add(wx.StaticText(self.pnl, label='Resolution (M/FWHM): '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox7.Add(self.resbox, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox9.Add(self.hbox7, 0, wx.ALIGN_RIGHT)

        self.hbox11.Add(self.hbox9, 1, wx.ALIGN_RIGHT)
        self.sbs.Add(self.hbox11, 0, wx.EXPAND)
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

        self.plotbutton.Bind(wx.EVT_BUTTON, self.OnPlot)
        self.fitbutton.Bind(wx.EVT_BUTTON, self.OnFit)
        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)
        self.CenterOnParent()

    def OnClose(self, e):
        self.config.psfun = self.ctlpsfun.GetSelection()
        self.config.mzsig = simp_string_to_value(self.ctlmzsig.GetValue())
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def OnReset(self, e):
        self.ctlmzsig.SetValue(str(self.mzsig))
        newmin = np.amin(self.data[:, 0])
        newmax = np.amax(self.data[:, 0])
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")

    def OnCenter(self, e):
        xlim = self.plot1.subplot1.get_xlim()
        newmin = xlim[0]
        newmax = xlim[1]
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")

    def OnPlot(self, e):
        self.OnCenter(e)
        self.psfun = self.ctlpsfun.GetSelection()
        self.mzsig2 = simp_string_to_value(self.ctlmzsig.GetValue())
        self.mid = np.argmax(self.centdat[:, 1])
        self.scale = np.amax(self.centdat[:, 1])
        self.psguess = make_peak_shape(self.centdat[:, 0], self.psfun, self.mzsig2, self.centdat[self.mid, 0])
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")
        self.plot1.plotadd(self.centdat[:, 0], self.psguess * self.scale, "blue", "Peak Shape Guess", nopaint=False)

    def OnFit(self, e):
        self.OnCenter(e)
        self.psfun = self.ctlpsfun.GetSelection()

        def f0(x, s, m):
            s2 = s / 2.35482
            return np.exp(-(x - m) * (x - m) / (2. * s2 * s2))

        def f1(x, s, m):
            return s / 2. * (s / 2.) / ((x - m) * (x - m) + (s / 2.) * (s / 2.))

        def f2(x, s, m):
            if x > 0:
                return s / 2. * (s / 2.) / ((x - m) * (x - m) + (s / 2.) * (s / 2.))
            else:
                return np.exp(-(x - m) * (x - m) / (2 * s * s))

        def fitfunc_vec_self(x, p, m):
            y = np.zeros(x.shape)
            for i in range(len(y)):
                y[i] = f2(x[i], p, m)
            return y

        def fit(x2, cent, psfun, mid, sig):
            guess = [sig, mid]
            fits = [curve_fit(f0, x2, cent, p0=guess)[0], curve_fit(f1, x2, cent, p0=guess)[0],
                    curve_fit(fitfunc_vec_self, x2, cent, p0=guess)[0]]
            print fits
            if psfun == 0:
                return [abs(fits[0])]
            if psfun == 1:
                return [abs(fits[1])]
            if psfun == 2:
                return [abs(fits[2])]

        self.scale = np.amax(self.centdat[:, 1])
        self.mid = np.argmax(self.centdat[:, 1])
        self.mzsig2 = simp_string_to_value(self.ctlmzsig.GetValue())
        fitout = \
            fit(self.centdat[:, 0], self.centdat[:, 1] / self.scale, self.psfun, self.centdat[self.mid, 0],
                self.mzsig2)[0]

        self.psguess = make_peak_shape(self.centdat[:, 0], self.psfun, fitout[0], fitout[1])
        self.res = self.centdat[self.mid, 0] / fitout[0]
        self.ctlmzsig.SetValue(str(fitout[0]))
        self.resbox.SetValue(str(self.res))
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1],
                                  title="Data", xlabel="m/z (Th)", ylabel="Normalize Intensity", zoom="span")
        self.plot1.plotadd(self.centdat[:, 0], self.psguess * self.scale, "blue", "Peak Shape Guess", nopaint=False)

        self.cent = np.array(self.centdat[:, 1])
        self.cent2 = np.array(self.psguess)
        error = np.sum((self.cent2 - self.cent) * (self.cent2 - self.cent))
        self.errorbox.SetValue(str(error))
        pass


class PeakTools2d(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((700, 700))
        self.SetTitle("Peak Fitting Tools")

    def InitUI(self, data3, data2, config):
        self.mz = np.unique(data3[:, 0])
        self.dt = np.unique(data3[:, 1])
        self.int = data2[:, 1]
        self.C = data3[:, 2].reshape((len(self.mz), len(self.dt)))
        self.data = data2
        self.topdata = self.data
        self.topmax = np.amax(self.data[:, 0])
        self.topmin = np.amin(self.data[:, 0])
        self.config = config
        self.mzsig = config.mzsig
        self.dtsig = config.dtsig
        self.length = len(self.data)
        self.psfun = config.psfun

        self.fitmid = 0
        self.fitsig = 0

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Peak Shape Tool')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.plot1 = plot1d.Plot1d(self.pnl)

        self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")
        self.sbs.Add(self.plot1, 1, wx.EXPAND)

        self.hbox11 = wx.BoxSizer(wx.HORIZONTAL)

        self.hbox10 = wx.BoxSizer(wx.VERTICAL)

        self.centerbutton = wx.Button(self.pnl, label='Reset Range')
        self.flipbutton = wx.Button(self.pnl, label='Flip m/z and AT')
        self.centerbutton.Bind(wx.EVT_BUTTON, self.OnReset)
        self.flipbutton.Bind(wx.EVT_BUTTON, self.OnFlip)

        self.hboxsigs = wx.BoxSizer(wx.HORIZONTAL)
        self.outmzsig = wx.TextCtrl(self.pnl, value="")
        self.hboxsigs.Add(wx.StaticText(self.pnl, label='m/z Width Fit: '), wx.ALIGN_CENTER_VERTICAL)
        self.hboxsigs.Add(self.outmzsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        self.hboxsigs2 = wx.BoxSizer(wx.HORIZONTAL)
        self.outdtsig = wx.TextCtrl(self.pnl, value="")
        self.hboxsigs2.Add(wx.StaticText(self.pnl, label='A.T. Width Fit: '), wx.ALIGN_CENTER_VERTICAL)
        self.hboxsigs2.Add(self.outdtsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        self.hbox10.Add(self.centerbutton, 0, wx.ALIGN_LEFT)
        self.hbox10.Add(self.flipbutton, 0, wx.ALIGN_LEFT)
        self.hbox10.Add(self.hboxsigs, 0, wx.ALIGN_LEFT)
        self.hbox10.Add(self.hboxsigs2, 0, wx.ALIGN_LEFT)
        self.hbox11.Add(self.hbox10, 0, wx.ALIGN_LEFT)

        self.hbox9 = wx.BoxSizer(wx.VERTICAL)
        self.ctlpsfun = wx.RadioBox(self.pnl, label="Peak Shape Fit Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        self.ctlpsfun.SetSelection(self.psfun)
        self.fitbutton = wx.Button(self.pnl, label='Fit Peak Shape')
        self.plotbutton = wx.Button(self.pnl, label='Plot Guess')
        self.hbox9.Add(self.ctlpsfun, 0, wx.ALIGN_RIGHT)
        self.hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmzsig = wx.TextCtrl(self.pnl, value=str(self.mzsig))
        self.hbox8.Add(wx.StaticText(self.pnl, label='Guess for Peak Width: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox8.Add(self.ctlmzsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox9.Add(self.hbox8, 0, wx.ALIGN_RIGHT)
        self.hbox9.Add(self.plotbutton, 0, wx.ALIGN_RIGHT)
        self.hbox9.Add(self.fitbutton, 0, wx.ALIGN_RIGHT)
        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.errorbox = wx.TextCtrl(self.pnl, value="", style=wx.TE_READONLY)
        self.hbox6.Add(wx.StaticText(self.pnl, label='Error in Fit: '), wx.ALIGN_CENTER_VERTICAL)
        self.hbox6.Add(self.errorbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox9.Add(self.hbox6, 0, wx.ALIGN_RIGHT)
        self.hbox11.Add(self.hbox9, 1, wx.ALIGN_RIGHT)
        self.sbs.Add(self.hbox11, 0, wx.EXPAND)
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

        self.plotbutton.Bind(wx.EVT_BUTTON, self.OnPlot)
        self.fitbutton.Bind(wx.EVT_BUTTON, self.OnFit)
        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)
        self.CenterOnParent()
        self.flag = 0

    def OnClose(self, e):
        self.config.mzsig = simp_string_to_value(self.outmzsig.GetValue())
        self.config.dtsig = simp_string_to_value(self.outdtsig.GetValue())
        self.config.psfun = self.psfun
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def OnFlip(self, e):
        if self.flag == 0:
            self.psfun = float(self.ctlpsfun.GetSelection())
            self.pos = np.argmin(np.abs(self.mz - self.fitmid))
            self.intdt = self.C[self.pos]
            self.data2 = np.column_stack((self.dt, self.intdt))
            self.data, self.data2 = self.data2, self.data
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel="Drift Time (ms)",
                                      ylabel="Normalize Intensity", zoom="span")
            self.outmzsig.SetValue(str(self.fitsig))
            self.flag = 1
            if self.outdtsig.GetValue() != "":
                self.ctlmzsig.SetValue(self.outdtsig.GetValue())
            else:
                self.ctlmzsig.SetValue(str(self.dtsig))
            self.ctlpsfun.SetSelection(0)
        else:
            self.ctlpsfun.SetSelection(self.psfun)
            self.data, self.data2 = self.data2, self.data
            self.outdtsig.SetValue(str(self.fitsig))
            self.OnReset(0)
            self.flag = 0
            self.ctlmzsig.SetValue(self.outmzsig.GetValue())
        pass

    def OnReset(self, e):
        if self.flag == 0:
            self.ctlmzsig.SetValue(str(self.mzsig))
        else:
            self.ctlmzsig.SetValue(str(self.dtsig))
        newmin = np.amin(self.data[:, 0])
        newmax = np.amax(self.data[:, 0])
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")

    def OnCenter(self, e):
        xlim = self.plot1.subplot1.get_xlim()
        newmin = xlim[0]
        newmax = xlim[1]
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")

    def OnPlot(self, e):
        self.OnCenter(e)
        self.mzsig2 = simp_string_to_value(self.ctlmzsig.GetValue())
        self.mid = np.argmax(self.centdat[:, 1])
        self.scale = np.amax(self.centdat[:, 1])
        if self.flag == 0:
            self.psfun = self.ctlpsfun.GetSelection()
            self.psguess = make_peak_shape(self.centdat[:, 0], self.psfun, self.mzsig2, self.centdat[self.mid, 0])
        else:
            self.psguess = make_peak_shape(self.centdat[:, 0], 0, self.mzsig2, self.centdat[self.mid, 0])

        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")

        self.plot1.plotadd(self.centdat[:, 0], self.psguess * self.scale, "blue", "Peak Shape Guess", nopaint=False)

    def OnFit(self, e):
        self.OnCenter(e)

        def f0(x, s, m, a):
            s /= 2.35482
            return a * np.exp(-(x - m) * (x - m) / (2. * s * s))

        def f1(x, s, m, a):
            return a * s / 2. * (s / 2.) / ((x - m) * (x - m) + (s / 2.) * (s / 2.))

        def f2(x, s, m, a):
            if x > 0:
                return a * s / 2. * (s / 2.) / ((x - m) * (x - m) + (s / 2.) * (s / 2.))
            else:
                return a * np.exp(-(x - m) * (x - m) / (2. * s * s))

        def fitfunc_vec_self(x, p, m, a):
            y = np.zeros(x.shape)
            for i in range(len(y)):
                y[i] = f2(x[i], p, m, a)
            return y

        def fit(x2, cent, psfun, mid, sig):
            guess = [sig, mid, 1]
            fits = [curve_fit(f0, x2, cent, p0=guess)[0], curve_fit(f1, x2, cent, p0=guess)[0],
                    curve_fit(fitfunc_vec_self, x2, cent, p0=guess)[0]]
            # print fits
            if psfun == 0:
                return [abs(fits[0])]
            if psfun == 1:
                return [abs(fits[1])]
            if psfun == 2:
                return [abs(fits[2])]

        self.scale = np.amax(self.centdat[:, 1])
        self.mid = np.argmax(self.centdat[:, 1])
        self.mzsig2 = simp_string_to_value(self.ctlmzsig.GetValue())
        if self.flag == 0:
            self.psfun = self.ctlpsfun.GetSelection()
            fitout = fit(self.centdat[:, 0], self.centdat[:, 1] / self.scale, self.psfun, self.centdat[self.mid, 0],
                         self.mzsig2)[0]
        else:
            fitout = \
                fit(self.centdat[:, 0], self.centdat[:, 1] / self.scale, 0, self.centdat[self.mid, 0], self.mzsig2)[0]
        self.fitmid = fitout[1]
        self.fitsig = fitout[0]
        print self.fitmid, self.fitsig
        if self.flag == 0:
            self.psguess = make_peak_shape(self.centdat[:, 0], self.psfun, fitout[0], fitout[1])
            self.outmzsig.SetValue(str(fitout[0]))
        else:
            self.psguess = make_peak_shape(self.centdat[:, 0], 0, fitout[0], fitout[1])
            self.outdtsig.SetValue(str(fitout[0]))

        self.ctlmzsig.SetValue(str(fitout[0]))

        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalize Intensity", zoom="span")
        self.plot1.plotadd(self.centdat[:, 0], self.psguess * self.scale, "blue", "Peak Shape Guess", nopaint=False)

        self.cent = np.array(self.centdat[:, 1])
        self.cent2 = np.array(self.psguess)
        error = np.sum((self.cent2 - self.cent) * (self.cent2 - self.cent))
        self.errorbox.SetValue(str(error))
        pass

import wx
import numpy as np
from unidec.modules import PlottingWindow
from unidec.tools import make_peak_shape, simp_string_to_value, isolated_peak_fit

__author__ = 'Michael.Marty'


class PeakTools1d(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create a dialog for fitting the peak shapes from the data.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((700, 740))
        self.SetTitle("Peak Fitting Tools")
        self.config = None
        self.data = None
        self.centdat = None
        self.psfun = 0

        self.plot1 = None
        self.resbox = None
        self.ctlpsfun = None
        self.ctlmzsig = None
        self.errorbox = None

    def initialize_interface(self, config, data):
        """
        Create the GUI, import the paramters, and plot the intial data.
        :param config: UniDecConfig Object
        :param data: Data array (N x 2)
        :return: None
        """
        # Initialize parameters
        self.config = config
        self.data = data

        self.psfun = config.psfun

        # Create the GUI
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Peak Shape Tool')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot1 = PlottingWindow.Plot1d(pnl)

        self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)
        sbs.Add(self.plot1, 0, wx.EXPAND)

        hbox11 = wx.BoxSizer(wx.HORIZONTAL)

        hbox10 = wx.BoxSizer(wx.VERTICAL)
        centerbutton = wx.Button(pnl, label='Reset Range')
        centerbutton.Bind(wx.EVT_BUTTON, self.on_reset)
        hbox10.Add(centerbutton, 0, wx.ALIGN_LEFT)
        hbox11.Add(hbox10, 0, wx.ALIGN_LEFT)

        hbox9 = wx.BoxSizer(wx.VERTICAL)
        self.ctlpsfun = wx.RadioBox(pnl, label="Peak Shape Fit Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        self.ctlpsfun.SetSelection(self.psfun)
        fitbutton = wx.Button(pnl, label='Fit Peak Shape')
        plotbutton = wx.Button(pnl, label='Plot Guess')
        hbox9.Add(self.ctlpsfun, 0, wx.ALIGN_RIGHT)
        hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlmzsig = wx.TextCtrl(pnl, value=str(self.config.mzsig))
        hbox8.Add(wx.StaticText(pnl, label='Guess for Peak FWHM: '), wx.ALIGN_CENTER_VERTICAL)
        hbox8.Add(self.ctlmzsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox9.Add(hbox8, 0, wx.ALIGN_RIGHT)
        hbox9.Add(plotbutton, 0, wx.ALIGN_RIGHT)
        hbox9.Add(fitbutton, 0, wx.ALIGN_RIGHT)
        hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.errorbox = wx.TextCtrl(pnl, value="", style=wx.TE_READONLY)
        hbox6.Add(wx.StaticText(pnl, label='Error in Fit: '), wx.ALIGN_CENTER_VERTICAL)
        hbox6.Add(self.errorbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox9.Add(hbox6, 0, wx.ALIGN_RIGHT)

        hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        self.resbox = wx.TextCtrl(pnl, value=str(""))
        hbox7.Add(wx.StaticText(pnl, label='Resolution (M/FWHM): '), wx.ALIGN_CENTER_VERTICAL)
        hbox7.Add(self.resbox, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox9.Add(hbox7, 0, wx.ALIGN_RIGHT)

        hbox11.Add(hbox9, 1, wx.EXPAND)
        sbs.Add(hbox11, 0, wx.EXPAND)
        pnl.SetSizer(sbs)
        # Create the bottom buttons
        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)

        # Bind some things
        plotbutton.Bind(wx.EVT_BUTTON, self.on_plot)
        fitbutton.Bind(wx.EVT_BUTTON, self.on_fit)
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        self.CenterOnParent()

    def on_close(self, e):
        """
        Close the window and update self.config.psfun and self.config.mzsig from the parameters in the GUI.
        :param e: Unused event
        :return: None
        """
        self.config.psfun = self.ctlpsfun.GetSelection()
        self.config.mzsig = simp_string_to_value(self.ctlmzsig.GetValue())
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but do not update the parameters.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def on_reset(self, e):
        """
        Reset the plot to the full range.
        :param e: Unused event
        :return: 
        """
        self.ctlmzsig.SetValue(str(self.config.mzsig))
        newmin = np.amin(self.data[:, 0])
        newmax = np.amax(self.data[:, 0])
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)

    def on_center(self, e):
        """
        Crop the range to the values from the plot and replot the cropped data.
        Note this is a bit clumsy. Hopefully someone will come along and make this more eligant someday...
        :param e: Unused event
        :return: None
        """
        xlim = self.plot1.subplot1.get_xlim()
        newmin = xlim[0]
        newmax = xlim[1]
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)

    def on_plot(self, e):
        """
        For the sigma value in the self.ctrlmzsig, plot a simulated peak.
        :param e: Unused event
        :return: None
        """
        self.on_center(e)
        self.psfun = self.ctlpsfun.GetSelection()
        sigguess = simp_string_to_value(self.ctlmzsig.GetValue())
        midguess = np.argmax(self.centdat[:, 1])
        aguess = np.amax(self.centdat[:, 1])
        fitdat = make_peak_shape(self.centdat[:, 0], self.psfun, sigguess, self.centdat[midguess, 0])
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)
        self.plot1.plotadd(self.centdat[:, 0], fitdat * aguess, "blue", "Peak Shape Guess", nopaint=False)

    def on_fit(self, e):
        """
        Perform the fit, update the GUI, and plot the results.
        :param e: Unused event
        :return: None
        """
        self.on_center(e)
        self.psfun = self.ctlpsfun.GetSelection()
        # Fit the results
        fitout, fitdat = isolated_peak_fit(self.centdat[:, 0], self.centdat[:, 1], self.psfun)
        fitout = fitout[:, 0]
        resolution = self.centdat[np.argmax(self.centdat[:, 1]), 0] / fitout[0]
        error = np.sum((fitdat - self.centdat[:, 1]) * (fitdat - self.centdat[:, 1]))

        # Update the GUI
        self.ctlmzsig.SetValue(str(fitout[0]))
        self.resbox.SetValue(str(resolution))
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1],
                                  title="Data", xlabel="m/z (Th)", ylabel="Normalized Intensity", zoom="span", zoomout=True)
        self.plot1.plotadd(self.centdat[:, 0], fitdat, "blue", "Peak Shape Guess", nopaint=False)
        self.errorbox.SetValue(str(error))
        pass


class PeakTools2d(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create a dialog for fitting peak shapes in 2D IM-MS data. Fits 1D first in MS then flips to other dimension.
        :param args:
        :param kwargs:
        :return:
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        self.SetSize((700, 750))
        self.SetTitle("Peak Fitting Tools")
        self.plot1 = None
        self.flipbutton = None
        self.outmzsig = None
        self.outdtsig = None
        self.ctlpsfun = None
        self.ctlsigguess = None
        self.errorbox = None
        self.resbox = None

        self.config = None
        self.psfun = 0
        self.flipflag = 0
        self.centdat = np.array([])
        self.data = np.array([])
        self.data2 = np.array([])
        self.mz = np.array([])
        self.dt = np.array([])
        self.C = np.array([])
        self.fitmid = 0
        self.fitsig = 0

    def initialize_interface(self, data3, data2, config):
        """
        Initialize the parameters, create the interface, plot the initial data.
        :param data3: 2D IM-MS data
        :param data2: 1D MS data
        :param config: UniDecConfig Object
        :return: None
        """
        # Initialize the parameters
        self.mz = np.unique(data3[:, 0])
        self.dt = np.unique(data3[:, 1])
        self.C = data3[:, 2].reshape((len(self.mz), len(self.dt)))
        self.data = data2

        self.config = config
        self.psfun = config.psfun

        if self.config.imflag ==1:
            self.label = "Drift Time"
        else:
            self.label = "Charge"

        # Create the interface
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Peak Shape Tool')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot1 = PlottingWindow.Plot1d(pnl)
        sbs.Add(self.plot1, 1, wx.EXPAND)

        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
        hbox10 = wx.BoxSizer(wx.VERTICAL)

        centerbutton = wx.Button(pnl, label='Reset Range')
        self.flipbutton = wx.Button(pnl, label='Flip m/z and '+self.label)
        centerbutton.Bind(wx.EVT_BUTTON, self.on_reset)
        self.flipbutton.Bind(wx.EVT_BUTTON, self.on_flip)

        hboxsigs = wx.BoxSizer(wx.HORIZONTAL)
        self.outmzsig = wx.TextCtrl(pnl, value="")
        hboxsigs.Add(wx.StaticText(pnl, label='m/z Width Fit: '), wx.ALIGN_CENTER_VERTICAL)
        hboxsigs.Add(self.outmzsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        hboxsigs2 = wx.BoxSizer(wx.HORIZONTAL)
        self.outdtsig = wx.TextCtrl(pnl, value="")
        hboxsigs2.Add(wx.StaticText(pnl, label=self.label+' Width Fit: '), wx.ALIGN_CENTER_VERTICAL)
        hboxsigs2.Add(self.outdtsig, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        hbox10.Add(centerbutton, 0, wx.ALIGN_LEFT)
        hbox10.Add(self.flipbutton, 0, wx.ALIGN_LEFT)
        hbox10.Add(hboxsigs, 0, wx.ALIGN_LEFT)
        hbox10.Add(hboxsigs2, 0, wx.ALIGN_LEFT)
        hbox11.Add(hbox10, 0, wx.ALIGN_LEFT)

        hbox9 = wx.BoxSizer(wx.VERTICAL)
        self.ctlpsfun = wx.RadioBox(pnl, label="m/z Peak Shape Fit Function",
                                    choices=["Gaussian", "Lorentzian", "Split G/L"])
        self.ctlpsfun.SetSelection(self.psfun)
        fitbutton = wx.Button(pnl, label='Fit Peak Shape')
        plotbutton = wx.Button(pnl, label='Plot Guess')
        hbox9.Add(self.ctlpsfun, 0, wx.ALIGN_RIGHT)
        hbox8 = wx.BoxSizer(wx.HORIZONTAL)
        self.ctlsigguess = wx.TextCtrl(pnl, value=str(self.config.mzsig))
        hbox8.Add(wx.StaticText(pnl, label='Guess for Peak Width: '), wx.ALIGN_CENTER_VERTICAL)
        hbox8.Add(self.ctlsigguess, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        hbox9.Add(hbox8, 0, wx.ALIGN_RIGHT)
        hbox9.Add(plotbutton, 0, wx.ALIGN_RIGHT)
        hbox9.Add(fitbutton, 0, wx.ALIGN_RIGHT)
        hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.errorbox = wx.TextCtrl(pnl, value="", style=wx.TE_READONLY)
        hbox6.Add(wx.StaticText(pnl, label='Error in Fit: '), wx.ALIGN_CENTER_VERTICAL)
        hbox6.Add(self.errorbox, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        self.resbox = wx.TextCtrl(pnl, value=str(""))
        hbox7.Add(wx.StaticText(pnl, label='Resolution (M/FWHM): '), wx.ALIGN_CENTER_VERTICAL)
        hbox7.Add(self.resbox, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        hbox9.Add(hbox6, 0, wx.ALIGN_RIGHT)
        hbox9.Add(hbox7, 0, wx.ALIGN_RIGHT)
        hbox11.Add(hbox9, 1)
        sbs.Add(hbox11, 0, wx.EXPAND)
        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)

        plotbutton.Bind(wx.EVT_BUTTON, self.on_plot)
        fitbutton.Bind(wx.EVT_BUTTON, self.on_fit)
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        self.CenterOnParent()

        # Plot the result
        self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)

    def on_close(self, e):
        """
        Close the window and update self.config from the new values.
        :param e: Unused event
        :return: None
        """
        self.config.mzsig = simp_string_to_value(self.outmzsig.GetValue())
        self.config.dtsig = simp_string_to_value(self.outdtsig.GetValue())
        self.config.psfun = self.ctlpsfun.GetSelection()
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but do not update self.config.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)

    def on_flip(self, e):
        """
        Flip from m/z to arrival time at the previously fit m/z value.
        :param e: Unused event
        :return: None
        """
        if self.flipflag == 0:
            self.psfun = float(self.ctlpsfun.GetSelection())
            pos = np.argmin(np.abs(self.mz - self.fitmid))
            intdt = self.C[pos]
            self.data2 = np.column_stack((self.dt, intdt))
            self.data, self.data2 = self.data2, self.data
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1], title="Data", xlabel=self.label,
                                      ylabel="Normalized Intensity", zoom="span", zoomout=True)
            self.outmzsig.SetValue(str(self.fitsig))
            self.flipflag = 1
            if self.outdtsig.GetValue() != "":
                self.ctlsigguess.SetValue(self.outdtsig.GetValue())
            else:
                self.ctlsigguess.SetValue(str(self.config.dtsig))
            self.ctlpsfun.SetSelection(0)
        else:
            self.ctlpsfun.SetSelection(self.psfun)
            self.data, self.data2 = self.data2, self.data
            self.outdtsig.SetValue(str(self.fitsig))
            self.on_reset(0)
            self.flipflag = 0
            self.ctlsigguess.SetValue(self.outmzsig.GetValue())
        pass

    def on_reset(self, e):
        """
        Resets the plot.
        :param e: Unused event
        :return: None
        """
        if self.flipflag == 0:
            self.ctlsigguess.SetValue(str(self.config.mzsig))
        else:
            self.ctlsigguess.SetValue(str(self.config.dtsig))
        newmin = np.amin(self.data[:, 0])
        newmax = np.amax(self.data[:, 0])
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                  ylabel="Normalized Intensity", zoom="span", zoomout=True)

    def on_center(self, e):
        """
        Gets the range from the current plot zoom, cuts the data, and plots the cut data.
        :param e: Unused event
        :return: None
        """
        xlim = self.plot1.subplot1.get_xlim()
        newmin = xlim[0]
        newmax = xlim[1]
        boo1 = np.logical_and(self.data[:, 0] < newmax, self.data[:, 0] > newmin)
        self.centdat = self.data[boo1]
        if self.flipflag == 0:
            self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel="m/z (Th)",
                                      ylabel="Normalized Intensity", zoom="span", zoomout=True)
        else:
            self.plot1.plotrefreshtop(self.centdat[:, 0], self.centdat[:, 1], title="Data", xlabel=self.label,
                                      ylabel="Normalized Intensity", zoom="span", zoomout=True)

    def on_plot(self, e):
        """
        Plot a predicted peak shape based on the values input into the control.
        :param e: Unused event
        :return: None
        """
        self.on_center(e)
        sigguess = simp_string_to_value(self.ctlsigguess.GetValue())
        midguess = np.argmax(self.centdat[:, 1])
        aguess = np.amax(self.centdat[:, 1])
        if self.flipflag == 0:
            self.psfun = self.ctlpsfun.GetSelection()
            fitdat = make_peak_shape(self.centdat[:, 0], self.psfun, sigguess, self.centdat[midguess, 0])
        else:
            fitdat = make_peak_shape(self.centdat[:, 0], 0, sigguess, self.centdat[midguess, 0])

        self.plot1.plotadd(self.centdat[:, 0], fitdat * aguess, "blue", "Peak Shape Guess", nopaint=False)

    def on_fit(self, e):
        """
        Fits the spectrum and outputs the results.
        :param e:
        :return:
        """
        self.on_center(e)
        if self.flipflag == 0:
            psfun = self.ctlpsfun.GetSelection()
        else:
            psfun = 0

        # Fit the results
        fitout, fitdat = isolated_peak_fit(self.centdat[:, 0], self.centdat[:, 1], psfun)
        fitout = fitout[:, 0]
        resolution = self.centdat[np.argmax(self.centdat[:, 1]), 0] / fitout[0]
        error = np.sum((fitdat - self.centdat[:, 1]) * (fitdat - self.centdat[:, 1]))
        self.fitmid = fitout[1]
        self.fitsig = fitout[0]

        # Show the results
        print(self.fitmid, self.fitsig)
        if self.flipflag == 0:
            self.outmzsig.SetValue(str(fitout[0]))
        else:
            self.outdtsig.SetValue(str(fitout[0]))

        self.ctlsigguess.SetValue(str(fitout[0]))
        self.resbox.SetValue(str(resolution))

        self.plot1.plotadd(self.centdat[:, 0], fitdat, "blue", "Peak Shape Guess", nopaint=False)
        self.errorbox.SetValue(str(error))
        pass

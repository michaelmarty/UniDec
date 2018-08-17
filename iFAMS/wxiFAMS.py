import os
from copy import deepcopy
import numpy as np
import wx

from unidec_modules import unidecstructure, plot1d, plot2d
import unidec_modules.unidectools as ud
import matplotlib._cm as cm

from iFAMS import iFAMSfun


class iFAMS_Window(wx.Frame):
    def __init__(self, parent, config=None, directory=""):
        wx.Frame.__init__(self, parent, title="iFAMS")  # ,size=(-1,-1))

        self.parent = parent


        self.directory=directory

        # Set up the config file
        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config

        # Make the menu
        filemenu = wx.Menu()
        menu_load_data = filemenu.Append(wx.ID_ANY, "Load Spectrum", "Loads the Mass Spectrum from a .txt file")
        menu_save_fig_png = filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                            "Save all figures as PNG in central directory")
        menu_save_fig_pdf = filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                            "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_load_spectrum, menu_load_data)
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig_png)
        self.Bind(wx.EVT_MENU, self.on_save_fig_pdf, menu_save_fig_pdf)

        self.plotmenu = wx.Menu()

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the Plots
        panel = wx.Panel(self)

        self.plot1 = plot1d.Plot1d(panel)
        self.plot2 = plot1d.Plot1d(panel)
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
        controlsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer3 = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer4 = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer5 = wx.BoxSizer(wx.HORIZONTAL)

        # Set up the controls
        self.ctlm0 = wx.TextCtrl(panel, value=str(0))
        self.ctlwindow = wx.TextCtrl(panel, value=str(0))
        self.ctlwindow1 = wx.TextCtrl(panel, value=str(0))
        self.zerodata = wx.TextCtrl(panel, value=str(0))
        self.overtone = wx.TextCtrl(panel, value=str(0))
        self.harmavg = wx.TextCtrl(panel, value=str(0))
        self.lowcharge = wx.TextCtrl(panel, value=str(0))
        self.highcharge = wx.TextCtrl(panel, value=str(0))
        self.mansub = wx.TextCtrl(panel, value=str(0))
        controlsizer2.Add(wx.StaticText(panel, label="Minimum Frequency"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlm0, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(wx.StaticText(panel, label="Minimum Height"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlwindow, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(wx.StaticText(panel, label="Delta"), 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer2.Add(self.ctlwindow1, 0, wx.ALIGN_CENTER_VERTICAL)
        controlsizer4.Add(wx.StaticText(panel, label="# of zero frequency data"), 0, wx.ALIGN_RIGHT)
        controlsizer4.Add(self.zerodata, 0, wx.ALIGN_RIGHT)
        controlsizer4.Add(wx.StaticText(panel, label="# of harmonics for filter"), 0, wx.ALIGN_RIGHT)
        controlsizer4.Add(self.overtone, 0, wx.ALIGN_RIGHT)
        controlsizer3.Add(wx.StaticText(panel, label="# of harmonics for average"), 0, wx.ALIGN_RIGHT)
        controlsizer3.Add(self.harmavg, 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(wx.StaticText(panel, label="lowest charge state"), 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(self.lowcharge, 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(wx.StaticText(panel, label="highest charge state"), 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(self.highcharge, 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(wx.StaticText(panel, label="subunit mass"), 0, wx.ALIGN_RIGHT)
        controlsizer5.Add(self.mansub, 0, wx.ALIGN_RIGHT)




        recalcbutton = wx.Button(panel, label="Recalc. Maxima Finder")
        controlsizer2.Add(recalcbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_calc, recalcbutton)


        iFAMSbutton = wx.Button(panel, label="Run iFAMS analysis")
        controlsizer.Add(iFAMSbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_SubAndCharCalc, iFAMSbutton)

        Filbutton = wx.Button(panel, label="Fourier Filter")
        controlsizer4.Add(Filbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_fourier_fil, Filbutton)

        havgbutton = wx.Button(panel, label="Harmonic Average")
        controlsizer3.Add(havgbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_harm_average, havgbutton)

        self.realdatasel = wx.CheckBox(panel,label="plot real data")
        controlsizer.Add(self.realdatasel, 0, wx.EXPAND)
        self.Bind(wx.EVT_CHECKBOX, self.onChecked)

        mancalcbutton = wx.Button(panel, label="Man Calc Subunit and Charge")
        controlsizer5.Add(mancalcbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_man_calc, mancalcbutton)




        sizer.Add(controlsizer, 0, wx.EXPAND)
        sizer.Add(controlsizer2, 0, wx.EXPAND)
        sizer.Add(controlsizer3, 0, wx.EXPAND)
        sizer.Add(controlsizer4, 0, wx.EXPAND)
        sizer.Add(controlsizer5, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()



    def on_close(self, e):
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
    def on_load_spectrum(self, e=None):
        """
        Loads a spectrum from a .txt file
        :return: None
        """
        openFileDialog = wx.FileDialog(frame, "Open", "", "",
                                       "*.txt",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        openFileDialog.ShowModal()
        name = openFileDialog.GetPath()
        namestr = str(name)
        self.data=np.loadtxt(namestr)
        self.makeplot()


    def makeplot(self):
        """
        Runs the kendrick analysis on a single of the mass distributions in self.datalist set by self.pos.
        Plots the results in 1D and 2D.
        :return: None
        """
        self.xdata = self.data[:, 0]
        self.ydata = self.data[:, 1]
        self.getfromgui()

        self.plot1.plotrefreshtop(self.xdata, self.ydata, "Mass Spectrum", "m/z (Th)",
                                  "Intensity", "mass spectrum", color="b", config=self.config)
        self.plot1.add_legend(anchor=(1, 1))
        self.on_fft()

    def on_fft(self, e=None):
        self.getfromgui()
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y',
                  'k','b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']

        self.yfull, self.expandedspan, self.xnew, self.ftspacing, self.paddedxnew, self.ynew = iFAMSfun.plot_function(
            self.xdata, self.ydata)

        self.maxfreq = iFAMSfun.maxfreq(self, self.yfull, self.expandedspan)
        self.ftx, self.ABFT, self.FT = iFAMSfun.Fourier(self, self.maxfreq, self.yfull)
        self.refmaxtab = iFAMSfun.findmax(self.expandedspan, self.ABFT, self.ftx, 0.001, 5, 10)
        self.plot2.plotrefreshtop(self.ftx, self.ABFT, "FFT", "Frequency", "Amplitude", color='k', config=self.config)
        self.plot2.plotadddot(np.array(self.refmaxtab)[:, 0], np.array(self.refmaxtab)[:, 1],"r","o")
        self.plot2.repaint()

    def on_calc(self, e=None):
        """
        This module will help determine local maxima in your
        Fourier Spectrum.  You'll notice dots turn red to blue
        to indicate the change.  It takes numbers from self.clm0,
        self.ctlwindow, and self.ctlwindow1.
        :param e:
        :return: None
        """
        lowend = float(self.ctlm0.GetLineText(lineNo=0))
        pctpkht = float(self.ctlwindow.GetLineText(lineNo=0))
        delta = int(self.ctlwindow1.GetLineText(lineNo=0))
        self.refmaxtab = iFAMSfun.findmax(self.expandedspan, self.ABFT, self.ftx, lowend, delta, pctpkht)

        self.plot2.plotrefreshtop(self.ftx, self.ABFT, "FFT", "Frequency", "Amplitude", color='k', config=self.config)
        self.plot2.plotadddot(np.array(self.refmaxtab)[:, 0], np.array(self.refmaxtab)[:, 1],"b","o")
        self.plot2.repaint()



    def on_SubAndCharCalc(self, e=None):
        """
        This function will calculate the charge states and subunit mass based on the
        parameters above.  It will turn the dots green to indicate which maxima iFAMS
        uses to calculate the charge states and subunit mass
        :param e:
        :return:
        """

        self.newcalcX = []
        self.newcalcY = []
        self.refmaxtabCalc = np.array(self.refmaxtab)
        self.numchar = iFAMSfun.spacing(self, self.refmaxtab)
        self.omega = iFAMSfun.omega(self.refmaxtab, self.numchar)
        self.chargestates, self.chargestatesr = iFAMSfun.charge(self.refmaxtab, self.numchar, self.omega)
        self.chargestateints = [int(self.chargestatesr[i]) for i in range(0, len(self.chargestatesr))]
        for i in range(0, len(self.chargestatesr)):
            self.newcalcX.append(self.refmaxtabCalc[i, 0])
            self.newcalcY.append(self.refmaxtabCalc[i, 1])
        self.submass, self.stdevmass  = iFAMSfun.subunit(self.refmaxtab,self.numchar,self.omega,self.chargestates,
                                                         self.chargestatesr,self.ftx,self.ABFT)
        self.submassr = round(self.submass,2)

        newtext = "the subunit mass is " + str(self.submassr) + " +/- " \
                  + str(round(self.stdevmass,2)) + "\n with charge states" + str(self.chargestatesr)
        self.plot2.plotrefreshtop(self.ftx, self.ABFT, "FFT", "Frequency", "Amplitude",label=newtext, color='k', config=self.config)
        self.plot2.plotadddot(np.array(self.newcalcX), np.array(self.newcalcY),"g","o")
        self.plot2.addtext(newtext,(max(self.ftx)/2),max(self.ABFT),vlines=False)
        self.plot2.repaint()




        ############################### things you may need!!!!! #############################################

        print(self.chargestates)# This is a list with the calculated charge states
        print(self.chargestatesr)# This is the same list as above, only with the numbers rounded
        print(self.submass)# This is the calculated subunit mass

        ############################### things you may need!!!!! #############################################

        self.on_envelope_fun()

    def on_envelope_fun(self,e=None):
        """
        This function will inverse Fourier transform specific peaks in
        the Fourier domain.  It works by using the charge states and
        subunit mass above.  It will place a square window around each
        peak, whose centroid is determine by the charge state, and
        width is determine by the subunit mass.  The window is
        placed, and the Fourier spectrum is inverse Fourier transformed
        thus giving you the contribution from the each charge state

        :param e:
        :return:
        """

        self.ABIFT = iFAMSfun.envelope_calc(self.chargestatesr,self.expandedspan,self.submass,self.ftx,self.ftspacing,
                                            self.FT,self.ydata)
        self.plot1.plotrefreshtop(self.xdata, self.ydata, "Data", "m/z (Th)",
                                  "Intensity",label= "mass spectrum", color="k", config=self.config)
        for i in range(0,len(self.chargestatesr)):
            self.plot1.plotadd(self.xnew, self.ABIFT[i][0:int(len(self.xnew))], self.colors[i],newlabel=str(self.chargestatesr[i]))
            self.plot1.repaint()
        self.plot1.add_legend(anchor=(1, 1))
        self.zero_charge()

    def zero_charge(self,e=None):
        """
        This will calculate a zero charge spectrum based on the charge states
        and envelope functions from the previous python function.
        :param e:
        :return:
        """
        self.xrange, self.yfinal, self.yrangespec = iFAMSfun.zerocharge(self.ABIFT, self.xnew, self.chargestatesr)
        self.plot3.plotrefreshtop(self.xrange,self.yfinal, "Zero-Charge", "mass (Da)",
                                  "Intensity", "zero charge", color="k", config=self.config)
        for i in range(0,len(self.yrangespec)):
            self.plot3.plotadd(self.xrange,self.yrangespec[i],self.colors[i],newlabel=str(self.chargestatesr[i]))
            self.plot3.repaint()
        self.plot3.add_legend(anchor=(1,1))
    def on_fourier_fil(self, e=None):

        """
        This will filter the spectrum.  The filter uses two inputs,
        a Zero frequency data, which is basically how many fundamental frequency
        lengths should the zero frequency go out.  So in other words, if the
        fundamental frequency of the the subunit mass is 0.5, (or,
        the 1+ frequency), then typing in "3" for the zero frequency data
        will use all of the data to frequency 1.5 as the zero frequency data.

        The other input is the harmonic number, which is simply how many
        harmonics you would like to use.  Generally, the more baseline resolved
        harmonics, the better.
        :param e:
        :return:
        """

        ZFreqData = int(self.zerodata.GetLineText(lineNo=0))
        OTnum = int(self.overtone.GetLineText(lineNo=0))
        zeropointlabel = "zero frequency data"
        Fourierfilterlab = "Fourier filtered data"
        baselinelabel = "subtracted zero freq data"
        self.reconstspec, self.reconstbaseline = iFAMSfun.FFTFilter(self.expandedspan,self.submass,ZFreqData,self.ftx,
                                                                  self.ftspacing,self.FT,self.chargestatesr,OTnum)
        self.plot4.plotrefreshtop(self.xdata,self.ydata,"Fourier Filter", "m/z (Th)",
                                  "Intensity", "", color="k", config=self.config)
        self.plot4.plotadd(self.xnew,self.reconstspec[0:int(len(self.xnew))],'r',newlabel=Fourierfilterlab)
        self.plot4.plotadd(self.xnew,self.reconstbaseline[0:int(len(self.xnew))],'g',newlabel=zeropointlabel)
        self.plot4.plotadd(self.xnew,self.ynew-self.reconstbaseline[0:int(len(self.xnew))],'b',newlabel=baselinelabel)
        self.plot4.add_legend(anchor=(1,1))
        self.plot4.repaint()



        ############################### things you may need!!!!! #############################################
        """ currently have the x and y part as separate lists for the fourier filtered spectrum
         Wouldn't be too hard to make this one list if you prefer"""
        self.filtspecY = self.reconstspec[0:int(len(self.xnew))]
        print(len(self.xnew)) # x data for the filtered spectrum
        print(len(self.filtspecY)) # y data for the filtered spectrum
        ############################### things you may need!!!!! #############################################



    def on_harm_average(self,e=None):
        """
        This function will calculate a "Harmonic Average" for the envelope functions.
        This means that iFAMS will calculate multiple envelope functions for each charge
        state (one at the fundamental, and one for every additional charge state specified)
        and "average" them.
        It is advised that if multiple sets of resolvable harmonics exist, then one
        should use the harmonic average over just one set.
        :param e:
        :return:
        """
        self.plot1.plotrefreshtop(self.xdata, self.ydata, "Data", "m/z (Th)",
                                  "Intensity", "", color="k", config=self.config)
        self.chargestateints = [int(self.chargestatesr[i]) for i in range(0, len(self.chargestatesr))]
        ov =  int(self.harmavg.GetLineText(lineNo=0))
        self.ABIFT = iFAMSfun.AverageHarmFun(self.chargestatesr, self.expandedspan,self.submass,self.ftx,self.ftspacing,
                                             self.FT, self.paddedxnew,self.xnew, self.ynew,self.ydata,ov)
        for i in range(0,len(self.chargestatesr)):
            self.plot1.plotadd(self.xnew, self.ABIFT[i][0:int(len(self.xnew))], self.colors[i])
            self.plot1.repaint()
        self.zero_charge()

    def onChecked(self, e=None):

        if self.realdatasel.GetValue():#here you check if it is true or not
            self.ABIFT = iFAMSfun.realdata(self.chargestatesr,self.expandedspan,self.submass,self.ftx,self.ftspacing,
                                                self.FT,self.ydata)
            self.plot1.plotrefreshtop(self.xdata, self.ydata, "Data", "m/z (Th)",
                                      "Intensity", "mass spectrum", color="k", config=self.config)
            for i in range(0,len(self.chargestatesr)):
                self.plot1.plotadd(self.xnew, self.ABIFT[i][0:int(len(self.xnew))],
                                   self.colors[i],newlabel=str(self.chargestatesr[i]))
                self.plot1.repaint()
            self.plot1.add_legend(anchor=(1,1))
            self.zero_charge()
        else: self.on_envelope_fun()

    def on_man_calc(self,e=None):
        """
        This function allows the user to manually put in the charge states
        and subunit mass, in case the maxima finder has trouble picking
        the right peaks
        :return:
        """
        lowcharge = int(self.lowcharge.GetLineText(lineNo=0))
        highcharge = int(self.highcharge.GetLineText(lineNo=0))
        self.submass = float(self.mansub.GetLineText(lineNo=0))
        self.chargestatesr = np.arange(lowcharge,highcharge+1)
        self.refmaxtab = iFAMSfun.inputmax(self.chargestatesr,self.submass,self.ABFT,self.ftx)
        self.on_SubAndCharCalc()




    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Figure1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure2.png")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "Figure3.png")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure4.png")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, "Figure1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure2.pdf")
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
            # print name2
        name1 = os.path.join(self.directory, "Figure3.pdf")
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name1)
            # print name1
        name2 = os.path.join(self.directory, "Figure4.pdf")
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name2)
            # print name2


# Main App Execution
if __name__ == "__main__":

    app = wx.App(False)
    frame = iFAMS_Window(None)
    app.MainLoop()

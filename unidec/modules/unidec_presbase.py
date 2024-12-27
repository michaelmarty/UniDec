import wx
from unidec.modules.isolated_packages import FileDialogs, biopolymer_calculator
from unidec.modules import peakwidthtools
import os
import numpy as np
import unidec.tools as ud
from unidec.modules import ManualSelectionWindow, AutocorrWindow, peakstructure
from unidec.modules import SubDiv, miscwindows
from unidec.metaunidec.mudstruct import MetaDataSet
import sys
import getopt


class UniDecPres(object):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, ignore_args=False, *args, **kwargs):
        self.timer = None
        self.wx_app = wx.App(redirect=False)
        self.eng = None
        self.view = None
        self.recent_files = []

        self.infile = None
        self.top_path = None
        opts = None
        if not ignore_args:
            try:
                opts, args = getopt.getopt(sys.argv[1:], "ucmdf:", ["file=", "unidec", "meta", "chrom", "ucd"])
            except getopt.GetoptError as e:
                print("Error in Argv. Likely unknown option: ", sys.argv, e)
                print("Known options: -u, -m, -c, -f, -d")

            # print("ARGS:", args)
            # print("KWARGS:", kwargs)
            # print("OPTS:", opts)
            if opts is not None:
                for opt, arg in opts:
                    if opt in ("-f", "--file"):
                        self.infile = arg
                        print("Opening File:", self.infile)
            if ud.isempty(opts):
                if len(args) > 0:
                    self.infile = args[0]
                    print("Opening File:", self.infile)
            pass

    def start(self):
        """
        Launch view and start MainLoop
        :return:
        """
        self.wx_app.SetTopWindow(self.view)
        self.wx_app.SetAppName("GUniDec")
        self.wx_app.SetVendorName("Michael Marty - University of Arizona")

        self.view.Show()
        self.wx_app.MainLoop()

    # def on_end_session(self):
    #    wx.CallAfter(self.quit_application, force=True)

    def quit_application(self):
        # self.wx_app.ProcessIdle()
        self.wx_app.ExitMainLoop()
        return True

    def import_config(self, file_name=None):
        """
        Import configuration from file to engine (if file_name is specified) and from engine to GUI.
        :param file_name: Path of file to import
        :return: None
        """
        # tstart = time.perf_counter()
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.read_hdf5(file_name)
            else:
                self.eng.config.config_import(file_name)
        self.view.import_config_to_gui()
        if self.eng.config.filetype == 1:
            self.eng.config.write_hdf5()
        self.eng.update_history()


    def export_config(self, file_name=None):
        """
        Get configuration from GUI and (if file_name is specified) write from engine to file_name
        :param file_name: Path of file to save config to
        :return: None
        """
        self.view.export_gui_to_config()
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.write_hdf5(file_name)
                pass
            else:
                self.eng.config.config_export(file_name)
        else:
            if self.eng.config.filetype == 1:
                self.eng.config.write_hdf5()
                pass
        self.eng.update_history()

    def check_badness(self):
        """
        Check for any bad parameters and warn the user if something is off.
        :return: code (see unidecstructure.check_badness)
        """
        badness, warning = self.eng.config.check_badness()
        if warning != "":
            self.warn(warning)
        self.eng.check_isomode()
        return badness

    def warn(self, message, caption='Warning!'):
        """
        Send the user a message box.
        :param message: Message string
        :param caption: Caption string
        :return: None
        """
        dlg = wx.MessageDialog(self.view, message, caption, wx.OK | wx.ICON_WARNING)
        dlg.ShowModal()
        dlg.Destroy()

    def on_dataprep_button(self, e=None):
        print("Empty Function")

    def on_open_file(self, file_name, file_directory, time_range=None, refresh=False):
        print("Empty Function")

    def on_unidec_button(self, e=None):
        print("Empty Function")

    def on_integrate(self, e=None):
        print("Empty Function")

    def on_replot(self, e=None, plotsums=True):
        print("Empty Function")

    def peak_plots(self, e=None):
        print("Empty Function")

    def on_mass_tools(self, e=None, show=True):
        print("Empty Function")

    def on_get_mzlimits(self):
        try:
            if not wx.GetKeyState(wx.WXK_CONTROL):
                limits = self.view.plot1.subplot1.get_xlim()
                self.view.controls.ctlminmz.SetValue(str(limits[0]))
                self.view.controls.ctlmaxmz.SetValue(str(limits[1]))
                print("New m/z limits:", limits)
                try:
                    self.on_dataprep_button()
                except Exception:
                    pass
        except Exception as e:
            print(e)

    def on_load_conf_file(self, e=None):
        """
        Opens a file dialog and then imports a new _conf.dat config file
        :param e: unused space for event
        :return: None
        """
        cfilename = FileDialogs.open_file_dialog("Open Configuration File (_conf.dat)", file_types="*.*")
        if cfilename is not None:
            self.import_config(cfilename)
        pass

    def on_save_default(self, e=None):
        """
        Saves the default config file to self.eng.config.defaultconfig
        :param e: unused space for event
        :return: None
        """
        print("Saved: ", self.eng.config.defaultconfig)
        self.export_config(self.eng.config.defaultconfig)
        pass

    def on_load_default(self, e=None):
        """
        Resets the configuration then loads the default config file from self.eng.config.defaultconfig.
        Ignores min and max data cutoffs.
        :param e: unused space for event
        :return: None
        """
        self.on_reset(e)
        if os.path.isfile(self.eng.config.defaultconfig):
            try:
                self.import_config(self.eng.config.defaultconfig)
                try:
                    if not ud.isempty(self.eng.data.rawdata):
                        self.view.controls.ctlminmz.SetValue(str(np.amin(self.eng.data.rawdata[:, 0])))
                        self.view.controls.ctlmaxmz.SetValue(str(np.amax(self.eng.data.rawdata[:, 0])))
                        if self.eng.config.imflag == 1:
                            self.view.controls.ctlmindt.SetValue(str(np.amin(self.eng.data.rawdata3[:, 1])))
                            self.view.controls.ctlmaxdt.SetValue(str(np.amax(self.eng.data.rawdata3[:, 1])))
                except:
                    pass
                print("Loaded: ", self.eng.config.defaultconfig)
            except (ValueError, IndexError, TypeError):
                print("Failed to Load: ", self.eng.config.defaultconfig)
        self.view.SetStatusText("Loaded Default", number=5)
        pass

    def on_reset(self, e=None):
        """
        Resets UniDecConfig to default and then loads to self.view.
        :param e: unused space for event
        :return:
        """
        self.eng.reset_config()
        self.import_config(None)
        self.view.SetStatusText("Reset", number=5)

    def makeplot5(self, e=None, plot=None):
        """
        Plot mass vs charge data in 2D in self.view.plot5
        :param e: unused event
        :param plot: Plot object. Default is None, which will set plot to self.view.plot5
        :return: None
        """
        if plot is None:
            plot = self.view.plot5
        self.eng.makeplot5(plot)

    def makeplot3(self, e=None, plot=None):
        """
        Plot m/z vs charge data in 2D in self.view.plot3
        :param e: unused event
        :param plot: Plot object. Default is None, which will set plot to self.view.plot3
        :return: None
        """
        if plot is None:
            plot = self.view.plot3
        self.eng.makeplot3(plot)

    def makeplot2(self, e=None, plot=None, data=None, pks=None):
        """
        Plot mass data and peaks if possible in self.view.plot2
        :param e: unused event
        :param plot: Plot object. Default is None, which will set plot to self.view.plot2
        :param data: Data to plot. Default is None, which will set data to self.eng.data.massdat
        :param pks: Peaks to plot. Default is None, which will set pks to self.eng.pks
        :return: None
        """
        if plot is None:
            plot = self.view.plot2
        if data is None:
            data = self.eng.data.massdat
        if pks is None:
            pks = self.eng.pks
        self.eng.makeplot2(plot, data, pks)

    def makeplot4(self, e=None, plot=None, data=None, pks=None):
        """
        Plots isolated peaks against the data in self.view.plot4.
        Will plot dots at peak positions.
        If possible, will plot full isolated spectra.
        :param e: unused event
        :param plot: Plot object. Default is None, which will set plot to self.view.plot4
        :return: None
        """
        if plot is None:
            plot = self.view.plot4
        if data is None:
            data = self.eng.data.data2
        if pks is None:
            pks = self.eng.pks
        self.eng.makeplot4(plot, data, pks)

    def on_gen_html_report(self, e=None):
        """
        Generates an HTML report of the current data and parameters.
        :param e: unused event
        :return: None
        """
        # plots = [[self.view.plot4, self.view.plot2]]
        self.eng.gen_html_report(plots=self.view.plots)
        pass

    def on_filter_peaks(self, e=None):
        defaultvalue = "40"
        try:
            defaultvalue = str(self.eng.fdrs[0, 1] * 100)
        except:
            pass
        dialog = miscwindows.SingleInputDialog(self.view)
        dialog.initialize_interface(title="Minimum DScore", message="Set Minimum DScore Value (%): ",
                                    defaultvalue=defaultvalue)
        dialog.ShowModal()

        try:
            minval = float(dialog.value) / 100.
        except:
            print("Error with Score Input:", dialog.value)
            minval = 0.4

        print("Using DScore Cutoff (%): ", minval * 100)
        self.eng.filter_peaks(minscore=minval)
        self.view.peakpanel.add_data(self.eng.pks, show="dscore")

    def plot_integral(self, limits, color=None, filled=True, plot=None, repaint=False):
        """
        Plot a filled integral under a peak.
        :param limits: Limits of integration
        :param color: Color of fill and plot
        :param filled: If True, will fill the area under the curve
        :param plot: Plot object. Default is None, which will set plot to self.view.plot2
        :param repaint: If True, will repaint the plot
        :return: None
        """
        if plot is None:
            plot = self.view.plot2
        integral, intdat = self.eng.integrate(limits=limits)
        if filled:
            plot.filledplot(intdat[:, 0], intdat[:, 1], color)

        # Get Line width from matplotlib plot
        try:
            linewidth = plot.get_linewidth()
            linewidth *= 1.75
        except:
            linewidth = 1.75

        plot.plotadd(intdat[:, 0], intdat[:, 1], color, linewidth=linewidth)

        if repaint:
            plot.repaint()

    def on_color_peaks(self, e=None, plot2=None, plot4=None, mass=None, peakpanel=None, pks=None, filled=True):
        """
        Color the peaks in plot2 and plot4. Will also plot the integral under the peak.
        :param e: Unused event
        :param plot2: Plot object. Default is None, which will set plot to self.view.plot2
        :param plot4: Plot object. Default is None, which will set plot to self.view.plot4
        :param mass: Mass list to color. Default is None, which will set mass to self.view.peakpanel.selection2. If set, it will override peakpanel.
        :param peakpanel: Peak panel object. Default is None, which will set peakpanel to self.view.peakpanel.
        :param pks: Peaks object. Default is None, which will set pks to self.eng.pks
        :param filled: Whether to fill the integral under the peak
        :return: None
        """
        if plot2 is None:
            plot2 = self.view.plot2
        if plot4 is None:
            try:
                plot4 = self.view.plot4
            except Exception:
                plot4 = None

        if mass is None:
            if peakpanel is None:
                peaksel = self.view.peakpanel.selection2
            else:
                peaksel = peakpanel.selection2
        else:
            peaksel = mass

        if pks is None:
            pks = self.eng.pks

        # Get Line width from matplotlib plot
        try:
            linewidth = plot4.get_linewidth()
            linewidth *= 1.75
        except:
            linewidth = 1.75

        # TODO: Take away black background on colored lines
        for p in self.eng.pks.peaks:
            if p.mass in peaksel:
                if not ud.isempty(p.intervalFWHM):
                    # limits = p.integralrange
                    limits = 2 * (np.array(p.intervalFWHM) - p.mass) + p.mass
                    # print(p.intervalFWHM, limits, p.mass)

                    color = p.color
                    self.plot_integral(limits, color=color, filled=filled, plot=plot2, repaint=False)
                    if plot4 is not None:
                        for i, z in enumerate(self.eng.data.ztab):
                            if p.mztab[i, 1] > self.eng.config.peakplotthresh * np.amax(p.mztab[:, 1]):
                                mzlimits = (np.array(limits) / float(z)) + self.eng.config.adductmass
                                boo1 = self.eng.data.data2[:, 0] < mzlimits[1]
                                boo2 = self.eng.data.data2[:, 0] > mzlimits[0]
                                intdat = self.eng.data.data2[np.all([boo1, boo2], axis=0)]
                                if filled:
                                    plot4.filledplot(intdat[:, 0], intdat[:, 1], color)

                                plot4.plotadd(intdat[:, 0], intdat[:, 1], color, linewidth=linewidth)
        plot2.repaint()
        if plot4 is not None:
            plot4.repaint()

    def on_charge_states(self, e=None, mass=None, plot=None, peakpanel=None, data=None):
        """
        Triggered by right click "plot charge states" on self.view.peakpanel.
        Plots a line with text listing the charge states of a specific peak.
        :param e: unused event
        :return: None
        """
        if plot is None:
            plot = self.view.plot4

        if self.eng.config.adductmass > 0:
            sign = "+"
        else:
            sign = "-"
        charges = np.arange(self.eng.config.startz, self.eng.config.endz + 1)
        if mass is None:
            if peakpanel is None:
                peaksel = self.view.peakpanel.selection2[0]
            else:
                peaksel = peakpanel.selection2[0]
        else:
            peaksel = mass
        peakpos = (peaksel + charges * self.eng.config.adductmass) / charges
        boo1 = np.all([peakpos < self.eng.config.maxmz, peakpos > self.eng.config.minmz], axis=0)
        peakpos = peakpos[boo1]
        charges = charges[boo1]
        index = 0
        plot.textremove()
        if data is None:
            data = self.eng.data.data2
        for i in charges:
            plot.addtext(sign + str(i), peakpos[index], np.amax(data[:, 1]) * 0.99)
            index += 1

    def on_differences(self, e=None, pks=None, plot=None, massdat=None):
        """
        Triggered by right click "Display Differences" on self.view.peakpanel.
        Plots a line with text listing the difference between each mass and a specific peak.
        Updates the peakpanel to show the differences.
        :param e: unused event
        :return: None
        """

        if pks is None:
            pks = self.eng.pks
        if plot is None:
            plot = self.view.plot2
        if massdat is None:
            massdat = self.eng.data.massdat

        b1 = pks.get_bool()

        pmasses = np.array([p.mass for p in pks.peaks])[b1]
        peakdiff = np.array([p.diff for p in pks.peaks])[b1]
        mval = np.amax(massdat[:, 1])
        # print peakdiff

        plot.textremove()
        for i, d in enumerate(peakdiff):
            if d != 0:
                label = ud.decimal_formatter(d, self.eng.config.massbins)
                plot.addtext(label, pmasses[i], mval * 0.99 - (i % 7) * 0.05 * mval, nopaint=True)
            else:
                plot.addtext("0", pmasses[i], mval * 0.99 - (i % 7) * 0.05 * mval, nopaint=True)
        plot.repaint()

    def on_label_masses(self, e=None, peakpanel=None, pks=None, plot=None, dataobj=None):
        """
        Triggered by right click "Label Masses" on self.view.peakpanel.
        Plots a line with text listing the mass of each specific peak.
        Updates the peakpanel to show the masses.
        :param e: unused event
        :return: None
        """
        if peakpanel is None:
            peakpanel = self.view.peakpanel
        if pks is None:
            pks = self.eng.pks
        if plot is None:
            plot = self.view.plot2
        if dataobj is None:
            dataobj = self.eng.data

        peaksel = peakpanel.selection2
        pmasses = np.array([p.mass for p in pks.peaks])
        pint = np.array([p.height for p in pks.peaks])
        mval = np.amax(dataobj.massdat[:, 1])

        if isinstance(dataobj, MetaDataSet):
            mval = (mval * 1.3 + self.eng.config.separation * dataobj.len)

        plot.textremove()
        for i, d in enumerate(pmasses):
            if d in peaksel:
                label = ud.decimal_formatter(d, self.eng.config.massbins)
                plot.addtext(label, pmasses[i], mval * 0.06 + pint[i], vlines=False, nopaint=True)
        plot.repaint()

    def on_label_integral(self, e=None, peakpanel=None, pks=None, plot=None, dataobj=None):
        """
        Triggered by right click "Label Masses" on self.view.peakpanel.
        Plots a line with text listing the mass of each specific peak.
        Updates the peakpanel to show the masses.
        :param e: unused event
        :return: None
        """

        try:
            self.on_integrate()
            print("Integrated Peaks")
            int_success = True
        except:
            print("Unable to integrate. Using Peak Heights")
            int_success = False
            pass
        if peakpanel is None:
            peakpanel = self.view.peakpanel
        if pks is None:
            pks = self.eng.pks
        if plot is None:
            plot = self.view.plot2
        if dataobj is None:
            dataobj = self.eng.data

        if int_success:
            pintegral = np.array([p.integral for p in pks.peaks])
        else:
            pintegral = np.array([p.height for p in pks.peaks])

        peaksel = peakpanel.selection2
        pmasses = np.array([p.mass for p in pks.peaks])
        if ud.isempty(peaksel):
            peaksel = pmasses
        pint = np.array([p.height for p in pks.peaks])
        mval = np.amax(dataobj.massdat[:, 1])

        if isinstance(dataobj, MetaDataSet):
            mval = (mval * 1.3 + self.eng.config.separation * dataobj.len)

        plot.textremove()
        for i, d in enumerate(pmasses):
            if d in peaksel:
                label = str(np.round(pintegral[i], 2))
                plot.addtext(label, pmasses[i], mval * 0.06 + pint[i], vlines=False)

    def on_label_names(self, e=None, peakpanel=None, pks=None, plot=None, dataobj=None):
        """
        Triggered by right click "Label Masses" on self.view.peakpanel.
        Plots a line with text listing the mass of each specific peak.
        Updates the peakpanel to show the masses.
        :param e: unused event
        :return: None
        """

        if peakpanel is None:
            peakpanel = self.view.peakpanel
        if pks is None:
            pks = self.eng.pks
        if plot is None:
            plot = self.view.plot2
        if dataobj is None:
            dataobj = self.eng.data

        pnames = np.array([p.label for p in pks.peaks])

        peaksel = peakpanel.selection2
        pmasses = np.array([p.mass for p in pks.peaks])
        if ud.isempty(peaksel):
            peaksel = pmasses
        pint = np.array([p.height for p in pks.peaks])
        mval = np.amax(dataobj.massdat[:, 1])

        if isinstance(dataobj, MetaDataSet):
            mval = (mval * 1.3 + self.eng.config.separation * dataobj.len)

        plot.textremove()
        for i, d in enumerate(pmasses):
            if d in peaksel:
                plot.addtext(pnames[i], pmasses[i], mval * 0.06 + pint[i], vlines=False)

    def on_auto_peak_width(self, e=None, set=True):
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.eng.get_auto_peak_width(set_it=set)
            self.import_config()
        else:
            print("Need to process data first")

    def on_peak_width_tool(self, e=None):
        """
        Open peak width tool window. After it has returned, update the GUI to reflect the new peak widths.
        :param e: unused event
        :return: None
        """
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.export_config(None)
            if self.eng.config.imflag == 0:
                dlg = peakwidthtools.PeakTools1d(self.view)
                dlg.initialize_interface(self.eng.config, self.eng.data.data2)
            else:
                dlg = peakwidthtools.PeakTools2d(self.view)
                dlg.initialize_interface(self.eng.data.data3, self.eng.data.data2, self.eng.config)
            dlg.ShowModal()
            self.import_config(None)
        else:
            print("Need to process data first")
        pass

    def on_manual(self, e=None):
        """
        Opens window for setting manual assignments. Window directly modifies self.eng.config.
        :param e: unused event
        :return: None
        """
        dlg = ManualSelectionWindow.ManualSelection(self.view)
        if self.eng.config.imflag == 0:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data2)
        else:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data3)
        dlg.ShowModal()

    def on_smash_window(self, e=None):
        """
        Opens window for selecting peaks to smash. Window directly modifies self.eng.config.
        :param e: unused event
        :return: None
        """
        dlg = ManualSelectionWindow.SmashSelection(self.view)
        if self.eng.config.imflag == 0 and self.eng.config.cdmsflag == 0:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data2)
        elif self.eng.config.cdmsflag == 1:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data3)
        else:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data3)
        dlg.ShowModal()

    def on_match(self, e=None):
        """
        Automatic matching to present oligomer list.

        Opens masstools window but doesn't make it visible.
        Matches results to any combination of oligomers.
        Hits ok and plots results.

        :param e:
        :return:
        """
        # TODO: Rewrite this so that it doesn't need to fake open the window
        self.on_mass_tools(0, show=False)

    def on_autocorr_window(self, e=None):
        """
        Opens the autocorrelation window. Feed the config and massdat from the engine.
        :param e: Unused event
        :return: None
        """
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, self.eng.data.massdat)
        dlg.ShowModal()

    def on_biopolymer(self, e=0):
        bp = biopolymer_calculator.BiopolymerFrame(self.view)
        result = bp.Show()
        '''
        if result == 0:
            print("Calculated Mass from Sequence:", bp.mass)
            calc_mass = 0
            try:
                calc_mass = float(bp.mass)
                if calc_mass > 0:
                    label = bp.mass
                    mval = np.amax(self.eng.data.massdat[:, 1])
                    self.view.plot2.addtext(label, calc_mass, mval * 0.96, vlines=False)
                    pass
            except:
                pass
        bp.Destroy()'''

    def plot_theo_mass(self, e=0):
        dialog = miscwindows.SingleInputDialog(self.view)
        dialog.initialize_interface(title="Theoretical Mass", message="Set Theoretical Mass to Plot (Da): ",
                                    defaultvalue="")
        dialog.ShowModal()

        try:
            mass = float(dialog.value)
        except:
            print("Error with theoretical mass input:", dialog.value)
            mass = 0

        if self.eng.config.massbins >= 1:
            label = str(int(round(mass)))
        else:
            label = str(mass)
        mval = np.amax(self.eng.data.massdat[:, 1])
        self.view.plot2.addtext(label, mass, mval * 0.96, vlines=True)
        self.on_charge_states(mass=mass)

    # def register(self, e=None):
    #     from unidec.modules.data_reader import register
    #     register()

    def on_linreg(self, e=0):
        fit, rsquared = self.eng.linear_regression_peaks()
        message = "Slope: " + str(fit[0]) + "\nIntercept: " + str(fit[1]) + "\nR-Squared: " + str(rsquared)
        self.warn(message, caption="Linear Regression Results")

    def fpop(self, e=0):
        output = self.eng.oxidation_analysis()
        outstring = ""
        for o in output:
            outstring += str(o) + "\n"
        # Create text data object
        clipboard = wx.TextDataObject()

        # Set data object value
        clipboard.SetText(outstring)

        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()

    def sub_div(self, e=0):
        try:
            masses = []
            for p in self.eng.pks.peaks:
                if p.ignore == 0:
                    masses.append(p.mass)
            defaultsub = np.amin(masses)
        except:
            defaultsub = 44088

        defaultdiv = self.eng.config.molig
        dlg = miscwindows.DoubleInputDialog(self.view)
        dlg.initialize_interface("Subtract and Divide", "Subtract:", str(defaultsub),
                                 "Divide:", str(defaultdiv))
        dlg.ShowModal()
        try:
            sub = float(dlg.value)
            div = float(dlg.value2)
        except:
            print("Error with Subtract and Divide Inputs:", dlg.value, dlg.value2)
            return 0

        try:
            sd_results = []
            message = "Average Masses:\n\n"
            outstring = ""
            for s in self.eng.data.spectra:
                pks = peakstructure.Peaks()
                peaks = ud.peakdetect(s.massdat, self.eng.config)
                pks.add_peaks(peaks, massbins=self.eng.config.massbins)
                sd_result = ud.subtract_and_divide(pks, sub, div)
                sd_results.append(sd_result)

                outstring += str(sd_result) + "\n"

            avg = np.mean(sd_results)
            message += outstring
            message += "\nOverall Average: " + str(avg)
            outstring += "\n" + str(avg)
            self.copy_to_clipboard(outstring)
            self.warn(message, caption="Subtract and Divide Results")

        except:
            sd_result = ud.subtract_and_divide(self.eng.pks, sub, div)
            outstring = str(sd_result)
            message = "Average Mass: " + outstring
            self.copy_to_clipboard(outstring)
            self.warn(message, caption="Subtract and Divide Results")

    def sub_div_window(self, e=None):
        sd = SubDiv.SubDivFrame(self.view, self.eng.data.massdat, self.eng.pks, self.eng.config)

    def copy_to_clipboard(self, outstring):
        # Create text data object
        clipboard = wx.TextDataObject()

        # Set data object value
        clipboard.SetText(outstring)

        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()

    def write_to_recent(self, path=None):
        if path is None:
            path = os.path.join(self.eng.config.dirname, self.eng.config.filename)

        with open(self.eng.config.recentfile, 'a') as file:
            file.write(path + "\n")

        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)

    def read_recent(self):
        if os.path.isfile(self.eng.config.recentfile):
            with open(self.eng.config.recentfile, 'r') as file:
                lines = []
                for l in file:
                    p = l.strip("\n")
                    try:
                        if os.path.isfile(p):
                            lines.append(p)
                    except:
                        pass
        else:
            self.eng.config.recentfile = os.path.join(self.eng.config.UniDecDir, "recent.txt")
            with open(self.eng.config.recentfile, 'w') as file:
                pass  # Create empty file
            return []

        if not ud.isempty(lines):
            lines = np.array(lines)
            lines = lines[::-1][:10]
            unique, indexes = np.unique(lines, return_index=True)
            lines = np.array(unique[indexes.argsort()])
        else:
            lines = []
        return lines

    def cleanup_recent_file(self, recent_files):
        with open(self.eng.config.recentfile, "r+") as f:
            for p in recent_files[::-1]:
                l = p + '\n'
                f.write(l)
            f.truncate()

    def open_most_recent(self, e=None):
        recent_files = self.read_recent()
        if not ud.isempty(recent_files):
            directory, file = os.path.split(recent_files[0])
            self.on_open_file(file, directory, refresh=True)
        else:
            print("No Recent Files")

    def auto_refresh_stop(self, e=None):
        self.timer.Stop()

    def create_timer(self, e=None):
        self.timer = wx.Timer(self.view)
        self.view.Bind(wx.EVT_TIMER, self.on_timer)
        self.timer.Start(10000)

    def create_timer2(self, e=None):
        self.timer = wx.Timer(self.view)
        self.view.Bind(wx.EVT_TIMER, self.on_timer2)
        self.timer.Start(10000)

    def on_timer(self, e=None):
        print("Event")
        print("Refreshing: ", self.top_path)
        file_directory = os.path.dirname(self.top_path)
        file_name = os.path.basename(self.top_path)
        self.on_open_file(file_name, file_directory, refresh=True)
        self.quick_auto()

    def on_timer2(self, e=None):
        print("Event")
        print("Refreshing: ", self.top_path)
        file_directory = os.path.dirname(self.top_path)
        file_name = os.path.basename(self.top_path)
        maxtime = ud.get_max_time(self.top_path)
        print("Latest Time:", maxtime)

        mintime = maxtime - 1
        if mintime < 0:
            mintime = 0
        timerange = [mintime, maxtime]
        print(timerange)

        self.on_open_file(file_name, file_directory, time_range=timerange, refresh=True)
        self.quick_auto()

    def quick_auto(self, e=None):
        self.on_dataprep_button()
        self.on_unidec_button()

    def on_ex(self, e=0, pos=1):
        print("Loading Example Data")
        # Load the example data from the event. If there is an error, grab the pos value and load that file.
        try:
            self.view.menu.on_example_data(e)
        except:
            self.view.menu.load_example_data(pos)

        # If you hold down control, it will load everything
        if wx.GetKeyState(wx.WXK_CONTROL):
            self.on_load_everything()
        pass

    def on_load_everything(self, e=None):
        pass

    def on_delete(self, evt=None):
        pass

    def on_full(self, evt=None):
        pass

    def on_plot_peaks(self, evt=None):
        pass

    def on_pick_peaks(self, evt=None):
        pass

import time
import os
import thread
import wx
import numpy as np
import unidec
from wx.lib.pubsub import setupkwargs
# from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub

import unidec_modules.unidectools as ud
import unidec_modules.IM_functions as IM_func
from unidec_modules import Extract2D, peakwidthtools, masstools, miscwindows, \
    MassDefects, mainwindow, nativez, ManualSelectionWindow, AutocorrWindow, fft_window, GridDecon
from unidec_modules.isolated_packages import FileDialogs, texmaker, twitter_interface
import datacollector
import import_wizard
import unidec_modules.IM_windows as IM_wind
# import UniMin
from copy import deepcopy
import wx.py as py
import platform
import multiprocessing

# import FileDialog  # Needed for pyinstaller

__author__ = 'Michael.Marty'


# noinspection PyUnusedLocal
class Shell(object):
    def __init__(self, *args, **kwargs):
        self.__wx_app = wx.App(redirect=True)

        self.shell = py.shell.Shell(wx.Frame(None))

        self.shellwindow = py.shell.ShellFrame(self.shell, title="UniDecShell").Show()

        self.shell.Execute('app=UniDecApp()')
        # self.shell.Execute('app.start()')
        # self.shellwindow.Center()
        # self.shell.setFocus()
        self.__wx_app.MainLoop()


# noinspection,PyBroadException,PyUnusedLocal,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyUnusedLocal,PyBroadException
class UniDecApp(object):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize App
        :param args:
        :param kwargs:
        :return: UniDecApp object
        """
        self.__wx_app = wx.App(redirect=False)
        self.eng = None
        self.view = None
        self.twittercodes = None
        self.init(*args, **kwargs)

    def start(self):
        """
        Launch view and start MainLoop
        :return:
        """
        self.view.Show()
        self.__wx_app.MainLoop()

    def on_end_session(self):
        wx.CallAfter(self.quit_application, force=True)

    def quit_application(self):
        self.__wx_app.ProcessIdle()
        self.__wx_app.ExitMainLoop()
        return True

    def init(self, *args, **kwargs):
        """
        Initialize Engine and View. Load defaults.
        :param args:
        :param kwargs:
        :return:
        """
        self.eng = unidec.UniDec()
        self.twittercodes = None

        self.view = mainwindow.Mainwindow(self, "UniDec", self.eng.config)
        self.__wx_app.SetTopWindow(self.view)
        self.__wx_app.SetAppName("GUniDec")
        self.__wx_app.SetVendorName("Michael Marty - University of Oxford")

        pub.subscribe(self.on_integrate, 'integrate')
        pub.subscribe(self.on_smash, 'smash')
        pub.subscribe(self.on_get_mzlimits, 'mzlimits')
        pub.subscribe(self.on_left_click, 'left_click')

        self.on_load_default(0)

        if "path" in kwargs:
            newdir, fname = os.path.split(kwargs["path"])
            self.on_open_file(fname, newdir)
            # self.on_dataprep_button(0)
            self.on_auto(0)

        # For testing, load up a spectrum at startup. Used only on MTM's computer.
        if False and platform.node() == "RobMike":
            # fname = "HSPCID.txt"
            fname = "0.txt"
            # fname = "250313_AQPZ_POPC_100_imraw_input.dat"
            newdir = "C:\\cprog\\UniDecDemo"
            self.on_open_file(fname, newdir)
            # self.on_dataprep_button(0)
            self.on_auto(0)
            # self.on_integrate()
            # self.on_grid_decon(0)
            # self.make_cube_plot(0)
            # self.on_plot_peaks(0)

    # ..............................
    #
    #  Main Utility Functions
    #
    # ........................................

    def import_config(self, file_name=None):
        """
        Import configuration from file to engine (if file_name is specified) and from engine to GUI.
        :param file_name: Path of file to import
        :return: None
        """
        if file_name is not None:
            self.eng.config.config_import(file_name)
        self.view.import_config_to_gui()

    def export_config(self, file_name=None):
        """
        Get configuration from GUI and (if file_name is specified) write from engine to file_name
        :param file_name: Path of file to save config to
        :return: None
        """
        self.view.export_gui_to_config()
        if file_name is not None:
            self.eng.config.config_export(file_name)

    def on_open(self, e=None):
        """
        Open dialog for file opening
        :param e: unused space for event
        :return: None
        """
        dlg = wx.FileDialog(self.view, "Choose a data file in x y list format", '', "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.view.SetStatusText("Opening", number=5)
            self.eng.config.filename = dlg.GetFilename()
            print "Openening: ", self.eng.config.filename
            if os.path.splitext(self.eng.config.filename)[1] == ".zip":
                print "Can't open zip, try Load State."
                return
            self.eng.config.dirname = dlg.GetDirectory()
            self.on_open_file(self.eng.config.filename, self.eng.config.dirname)
        dlg.Destroy()

    def on_open_file(self, filename, directory, skipengine=False):
        """
        Opens a file. Run self.eng.open_file.
        :param filename: File name
        :param directory: Directory containing file
        :param skipengine: Boolean, Whether to skip running the engine (used when loading state)
        :return: None
        """
        # Clear other plots and panels
        self.view.peakpanel.clear_list()
        self.view.clear_all_plots()
        if not skipengine:
            # Open File in Engine
            self.eng.open_file(filename, directory)

        # Set Status Bar Text Values
        self.view.SetStatusText("File: " + filename, number=1)
        # print self.view.imflag, self.eng.config.imflag
        if self.view.imflag != self.eng.config.imflag:
            print "Changing Modes"
            self.on_flip_mode(0)
        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText(u"R\u00B2 ", number=3)
        # Update view with data limits
        if self.eng.config.batchflag != 1:
            self.view.ctlminmz.SetValue(str(np.amin(self.eng.data.data2[:, 0])))
            self.view.ctlmaxmz.SetValue(str(np.amax(self.eng.data.data2[:, 0])))
        # Plot 1D
        if self.eng.config.batchflag == 0:
            self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], "Data", "m/z",
                                           "Intensity", "Data", self.eng.config)
        # IM Loading and Plotting
        if self.eng.config.imflag == 1 and self.eng.config.batchflag != 1:
            self.view.ctlmindt.SetValue(str(np.amin(self.eng.data.data3[:, 1])))
            self.view.ctlmaxdt.SetValue(str(np.amax(self.eng.data.data3[:, 1])))
            if self.eng.config.batchflag == 0:
                self.view.plot1im.contourplot(self.eng.data.rawdata3, self.eng.config, xlab="m/z (Th)",
                                              ylab="Arrival Time (ms)", title="IM-MS Data")
        # Load Config to GUI
        self.import_config()
        self.view.SetStatusText("Ready", number=5)

    def on_save_state(self, e=None, filenew=None):
        """
        Saves the state by running self.eng.save_state. If no file is specified, opens a file dialog.
        :param e: unused space for event
        :param filenew: Save state zip file name
        :return: None
        """
        if filenew is None:
            filenew = FileDialogs.save_file_dialog(message="Save to UniDec Zip File", file_types="*.zip",
                                                   default_file=self.eng.config.outfname + ".zip")
        if filenew is not None:
            self.view.SetStatusText("Saving", number=5)
            self.eng.save_state(filenew)
            self.view.SetStatusText("Saved", number=5)
        pass

    def on_load_state(self, e=None):
        """
        Loads a saved file state by running self.eng.load_state and then updating self.view.
        Runs several engine components rather than importing.
        :param e: unused space for event
        :return: None
        """
        filenew = FileDialogs.open_file_dialog(message="Select UniDec Zip File to Open", file_types="*.zip")
        print filenew
        if filenew is not None:
            # Reset GUI
            self.on_reset(0)
            self.view.SetStatusText("Loading", number=5)
            tstart = time.clock()
            dirname, filename = os.path.split(filenew)
            self.view.SetStatusText("File: " + filename, number=1)
            self.view.SetStatusText(u"R\u00B2 ", number=3)
            self.view.SetStatusText("Data Length: ", number=2)
            # Load Into Engine and Presenter
            self.eng.load_state(filenew)
            self.on_open_file(filename, dirname, skipengine=True)

            # Load Processed Data
            if self.eng.config.procflag == 1:
                self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], "Data", "m/z",
                                               "Intensity", "Data", self.eng.config)
                if self.eng.config.imflag == 1:
                    self.view.plot1im.contourplot(self.eng.data.data3, self.eng.config, xlab="m/z (Th)",
                                                  ylab="Arrival Time (ms)", title="IM-MS Data")
                self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
            # Load UniDec Plots
            if os.path.isfile(self.eng.config.outfname + "_error.txt"):
                self.after_unidec_run()
            # RePick Peaks
            if os.path.isfile(self.eng.config.peaksfile):
                self.on_pick_peaks(e)
            # Match Again
            # TODO: Figure out what is happening with matchlist
            if os.path.isfile(self.eng.config.matchfile):
                self.eng.config.matchlist = np.transpose(
                    np.genfromtxt(self.eng.config.matchfile, dtype='str', delimiter=","))

            tend = time.clock()
            print "Loading Time: %.2gs" % (tend - tstart)
            self.view.SetStatusText("Ready", number=5)
        pass

    def on_raw_open(self, e=None, dirname=None):
        """
        Opens a file dialog for opening a raw file. Runs self.eng.raw_process and then opens the processed file.
        :param e: unused space for event
        :return:
        """
        if dirname is None:
            self.eng.config.dirname = FileDialogs.open_single_dir_dialog("Choose a raw file", '')
        else:
            self.eng.config.dirname = dirname

        if self.eng.config.imflag == 1:
            if int(self.view.ctlconvertflag.GetValue()) == 1:
                binsize = str(ud.string_to_value(self.view.ctlbinsize.GetValue()))
                print "Converting at resolution of: " + binsize
            else:
                binsize = "0"
                print "Converting using full resolution"
        else:
            binsize = None

        if self.eng.config.dirname is not None:
            self.view.SetStatusText("Converting", number=5)
            self.eng.config.dirname = os.path.abspath(self.eng.config.dirname)
            print "Loading Raw File: ", self.eng.config.dirname
            self.eng.config.filename, self.eng.config.dirname = self.eng.raw_process(self.eng.config.dirname, True,
                                                                                     binsize=binsize)
            if self.eng.config.filename is not None:
                self.on_open_file(self.eng.config.filename, self.eng.config.dirname)
        pass

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
        print "Saved: ", self.eng.config.defaultconfig
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
                if not ud.isempty(self.eng.data.rawdata):
                    self.view.ctlminmz.SetValue(str(np.amin(self.eng.data.rawdata[:, 0])))
                    self.view.ctlmaxmz.SetValue(str(np.amax(self.eng.data.rawdata[:, 0])))
                    if self.eng.config.imflag == 1:
                        self.view.ctlmindt.SetValue(str(np.amin(self.eng.data.rawdata3[:, 1])))
                        self.view.ctlmaxdt.SetValue(str(np.amax(self.eng.data.rawdata3[:, 1])))
                print "Loaded: ", self.eng.config.defaultconfig
            except (ValueError, IndexError, TypeError):
                print "Failed to Load: ", self.eng.config.defaultconfig
        self.view.SetStatusText("Loaded Default", number=5)
        pass

    def on_paste_spectrum(self, e=None):
        """
        Gets spectum from the clipboard, writes it to a new file, and then opens that new file.
        :param e: unused space for event
        :return: None
        """
        try:
            wx.TheClipboard.Open()
            do = wx.TextDataObject()
            wx.TheClipboard.GetData(do)
            wx.TheClipboard.Close()
            text = do.GetText()
            text = text.splitlines()
            data = []
            for t in text:
                line = t.split()
                if len(line) == 2:
                    try:
                        mz = float(line[0])
                        i = float(line[1])
                        data.append([mz, i])
                    except (ValueError, TypeError):
                        pass
            data = np.array(data)
            if len(data) > 0:
                if os.getcwd() == self.eng.config.UniDecDir:
                    topdir = os.path.dirname(os.getcwd())
                else:
                    topdir = os.path.dirname(os.path.dirname(os.getcwd()))
                newdir = os.path.join(topdir, "UniDecPastedSpectra")
                if not os.path.isdir(newdir):
                    os.mkdir(newdir)
                os.chdir(newdir)
                fname = "PastedSpectrum_" + str(time.strftime("%Y_%b_%d_%H_%M_%S")) + ".txt"
                np.savetxt(fname, data)
                print "Saved Pasted Spectrum as File:", fname, " in directory:", newdir
                self.on_open_file(fname, newdir)
            else:
                print "Paste failed, got: ", data
        except Exception, e:
            print e
            wx.MessageBox("Unable to open the clipboard", "Error")

    def on_reset(self, e=None):
        """
        Resets UniDecConfig to default and then loads to self.view.
        :param e: unused space for event
        :return:
        """
        self.eng.reset_config()
        self.import_config(None)
        self.view.SetStatusText("Reset", number=5)

    # ..........................
    #
    #  Core engine functions
    #
    # ..........................................

    def on_dataprep_button(self, e=None):
        """
        Prepares data. Runs eng.process_data() and then updates and plots in view
        :param e: unused space for event
        :return: None
        """
        tstart = time.clock()
        self.view.SetStatusText("Data Prep", number=5)
        self.export_config(self.eng.config.confname)
        self.eng.process_data()
        self.import_config()
        self.view.clear_all_plots()
        self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], "Data Sent to UniDec",
                                       "m/z (Th)", "Normalized Intensity", "Data", self.eng.config)
        if self.eng.config.intthresh != 0 and self.eng.config.imflag == 0:
            self.view.plot1.plotadd(self.eng.data.data2[:, 0],
                                    np.zeros_like(self.eng.data.data2[:, 1]) + self.eng.config.intthresh, "red",
                                    "Noise Threshold")
            self.view.plot1.add_legend()

        if self.eng.config.imflag == 1:
            self.view.plot1im.contourplot(self.eng.data.data3, self.eng.config, xlab="m/z (Th)",
                                          ylab="Arrival Time (ms)", title="IM-MS Data")
        self.view.plot1.repaint()

        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText(u"R\u00B2 ", number=3)
        self.view.SetStatusText("Data Prep Done", number=5)
        tend = time.clock()
        # print "Data Prep Done. Time: %.2gs" % (tend - tstart)
        pass

    def on_unidec_button(self, e=None):
        """
        Runs UniDec, imports the results, and plots the output.
        :param e: unused space for event
        :return: None
        """
        if self.eng.config.rawflag == 2:
            self.eng.config.rawflag = 0
        if self.eng.config.rawflag == 3:
            self.eng.config.rawflag = 1
        if self.eng.config.procflag == 0:
            # self.warn("Need to process data first.")
            self.on_dataprep_button()
        self.export_config(self.eng.config.confname)
        self.view.SetStatusText("UniDec Run", number=5)

        out = self.eng.run_unidec()
        if out == 0:
            self.after_unidec_run()
            self.view.SetStatusText("UniDec Done %.2gs" % self.eng.config.runtime, number=5)
        else:
            self.view.SetStatusText("Error %.0g" % out, number=5)
            print "Error ", out
            self.warn("Error %.0g" % out)
        pass

    def after_unidec_run(self):
        """
        Plots to make after running UniDec.
        :return: None
        """
        if self.eng.config.imflag == 0:
            self.view.SetStatusText("UniDec Plot", number=5)
            if self.view.system == "Linux":
                self.makeplot1(1)
                self.makeplot2(1)
                self.makeplot3(1)
                self.makeplot5(1)
            else:
                thread.start_new_thread(self.makeplot3, (1,))
                thread.start_new_thread(self.makeplot5, (1,))
                self.makeplot1(1)
                self.makeplot2(1)
        else:
            self.view.SetStatusText("UniDec Plot", number=5)
            self.make_im_plots()

        self.view.SetStatusText(u"R\u00B2: " + str(self.eng.config.error), number=3)

        self.view.plot4.clear_plot()
        self.view.plot6.clear_plot()
        self.view.peakpanel.clear_list()

    def on_pick_peaks(self, e=None):
        """
        Pick peaks and perform initial plots on them.
        :param e: unused space for event
        :return: None
        """
        self.view.SetStatusText("Detecting Peaks", number=5)
        self.export_config(self.eng.config.confname)
        self.eng.pick_peaks()
        self.view.SetStatusText("Plotting Peaks", number=5)

        if self.eng.config.batchflag == 0:
            self.view.peakpanel.add_data(self.eng.pks)
            self.makeplot2(1)
            self.makeplot6(1)
            self.makeplot4(1)

        self.view.SetStatusText("Peak Pick Done", number=5)
        print "Peak Picking"
        pass

    def on_plot_peaks(self, e=None):
        """
        Plots individual peaks as simulated spectra in plot 4.
        :param e: unused space for event
        :return: None
        """
        self.view.SetStatusText("Convolving", number=5)
        self.export_config(None)
        tstart = time.clock()
        self.eng.convolve_peaks()
        tend = time.clock()
        print "Convolving: %.2gs" % (tend - tstart)
        self.view.SetStatusText("Plotting", number=5)
        self.makeplot4(1)
        self.export_config(self.eng.config.confname)
        self.view.SetStatusText("Peak Plot Done", number=5)
        print "peak plotting"
        pass

    def on_peak_errors(self, e=None):
        """
        Get peak errors and plot on Plot2.
        :param e: unused space for event
        :return: None
        """
        print "Getting Errors"
        self.eng.get_errors()
        for p in self.eng.pks.peaks:

            masserr = p.masserr
            if masserr == 0:
                masserr = np.array([[np.amin(self.eng.data.massdat[:, 0]), np.amax(self.eng.data.massdat[:, 0])]])

            mass = p.massavg / self.view.plot2.kdnorm
            masserr /= self.view.plot2.kdnorm

            self.view.plot2.subplot1.errorbar(mass, p.corrint, xerr=masserr, yerr=p.correrr, color=p.color)

        self.view.plot2.repaint()

    def check_badness(self):
        """
        Check for any bad parameters and warn the user if something is off.
        :return: code (see unidecstructure.check_badness)
        """
        badness, warning = self.eng.config.check_badness()
        if warning is not "":
            self.warn(warning)
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

    def on_auto(self, e=None):
        """
        Prepare data, run UniDec, and pick peaks.
        :param e: unused space for event
        :return: None
        """
        self.on_dataprep_button(e)
        self.on_unidec_button(e)
        self.on_pick_peaks(e)

    # ..............................
    #
    #  Plotting Functions
    #
    # ...........................................

    def makeplot1(self, e=None):
        """
        Plot data and fit in self.view.plot1 and optionally in plot1fit
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.clock()
            if self.eng.config.imflag == 1:
                self.view.plot1fit.contourplot(self.eng.data.fitdat2d, self.eng.config, xlab="m/z (Th)",
                                               ylab="Arrival Time (ms)", title="IM-MS Fit")
            self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], "Data and UniDec Fit",
                                           "m/z (Th)", "Normalized Intensity", "Data", self.eng.config, nopaint=True)
            self.view.plot1.plotadd(self.eng.data.data2[:, 0], self.eng.data.fitdat, 'red', "Fit Data")
            self.view.plot1.add_legend()
            tend = time.clock()
            print "Plot 1: %.2gs" % (tend - tstart)

    def makeplot2(self, e=None):
        """
        Plot mass data and peaks if possible in self.view.plot2
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.clock()
            self.view.plot2.plotrefreshtop(self.eng.data.massdat[:, 0], self.eng.data.massdat[:, 1],
                                           "Zero-charge Mass Spectrum", "Mass (Da)",
                                           "Intensity", "Mass Distribution", self.eng.config, test_kda=True)
            if self.eng.pks.plen > 0:
                for p in self.eng.pks.peaks:
                    if p.ignore == 0:
                        self.view.plot2.plotadddot(p.mass, p.height, p.color, p.marker)
            self.view.plot2.repaint()
            tend = time.clock()
            print "Plot 2: %.2gs" % (tend - tstart)

    def makeplot3(self, e=None):
        """
        Plot m/z vs charge grid.
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.clock()
            self.view.plot3.contourplot(self.eng.data.mzgrid, self.eng.config)
            tend = time.clock()
            print "Plot 3: %.2gs" % (tend - tstart)

    def makeplot4(self, e=None):
        """
        Plots isolated peaks against the data in self.view.plot4.
        Will plot dots at peak positions.
        If possible, will plot full isolated spectra.
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.clock()
            self.view.plot4.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1],
                                           "Data with Offset Isolated Species", "m/z (Th)",
                                           "Normalized and Offset Intensity", "Data", self.eng.config, nopaint=True)
            num = 0
            if self.eng.config.isotopemode == 1:
                try:
                    stickmax = np.amax(np.array([p.stickdat for p in self.eng.pks.peaks]))
                except AttributeError:
                    stickmax = 1.0
            else:
                stickmax = 1.0
            for i in range(0, self.eng.pks.plen):
                p = self.eng.pks.peaks[i]
                if p.ignore == 0:
                    list1 = []
                    list2 = []
                    if (not ud.isempty(p.mztab)) and (not ud.isempty(p.mztab2)):
                        mztab = np.array(p.mztab)
                        mztab2 = np.array(p.mztab2)
                        maxval = np.amax(mztab[:, 1])
                        for k in range(0, len(mztab)):
                            if mztab[k, 1] > self.eng.config.peakplotthresh * maxval:
                                list1.append(mztab2[k, 0])
                                list2.append(mztab2[k, 1])
                                # print mztab[k]
                        self.view.plot4.plotadddot(np.array(list1), np.array(list2), p.color, p.marker)
                    if not ud.isempty(p.stickdat):
                        self.view.plot4.plotadd(self.eng.data.data2[:, 0], np.array(p.stickdat) / stickmax - (
                            num + 1) * self.eng.config.separation, p.color, "useless label")
                    num += 1
            self.view.plot4.repaint()
            tend = time.clock()
            print "Plot 4: %.2gs" % (tend - tstart)

    def makeplot5(self, e=None):
        """
        Plot mass vs charge grid
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.clock()
            self.view.plot5.contourplot(
                xvals=self.eng.data.massdat[:, 0], yvals=self.eng.data.ztab, zgrid=self.eng.data.massgrid,
                config=self.eng.config, title="Mass vs. Charge", test_kda=True)
            tend = time.clock()
            print "Plot 5: %.2gs" % (tend - tstart)

    def makeplot6(self, e=None, show="height"):
        """
        Plots bar chart of peak heights or areas in self.view.plot6.
        :param e: unused event
        :param show: What parameter to plot
        "height" will plot p.height for p in self.eng.pks.peaks
        "integral" will plot p.integral
        :return: None
        """
        if self.eng.config.batchflag == 0:
            if self.eng.pks.plen > 0:
                num = 0
                ints = []
                cols = []
                labs = []
                for p in self.eng.pks.peaks:
                    if p.ignore == 0:
                        num += 1
                        if show == "height":
                            ints.append(p.height)
                        elif show == "integral":
                            ints.append(p.integral)
                        cols.append(p.color)
                        labs.append(p.label)
                self.view.plot6.barplottop(range(0, num), ints, labs, cols, "Species", "Intensity", "Peak Intensities",
                                           repaint=False)
            for i in range(0, self.eng.pks.plen):
                p = self.eng.pks.peaks[i]
                if p.ignore == 0:
                    if show == "height":
                        y = p.height
                    elif show == "integral":
                        y = p.integral
                    else:
                        y = 0
                    self.view.plot6.plotadddot(i + 0.5, y, p.color, p.marker)
            self.view.plot6.repaint()

    def on_plot_composite(self, e=None):
        """
        Plots a composite of isolated peaks that are selected.
        :param e: unused event
        :return: None
        """
        try:
            self.eng.pks.composite = np.zeros(len(self.eng.data.data2))
            for i in xrange(0, self.eng.pks.plen):
                if self.eng.pks.peaks[i].ignore == 0:
                    self.eng.pks.composite += self.eng.pks.peaks[i].stickdat

            self.view.plot4.plotadd(self.eng.data.data2[:, 0], self.eng.pks.composite, "b", "useless label")
            self.view.plot4.repaint()
        except ValueError:
            print "Need to hit Plot Species button first"

    def make_im_plots(self):
        """
        Make Ion Mobility plots (but not cube plots)
        :return: None
        """
        if self.eng.config.batchflag == 0:
            self.makeplot1(1)
            self.makeplot2(1)
            self.makeplot5(1)
            self.makeplot3(1)
            self.view.plot2ccs.plotrefreshtop(self.eng.data.ccsdata[:, 0], self.eng.data.ccsdata[:, 1],
                                              title="CCS Distribution", xlabel="CCS (${\AA}$$^2$)", ylabel="Intensity",
                                              label="CCS Summation", config=self.eng.config, nopaint=False)
            self.view.plot5mccs.contourplot(xvals=self.eng.data.massdat[:, 0], yvals=self.eng.data.ccsdata[:, 0],
                                            zgrid=np.ravel(self.eng.data.massccs), config=self.eng.config,
                                            ylab="CCS (${\AA}$$^2$)", title="Mass vs. CCS", test_kda=True)

            ccsgrid2, zgrid2 = np.meshgrid(self.eng.data.ztab, self.eng.data.ccsdata[:, 0], sparse=False, indexing='ij')
            self.view.plot5ccsz.contourplot(
                np.transpose([np.ravel(ccsgrid2), np.ravel(zgrid2), np.ravel(self.eng.data.ccsz)]), self.eng.config,
                xlab="Charge", ylab="CCS (${\AA}$$^2$)", title="CCS vs. Charge")

            try:
                self.view.plot3color.make_color_plot(self.eng.data.mztgrid, np.unique(self.eng.data.data3[:, 0]),
                                                     np.unique(self.eng.data.data3[:, 1]), self.eng.data.ztab)
            except Exception, e:
                print "Color Plot Error", e

    def on_plot_nativeccs(self, e=None):
        """
        Plot native CCS as a red line on self.view.plot5mccs (mass vs ccs) plot.
        :param e: unused event
        :return: None
        """
        if not ud.isempty(self.eng.data.massdat) and self.eng.config.imflag == 1:
            ccses = [IM_func.calc_native_ccs(mass, self.eng.config.gasmass) for mass in self.eng.data.massdat[:, 0]]
            self.view.plot5mccs.subplot1.plot(self.eng.data.massdat[:, 0] / self.view.plot5mccs.kdnorm, ccses,
                                              color="r")
            self.view.plot5mccs.repaint()
            print "Plotted predicted native CCS values"

    def on_replot(self, e=None):
        """
        Refresh the parameters from the GUI and plot everything again with the new parameters.
        :param e: unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        if self.eng.config.imflag == 1:
            self.make_im_plots()
            self.makeplot4(1)
            self.makeplot6(1)
            if self.view.plot9.flag and self.view.plot10.flag:
                self.make_cube_plot(0)
        else:
            self.makeplot1(1)
            self.makeplot2(1)
            self.makeplot3(1)
            self.makeplot4(1)
            self.makeplot5(1)
            self.makeplot6(1)

    def make_cube_plot(self, e=None):
        """
        Make cube plots for IMMS.
        :param e: unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        try:
            starttime = time.clock()
            self.view.plot9.cubeplot(np.unique(self.eng.data.data3[:, 0]), np.unique(self.eng.data.data3[:, 1]),
                                     self.eng.data.ztab, np.sum(self.eng.data.mztgrid, axis=2),
                                     np.sum(self.eng.data.mztgrid, axis=1), np.sum(self.eng.data.mztgrid, axis=0),
                                     xlab="m/z (Th)", ylab="Arrival Time (ms)", zlab="Charge",
                                     cmap=self.eng.config.cmap)
            endtime = time.clock()
            print "Finished m/z Cube in: ", (endtime - starttime), " s"
        except Exception, ex:
            print "Failed m/z cube", ex
            pass
        try:
            starttime = time.clock()
            self.view.plot10.cubeplot(self.eng.data.massdat[:, 0], self.eng.data.ccsdata[:, 0], self.eng.data.ztab,
                                      self.eng.data.massccs, self.eng.data.massgrid.reshape(
                    (len(self.eng.data.massdat), len(self.eng.data.ztab))), self.eng.data.ccsz.transpose(),
                                      xlab="Mass (Da)",
                                      ylab="CCS (${\AA}$$^2$)", zlab="Charge", cmap=self.eng.config.cmap)
            endtime = time.clock()
            print "Finished Final Cube in: ", (endtime - starttime), " s"
        except Exception, ex:
            print "Failed final cube", ex
            pass

    def on_delete(self, e=None):
        """
        Triggered on deletion or modification of selection in self.view.peakpanel.
        Replots everything with peaks involved to reflect ignored peaks.
        :param e: unused event
        :return: None
        """
        self.view.SetStatusText("Changing Peaks", number=5)
        self.makeplot2(1)
        self.makeplot4(1)
        self.makeplot6(1)
        self.view.SetStatusText("Peak Change Done", number=5)

    def on_charge_states(self, e=None):
        """
        Triggered by right click "plot charge states" on self.view.peakpanel.
        Plots a line with text listing the charge states of a specific peak.
        :param e: unused event
        :return: None
        """
        charges = np.arange(self.eng.config.startz, self.eng.config.endz + 1)
        peaksel = self.view.peakpanel.selection2[0]
        peakpos = (peaksel + charges * self.eng.config.adductmass) / charges
        boo1 = np.all([peakpos < self.eng.config.maxmz, peakpos > self.eng.config.minmz], axis=0)
        peakpos = peakpos[boo1]
        charges = charges[boo1]
        index = 0
        self.view.plot4.textremove()
        for i in charges:
            self.view.plot4.addtext(str(i), peakpos[index], np.amax(self.eng.data.data2[:, 1]) * 0.99)
            index += 1

    def on_differences(self, e=None):
        """
        Triggered by right click "Display Differences" on self.view.peakpanel.
        Plots a line with text listing the difference between each mass and a specific peak.
        Updates the peakpanel to show the differences.
        :param e: unused event
        :return: None
        """
        peaksel = self.view.peakpanel.selection2
        pmasses = np.array([p.mass for p in self.eng.pks.peaks])
        peakdiff = pmasses - peaksel
        print peakdiff

        self.view.plot2.textremove()
        for i, d in enumerate(peakdiff):
            if d != 0:
                self.view.plot2.addtext(str(d), pmasses[i], np.amax(self.eng.data.massdat[:, 1]) * 0.99)
            else:
                self.view.plot2.addtext("0", pmasses[i], np.amax(self.eng.data.massdat[:, 1]) * 0.99)

    def on_plot_offsets(self, e=None):
        """
        Transform the mass vs charge grid to a mass vs charge offset grid. Plot in self.view.plot5.
        :param e: unused event
        :return: None
        """
        print "Tranforming charge to charge offset..."
        oaxis, outgrid = self.eng.mass_grid_to_f_grid()
        print "Plotting offsets"
        self.view.plot5.contourplot(xvals=self.eng.data.massdat[:, 0], yvals=oaxis, zgrid=outgrid,
                                    config=self.eng.config,
                                    ylab="Charge Offset", title="Mass vs. Charge Offset", test_kda=True)

    def on_color_plot1d(self, e=None, filled=False):
        self.on_integrate(plot=False)
        self.eng.get_peaks_scores()
        self.makeplot2(1)
        # TODO: Take away black background on colored lines
        for p in self.eng.pks.peaks:
            if not ud.isempty(p.integralrange):
                limits = p.integralrange
                color = p.color
                self.plot_integral(limits, color=color, filled=filled)
                self.view.plot2.repaint()

                for i, z in enumerate(self.eng.data.ztab):
                    if p.mztab[i, 1] > self.eng.config.peakplotthresh * np.amax(p.mztab[:, 1]):
                        mzlimits = np.array(limits) / float(z)
                        boo1 = self.eng.data.data2[:, 0] < mzlimits[1]
                        boo2 = self.eng.data.data2[:, 0] > mzlimits[0]
                        intdat = self.eng.data.data2[np.all([boo1, boo2], axis=0)]
                        if filled:
                            self.view.plot4.filledplot(intdat[:, 0], intdat[:, 1], color)
                        else:
                            self.view.plot4.plotadd(intdat[:, 0], intdat[:, 1], color)
                        self.view.plot4.repaint()

    def plot_integral(self, limits, color=None, filled=True):
        """
        Plot a filled integral under a peak.
        :param limits: Limits of integration
        :param color: Color of fill and plot
        :return: None
        """
        integral, intdat = self.eng.integrate(limits=limits)
        if filled:
            self.view.plot2.filledplot(intdat[:, 0], intdat[:, 1], color)
        else:
            self.view.plot2.plotadd(intdat[:, 0], intdat[:, 1], color)

    def on_integrate(self, plot=True, filled=True):
        """
        Triggered by right click on plot2. Integrates peaks.

        If plot2 is zoomed out, it will use self.eng.autointegrate() to integrate the peaks.
        If plot2 is zoomed in to a single peak, the integral for that peak is recalculated from the plot limits.
        If plot2 is zoomed in to more than one peak, the integral is not set for a single peak but simply printed on
            the plot

        :param plot: Boolean, whether to add filled areas to plot.
        :return: None
        """
        # Get Limits
        self.export_config(self.eng.config.confname)
        limits = self.view.plot2.subplot1.get_xlim()
        olimits = deepcopy(limits)

        limits = np.array(limits) * self.view.plot2.kdnorm

        # Run Integration
        if limits[0] <= np.amin(self.eng.data.massdat[:, 0]) and limits[1] >= np.amax(self.eng.data.massdat[:, 0]):
            print "Auto Integrating"
            self.eng.autointegrate()
        else:
            integral = self.eng.integrate(limits)
            if self.eng.pks.plen > 0:
                boo1 = self.eng.pks.masses < limits[1]
                boo2 = self.eng.pks.masses > limits[0]
                peaksinrange = self.eng.pks.masses[np.all([boo1, boo2], axis=0)]
            else:
                peaksinrange = []
            if len(peaksinrange) == 1:
                peak = peaksinrange[0]
                # print peak
                i = np.argmin((self.eng.pks.masses - peak) ** 2)
                self.eng.pks.peaks[i].integral = integral
                self.eng.pks.peaks[i].integralrange = limits
                print "Differences: ", limits - self.eng.pks.peaks[i].mass

            else:
                boo1 = self.eng.data.massdat[:, 0] < limits[1]
                boo2 = self.eng.data.massdat[:, 0] > limits[0]
                intdat = self.eng.data.massdat[np.all([boo1, boo2], axis=0)]
                self.view.plot2.addtext(str(integral), np.mean(np.array(olimits)),
                                        np.amax(intdat[:, 1]) + 0.05 * np.amax(self.eng.data.massdat[:, 1]),
                                        range=olimits)
                return 0

        # Normalize and write
        self.eng.normalize_peaks()
        areas = [[p.mass, p.integral] for p in self.eng.pks.peaks]
        np.savetxt(self.eng.config.outfname + "_peak_areas.dat", areas)

        # Plot
        if plot:
            self.makeplot2(1)
            for p in self.eng.pks.peaks:
                if not ud.isempty(p.integralrange):
                    limits = p.integralrange
                    color = p.color
                    self.plot_integral(limits, color=color, filled=filled)
            self.view.plot2.repaint()
            self.view.peakpanel.add_data(self.eng.pks, show="integral")
            try:
                self.makeplot6(1, show="integral")
            except Exception, ex:
                print "Didn't update bar chart", ex
                pass
        pass

    def on_smash(self):
        """
        Triggered by right click on plot1.
        Smashes the zoomed region to 0. Used to eliminate unwanted peaks.
        :return: None
        """
        print "Smashing!"
        self.export_config(None)
        limits = self.view.plot1.subplot1.get_xlim()
        bool1 = self.eng.data.data2[:, 0] > limits[0]
        bool2 = self.eng.data.data2[:, 0] < limits[1]
        bool3 = np.all(np.array([bool1, bool2]), axis=0)
        self.eng.data.data2[bool3, 1] = 0
        ud.dataexport(self.eng.data.data2, self.eng.config.infname)

        self.view.clear_all_plots()
        self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], "Data Sent to UniDec",
                                       "m/z (Th)",
                                       "Normalized Intensity", "Data", self.eng.config)
        if self.eng.config.intthresh != 0:
            self.view.plot1.plotadd(self.eng.data.data2[:, 0],
                                    np.zeros_like(self.eng.data.data2[:, 1]) + self.eng.config.intthresh, "red",
                                    "Noise Threshold")
            self.view.plot1.add_legend()
        self.view.plot1.repaint()

        self.view.SetStatusText(u"R\u00B2 ", number=3)
        message = "Smashing peaks from " + str(limits[0]) + " to " + str(limits[1]) + "\nReprocess data to undo"
        self.warn(message)
        pass

    def on_get_mzlimits(self):
        limits = self.view.plot1.subplot1.get_xlim()
        self.view.ctlminmz.SetValue(str(limits[0]))
        self.view.ctlmaxmz.SetValue(str(limits[1]))
        print "New m/z limits:", limits

    def on_charge_plot(self, e=None):
        """
        Change peaks and plots so that the total charge distribution is plotted in plot2.
        Picks new peaks as each charge state and updates peakpanel.
        :param e: unused event
        :return: None
        """
        self.export_config(None)
        cdat = self.eng.get_charge_peaks()
        self.view.peakpanel.add_data(self.eng.pks, collab1="Charge")

        self.view.plot2.plotrefreshtop(cdat[:, 0], cdat[:, 1], "Charge Distribution", "Charge", "Total Intensity",
                                       "useless label", self.eng.config)
        if self.eng.pks.plen > 0:
            for p in self.eng.pks.peaks:
                if p.ignore == 0:
                    self.view.plot2.plotadddot(p.mass, p.height, p.color, p.marker)
        self.view.plot2.repaint()
        self.makeplot4(0)
        self.makeplot6(0)
        cavg, cstd = ud.center_of_mass(cdat)
        print "Weighted average charge state:", cavg
        print "Weighted standard deviation:", cstd

    # ..................................
    #
    # Tools
    #
    # .................................................

    def on_batch_raw(self, e=None, dirs=None, clip=True):
        """
        Batch process Waters .Raw files into .txt files.
        Opens a multiple directory dialog and then sends the results to self.eng.raw_process
        :param e: unused event
        :return: None
        """
        if not (self.eng.config.system == "Windows"):
            print "Sorry. Waters Raw file converter only works on windows. Convert to txt file on a seperate machine."
            return None
            # self.eng.config.dirname="C:\\MassLynx\\"
            # if(not os.path.isdir(self.eng.config.dirname)):
        self.eng.config.dirname = ''
        if dirs is None:
            dirs = FileDialogs.open_multiple_dir_dialog("Select Raw Folders", self.eng.config.dirname)
        print dirs
        if dirs is not None:
            for d in dirs:
                if clip:
                    d = 'C:\\' + d[15:]
                self.eng.raw_process(d, False)

        print "Batch Converted"
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

    def on_mass_tools(self, e=None, show=True):
        """
        Opens masstools window.

        If a match was performed, it will update plot 6 and the peak panel.
        If the user decides to use simulated masses, it will make plot these new peaks.
        :param e: unused event
        :param show: Whether to thow the window (True) or simply match and return (False)
        :return: None
        """
        dlg = masstools.MassSelection(self.view)
        dlg.init_dialog(self.eng.config, self.eng.pks, massdat=self.eng.data.massdat)
        if show:
            result = dlg.ShowModal()
        else:
            result = 0
            dlg.on_match_all(0)
            dlg.on_close(0)
        # TODO: Rewrite so p.match isn't overwritten somehow if cancel is selected
        if self.eng.config.matchlist != [] and result == 0:
            if len(self.eng.config.matchlist[3]) == self.eng.pks.plen:
                self.view.SetStatusText("Matching", number=5)
                np.savetxt(self.eng.config.matchfile, np.transpose(self.eng.config.matchlist), fmt='%s', delimiter=",")
                self.view.peakpanel.add_data(self.eng.pks)
                self.makeplot6(1)
            else:
                self.eng.config.matchlist = []

        if self.eng.pks.changed == 1:
            print "Simulating Peaks"
            mztab = ud.make_peaks_mztab(self.eng.data.mzgrid, self.eng.pks, self.eng.config.adductmass)
            ud.make_peaks_mztab_spectrum(self.eng.data.mzgrid, self.eng.pks, self.eng.data.data2, mztab)
            self.view.peakpanel.add_data(self.eng.pks)
            self.makeplot2(1)
            self.makeplot6(1)
            self.makeplot4(1)
            self.eng.pks.changed = 0
            self.on_plot_peaks(0)
        self.export_config(self.eng.config.confname)
        self.view.SetStatusText("Match Done", number=5)
        pass

    def on_auto_peak_width(self, e=None):
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.eng.get_auto_peak_width()
            self.import_config()
        else:
            print "Need to process data first"

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
            print "Need to process data first"
        pass

    def on_additional_parameters(self, e=None):
        """
        Open additional data parameter window, which hads some of the experimental and obscure variables.
        Window directly modifies self.eng.config.
        Saves config to file after window is closed.
        :param e: unused event
        :return: None
        """
        dlg = miscwindows.AdditionalParameters(self.view)
        dlg.initialize_interface(self.eng.config)
        dlg.ShowModal()
        self.export_config(self.eng.config.confname)

    def on_autocorr_window(self, e=None):
        """
        Opens the autocorrelation window. Feed the config and massdat from the engine.
        :param e: Unused event
        :return: None
        """
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, self.eng.data.massdat)
        dlg.ShowModal()

    def on_unidec_path(self, e=None):
        """
        Opens file dialog to specify new UniDec binary.
        Only useful for advanced users or if something has gone horribly wrong.
        Modifies self.eng.config.UniDecPath and others.
        :param e: unused event
        :return: None
        """
        dlg = wx.FileDialog(self.view, "Locate the UniDec Executable", self.eng.config.UniDecDir, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.eng.config.UniDecName = dlg.GetFilename()
            self.eng.config.UniDecDir = dlg.GetDirectory()
            self.eng.config.UniDecPath = os.path.join(self.eng.config.UniDecDir, self.eng.config.UniDecName)
            print "New Path:", self.eng.config.UnIDecPath
        dlg.Destroy()

    def on_file_name(self, e=None):
        """
        Opens window with default file names in it.
        Possible to tweak default paths here, but probably dangerous.
        :param e: unused event
        :return: None
        """
        dlg = miscwindows.FileNameDialog(self.view)
        dlg.initialize_interface(self.eng.config)
        dlg.ShowModal()
        pass

    def on_data_collector(self, e=None):
        """
        Spawns separate DataCollector window.
        :param e: unused event
        :return: None
        """
        datacollector.DataCollector(None, "Data Collector", config=self.eng.config, pks=self.eng.pks,
                                    directory=self.eng.config.dirname)

    def on_import_wizard(self, e=None):
        """
        Opens Import Wizard for extreme file conversion.
        :param e: unused event
        :return: None
        """
        dlg = import_wizard.ImportWizard(self.view, dir=self.eng.config.UniDecDir)
        dlg.Show()

    def on_im_tools(self, e=None):
        """
        Open IM Parameters tool window for guessing at IM parameters.
        :param e: unused event
        :return: None
        """
        self.export_config()
        if not ud.isempty(self.eng.data.data3):
            self.export_config(None)
            if self.eng.config.imflag == 1:
                dlg = IM_wind.IMTools(self.view)
                dlg.initialize_interface(self.eng.data.data3, self.eng.config)
                dlg.ShowModal()
                self.eng.config = dlg.defaultconfig
            self.import_config(None)
        else:
            print "Load Data First"
        pass

    def on_im_extract(self, e=None):
        """
        Open IM extraction window for extracting data from IM plots.
        Has to run UniDec again with self.eng.config.zout=-1,
            which causes the export of all charge state slices of mass vs ccs plots.
            Normally only the faces of the mass vs ccs vs charge cube are output.
        :param e: unused event
        :return: None
        """
        if not ud.isempty(self.eng.data.ccsdata):
            print "Running UniDec to Generate Outputs"
            self.eng.config.zout = -1
            self.export_config(self.eng.config.confname)
            ud.unidec_call(self.eng.config.UniDecIMPath, self.eng.config.confname)
            dlg = IM_wind.IMToolExtract(self.view)
            dlg.initialize_interface(self.eng.data.massdat, self.eng.data.ccsdata, self.eng.data.massccs,
                                     self.eng.config,
                                     self.eng.pks)
            dlg.ShowModal()
            self.eng.config.zout = 0
        pass

    def on_tweet(self, e=None):
        """
        Opens Twitter Extension Window.

        First makes PNG files of all the figures it can. Those are fed to the window.

        self.twittercodes is modified if the person has logged in, so that the person only has to log in once.

        Note: don't mess with self.twittercodes. It DOES NOT contain log in information such as user name or password,
        but it shouldn't be printed, which is why it lives in the memory only.

        :param e: event
        :return: None
        """
        try:
            self.view.on_save_figure_png(e, transparent=False)
        except Exception, ex:
            print "Couldn't make figures for Twitter", ex
        os.environ['REQUESTS_CA_BUNDLE'] = os.path.join(self.eng.config.UniDecDir, 'cacert.pem')
        print "Will look for file: ", os.path.join(self.eng.config.UniDecDir, 'cacert.pem')
        tweetwindow = twitter_interface.TwitterWindow(self.view, pngs=self.view.pngs, codes=self.twittercodes,
                                                      imflag=self.eng.config.imflag)
        tweetwindow.ShowModal()
        self.twittercodes = tweetwindow.codes
        # print "Tweet"

    def on_kendrick(self, e=None):
        """
        Opens Kendrick Mass Defect Analysis window.
        :param e: unused event
        :return: None
        """
        MassDefects.MassDefectWindow(self.view, [self.eng.data.massdat], config=self.eng.config,
                                     pks=self.eng.pks, value=self.eng.config.molig)

    def on_2d_grid(self, e=None):
        """
        Opens 2D grid extraction window. Grid extraction parameters are stored at self.eng.config.gridparams.
        :param e: unused event
        :return: None
        """
        extract2d = Extract2D.Extract2DPlot(self.view, [self.eng.data.massdat], config=self.eng.config,
                                            params=self.eng.config.gridparams)
        self.eng.config.gridparams = extract2d.params
        pass

    def on_nativez_tools(self, e=None):
        """
        Opens Native Charge State analysis tools window. Parameters are stored in self.eng.config.zoffs.
        :param e: unused event
        :return: None
        """
        dlg = nativez.NativeZ(self.view)
        dlg.initialize_interface(self.eng.data.massdat[:, 0], np.unique(self.eng.data.mzgrid[:, 1]),
                                 self.eng.data.massgrid,
                                 self.eng.config, self.eng.pks)
        dlg.ShowModal()

    def on_export_params(self, e=None):
        """
        Runs self.eng.export_params(), which gets critical peak parameters and writes them to a file.
        :param e: event or arguments passed to self.eng.export_params()
        :return: None
        """
        self.eng.export_params(e)

    def on_mass_process(self, e=None):
        """
        Updates config from the GUI, processes the zero-charge mass spectrum with the same parameters used for the
        raw data, then replots the zero-charge mass spectrum.
        :param e: unused event
        :return: None
        """
        self.export_config(None)
        self.eng.process_mass_data()
        self.makeplot2(0)

    def on_calibrate(self, e=None):
        # Launch window to input calibration parameters
        dialog = miscwindows.SingleInputDialog(self.view)
        dialog.initialize_interface(title="Calibration Parameters",
                                    message="Polynomial Coefficents (order=n n-1 ... 0): ")
        dialog.ShowModal()

        try:
            result = dialog.value
            if result != "None":
                coeff = np.array(result.split())
                coeff = coeff.astype(np.float)
            else:
                coeff = None
            self.eng.data.rawdata = ud.cal_data(self.eng.data.rawdata, coeff)
            print "Calibration Success! Polynomial Coefficients (order = n to 0):", coeff
        except Exception, e:
            print "Calibration failed:", e

    def on_center_of_mass(self, e=None):
        """
        Determines the center of mass for the region zoomed on plot2.

        Get limits from plot2.
        Get Center of Mass from self.eng.center_of_mass
        Reports these as text on plot2.
        :param e:
        :return:
        """
        if not ud.isempty(self.eng.data.massdat):
            limits = self.view.plot2.subplot1.get_xlim()
            limits = np.array(limits) * self.view.plot2.kdnorm
            # print "limits", limits
            com, std = self.eng.center_of_mass(limits=limits)
            print "Center of Mass over region", limits, ":", com
            self.view.plot2.addtext(str(com), com, 0.99 * np.amax(self.eng.data.massdat[:, 1]))
            self.view.plot2.repaint()
        else:
            print "Need to get zero-charge mass spectrum first."

    def on_zerocharge_mass(self, e=None):
        """
        The old switcheroo. Takes the mass output and makes it the m/z input.
        Then changes some of the deconvolution parameters to reflect this.

        The reason for this type of function is to separate deconvolution from deisotoping when isotope mode is on.
        Deconvolution can be performed to generate a zero-charge mass spectrum, which can then be loaded as a new
        spectrum and deisotoped with isotope mode. Still a bit experimental...
        :param e: unused event
        :return: None
        """
        self.export_config()
        massoutputfile = self.eng.config.outfname + "_mass.txt"
        self.eng.config.default_zero_charge()
        self.eng.config.minmz = np.amin(self.eng.data.massdat[:, 0])
        self.eng.config.maxmz = np.amax(self.eng.data.massdat[:, 0])
        self.import_config(None)
        self.on_open_file(massoutputfile, os.getcwd())
        pass

    def on_fit_masses(self, e=None):
        """
        Using the peaks as a starting guess, it will fit the zero-charge mass spectrum to a set of overlapping peaks.

        Gets guesses using self.on_export_params()
        Calls self.eng.fit_all_masses()
        Then, plots resulting fit.
        Updates peakpanel to reflect normalized peak area based on fits.
        Export parameters from fit using self.on_export_params() and sending the "PostFit" argument
        :param e: event or argument passed to self.on_export_params
        :return: None
        """
        self.on_export_params(e)
        print "Fitting Zero-charge Mass Spectra to Peak Shapes"
        massfitdat, massfit = self.eng.fit_all_masses()
        print "Fit: ", massfit
        self.makeplot2(1)
        self.view.plot2.plotadd(self.eng.data.massdat[:, 0],
                                massfitdat / np.amax(massfitdat) * np.amax(self.eng.data.massdat[:, 1]), "green",
                                "Minimization")
        self.view.plot2.repaint()
        for i in xrange(0, self.eng.pks.plen):
            p = self.eng.pks.peaks[i]
            p.area = "%.2f" % float(massfit[i, 2] / np.sum(massfit[:, 2]))
        self.view.peakpanel.add_data(self.eng.pks)
        self.on_export_params("PostFit")

    def on_batch(self, e=None, flag=0, batchfiles=None):
        """
        Batch processing!

        Spawns a directory to select multiple files.
        Feeds the file names into a loop that:

        Opens the file (self.on_open_file)
        Prepares the data (self.on_dataprep_button)
        Runs UniDec (self.unidec_button)
        Picks Peaks (self.on_pick_peaks)
        Exports Peak Parameters (self.on_export_params)
        Saves State (self.on_save_state)

        If uses self.eng.config.batchflag to prevent certain things from plotting and all key parameters from changing.
        If batchflag is 2 (flag=1),
            all key paramters are kept, but the data ranges are refreshed from the individual files.
        :param e: event passed to some function (unused)
        :param flag: flag added to self.eng.config.batchflag
        :param batchfiles: List of files to run in batch. If None, will open dialog.
        :return:
        """
        if batchfiles is None:
            batchfiles = FileDialogs.open_multiple_files_dialog(
                message="Select Files To Process With Current Parameters",
                file_type="Text (.txt)|*.txt|Any Type|*.*")

        self.eng.config.batchflag = 1 + flag
        tstarttop = time.clock()
        if batchfiles is not None:
            self.view.clear_all_plots()
            for i, path in enumerate(batchfiles):
                print path
                tstart = time.clock()
                dirname, filename = os.path.split(path)
                self.on_open_file(filename, dirname)
                self.on_dataprep_button(e)
                self.on_unidec_button(e)
                self.on_pick_peaks(e)
                self.on_export_params(e)
                outfile = os.path.join(dirname, self.eng.config.outfname + ".zip")
                self.on_save_state(0, outfile)
                tend = time.clock()
                print "File saved to: " + str(outfile)
                print "Run Time: %.2gs" % (tend - tstart)
                print "\n"
        self.eng.config.batchflag = 0
        tend = time.clock()
        print "\nTotal Batch Run Time: %.2gs" % (tend - tstarttop)

    def on_batch2(self, e=None):
        """
        Runs batch processing without restricting the data ranges.

        The date range can be set ahead of time and will update independently for each file.
        Other parameters are fixed for the whole batch.
        :param e: unused event passed to self.on_batch
        :return: None
        """
        self.on_batch(e, flag=1)

    def on_super_batch(self, e=None):
        """
        Speedy, minimal batch procssing. Data must be processed in advance.

        Does not load file or process data.
        Simply writes the configuration with updated defaults to a conf.dat file.
        Then runs the conf.dat file with UniDec binary.
        Does not do any post processing.

        This function is useful for largee data sets where you are using DataCollector or some other program to
        analyze the results outside of UniDec. It allows new parameters to be run with the maximum speed.

        :param e: unused event
        :return: None
        """
        batchfiles = FileDialogs.open_multiple_files_dialog(message="Select Files To Process With Current Parameters",
                                                            file_type="Text (.txt)|*.txt|Any Type|*.*")
        if batchfiles is not None:
            self.view.SetStatusText("Speedy Batch Run", number=5)
            tstarttop = time.clock()
            for i, path in enumerate(batchfiles):
                print path
                tstart = time.clock()
                # Open File Stripped
                self.eng.config.dirname, self.eng.config.filename = os.path.split(path)
                self.view.SetStatusText("File: " + self.eng.config.filename, number=1)
                self.eng.config.outfname = os.path.splitext(self.eng.config.filename)[0]
                self.eng.config.default_file_names()
                os.chdir(self.eng.config.dirname)
                dirnew = self.eng.config.outfname + "_unidecfiles"
                if not os.path.isdir(dirnew):
                    print "Error: Need to process data in advance for Speed Batch Mode"
                    return False
                os.chdir(dirnew)

                # Write New Config and Run It
                self.export_config(self.eng.config.confname)
                if self.eng.config.imflag == 0:
                    ud.unidec_call(self.eng.config.UniDecPath, self.eng.config.confname)
                else:
                    ud.unidec_call(self.eng.config.UniDecIMPath, self.eng.config.confname)
                tend = time.clock()
                print "Run Time: %.2gs" % (tend - tstart)
                print "\n"
            tend = time.clock()
            print "\nTotal Speedy Batch Run Time: %.2gs" % (tend - tstarttop)

    def on_cross_validate(self, e=None):
        """
        Experimental...

        Run cross validation on spectrum by splitting the spectrum into subspectra and fitting each independently.

        Runs self.eng.cross_validate()
        Adds plots of mean and std deviation to plot2
        :param e: unused event
        :return: None
        """
        print "Cross Validation"
        self.export_config(self.eng.config.confname)
        mean, stddev = self.eng.cross_validate()
        norm = np.amax(self.eng.data.massdat[:, 1]) / np.amax(mean)
        self.view.plot2.plotadd(self.eng.data.massdat[:, 0], mean * norm, 'r', 'Mean')
        self.view.plot2.plotadd(self.eng.data.massdat[:, 0], (mean + stddev) * norm, 'y', 'Mean+STD')
        self.view.plot2.plotadd(self.eng.data.massdat[:, 0], (mean - stddev) * norm, 'y', 'Mean-STD')
        self.view.plot2.repaint()
        pass

    def on_pdf_report(self, e=None):
        """
        Creates PDF report.

        First, writes figures to PDF.
        Then sends results to texmaker, which creates a .tex file.
        Finally, runs pdflatex as commandline subprocess to convert .tex to PDF.

        :param e: event passed to self.view.on_save_figur_pdf
        :return: None
        """
        figureflags, files = self.view.on_save_figure_pdf(e)
        textmarkertab = [p.textmarker for p in self.eng.pks.peaks]
        peaklabels = [p.label for p in self.eng.pks.peaks]
        peakcolors = [p.color for p in self.eng.pks.peaks]
        peaks = np.array([[p.mass, p.height] for p in self.eng.pks.peaks])
        if self.eng.config.imflag == 0:
            texmaker.MakeTexReport(self.eng.config.outfname + '_report.tex', self.eng.config, self.eng.config.dirname,
                                   peaks, textmarkertab, peaklabels, peakcolors, figureflags)
            self.view.SetStatusText("TeX file Written", number=5)
            try:
                texmaker.PDFTexReport(self.eng.config.outfname + '_report.tex')
                self.view.SetStatusText("PDF Report Finished", number=5)
            except Exception, ex:
                self.view.SetStatusText("PDF Report Failed", number=5)
                print "PDF Report Failed to Generate. Check LaTeX installation.Need pdflatex in path.", ex
        else:
            print "PDF Figures written."
        pass

    def on_fft_window(self, e):
        print "FFT window..."
        fft_window.FFTWindow(self.view, self.eng.data.rawdata, self.eng.config)
        pass

    def on_left_click(self, xpos, ypos):
        """
        Triggered by pubsub from plot windows.
        Gets a m/z peak near the click, stores it, and waits for another click.
        When two clicks has been performed, it tries to calculate the mass from their m/z value.
        :param xpos: x position fed from event
        :param ypos: y position fed from event
        :return: None
        """
        plot = True
        if xpos is not None and ypos is not None:
            # print "x=%.2f y=%.2f" % (xpos, ypos)
            # Determine the limits for local max determination
            xlimits = self.view.plot1.subplot1.get_xlim()
            limdiff = abs(xlimits[1] - xlimits[0])
            window = limdiff * 0.01

            # Find the local max near the clicked position
            newxpos = ud.localmaxpos(self.eng.data.data2, xpos - window, xpos + window)
            if newxpos > 0:
                # If a suitable local max was found, use it.
                xpos = newxpos

            if self.view.plot1.x1 is None or xpos == self.view.plot1.x1:
                # Store the first value
                self.view.plot1.x1 = xpos
            else:
                # Store the second value
                self.view.plot1.x2 = xpos
                # Switch them if mixed up
                if self.view.plot1.x2 < self.view.plot1.x1:
                    self.view.plot1.x1, self.view.plot1.x2 = self.view.plot1.x2, self.view.plot1.x1
                print "m/z values:", self.view.plot1.x1, self.view.plot1.x2
                # Solve for the mass and charges
                mass, z1, z2 = ud.solve_for_mass(self.view.plot1.x1, self.view.plot1.x2)
                outstring = "Mass=%.2f z=%d, %d" % (mass, z1, z2)

                if np.all(np.abs(np.array(self.view.plot1.mlist) - mass) > window * z1 * 0.0) and plot:
                    self.view.plot1.mlist.append(mass)

                    newcolor = 'ybgrcmk'[len(self.view.plot1.mlist) % 6]
                    self.view.plot1.colors.append(newcolor)

                    try:
                        self.view.plot1.subplot1.legend_.remove()
                    except AttributeError:
                        pass
                    # Add new things
                    maxy = np.amax(self.eng.data.data2[:, 1])
                    self.view.plot1.addtext(str(mass), np.amax(self.eng.data.data2[:, 0]) * 0.97,
                                            maxy - 0.05 * len(self.view.plot1.mlist) * maxy, vlines=False,
                                            color=newcolor)
                elif plot:
                    index = ud.nearestunsorted(np.array(self.view.plot1.mlist), mass)
                    newcolor = self.view.plot1.colors[index]

                if plot:
                    # Add the charge state assignments to the plot
                    pad = 0.05 * np.amax(self.eng.data.data2[:, 1])
                    y1 = ud.interp_val(self.eng.data.data2, self.view.plot1.x1) + pad
                    y2 = ud.interp_val(self.eng.data.data2, self.view.plot1.x2) + pad
                    self.view.plot1.addtext(str(int(z1)), self.view.plot1.x1, y1, color=newcolor)
                    self.view.plot1.addtext(str(int(z2)), self.view.plot1.x2, y2, color=newcolor)
                    # Remove the legend

                # Reset and write out values
                self.view.SetStatusText(outstring, number=5)
                self.view.plot1.x1, self.view.plot1.x2 = None, None
        pass

    def on_grid_decon(self, e):
        GridDecon.GridDeconWindow(self.view, self.eng.data.data2, config=self.eng.config)

    def on_flip_mode(self, e):
        """
        Flips between MS and IM-MS mode
        :param e: wx event or anything (will flip if not 0)
        :return: None
        """
        if e is not 0:
            self.eng.config.imflag = (self.eng.config.imflag + 1) % 2
        self.remake_mainwindow(self.view.tabbed)

        self.view.import_config_to_gui()
        if self.eng.config.imflag == 1:
            print "Ion Mobility Mode"
        elif self.eng.config.imflag == 0:
            print "Mass Spec Mode"

    def on_flip_tabbed(self, e):
        """
        Flips between tabbed plots and a single window of plots
        :param e: wx Event or anything (will flip if not 0)
        :return: None
        """
        if e is not 0:
            tabbed = (self.view.tabbed + 1) % 2
        else:
            tabbed = self.view.tabbed
        self.remake_mainwindow(tabbed=tabbed)
        try:
            self.on_replot(e)
            self.view.peakpanel.add_data(self.eng.pks)
        except Exception, exc:
            print "Failed to replot when making window:", exc
        if self.view.tabbed == 1:
            print "Tabbed Mode"
        elif self.view.tabbed == 0:
            print "Single Plot Window Mode"

    def on_flip_twave(self, e):
        """
        Flips between T-Wave and Linear IM-MS
        :param e: wx Event or anything (will get value from Selection if not 0)
        :return: None
        """
        if e is not 0:
            self.eng.config.twaveflag = self.view.ctltwave.GetSelection()

        if self.eng.config.twaveflag == 0:
            self.eng.config.gasmass = 4.002602
            print "Using Linear Cell"
        elif self.eng.config.twaveflag > 0:
            self.eng.config.gasmass = 28.0134
            print "Using Travelling Wave"
        else:
            print "Error: Unsupported twaveflag.", self.eng.config.twaveflag
        self.remake_mainwindow(self.view.tabbed)
        # self.view.ctltwave.SetSelection(self.eng.config.twaveflag)
        # self.view.import_config_to_gui()

    def remake_mainwindow(self, tabbed=None):
        iconfile = self.view.icon_path
        self.view.Destroy()
        self.view = []
        self.view = mainwindow.Mainwindow(self, "UniDec", self.eng.config, iconfile=iconfile, tabbed=tabbed)
        self.view.Show()
        self.view.import_config_to_gui()


if __name__ == '__main__':
    # app2 = Shell()
    multiprocessing.freeze_support()
    app = UniDecApp()
    app.start()

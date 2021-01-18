import time
import os
import _thread
import wx
import numpy as np
import unidec

#
from pubsub import pub

import unidec_modules.unidectools as ud
import unidec_modules.IM_functions as IM_func
from unidec_modules import Extract2D, peakwidthtools, masstools, miscwindows, \
    MassDefects, mainwindow, nativez, ManualSelectionWindow, AutocorrWindow, fft_window, GridDecon, isotopetools
from unidec_modules.isolated_packages import FileDialogs, texmaker, score_window, texmaker_nmsgsb, navia_importer
import datacollector
import import_wizard
import unidec_modules.IM_windows as IM_wind
# import UniMin
from copy import deepcopy
import platform
import multiprocessing
from unidec_modules.unidec_presbase import UniDecPres
import Launcher
from iFAMS.wxiFAMS import iFAMS_Window

try:
    import unidec_modules.thermo_reader.rawreader as rawreader
except Exception as e:
    print("Error importing Thermo Raw Reader, try installing MSFileReader from Thermo and pymsfilereader")
    print(e)

# import FileDialog  # Needed for pyinstaller

__author__ = 'Michael.Marty'


# noinspection,PyBroadException,PyUnusedLocal,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyBroadException,PyUnusedLocal,PyBroadException
class UniDecApp(UniDecPres):
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
        UniDecPres.__init__(self, *args, **kwargs)
        self.init(*args, **kwargs)

    def init(self, *args, **kwargs):
        """
        Initialize Engine and View. Load defaults.
        :param args:
        :param kwargs:
        :return:
        """
        self.eng = unidec.UniDec()

        self.view = mainwindow.Mainwindow(self, "UniDec", self.eng.config)

        pub.subscribe(self.on_integrate, 'integrate')
        pub.subscribe(self.on_smash, 'smash')
        pub.subscribe(self.on_get_mzlimits, 'mzlimits')
        pub.subscribe(self.on_left_click, 'left_click')

        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)
        self.view.menu.update_recent()

        self.on_load_default(0)

        if "path" in kwargs:
            newdir, fname = os.path.split(kwargs["path"])
            self.on_open_file(fname, newdir)
            # self.on_dataprep_button(0)
            self.on_auto(0)

        # For testing, load up a spectrum at startup. Used only on MTM's computer.
        if False and platform.node() == "DESKTOP-R236BN2":
            # fname = "HSPCID.txt"
            fname = "0.txt"
            # fname = "test.raw"
            # fname = "250313_AQPZ_POPC_100_imraw_input.dat"
            newdir = os.path.join(os.getcwd(), "TestSpectra")
            # newdir = "C:\\cprog\\UniDecDemo"
            self.on_open_file(fname, newdir)
            # self.view.on_save_figure_eps(0)
            # self.on_dataprep_button(0)
            # self.on_auto(0)
            # self.on_integrate()
            # self.on_grid_decon(0)
            # self.make_cube_plot(0)
            # self.on_plot_peaks(0)
            # self.on_flip_tabbed(None)
            # self.on_label_max_charge_states(0)
            # self.view.plot1.copy_to_clipboard()

    # ..............................
    #
    #  Main Utility Functions
    #
    # ........................................

    def on_open(self, e=None):
        """
        Open dialog for file opening
        :param e: unused space for event
        :return: None
        """
        dlg = wx.FileDialog(self.view, "Choose a data file in x y list, mzML, or Thermo Raw format", '', "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            self.view.SetStatusText("Opening", number=5)
            filename = dlg.GetFilename()
            print("Opening: ", filename)
            if os.path.splitext(filename)[1] == ".zip":
                print("Can't open zip, try Load State.")
                return
            dirname = dlg.GetDirectory()
            self.on_open_file(filename, dirname)
        dlg.Destroy()

    def on_open_file(self, filename, directory, skipengine=False, **kwargs):
        """
        Opens a file. Run self.eng.open_file.
        :param filename: File name
        :param directory: Directory containing file
        :param skipengine: Boolean, Whether to skip running the engine (used when loading state)
        :return: None
        """
        # tstart =time.perf_counter()

        # Clear other plots and panels
        self.view.peakpanel.clear_list()
        self.view.clear_all_plots()
        if not skipengine:
            # Open File in Engine
            self.eng.open_file(filename, directory, **kwargs)

        # Set Status Bar Text Values
        self.view.SetStatusText("File: " + filename, number=1)
        # print self.view.imflag, self.eng.config.imflag
        if self.view.imflag != self.eng.config.imflag:
            print("Changing Modes")
            self.on_flip_mode(0)
        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText("R\u00B2 ", number=3)
        # Update view with data limits
        if self.eng.config.batchflag != 1:
            self.view.controls.ctlminmz.SetValue(str(np.amin(self.eng.data.data2[:, 0])))
            self.view.controls.ctlmaxmz.SetValue(str(np.amax(self.eng.data.data2[:, 0])))

        # Plot 1D
        if self.eng.config.batchflag == 0:
            self.makeplot1(imfit=False)
        # IM Loading and Plotting
        if self.eng.config.imflag == 1 and self.eng.config.batchflag != 1:
            self.view.controls.ctlmindt.SetValue(str(np.amin(self.eng.data.data3[:, 1])))
            self.view.controls.ctlmaxdt.SetValue(str(np.amax(self.eng.data.data3[:, 1])))
            # if self.eng.config.batchflag == 0:
            #    self.view.plot1im.contourplot(self.eng.data.rawdata3, self.eng.config, xlab="m/z (Th)",
            #                                  ylab="Arrival Time (ms)", title="IM-MS Data")
        # tstart = time.perf_counter()
        # Load Config to GUI
        self.import_config()

        # self.write_to_recent()

        self.view.SetStatusText("Ready", number=5)

        self.write_to_recent()
        self.view.menu.update_recent()

        # print("ImportConfig: %.2gs" % (time.perf_counter() - tstart))
        # if False:
        #    try:
        #        self.eng.unidec_imports(everything=False)
        #        self.after_unidec_run()
        #    except:
        #        pass

        # print("ImportData: %.2gs" % (time.perf_counter() - tstart))

    def on_open_kernel(self):
        """
        Opens a kernel file for use in DoubleDec deconvolution.
        :return: A string indicating the chosen file path
        """
        dlg = wx.FileDialog(self.view, "Open a kernel file")
        if dlg.ShowModal() == wx.ID_CANCEL:
            return
        file_path = dlg.GetPath()
        if os.path.splitext(os.path.basename(file_path))[1] == ".zip":
            print("Can't open zip")
            return
        dlg.Destroy()
        return file_path

    def on_load_everything(self, e=None):
        self.eng.unidec_imports(everything=False)
        self.after_unidec_run()
        self.on_pick_peaks()

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

    def on_load_state(self, e=None, filenew=None):
        """
        Loads a saved file state by running self.eng.load_state and then updating self.view.
        Runs several engine components rather than importing.
        :param e: unused space for event
        :return: None
        """
        if filenew is None:
            filenew = FileDialogs.open_file_dialog(message="Select UniDec Zip File to Open", file_types="*.zip")
            print(filenew)
        if filenew is not None:
            # Reset GUI
            self.on_reset(0)
            self.view.SetStatusText("Loading", number=5)
            tstart = time.perf_counter()
            dirname, filename = os.path.split(filenew)
            self.view.SetStatusText("File: " + filename, number=1)
            self.view.SetStatusText("R\u00B2 ", number=3)
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
            if os.path.isfile(self.eng.config.errorfile):
                self.after_unidec_run()
            # RePick Peaks
            if os.path.isfile(self.eng.config.peaksfile):
                self.on_pick_peaks(e)
            # Match Again
            # TODO: Figure out what is happening with matchlist
            if os.path.isfile(self.eng.config.matchfile):
                self.eng.config.matchlist = np.transpose(
                    np.genfromtxt(self.eng.config.matchfile, dtype='str', delimiter=","))

            tend = time.perf_counter()
            print("Loading Time: %.2gs" % (tend - tstart))
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
            if int(self.view.controls.ctlconvertflag.GetValue()) == 1:
                binsize = str(ud.string_to_value(self.view.controls.ctlbinsize.GetValue()))
                print("Converting at resolution of: " + binsize)
            else:
                binsize = "0"
                print("Converting using full resolution")
        else:
            binsize = None

        if self.eng.config.dirname is not None:
            self.view.SetStatusText("Converting", number=5)
            self.eng.config.dirname = os.path.abspath(self.eng.config.dirname)
            print("Loading Raw File: ", self.eng.config.dirname)
            self.eng.config.filename, self.eng.config.dirname = self.eng.raw_process(self.eng.config.dirname, True,
                                                                                     binsize=binsize)
            if self.eng.config.filename is not None:
                self.on_open_file(self.eng.config.filename, self.eng.config.dirname)
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
            fname = "PastedSpectrum_" + str(time.strftime("%Y_%b_%d_%H_%M_%S")) + ".txt"
            for t in text:
                if ".RAW" in t:
                    print(t)
                    fname = os.path.splitext(t)[0] + ".txt"
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
                np.savetxt(os.path.join(newdir, fname), data)
                print("Saved Pasted Spectrum as File:", fname, " in directory:", newdir)
                self.on_open_file(fname, newdir, pasted=True)
            else:
                print("Paste failed, got: ", data)
        except Exception as e:
            print(e)
            wx.MessageBox("Unable to open the clipboard", "Error")

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
        tstart = time.perf_counter()
        self.view.SetStatusText("Data Prep", number=5)
        self.export_config(self.eng.config.confname)
        self.eng.process_data()
        self.eng.get_auto_peak_width(set=False)
        self.import_config()
        # print("Data Prep Time1: %.2gs" % (time.perf_counter() - tstart))
        self.view.clear_all_plots()
        self.makeplot1(imfit=False)
        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText("R\u00B2 ", number=3)
        self.view.SetStatusText("Data Prep Done", number=5)
        tend = time.perf_counter()
        print("Data Prep Done. Time: %.2gs" % (tend - tstart))
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
            print("Error ", out)
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
                _thread.start_new_thread(self.makeplot3, (1,))
                _thread.start_new_thread(self.makeplot5, (1,))
                self.makeplot1(1)
                self.makeplot2(1)
        else:
            self.view.SetStatusText("UniDec Plot", number=5)
            self.make_im_plots()

        self.view.SetStatusText("R\u00B2: " + str(self.eng.config.error), number=3)

        self.view.plot4.clear_plot()
        self.view.plot6.clear_plot()
        self.view.peakpanel.clear_list()

    def on_pick_peaks(self, e=None):
        """
        Pick peaks and perform initial plots on them.
        :param e: unused space for event
        :return: None
        """
        print("Peak Picking")
        self.view.SetStatusText("Detecting Peaks", number=5)
        tstart = time.perf_counter()
        self.export_config(self.eng.config.confname)
        self.eng.pick_peaks()
        self.view.SetStatusText("Plotting Peaks", number=5)
        if self.eng.config.batchflag == 0:
            self.view.peakpanel.add_data(self.eng.pks)
            self.makeplot2(1)
            self.makeplot6(1)
            self.makeplot4(1)
        self.view.SetStatusText("Peak Pick Done", number=5)
        self.on_score()
        pass

    def on_plot_peaks(self, e=None):
        """
        Plots individual peaks as simulated spectra in plot 4.
        :param e: unused space for event
        :return: None
        """
        self.view.SetStatusText("Convolving", number=5)
        self.export_config(None)
        tstart = time.perf_counter()
        self.eng.convolve_peaks()
        tend = time.perf_counter()
        print("Convolving: %.2gs" % (tend - tstart))
        self.view.SetStatusText("Plotting", number=5)
        self.makeplot4(1)
        self.export_config(self.eng.config.confname)
        self.view.SetStatusText("Peak Plot Done", number=5)
        print("peak plotting")
        pass

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

    def makeplot1(self, e=None, intthresh=False, imfit=True):
        """
        Plot data and fit in self.view.plot1 and optionally in plot1fit
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.perf_counter()
            leg = False

            if self.eng.config.imflag == 1:
                try:
                    self.view.plot1im.contourplot(self.eng.data.data3, self.eng.config, xlab="m/z (Th)",
                                                  ylab="Arrival Time (ms)", title="IM-MS Data")
                except:
                    pass
                if imfit:
                    try:
                        self.view.plot1fit.contourplot(self.eng.data.fitdat2d, self.eng.config, xlab="m/z (Th)",
                                                       ylab="Arrival Time (ms)", title="IM-MS Fit")
                    except:
                        pass

            if self.eng.config.reductionpercent < 0:
                print("Making Dot Plot")
                data2 = ud.dataprep(self.eng.data.rawdata, self.eng.config, peaks=False, intthresh=False)
                self.view.plot1.plotrefreshtop(data2[:, 0], data2[:, 1],
                                               "Data and UniDec Fit",
                                               "m/z (Th)", "Normalized Intensity", "Data", self.eng.config,
                                               nopaint=True)
                self.view.plot1.plotadddot(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1], 'blue', "o", "Peaks")

                try:
                    if len(self.eng.data.fitdat) > 0 and imfit:
                        self.view.plot1.plotadddot(self.eng.data.data2[:, 0], self.eng.data.fitdat, 'red', "s",
                                                   "Fit Data")
                        leg = True
                    pass
                except:
                    pass

            else:
                self.view.plot1.plotrefreshtop(self.eng.data.data2[:, 0], self.eng.data.data2[:, 1],
                                               "Data and UniDec Fit",
                                               "m/z (Th)", "Normalized Intensity", "Data", self.eng.config,
                                               nopaint=True)

                if self.eng.config.intthresh != 0 and self.eng.config.imflag == 0 and intthresh:
                    self.view.plot1.plotadd(self.eng.data.data2[:, 0],
                                            np.zeros_like(self.eng.data.data2[:, 1]) + self.eng.config.intthresh, "red",
                                            "Noise Threshold")
                    leg = True

                try:
                    if len(self.eng.data.fitdat) > 0 and imfit:
                        self.view.plot1.plotadd(self.eng.data.data2[:, 0], self.eng.data.fitdat, 'red', "Fit Data")
                        leg = True
                    pass
                except:
                    pass
            if self.eng.config.aggressiveflag != 0 and len(self.eng.data.baseline) == len(self.eng.data.fitdat):
                self.view.plot1.plotadd(self.eng.data.data2[:, 0], self.eng.data.baseline, 'blue', "Baseline")
            if leg:
                self.view.plot1.add_legend()
            self.view.plot1.repaint()
            tend = time.perf_counter()
            print("Plot 1: %.2gs" % (tend - tstart))

    def makeplot3(self, e=None):
        """
        Plot m/z vs charge grid.
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.perf_counter()
            self.view.plot3.contourplot(self.eng.data.mzgrid, self.eng.config)
            tend = time.perf_counter()
            print("Plot 3: %.2gs" % (tend - tstart))

    def makeplot5(self, e=None):
        """
        Plot mass vs charge grid
        :param e: unused event
        :return: None
        """
        if self.eng.config.batchflag == 0:
            tstart = time.perf_counter()
            self.view.plot5.contourplot(
                xvals=self.eng.data.massdat[:, 0], yvals=self.eng.data.ztab, zgrid=self.eng.data.massgrid,
                config=self.eng.config, title="Mass vs. Charge", test_kda=True)
            tend = time.perf_counter()
            print("Plot 5: %.2gs" % (tend - tstart))

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
                marks = []
                for i, p in enumerate(self.eng.pks.peaks):
                    if p.ignore == 0:
                        num += 1
                        if show == "height":
                            ints.append(p.height)
                        elif show == "integral":
                            ints.append(p.integral)
                        else:
                            ints.append(0)
                        cols.append(p.color)
                        labs.append(p.label)
                        marks.append(p.marker)
                indexes = list(range(0, num))
                self.view.plot6.barplottop(indexes, ints, labs, cols, "Species", "Intensity",
                                           "Peak Intensities", repaint=False)
                for i in indexes:
                    self.view.plot6.plotadddot(i, ints[i], cols[i], marks[i])
            self.view.plot6.repaint()

    def on_plot_composite(self, e=None):
        """
        Plots a composite of isolated peaks that are selected.
        :param e: unused event
        :return: None
        """
        try:
            self.eng.pks.composite = np.zeros(len(self.eng.data.data2))
            for i in range(0, self.eng.pks.plen):
                if self.eng.pks.peaks[i].ignore == 0:
                    self.eng.pks.composite += self.eng.pks.peaks[i].stickdat

            self.view.plot4.plotadd(self.eng.data.data2[:, 0], self.eng.pks.composite, "b", "useless label")
            self.view.plot4.repaint()
        except ValueError:
            print("Need to hit Plot Species button first")

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
            print("Made IM Plots")
            try:
                self.view.plot3color.make_color_plot(self.eng.data.mztgrid, np.unique(self.eng.data.data3[:, 0]),
                                                     np.unique(self.eng.data.data3[:, 1]), self.eng.data.ztab)
            except Exception as e:
                print("Color Plot Error", e)

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
            print("Plotted predicted native CCS values")

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
            starttime = time.perf_counter()
            self.view.plot9.cubeplot(np.unique(self.eng.data.data3[:, 0]), np.unique(self.eng.data.data3[:, 1]),
                                     self.eng.data.ztab, np.sum(self.eng.data.mztgrid, axis=2),
                                     np.sum(self.eng.data.mztgrid, axis=1), np.sum(self.eng.data.mztgrid, axis=0),
                                     xlab="m/z (Th)", ylab="Arrival Time (ms)", zlab="Charge",
                                     cmap=self.eng.config.cmap)
            endtime = time.perf_counter()
            print("Finished m/z Cube in: ", (endtime - starttime), " s")
        except Exception as ex:
            print("Failed m/z cube", ex)
            pass
        try:
            starttime = time.perf_counter()
            self.view.plot10.cubeplot(self.eng.data.massdat[:, 0], self.eng.data.ccsdata[:, 0], self.eng.data.ztab,
                                      self.eng.data.massccs, self.eng.data.massgrid.reshape(
                    (len(self.eng.data.massdat), len(self.eng.data.ztab))), self.eng.data.ccsz.transpose(),
                                      xlab="Mass (Da)",
                                      ylab="CCS (${\AA}$$^2$)", zlab="Charge", cmap=self.eng.config.cmap)
            endtime = time.perf_counter()
            print("Finished Final Cube in: ", (endtime - starttime), " s")
        except Exception as ex:
            print("Failed final cube", ex)
            pass

    def on_autoformat(self, e=None):
        self.eng.pks.auto_format()
        self.on_delete()
        self.view.peakpanel.add_data(self.eng.pks)

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

    def on_plot_offsets(self, e=None):
        """
        Transform the mass vs charge grid to a mass vs charge offset grid. Plot in self.view.plot5.
        :param e: unused event
        :return: None
        """
        print("Tranforming charge to charge offset...")
        oaxis, outgrid = self.eng.mass_grid_to_f_grid()
        print("Plotting offsets")
        self.view.plot5.contourplot(xvals=self.eng.data.massdat[:, 0], yvals=oaxis, zgrid=outgrid,
                                    config=self.eng.config,
                                    ylab="Charge Offset", title="Mass vs. Charge Offset", test_kda=True)

    def on_color_plot1d(self, e=None, filled=False):
        self.on_integrate(plot=False)
        # self.eng.get_peaks_scores()
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
        if not ud.isempty(self.eng.pks.peaks) and limits[0] <= np.amin(self.eng.data.massdat[:, 0]) and limits[
            1] >= np.amax(self.eng.data.massdat[:, 0]):
            print("Auto Integrating")
            self.eng.autointegrate()
        else:
            integral, intdat = self.eng.integrate(limits)
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
                print("Differences: ", limits - self.eng.pks.peaks[i].mass)

            else:
                self.view.plot2.addtext(str(integral), np.mean(np.array(limits)),
                                        np.amax(intdat[:, 1]) + 0.05 * np.amax(self.eng.data.massdat[:, 1]),
                                        range=limits)
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
            except Exception as ex:
                print("Didn't update bar chart", ex)
                pass
        pass

    def on_smash(self):
        """
        Triggered by right click on plot1.
        Smashes the zoomed region to 0. Used to eliminate unwanted peaks.
        :return: None
        """
        if wx.GetKeyState(wx.WXK_CONTROL):
            print("Smashing!")
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

            self.view.SetStatusText("R\u00B2 ", number=3)
            message = "Smashing peaks from " + str(limits[0]) + " to " + str(limits[1]) + "\nReprocess data to undo"
            self.warn(message)
            pass

    def on_full(self, e=None):
        maxmz = np.amax(self.eng.data.rawdata[:, 0])
        minmz = np.amin(self.eng.data.rawdata[:, 0])
        self.view.controls.ctlminmz.SetValue(str(minmz))
        self.view.controls.ctlmaxmz.SetValue(str(maxmz))
        self.on_dataprep_button()

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
        print("Weighted average charge state:", cavg)
        print("Weighted standard deviation:", cstd)

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
            print("Sorry. Waters Raw file converter only works on windows. Convert to txt file on a seperate machine.")
            return None
            # self.eng.config.dirname="C:\\MassLynx\\"
            # if(not os.path.isdir(self.eng.config.dirname)):
        self.eng.config.dirname = ''
        if dirs is None:
            dirs = FileDialogs.open_multiple_dir_dialog("Select Raw Folders", self.eng.config.dirname)
        print(dirs)
        if dirs is not None:
            for d in dirs:
                if clip:
                    if "(C:)\\" in d:
                        d = "C:\\" + d.split("(C:)\\")[1]
                    elif "C:\\" in d:
                        d = "C:\\" + d.split("C:\\")[1]
                    else:
                        d = 'C:\\' + d[15:]
                self.eng.raw_process(d, False)

        print("Batch Converted")
        pass

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
            print("Simulating Peaks")
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

    def on_unidec_path(self, e=None):
        """
        Opens file dialog to specify new UniDec binary.
        Only useful for advanced users or if something has gone horribly wrong.
        Modifies self.eng.config.UniDecPath and others.
        :param e: unused event
        :return: None
        """
        dlg = wx.FileDialog(self.view, "Locate the UniDec Executable", self.eng.config.UniDecDir, "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            self.eng.config.UniDecName = dlg.GetFilename()
            self.eng.config.UniDecDir = dlg.GetDirectory()
            self.eng.config.UniDecPath = os.path.join(self.eng.config.UniDecDir, self.eng.config.UniDecName)
            print("New Path:", self.eng.config.UnIDecPath)
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
                out = dlg.ShowModal()
                if out == 0:
                    self.import_config(None)
                else:
                    self.export_config()
        else:
            print("Load Data First")
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
            print("Running UniDec to Generate Outputs")
            self.eng.config.zout = -1
            self.export_config(self.eng.config.confname)
            ud.unidec_call(self.eng.config)
            dlg = IM_wind.IMToolExtract(self.view)
            dlg.initialize_interface(self.eng.data.massdat, self.eng.data.ccsdata, self.eng.data.massccs,
                                     self.eng.config, self.eng.pks)
            dlg.ShowModal()
            self.eng.config.zout = 0
        pass

    def on_kendrick(self, e=None):
        """
        Opens Kendrick Mass Defect Analysis window.
        :param e: unused event
        :return: None
        """
        MassDefects.MassDefectWindow(self.view, [self.eng.data.massdat], config=self.eng.config,
                                     pks=self.eng.pks, value=self.eng.config.molig, directory=self.eng.config.udir)

    def on_iFAMS(self, e=None):
        iFAMS_Window(self.view, self.eng.data.data2, config=self.eng.config, directory=os.getcwd())

    def on_navia(self, e=None):
        with wx.FileDialog(self.view, "Open NaViA session", wildcard="XYZ files (*.navia)|*.navia",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # the user changed their mind

            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            newpath = navia_importer.navia_import(pathname)
            newdir, newfile = os.path.split(newpath)
            self.on_open_file(newfile, newdir)

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
            print("Calibration Success! Polynomial Coefficients (order = n to 0):", coeff)
        except Exception as e:
            print("Calibration failed:", e)

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
            print("Center of Mass over region", limits, ":", com)
            self.view.plot2.addtext(str(com), com, 0.99 * np.amax(self.eng.data.massdat[:, 1]))
            self.view.plot2.repaint()
        else:
            print("Need to get zero-charge mass spectrum first.")

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
        massoutputfile = self.eng.config.massdatfile
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
        print("Fitting Zero-charge Mass Spectra to Peak Shapes")
        massfitdat, massfit = self.eng.fit_all_masses()
        print("Fit: ", massfit)
        self.makeplot2(1)
        self.view.plot2.plotadd(self.eng.data.massdat[:, 0],
                                massfitdat / np.amax(massfitdat) * np.amax(self.eng.data.massdat[:, 1]), "green",
                                "Minimization")
        self.view.plot2.repaint()
        for i in range(0, self.eng.pks.plen):
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
                file_type="Any Type|*.*|Text (.txt)|*.txt|Thermo (.raw)|*.raw|mzML (.mzML)|*.mzML")

        self.eng.config.batchflag = 1 + flag
        tstarttop = time.perf_counter()
        print(batchfiles)
        if batchfiles is not None:
            self.view.clear_all_plots()
            for i, path in enumerate(batchfiles):
                print(path)
                tstart = time.perf_counter()
                dirname, filename = os.path.split(path)
                self.on_open_file(filename, dirname)
                self.on_dataprep_button(e)
                self.on_unidec_button(e)
                self.on_pick_peaks(e)
                self.on_export_params(e)
                # outfile = os.path.join(dirname, self.eng.config.outfname + ".zip")
                # self.on_save_state(0, outfile)
                # print "File saved to: " + str(outfile)
                print("Completed: " + path)
                tend = time.perf_counter()
                print("Run Time: %.2gs" % (tend - tstart))
                print("\n")
        self.eng.config.batchflag = 0
        tend = time.perf_counter()
        print("\nTotal Batch Run Time: %.3gs" % (tend - tstarttop))

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
        This function is useful for large data sets where you are using DataCollector or some other program to
        analyze the results outside of UniDec. It allows new parameters to be run with the maximum speed.
        :param e: unused event
        :return: None
        """
        batchfiles = FileDialogs.open_multiple_files_dialog(message="Select Files To Process With Current Parameters",
                                                            file_type="Text (.txt)|*.txt|Any Type|*.*")
        if batchfiles is not None:
            self.view.SetStatusText("Speedy Batch Run", number=5)
            tstarttop = time.perf_counter()
            for i, path in enumerate(batchfiles):
                print(path)
                tstart = time.perf_counter()
                # Open File Stripped
                self.eng.config.dirname, self.eng.config.filename = os.path.split(path)
                self.view.SetStatusText("File: " + self.eng.config.filename, number=1)

                basename = os.path.splitext(self.eng.config.filename)[0]

                dirnew = os.path.join(self.eng.config.dirname, basename + "_unidecfiles")
                if not os.path.isdir(dirnew):
                    print("Error: Need to process data in advance for Speed Batch Mode")
                    return False
                self.eng.config.outfname = os.path.join(dirnew, basename)
                self.eng.config.default_file_names()

                # Write New Config and Run It
                self.export_config(self.eng.config.confname)

                ud.unidec_call(self.eng.config)

                tend = time.perf_counter()
                print("Run Time: %.2gs" % (tend - tstart))
                print("\n")
            tend = time.perf_counter()
            print("\nTotal Speedy Batch Run Time: %.2gs" % (tend - tstarttop))

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
            texmaker.MakeTexReport(self.eng.config.outfname + '_report.tex', self.eng.config, self.eng.config.udir,
                                   peaks, textmarkertab, peaklabels, peakcolors, figureflags)
            self.view.SetStatusText("TeX file Written", number=5)
            try:
                texmaker.PDFTexReport(self.eng.config.outfname + '_report.tex')
                self.view.SetStatusText("PDF Report Finished", number=5)
            except Exception as ex:
                self.view.SetStatusText("PDF Report Failed", number=5)
                print("PDF Report Failed to Generate. Check LaTeX installation.Need pdflatex in path.", ex)
        else:
            print("PDF Figures written.")
        pass

    def on_nmsgsb_report(self, e=0):

        """
        Creates PDF report for the Native MS Guided Structural Biology format.
        First, writes figures to PDF.
        Then sends results to texmaker, which creates a .tex file.
        Finally, runs pdflatex as commandline subprocess to convert .tex to PDF.
        :param e: event passed to self.view.on_save_figur_pdf
        :return: None
        """
        path = os.path.join(self.eng.config.dirname, self.eng.config.filename)
        print(path)
        rawsamplename = ""
        defaultvalue = ""
        if os.path.splitext(path)[1].lower() == ".raw":
            print("Getting Raw Data")
            defaultvalue = rawreader.get_raw_metadata(path)
            # try:
            #    rawoutput = rawreader.get_raw_metadata(path)
            # except:
            #    rawoutput = None
            rawsamplename = rawreader.get_raw_samplename(path)
        # Andrew - edit

        dialog = miscwindows.SingleInputDialog(self.view)
        dialog.initialize_interface(title="Report Info: Input1;Input2;...;InputN", message="Set Inputs Here: ",
                                    defaultvalue=defaultvalue)
        dialog.ShowModal()
        output = dialog.value

        self.view.shrink_all_figures(figsize=(6, 5))
        figureflags, files = self.view.on_save_figure_eps(e)
        figureflags, files = self.view.on_save_figure_pdf(e)
        textmarkertab = [p.textmarker for p in self.eng.pks.peaks]
        peaklabels = [p.label for p in self.eng.pks.peaks]
        peakcolors = [p.color for p in self.eng.pks.peaks]
        peaks = np.array([[p.mass, p.height] for p in self.eng.pks.peaks])
        uniscore = self.eng.pks.uniscore
        # str(round(self.eng.pks.uniscore * 100, 2))
        # oligos = np.array(oligos)
        # match = np.array([[p.peaks, p.matches, p.errors, p.names] for p in self.eng.config.matchlist])
        # match = np.transpose(self.eng.config.matchlist)
        # match = self.eng.config.matchlist
        # self.eng.config.matchlist = np.transpose(
        #    np.genfromtxt(self.eng.config.matchfile, dtype='str', delimiter=","))

        if os.path.isfile(self.eng.config.matchfile):
            match = np.transpose(self.eng.config.matchlist)
        else:
            match = "-"
        if self.eng.config.imflag == 0:
            texmaker_nmsgsb.MakeTexReport(self.eng.config.outfname + '_report.tex', self.eng.config,
                                          self.eng.config.udir,
                                          peaks, textmarkertab, peaklabels, peakcolors, figureflags, output,
                                          rawsamplename, match, uniscore)
            self.view.SetStatusText("TeX file Written", number=5)
            try:
                texmaker_nmsgsb.PDFTexReport(self.eng.config.outfname + '_report.tex')
                self.view.SetStatusText("PDF Report Finished", number=5)
            except Exception as ex:
                self.view.SetStatusText("PDF Report Failed", number=5)
                print("PDF Report Failed to Generate. Check LaTeX installation.Need pdflatex in path.", ex)
        else:
            print("PDF Figures written.")
        # self.on_replot()
        # self.view.shrink_all_figures(figsize=self.eng.config.figsize)
        # print("Resetting Figure Sizes", self.eng.config.figsize)
        # self.on_replot()
        self.on_flip_tabbed(e=0)
        pass

    def on_fft_window(self, e):
        print("FFT window...")
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
                print("m/z values:", self.view.plot1.x1, self.view.plot1.x2)
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

    def on_label_max_charge_states(self, e):
        charges = np.arange(self.eng.config.startz, self.eng.config.endz + 1)
        if self.eng.config.adductmass > 0:
            sign = "+"
        else:
            sign = "-"

        self.view.plot4.textremove()
        for i in range(0, self.eng.pks.plen):
            p = self.eng.pks.peaks[i]
            if p.ignore == 0:
                if (not ud.isempty(p.mztab)) and (not ud.isempty(p.mztab2)):
                    mztab = np.array(p.mztab)
                    mztab2 = np.array(p.mztab2)
                    maxval = np.amax(mztab[:, 1])
                    for k in range(0, len(mztab)):
                        if mztab[k, 1] == maxval:
                            self.view.plot4.addtext(sign + str(charges[k]), mztab2[k, 0],
                                                    mztab2[k, 1] + 0.075 * np.amax(self.eng.data.data2[:, 1]),
                                                    vlines=False, color=p.color)
        self.view.plot4.repaint()

    def on_label_avg_charge_states(self, e):
        charges = np.arange(self.eng.config.startz, self.eng.config.endz + 1)
        if self.eng.config.adductmass > 0:
            sign = "+"
        else:
            sign = "-"

        self.view.plot4.textremove()
        for i in range(0, self.eng.pks.plen):
            p = self.eng.pks.peaks[i]
            if p.ignore == 0:
                if (not ud.isempty(p.mztab)) and (not ud.isempty(p.mztab2)):
                    mztab = np.array(p.mztab)
                    avgcharge = ud.weighted_avg(charges, mztab[:, 1])
                    p.avgcharge = avgcharge
                    pos = (p.mass + self.eng.config.adductmass * avgcharge) / avgcharge

                    print("Mass:", p.mass, "Average Charge:", avgcharge)
                    mztab = ud.datachop(mztab, np.amin(self.eng.data.data2[:, 0]), np.amax(self.eng.data.data2[:, 0]))
                    self.view.plot4.plotadd(mztab[:, 0], mztab[:, 1], p.color)
                    self.view.plot4.addtext(sign + str(np.round(avgcharge, 2)), pos, np.amax(mztab[:, 1]) + 0.07,
                                            vlines=True, color=p.color)
        self.view.plot4.repaint()
        self.view.peakpanel.add_data(self.eng.pks, show="avgcharge")

    def on_plot_isotope_distribution(self, e=0):
        for i in range(0, self.eng.pks.plen):
            p = self.eng.pks.peaks[i]
            if p.ignore == 0:
                # print(p.mass)
                dist = isotopetools.calc_averagine_isotope_dist(p.mass)
                dist[:, 1] *= p.height / np.amax(dist[:, 1])
                self.view.plot2.plotadd(dist[:, 0], dist[:, 1], colval=p.color)
        self.view.plot2.repaint()

    def on_score(self, e=0):
        self.eng.dscore()
        self.view.peakpanel.add_data(self.eng.pks, show="dscore")
        self.view.SetStatusText("UniScore: " + str(round(self.eng.pks.uniscore * 100, 2)), number=3)

    def on_score2(self, e=0):
        self.on_filter_peaks(e)
        self.view.SetStatusText("UniScore: " + str(round(self.eng.pks.uniscore * 100, 2)), number=3)
        self.makeplot2()
        self.makeplot4()
        self.makeplot6()
        self.on_score_window()

    def on_score_window(self, e=0):
        self.on_score()
        sw = score_window.ScoreFrame(self.view)
        sw.populate(self.eng.pks)
        pass

    def on_score_label(self, e=0):
        self.on_score()
        offset = 0.08 * np.amax(self.eng.data.massdat[:, 1])
        for p in self.eng.pks.peaks:
            text = str(int(round(p.dscore * 100)))
            self.view.plot2.addtext(text, p.mass, p.height + offset, vlines=False)

    # def on_remove_noise_points(self, e=0):
    #    self.on_dataprep_button(removenoise=True)

    def on_score_FDR(self, e=0):
        self.eng.estimate_FDR()

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

    def on_flip_mode(self, e=None):
        """
        Flips between MS and IM-MS mode
        :param e: wx event or anything (will flip if not 0)
        :return: None
        """
        if e != 0:
            self.eng.config.imflag = (self.eng.config.imflag + 1) % 2
        self.remake_mainwindow(self.view.tabbed)
        if self.eng.config.imflag == 1:
            print("Ion Mobility Mode")
            if self.eng.config.mzbins == 0:
                self.eng.config.mzbins = 1
        elif self.eng.config.imflag == 0:
            print("Mass Spec Mode")
        self.view.import_config_to_gui()

    def on_flip_tabbed(self, e):
        """
        Flips between tabbed plots and a single window of plots
        :param e: wx Event or anything (will flip if not 0)
        :return: None
        """
        if e != 0:
            tabbed = (self.view.tabbed + 1) % 2
        else:
            tabbed = self.view.tabbed
        self.remake_mainwindow(tabbed=tabbed)
        try:
            self.on_replot(e)
            self.view.peakpanel.add_data(self.eng.pks)
        except Exception as exc:
            print("Failed to replot when making window:", exc)
        if self.view.tabbed == 1:
            print("Tabbed Mode")
        elif self.view.tabbed == 0:
            print("Single Plot Window Mode")

    def on_flip_twave(self, e):
        """
        Flips between T-Wave and Linear IM-MS
        :param e: wx Event or anything (will get value from Selection if not 0)
        :return: None
        """
        if e != 0:
            self.eng.config.twaveflag = self.view.controls.ctltwave.GetSelection()

        if self.eng.config.twaveflag == 0:
            self.eng.config.gasmass = 4.002602
            print("Using Linear Cell")
        elif self.eng.config.twaveflag > 0:
            self.eng.config.gasmass = 28.0134
            print("Using Travelling Wave")
        else:
            print("Error: Unsupported twaveflag.", self.eng.config.twaveflag)
        self.remake_mainwindow(self.view.tabbed)
        # self.view.ctltwave.SetSelection(self.eng.config.twaveflag)
        # self.view.import_config_to_gui()

    def remake_mainwindow(self, tabbed=None):
        iconfile = self.view.icon_path
        # evt=EventManager()
        # print evt.GetStats()
        wx.GetApp().Yield()
        self.view.on_exit()
        self.view = []
        self.view = mainwindow.Mainwindow(self, "UniDec", self.eng.config, iconfile=iconfile, tabbed=tabbed)
        self.view.Show()
        self.view.import_config_to_gui()

    def on_undo(self, e=None):
        self.view.export_gui_to_config()
        self.eng.undo()
        self.view.import_config_to_gui()
        # print "Undo"

    def on_redo(self, e=None):
        self.view.export_gui_to_config()
        self.eng.redo()
        self.view.import_config_to_gui()
        # print "Redo"

    def on_write_hdf5(self, e=None):
        self.eng.write_hdf5()
        print("Wrote: ", self.eng.config.hdf_file)

    def on_launcher(self, e=None):
        # self.view.Destroy()
        launcher = Launcher.UniDecLauncher()
        launcher.start()


# TODO: Charge state distributions of each peak

if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = UniDecApp()
    app.start()

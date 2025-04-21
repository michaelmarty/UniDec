import unidec.tools as ud
from unidec.modules.unidec_presbase import UniDecPres
from unidec.IsoDec.runtime import IsoDecRuntime
from unidec.IsoDec.IDGUI.IsoDecView import IsoDecView
from unidec.engine import UniDec
import os
import wx
import time
from unidec.modules.peakstructure import Peaks
from unidec.IsoDec.datatools import get_all_centroids
import numpy as np
from unidec.modules.isolated_packages import FileDialogs
import platform
from copy import deepcopy


class IsoDecPres(UniDecPres):
    def __init__(self):
        super().__init__()
        self.isodeceng = IsoDecRuntime()
        self.eng = UniDec()
        self.init_config()

        self.view = IsoDecView(self, "IsoDec", self.eng.config, iconfile=None)
        try:
            if platform.node() == 'MTM-VOSTRO':
                self.on_ex()
                self.on_unidec_button()
                self.on_plot_peaks()
                self.on_plot_dists()
        except:
            pass

    def init_config(self):
        self.translate_id_config()
        self.eng.config.massbins = 1
        self.eng.config.aggressiveflag = 8
        self.eng.config.numit = self.isodeceng.config.knockdown_rounds

    def on_init_config(self, e=None):
        self.isodeceng.config.__init__()
        self.eng.config.initialize()
        self.init_config()
        self.view.import_config_to_gui()

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

    def on_open_file(self, filename, directory, skipengine=False, refresh=False, **kwargs):
        """
        Opens a file. Run self.eng.open_file.
        :param filename: File name
        :param directory: Directory containing file
        :param skipengine: Boolean, Whether to skip running the engine (used when loading state)
        :return: None
        """
        # tstart =time.perf_counter()
        self.export_config()
        # Clear other plots and panels
        self.view.peakpanel.clear_list()
        self.view.clear_all_plots()
        if not skipengine:
            # Open File in Engine
            self.top_path = os.path.join(directory, filename)
            self.eng.open_file(filename, directory, refresh=refresh, prefer_centroid=True, isodeceng=self.isodeceng,
                               **kwargs)

        # Set Status Bar Text Values
        self.view.SetStatusText("File: " + filename, number=1)
        # print self.view.imflag, self.eng.config.imflag

        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        # self.view.SetStatusText("R\u00B2 ", number=3)

        # Plot 1D
        if self.eng.config.batchflag == 0:
            self.makeplot1(imfit=False)

        self.import_config()
        self.view.SetStatusText("Ready", number=5)

        self.write_to_recent()
        self.view.menu.update_recent()

    def on_raw_open(self, evt=None, dirname=None):
        # Export configuration file
        self.export_config(self.eng.config.confname)

        # If dirname is not passed, open the dialog and get the path
        if dirname is None:
            # This will open the dialog and get the directory path as a string
            self.eng.config.dirname = FileDialogs.open_single_dir_dialog("Choose a raw file", '')
        else:
            # If dirname is passed directly, use that
            self.eng.config.dirname = dirname

        # Ensure the directory is a valid string
        if self.eng.config.dirname is None:
            print("No directory selected.")
            return

        # Check if the directory exists and handle accordingly
        binsize = None

        # Process the raw file if a valid directory is chosen
        if self.eng.config.dirname:
            self.view.SetStatusText("Converting", number=5)

            # Ensure that the dirname is a valid path (avoid passing CommandEvent object)
            self.eng.config.dirname = os.path.abspath(self.eng.config.dirname)  # This should now be a valid string

            print("Loading Raw File: ", self.eng.config.dirname)

            # Call the raw processing function
            self.eng.config.filename, self.eng.config.dirname = self.eng.raw_process(self.eng.config.dirname, True,
                                                                                     binsize=binsize)

            if self.eng.config.filename is not None:
                self.on_open_file(self.eng.config.filename, self.eng.config.dirname)

    def makeplot1(self, e=None, intthresh=False, imfit=True):
        """
        Plot data and fit in self.view.plot1 and optionally in plot1fit
        :param e: unused event
        :return: None
        """
        self.view.plot1.centroid_plot(self.eng.data.data2, xlabel="m/z", ylabel="Intensity", color="k")
        # self.eng.makeplot1(plot=self.view.plot1, intthresh=intthresh, imfit=imfit)

    def makeplot2(self, e=None):
        """
        Plot mass spectrum in
        :param e: unused event
        :return: None
        """
        self.view.plot2.centroid_plot(self.eng.data.massdat, xlabel="Mass", ylabel="Intensity", color="k")

    def on_dataprep_button(self, e=None):
        """
        Run data preparation
        :param e: unused event
        :return: None
        """
        self.view.SetStatusText("Data Prep", number=5)
        self.fix_parameters()
        self.export_config(self.eng.config.confname)

        if not self.eng.config.centroided:
            path = self.top_path.split(".")
            #Everything going through IsoDec should be centroided
            if path[-1].lower() == 'raw' and not os.path.isdir(self.top_path):
                self.eng.data.data2 = self.isodeceng.reader.grab_all_centroid_dat()
            else:
                self.eng.data.data2 = get_all_centroids(self.eng.data.data2)
            self.eng.config.centroided = True
        self.eng.process_data()


        self.import_config()
        self.view.clear_all_plots()
        self.view.peakpanel.clear_list()
        self.makeplot1(imfit=False, intthresh=True)
        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText("Data Prep Done", number=5)
        pass

    def translate_config(self):
        self.isodeceng.config.peakwindow = self.eng.config.peakwindow
        self.isodeceng.config.peakthresh = self.eng.config.peakthresh
        self.isodeceng.config.knockdown_rounds = self.eng.config.numit
        self.isodeceng.config.phaseres = self.eng.config.aggressiveflag
        self.isodeceng.config.matchtol = self.eng.config.filterwidth
        self.isodeceng.config.maxshift = self.eng.config.msig
        self.isodeceng.config.css_thresh = self.eng.config.csig
        self.isodeceng.config.adductmass = self.eng.config.adductmass
        self.isodeceng.config.datathreshold = self.eng.config.exthresh
        self.isodeceng.config.mzwindow = [self.eng.config.integratelb, self.eng.config.integrateub]

    def translate_id_config(self):
        self.eng.config.peakwindow = self.isodeceng.config.peakwindow
        self.eng.config.peakthresh = self.isodeceng.config.peakthresh
        self.eng.config.numit = self.isodeceng.config.knockdown_rounds
        self.eng.config.aggressiveflag = self.isodeceng.config.phaseres
        self.eng.config.filterwidth = self.isodeceng.config.matchtol
        self.eng.config.msig = self.isodeceng.config.maxshift
        self.eng.config.csig = self.isodeceng.config.css_thresh
        self.eng.config.adductmass = self.isodeceng.config.adductmass
        self.eng.config.exthresh = self.isodeceng.config.datathreshold
        self.eng.config.integratelb = self.isodeceng.config.mzwindow[0]
        self.eng.config.integrateub = self.isodeceng.config.mzwindow[1]

    def on_unidec_button(self, e=None):
        """
        Run IsoDec
        :param e: unused event
        :return: None
        """
        self.view.export_gui_to_config()

        if self.eng.config.selection_type == "Peptide" or self.eng.config.selection_type == 'Rna':
            self.isodeceng.selection_type = self.eng.config.selection_type

        tstart = time.perf_counter()
        self.view.SetStatusText("Running IsoDec...", number=5)
        self.fix_parameters()
        self.export_config(self.eng.config.confname)
        self.translate_config()

        self.isodeceng.config.activescan = 1
        self.isodeceng.config.activescanorder = 2
        self.isodeceng.config.activescanrt = 0

        data = self.eng.data.data2
        # Run IsoDec
        self.isodeceng.batch_process_spectrum(data, centroided=True, refresh=True, type=self.isodeceng.selection_type)
        # Convert to mass
        self.eng.data.massdat = self.isodeceng.pks_to_mass(self.eng.config.massbins)

        if len(self.isodeceng.pks.masses) == 0:
            print("No Peaks Found")
            self.view.SetStatusText("No Peaks Found", number=5)
            return
        # Translate Pks
        self.translate_pks()

        self.import_config()
        self.view.clear_all_plots()
        self.makeplot1()
        self.makeplot2()
        self.view.SetStatusText("IsoDec Done: " + str(len(self.isodeceng.pks.peaks)), number=5)
        tend = time.perf_counter()
        print("IsoDec Done. Time: %.2gs" % (tend - tstart))
        self.export_results()
        pass

    def translate_pks(self):
        idpks = self.isodeceng.pks
        if len(idpks.masses) == 0:
            print("No Peaks Found")
            self.view.SetStatusText("No Peaks Found", number=5)
            return
        udpks = Peaks()
        udpks.merge_isodec_pks(idpks, self.eng.config.massbins, self.eng.config)
        self.eng.pks = udpks
        # Send pks to structure

    def on_full(self, evt=None):
        maxmz = np.amax(self.eng.data.rawdata[:, 0])
        minmz = np.amin(self.eng.data.rawdata[:, 0])
        self.view.controls.ctlminmz.SetValue(str(minmz))
        self.view.controls.ctlmaxmz.SetValue(str(maxmz))
        self.on_dataprep_button()

    def on_plot_peaks(self, evt=None):
        self.view.export_gui_to_config()
        print("Plotting Peaks")
        self.translate_pks()
        self.plot_mass_peaks()
        self.plot_mz_peaks()
        if self.eng.config.avgpeakmasses == 1:
            self.view.peakpanel.add_data(self.eng.pks, show="mass", collab1="Avg Mass")
            self.isodeceng.showavg = True
        else:
            self.view.peakpanel.add_data(self.eng.pks, show="zs")

    def plot_mass_peaks(self, evt=None):
        print("Plotting Mass Peaks")
        self.makeplot2()
        for p in self.eng.pks.peaks:
            norm = self.eng.pks.norm
            if not p.ignore:
                newmass = ud.round_to_nearest(p.mass, self.eng.config.massbins)
                self.view.plot2.plotadddot(newmass, p.height * norm, p.color, p.marker)
        self.view.plot2.repaint()

    def plot_mz_peaks(self, evt=None):
        print("Plotting m/z Peaks")
        self.makeplot1()
        for p in self.eng.pks.peaks:
            if not p.ignore:
                for d in p.mztab:
                    self.view.plot1.plotadddot(d[0], d[1], p.color, p.marker)
        self.view.plot1.repaint()

    def on_plot_dists(self, evt=None):
        print("Plotting Distributions")
        for p in self.eng.pks.peaks:
            if not p.ignore:
                isodist = deepcopy(p.stickdat)
                if len(isodist) == 0:
                    continue
                # elif len(isodist) == 1:
                #     isodist = isodist[0]
                # else:
                #     isodist = np.concatenate(isodist, axis=0)
                #     pass
                isodist[:, 1] = isodist[:, 1] * -1
                self.view.plot1.add_centroid(isodist, color=p.color, repaint=False)
        self.view.plot1.repaint()
        pass

    def on_delete(self, evt=None):
        self.plot_mass_peaks()
        self.plot_mz_peaks()

    def on_replot(self, evt=None):
        self.view.export_gui_to_config()
        self.on_plot_peaks()

    def on_charge_states(self, e=None, mass=None, plot=None, peakpanel=None, data=None):
        """
        Triggered by right click "plot charge states" on self.view.peakpanel.
        Plots a line with text listing the charge states of a specific peak.
        :param e: unused event
        :return: None
        """
        if plot is None:
            plot = self.view.plot1

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

        for p in self.isodeceng.pks.masses:
            if p.monoiso == peaksel:
                mzs = p.mzs
                zs = p.zs
                ints = p.mzints

                print("Showing Charge States for Mass: ", peaksel)
                print(zs)

                plot.textremove()
                if data is None:
                    data = self.eng.data.data2
                for i, mz in enumerate(mzs):
                    z = zs[i]
                    int = ints[i]
                    plot.addtext(sign + str(z), mz, int + np.amax(data[:, 1]) * 0.06, nopaint=True)
                plot.repaint()

                break

    def on_label_integral(self, e=None, peakpanel=None, pks=None, plot=None, dataobj=None):
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
        if ud.isempty(peaksel):
            peaksel = pmasses
        pint = np.array([p.height for p in pks.peaks])
        mval = np.amax(dataobj.massdat[:, 1])

        plot.textremove()
        for i, d in enumerate(pmasses):
            if d in peaksel:
                label = str(np.round(pint[i], 2))
                plot.addtext(label, pmasses[i], mval * 0.06 + pint[i], vlines=False, nopaint=True)
        plot.repaint()

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
                self.on_open_file(fname, newdir, refresh=True)
            else:
                print("Paste failed, got: ", data)
        except Exception as e:
            print(e)
            wx.MessageBox("Unable to open the clipboard", "Error")

    def on_batch(self, e=None):
        print("Batch Processing")
        dlg = wx.FileDialog(self.view, "Choose a data file mzML or Thermo Raw format", '', "", "*.*")
        if dlg.ShowModal() == wx.ID_OK:
            self.view.SetStatusText("Opening", number=5)
            path = dlg.GetPath()
            self.batch_process(path)
        dlg.Destroy()

    def batch_process(self, path):
        """
        Batch process a file
        :param path: File name
        :return: None
        """
        print("Batch Processing:", path)
        self.isodeceng.pks = None
        self.view.export_gui_to_config()
        self.isodeceng.process_file(path)
        # The output directory should be a directory with the same name as the input file + _unidecfiles
        outdir = os.path.dirname(path) + "\\" + os.path.splitext(os.path.basename(path))[0] + "_unidecfiles\\"
        # Check of the outdirectory exists
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Get the filename wihtout the path or extension
        filename = os.path.splitext(os.path.basename(path))[0]
        os.chdir(os.path.dirname(outdir))

        if self.eng.config.poolflag == 1:
            self.isodeceng.export_peaks("msalign", filename, self.isodeceng.config)

        if self.eng.config.compressflag == 1:
            self.isodeceng.export_peaks("tsv", filename + ".tsv")
        pass

    def export_results(self):
        if self.eng.config.poolflag == 1:
            self.isodeceng.export_peaks("msalign", filename=self.eng.config.outfname,
                                        reader=self.isodeceng.reader, max_precursors=1)
            print("Exported MSAlign File: ", self.eng.config.outfname)

        if self.eng.config.compressflag == 1:
            self.isodeceng.export_peaks("tsv", filename=self.eng.config.outfname + ".tsv")
            print("Exported TSV File: ", self.eng.config.outfname + ".tsv")

    def fix_parameters(self):
        print("Checking for Bad Parameters")
        self.view.export_gui_to_config()
        found_problem = False
        found_warning = False

        # min and maxmz cannot be negative
        try:
            if self.eng.config.minmz < 0:
                print("Min m/z cannot be negative")
                self.eng.config.minmz = ""
                found_problem = True
        except:
            pass
        try:
            if self.eng.config.maxmz < 0:
                print("Max m/z cannot be negative")
                self.eng.config.maxmz = ""
                found_problem = True
        except:
            pass

        # Data reduction should be between 0 and 100
        if self.eng.config.reductionpercent < 0:
            print("Data Reduction cannot be negative")
            self.eng.config.reductionpercent = 0
            found_problem = True
        elif self.eng.config.reductionpercent > 100:
            print("Data Reduction cannot be greater than 100")
            self.eng.config.reductionpercent = 80
            found_problem = True

        # Baseline subtraction should not be less than 0
        if self.eng.config.subbuff < 0:
            print("Baseline Subtraction cannot be negative")
            self.eng.config.subbuff = 0
            found_problem = True

        # Mass bins size should not be 0 or negative
        if self.eng.config.massbins <= 0:
            print("Mass Bins cannot be negative or 0")
            self.eng.config.massbins = 1
            found_problem = True

        # CSS should be between 0 and 1 inclusive
        if self.eng.config.csig < 0:
            print("CSS Threshold cannot be negative")
            self.eng.config.csig = 0
            found_problem = True
        elif self.eng.config.csig > 1:
            print("CSS Threshold cannot be greater than 1")
            self.eng.config.csig = 0.99
            found_problem = True
        elif self.eng.config.csig < 0.25:
            print("CSS Threshold is pretty low. Are you sure about this?")
            found_warning = True

        # Match Tolerance should be greater than 0
        if self.eng.config.filterwidth <= 0:
            print("Match Tolerance cannot be negative")
            self.eng.config.filterwidth = 5.0
            found_problem = True
        # Warn if it is greater than 1000
        if self.eng.config.filterwidth > 1000:
            print("Match Tolerance is very high. Are you sure? This is ppm we are talking about...")
            found_warning = True

        # Max shift should not be negative
        if self.eng.config.msig < 0:
            print("Max Shift cannot be negative")
            self.eng.config.msig = 0
            found_problem = True
        elif 10 > self.eng.config.msig > 5:
            print("Max Shift is very high. Are you sure?")
            found_warning = True
        elif self.eng.config.msig > 8:
            print("Max Shift is too high. We're concerned about you. "
                  "If you need us to let it go this high, send an email to the developers.")
            self.eng.config.msig = 3
            found_problem = True

        # Knockdown rounds should not be 0 or lower
        if self.eng.config.numit <= 0:
            print("Knockdown Rounds cannot be 0 or negative")
            self.eng.config.numit = 1
            found_problem = True
        elif 20 > self.eng.config.numit > 10:
            print("Knockdown Rounds is very high. Are you sure?")
            found_warning = True
        elif self.eng.config.numit > 20:
            print("Knockdown Rounds is too high. We're concerned about you. "
                  "If you need us to let it go this high, send an email to the developers.")
            self.eng.config.numit = 5
            found_problem = True

        # Adduct mass should be close to 1.007276 or its negative
        if abs(abs(self.eng.config.adductmass) - 1.007276) > 0.01:
            print("Adduct Mass is not close to 1.007276. Are you sure?")
            found_warning = True

        # Data encoding threshold should be between 0 and 1 inclusive
        if self.eng.config.exthresh < 0:
            print("Data Encoding Threshold cannot be negative")
            self.eng.config.exthresh = 0
            found_problem = True
        elif self.eng.config.exthresh > 1:
            print("Data Encoding Threshold cannot be greater than 1")
            self.eng.config.exthresh = 0.05
            found_problem = True
        elif self.eng.config.exthresh > 0.5:
            print("Data Encoding Threshold may be too high. Are you sure about this?")
            found_warning = True

        # Warn if either the mz window values are greater than 4
        if np.abs(self.eng.config.integratelb) > 4:
            print("Lower m/z window is very high. Are you sure?")
            found_warning = True
        if np.abs(self.eng.config.integrateub) > 4:
            print("Upper m/z window is very high. Are you sure?")
            found_warning = True

        # Peak detection window needs to be above 3
        if self.eng.config.peakwindow < 3:
            print("Peak Detection Window needs to be at least 3")
            self.eng.config.peakwindow = 3
            found_problem = True

        # Peak detection threshold needs to be above 0
        if self.eng.config.peakthresh < 0:
            print("Peak Detection Threshold needs to be 0 or above")
            self.eng.config.peakthresh = 0.0000001
            found_problem = True

        number = 1
        self.view.SetStatusText("", number=number)
        if found_problem:
            print("Problems were found in the parameters. Values have been automatically reset. "
                  "Please check the output above.")
            self.view.SetStatusText("Problems found in parameters. Values have been reset.", number=number)
            self.view.import_config_to_gui()

        if found_warning:
            print("Warnings were found in the parameters. Please check the output above.")
            self.view.SetStatusText("Warnings found in parameters. Check the console output.", number=number)




if __name__ == "__main__":
    app = IsoDecPres()
    app.start()

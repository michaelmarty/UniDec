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
import platform
from copy import deepcopy


class IsoDecPres(UniDecPres):
    def __init__(self):
        super().__init__()
        self.isodeceng = IsoDecRuntime()
        self.eng = UniDec()
        self.eng.config.peakwindow = self.isodeceng.config.peakwindow
        self.eng.config.peakthresh = self.isodeceng.config.peakthresh
        self.eng.config.massbins = 1
        self.eng.config.numit = self.isodeceng.config.knockdown_rounds

        self.view = IsoDecView(self, "IsoDec", self.eng.config, iconfile=None)
        if True:
            if platform.node() == 'MTM-VOSTRO':
                self.on_ex()
                self.on_unidec_button()
                self.on_plot_peaks()
                self.on_plot_dists()
        # except:
        #     pass

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
            self.eng.open_file(filename, directory, refresh=refresh, **kwargs)

        # Set Status Bar Text Values
        self.view.SetStatusText("File: " + filename, number=1)
        # print self.view.imflag, self.eng.config.imflag

        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        #self.view.SetStatusText("R\u00B2 ", number=3)

        # Plot 1D
        if self.eng.config.batchflag == 0:
            self.makeplot1(imfit=False)

        self.import_config()
        self.view.SetStatusText("Ready", number=5)

        self.write_to_recent()
        self.view.menu.update_recent()

    def makeplot1(self, e=None, intthresh=False, imfit=True):
        """
        Plot data and fit in self.view.plot1 and optionally in plot1fit
        :param e: unused event
        :return: None
        """
        self.view.plot1.centroid_plot(self.eng.data.data2, xlabel="m/z", ylabel="Intensity", color="k")
        #self.eng.makeplot1(plot=self.view.plot1, intthresh=intthresh, imfit=imfit)

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
        tstart = time.perf_counter()
        self.view.SetStatusText("Data Prep", number=5)
        self.export_config(self.eng.config.confname)
        self.eng.process_data()
        if not self.eng.config.centroided:
            self.eng.data.data2 = get_all_centroids(self.eng.data.data2)

        self.import_config()
        self.view.clear_all_plots()
        self.makeplot1(imfit=False, intthresh=True)
        self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        self.view.SetStatusText("Data Prep Done", number=5)
        tend = time.perf_counter()
        print("Data Prep Done. Time: %.2gs" % (tend - tstart))
        pass

    def on_unidec_button(self, e=None):
        """
        Run IsoDec
        :param e: unused event
        :return: None
        """
        tstart = time.perf_counter()
        self.view.SetStatusText("Running IsoDec...", number=5)
        self.export_config(self.eng.config.confname)
        self.isodeceng.config.peakwindow = self.eng.config.peakwindow
        self.isodeceng.config.peakthresh = self.eng.config.peakthresh
        self.isodeceng.config.knockdown_rounds = self.eng.config.numit

        data = self.eng.data.data2
        # Run IsoDec
        self.isodeceng.batch_process_spectrum(data, centroided=True)
        # Convert to mass
        self.eng.data.massdat = self.isodeceng.pks_to_mass(self.eng.config.massbins)
        # Translate Pks
        self.translate_pks()


        self.import_config()
        self.view.clear_all_plots()
        self.makeplot1()
        self.makeplot2()
        self.view.SetStatusText("IsoDec Done: " + str(len(self.isodeceng.pks.masses)), number=5)
        tend = time.perf_counter()
        print("IsoDec Done. Time: %.2gs" % (tend - tstart))
        pass

    def translate_pks(self):
        idpks = self.isodeceng.pks
        udpks = Peaks()
        udpks.merge_isodec_pks(idpks, self.eng.config.massbins, self.eng.config)
        self.eng.pks = udpks
        # Send pks to structure
        self.view.peakpanel.add_data(self.eng.pks, show="zs")


    def on_raw_open(self, evt=None):
        pass

    def on_paste_spectrum(self, evt=None):
        pass

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
                elif len(isodist) == 1:
                    isodist = isodist[0]
                else:
                    isodist = np.concatenate(isodist, axis=0)
                    pass
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

if __name__ == "__main__":
    app = IsoDecPres()
    app.start()
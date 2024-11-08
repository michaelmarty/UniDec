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

        self.on_ex()


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
        self.view.SetStatusText("R\u00B2 ", number=3)

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
        # Send pks to structure
        self.view.peakpanel.add_data(self.eng.pks)

        self.import_config()
        self.view.clear_all_plots()
        self.makeplot1()
        self.makeplot2()
        self.view.SetStatusText("IsoDec Done: " + str(len(self.isodeceng.pks)), number=5)
        tend = time.perf_counter()
        print("IsoDec Done. Time: %.2gs" % (tend - tstart))
        pass

    def translate_pks(self):
        idpks = self.isodeceng.pks
        udpks = Peaks()
        udpks.merge_isodec_pks(idpks, self.eng.config.massbins, self.eng.config)
        self.eng.pks = udpks

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




if __name__ == "__main__":
    app = IsoDecPres()
    app.start()
from UniDecCD import UniDecCDApp
import multiprocessing
import modules.HTEng as HTEng
from unidec.modules.gui_elements import CDWindow
from pubsub import pub
import wx
import unidec.tools as ud
import numpy as np
from unidec.modules import AutocorrWindow
from unidec.modules.unidecstructure import ChromatogramContainer
import os


class UniDecHTCDApp(UniDecCDApp):
    """
    Main UniDec GUI Application.
    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def init(self, *args, **kwargs):
        """
        Initialize Engine and View. Load defaults.
        :param args:
        :param kwargs:
        :return:
        """
        self.eng = HTEng.UniDecCDHT()
        self.cc = ChromatogramContainer()
        self.showht = False
        self.cycol = ud.create_color_cycle("bgcmy")

        self.view = CDWindow.CDMainwindow(self, "UniDecHT for HT-CD-MS Data",
                                          self.eng.config, htmode=True)
        self.comparedata = None

        pub.subscribe(self.on_select_mzz_region, 'mzlimits')
        pub.subscribe(self.on_smash, 'smash')

        self.eng.config.recentfile = self.eng.config.recentfileCD
        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)
        self.view.menu.update_recent()

        self.on_load_default(0)

        if "path" in kwargs:
            newdir, fname = os.path.split(kwargs["path"])
            self.on_open_file(fname, newdir)
            # self.on_dataprep_button(0)
            # self.on_auto(0)

        if self.infile is not None:
            newdir, fname = os.path.split(self.infile)
            self.on_open_file(fname, newdir)
            # self.on_dataprep_button(0)
            # self.on_auto(0)

        if True:  # and platform.node() == 'DESKTOP-08TGCJO':
            print("Opening Test File")
            path = ("Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM\\"
                    "20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56.dmt")

            self.on_open_file(None, None, path=path)

    def on_open_file(self, filename, directory, path=None, refresh=False):
        """
        Opens a file. Run self.eng.open_file.
        :param filename: File name
        :param directory: Directory containing file
        :param path: Full path to file
        :param refresh: Refresh the data ranges from the file. Default False.
        :return: None
        """
        self.cc.clear()
        self.on_open_cdms_file(filename, directory, path=path, refresh=refresh)
        self.load_chroms(self.eng.config.cdchrom)

    def on_dataprep_button(self, e=None):
        """
        Run data preparation. Run self.eng.process_data_scans.
        :param e: Unused event
        :return: None
        """
        self.cc.clear()
        self.dataprep()
        self.make_tic_plot()

    def on_select_mzz_region(self):
        self.export_config(self.eng.config.confname)
        if not wx.GetKeyState(wx.WXK_CONTROL):
            xlimits = self.view.plot1.subplot1.get_xlim()
            ylimits = self.view.plot1.subplot1.get_ylim()
            print("New limits:", xlimits, ylimits)
            self.view.plot1.reset_zoom()
            color = next(self.cycol)
            if self.showht:
                self.run_eic_ht(xlimits, ylimits, color=color)
            else:
                self.add_eic(xlimits, ylimits, color=color)

    def make_tic_plot(self):
        data = self.eng.get_tic(normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(data, color="black", label="TIC")

        self.showht = False
        self.plot_chromatograms(save=False)

    def add_eic(self, mzrange, zrange, color='b'):
        eicdata = self.eng.get_eic(mzrange, zrange, normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(eicdata, color=color, zrange=zrange, mzrange=mzrange)

        self.plot_chromatograms()

    def on_run_tic_ht(self, e=None):
        self.export_config(self.eng.config.confname)
        self.run_tic_ht()

    def on_run_eic_ht(self, e=None):
        self.export_config(self.eng.config.confname)
        for c in self.cc.chromatograms:
            if "TIC" not in c.label:
                if "HT" not in c.label:
                    self.run_eic_ht(c.mzrange, c.zrange, color=c.color, add_eic=False, plot=False)
        self.plot_chromatograms()

    def run_tic_ht(self):
        self.cc.clear()
        data = self.eng.get_tic(normalize=self.eng.config.datanorm)
        htdata = self.eng.tic_ht(normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(data, color="black", label="TIC")
        self.cc.add_chromatogram(htdata, color="red", label="TIC HT", ht=True)

        self.showht = True
        self.plot_chromatograms(save=False)

    def run_eic_ht(self, mzrange, zrange, color='b', add_eic=True, plot=True):
        htdata, eicdata = self.eng.eic_ht(mzrange, zrange, normalize=self.eng.config.datanorm)
        if add_eic:
            self.cc.add_chromatogram(eicdata, color=color, zrange=zrange, mzrange=mzrange)
        self.cc.add_chromatogram(htdata, color=color, zrange=zrange, mzrange=mzrange, ht=True)
        if plot:
            self.plot_chromatograms()

    def plot_chromatograms(self, e=None, save=True):
        self.makeplot1()
        self.view.plottic.clear_plot()
        self.view.chrompanel.list.populate(self.cc)
        for c in self.cc.chromatograms:
            if c.ignore:
                continue

            if self.eng.config.HTxaxis == "Scans":
                xdat = self.eng.fullscans
                xlab = "Scan Number"
            else:
                xdat = c.chromdat[:, 0]
                xlab = "Time"

            if not self.view.plottic.flag:
                self.view.plottic.plotrefreshtop(xdat, c.chromdat[:, 1], config=self.eng.config,
                                                 zoomout=True, label=c.label, xlabel=xlab,
                                                 color=c.color, nopaint=True)
            else:
                self.view.plottic.plotadd(xdat, c.chromdat[:, 1], colval=c.color, nopaint=True,
                                          newlabel=c.label)

            xlimits = c.mzrange
            ylimits = c.zrange
            if xlimits[0] != -1 and ylimits[0] != -1 and xlimits[1] != -1 and ylimits[1] != -1:
                self.view.plot1.add_rect(xlimits[0], ylimits[0], xlimits[1] - xlimits[0], ylimits[1] - ylimits[0],
                                         edgecolor=c.color, facecolor=c.color, nopaint=True)

        self.view.plottic.add_legend()
        self.view.plottic.repaint()
        self.view.plot1.repaint()
        if save:
            self.save_chroms()

    def save_chroms(self):
        chromarray = self.cc.to_array()
        np.savetxt(self.eng.config.cdchrom, chromarray, delimiter="\t", fmt="%s")

    def load_chroms(self, fname=None, array=None):
        if fname is not None:
            if not os.path.isfile(fname):
                print("Chrom file not found:", fname)
                return
            array = np.loadtxt(fname, delimiter="\t", dtype=str)

        if ud.isempty(array):
            return
        if len(array.shape) == 1:
            array = np.reshape(array, (1, len(array)))  # Make sure it is 2D
        print(array)
        for a in array:
            print(a)
            label = a[0]
            color = a[1]
            index = a[2]
            ht = a[8]
            mzrange = [float(a[3]), float(a[4])]
            zrange = [float(a[5]), float(a[6])]
            print(mzrange, zrange)
            if label == "TIC" or "HT" in label:
                continue
            chromdat = self.eng.get_eic(mzrange, zrange, normalize=self.eng.config.datanorm)

            self.cc.add_chromatogram(chromdat, color=color, label=label, ht=ht, zrange=zrange, mzrange=mzrange)

        self.plot_chromatograms()

    def on_ignore_repopulate(self, e=None):
        self.cc = self.view.chrompanel.list.get_data()
        self.plot_chromatograms()

    def on_run_all_ht(self, e=None):
        self.run_all_ht()

    def run_all_ht(self):
        self.eng.process_data_scans()
        self.eng.run_ht()

    def on_auto_set_ct(self, e=None):
        self.eng.get_cycle_time()
        self.eng.config.HTcycleindex = self.eng.cycleindex
        self.import_config()

    def on_autocorr2(self, index):
        """
        Manual Test - Passed
        :param index:
        :return:
        """
        data = self.cc.chromatograms[index].chromdat
        data[:, 0] = np.arange(0, len(data))
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, data, window=10)
        dlg.ShowModal()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = UniDecHTCDApp()
    app.start()

from UniDecCD import UniDecCDApp
import multiprocessing
import modules.HTEng as HTEng
from unidec.modules.gui_elements import CDWindow
from pubsub import pub
import wx
import unidec.tools as ud
import numpy as np

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
        self.on_open_cdms_file(filename, directory, path=path, refresh=refresh)
        self.eng.process_data_scans()
        self.make_tic_plot()

    def on_select_mzz_region(self):
        try:
            if not wx.GetKeyState(wx.WXK_CONTROL):
                xlimits = self.view.plot1.subplot1.get_xlim()
                ylimits = self.view.plot1.subplot1.get_ylim()
                print("New limits:", xlimits, ylimits)
                color = next(self.cycol)
                self.run_eic_ht(xlimits, ylimits, color=color)

                self.view.plot1.add_rect(xlimits[0], ylimits[0], xlimits[1] - xlimits[0], ylimits[1] - ylimits[0],
                                         edgecolor=color, facecolor=color)

        except Exception as e:
            print(e)

    def make_tic_plot(self):
        data = self.eng.get_tic()
        self.view.plot3.plotrefreshtop(data[:, 0], data[:, 1], config=self.eng.config, zoomout=True)

    def on_run_tic_ht(self, e=None):
        self.run_tic_ht()

    def run_tic_ht(self):
        self.make_tic_plot()
        data = self.eng.tic_ht(normalize=False)
        self.view.plot3.plotadd(data[:, 0], data[:, 1], colval="red", nopaint=False)

    def run_eic_ht(self, mzrange, zrange, color='b'):
        label = "m/z: " + str(round(mzrange[0])) + "-" + str(round(mzrange[1])) + \
                " z: " + str(round(zrange[0])) + "-" + str(round(zrange[1]))

        htdata, eicdata = self.eng.eic_ht(mzrange, zrange, normalize=False)
        self.view.plot3.plotadd(htdata[:, 0], htdata[:, 1], colval=color, nopaint=False)
        self.view.plot3.addtext(label, np.amax(self.eng.fulltime) * 0.8, np.amax(htdata[:, 1]), color=color,
                                vlines=False, hlines=False)
        pass

    def on_run_all_ht(self, e=None):
        self.run_all_ht()

    def run_all_ht(self):
        self.eng.run_ht()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = UniDecHTCDApp()
    app.start()

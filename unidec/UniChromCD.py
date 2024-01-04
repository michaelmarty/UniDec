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
import os, time


class UniChromCDApp(UniDecCDApp):
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

        self.view = CDWindow.CDMainwindow(self, "UniChrom for CD-MS Data",
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

        if False:  # and platform.node() == 'DESKTOP-08TGCJO':
            print("Opening Test File")
            path = ("Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM\\"
                    "20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56.dmt")

            self.on_open_file(None, None, path=path)
            # self.eng.process_data_scans()
            # self.make_cube_plot()
            # self.make_mass_time_2dplot()
            # self.run_all_ht()
            # self.run_all_mass_transform()
            # self.make_mass_cube_plot()

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
        Run data preparation. Makes the TIC and clears the chromatogram list.
        :param e: Unused event
        :return: None
        """
        self.cc.clear()
        self.dataprep()
        self.make_tic_plot()

    def on_pick_peaks(self, e=None):
        """
        Pick peaks. Runs the peak picking algorithm. Adds Mass EIC beyond conventional UCD.
        :param e: event (unused)
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.pick_peaks()
        if self.eng.fullmstack is not None:
            for p in self.eng.pks.peaks:
                mass = p.mass
                try:
                    ub = float(self.eng.config.integrateub)
                except:
                    ub = float(self.eng.config.peakwindow)
                try:
                    lb = float(self.eng.config.integratelb)
                except:
                    lb = -float(self.eng.config.peakwindow)

                massrange = [mass + lb, mass + ub]

                chromdat = self.eng.get_mass_eic(massrange)
                label = "Mass EIC: " + str(mass)
                self.cc.add_chromatogram(chromdat, color=p.color, label=label, ht=False)
                if self.eng.fullmstack_ht is not None:
                    htdata = self.eng.get_mass_eic(massrange, ht=True)
                    self.cc.add_chromatogram(htdata, color=p.color, label=label + " HT", ht=True)
        self.plot_chromatograms()

    def on_select_mzz_region(self):
        """
        Trigged by right click of m/z vs z 2D plot. Triggers EIC creation and HT if already setup.
        :return: None
        """
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
        """
        Make the TIC Plot
        :return: None
        """
        data = self.eng.get_tic(normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(data, color="black", label="TIC")

        self.showht = False
        self.plot_chromatograms(save=False)

    def add_eic(self, mzrange, zrange, color='b'):
        """
        Add an EIC to the list of chromatograms. Plot the chromatograms.
        :param mzrange: M/z range
        :param zrange: charge range
        :param color: Color of the plot. Default blue.
        :return: None
        """
        eicdata = self.eng.get_eic(mzrange, zrange, normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(eicdata, color=color, zrange=zrange, mzrange=mzrange)

        self.plot_chromatograms()

    def on_run_tic_ht(self, e=None):
        """
        Button to run TIC HT. Runs run_tic_ht.
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.run_tic_ht()

    def on_run_eic_ht(self, e=None):
        """
        Button to run EIC HT. Runs run_eic_ht on each EIC that is extracted.
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        for c in self.cc.chromatograms:
            if "TIC" not in c.label:
                if "HT" not in c.label:
                    self.run_eic_ht(c.mzrange, c.zrange, color=c.color, add_eic=False, plot=False)
        self.plot_chromatograms()

    def run_tic_ht(self):
        """
        Runs HT on the TIC. Clears the chromatogram list. Runs HT on TIC. Adds both to list. Plots chromatograms.
        :return: None
        """
        self.cc.clear()
        data = self.eng.get_tic(normalize=self.eng.config.datanorm)
        htdata = self.eng.tic_ht(normalize=self.eng.config.datanorm)
        self.cc.add_chromatogram(data, color="black", label="TIC")
        self.cc.add_chromatogram(htdata, color="red", label="TIC HT", ht=True)

        self.showht = True
        self.plot_chromatograms(save=False)

    def run_eic_ht(self, mzrange, zrange, color='b', add_eic=True, plot=True):
        """
        Function to generate and HT an EIC. Adds both to list. Plots chromatograms.
        :param mzrange: M/z range
        :param zrange: Charge range
        :param color: Color of plot. Default blue.
        :param add_eic: Add the EIC to the list. Default True.
        :param plot: Plot the chromatograms. Default True.
        :return: None
        """

        htdata, eicdata = self.eng.eic_ht(mzrange, zrange, normalize=self.eng.config.datanorm)
        if add_eic:
            self.cc.add_chromatogram(eicdata, color=color, zrange=zrange, mzrange=mzrange)
        self.cc.add_chromatogram(htdata, color=color, zrange=zrange, mzrange=mzrange, ht=True)
        if plot:
            self.plot_chromatograms()

    def plot_chromatograms(self, e=None, save=True):
        """
        Plot all chromatograms in the list. Add a box to the 2D plot if the chromatogram is an EIC.
        :param e: Unused event
        :param save: Whether to save the output to a text file. Default True
        :return: None
        """
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
        """
        Save the chromatograms to a text file.
        :return: None
        """
        chromarray = self.cc.to_array()
        np.savetxt(self.eng.config.cdchrom, chromarray, delimiter="\t", fmt="%s")

    def load_chroms(self, fname=None, array=None):
        """
        Load chromatograms from a text file. Add them to the list. Plot the chromatograms.
        :param fname: Input file name
        :param array: Input array. Not currently used, but could be used to load from a numpy array rather than file.
        :return: None
        """
        if fname is not None:
            if not os.path.isfile(fname):
                print("Chrom file not found:", fname)
                return
            else:
                print("Loading Chroms:", fname)
            array = np.loadtxt(fname, delimiter="\t", dtype=str)
        elif array is None:
            print("No file name or array supplied to load_chroms")
            return

        if ud.isempty(array):
            return
        if len(array.shape) == 1:
            array = np.reshape(array, (1, len(array)))  # Make sure it is 2D

        for a in array:
            label = a[0]
            color = a[1]
            index = a[2]
            ht = a[8]
            mzrange = [float(a[3]), float(a[4])]
            zrange = [float(a[5]), float(a[6])]

            if "TIC" in label or "HT" in label or "Mass EIC" in label:
                continue
            chromdat = self.eng.get_eic(mzrange, zrange, normalize=self.eng.config.datanorm)

            self.cc.add_chromatogram(chromdat, color=color, label=label, ht=ht, zrange=zrange, mzrange=mzrange)

        self.plot_chromatograms()

    def on_load_chroms(self, e=None):
        """
        Load chromatograms from a file. Triggered by a button.
        :param e: Unused event
        :return: None
        """
        # Load File dialog
        dlg = wx.FileDialog(self.view, "Choose a file")
        if dlg.ShowModal() == wx.ID_OK:
            fname = dlg.GetPath()
            self.load_chroms(fname=fname)
        dlg.Destroy()

    def on_ignore_repopulate(self, e=None):
        """
        Event triggered by updates to the chromatogram list. Repopulate the list based on adjusted ignore flags.
        :param e: Unused event
        :return: None
        """
        self.cc = self.view.chrompanel.list.get_data()
        self.plot_chromatograms()

    def on_select_time_range(self, e=None):
        """
        Event triggered by a click on the TIC plot. Select a time range and create a 2D m/z vs z sum from
        that time range post-HT.
        :param e: Unused event
        :return: None
        """
        if not self.showht or self.eng.fullhstack_ht is None:
            return
        if not wx.GetKeyState(wx.WXK_CONTROL):
            self.plot_chromatograms()
            xlimits = self.view.plottic.subplot1.get_xlim()
            print("New limits:", xlimits)
            self.view.plottic.reset_zoom()
            ylimits = self.view.plottic.subplot1.get_ylim()
            # Plot Red box on plottic
            self.view.plottic.add_rect(xlimits[0], ylimits[0], xlimits[1] - xlimits[0], ylimits[1] - ylimits[0],
                                       edgecolor="red", facecolor="red", nopaint=False)
            self.select_ht_range(range=xlimits)

    def on_run_all_ht(self, e=None):
        """
        Button to trigger running HT on all data. Runs run_all_ht.
        :param e: Unused event.
        :return: None
        """
        self.run_all_ht()

    def run_all_ht(self):
        """
        Run HT on all data. Clears the chromatogram list. Runs HT on each data point. Plots chromatograms after
        adding extracted TIC*_HT.
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.eng.process_data_scans()
        ticdat = self.eng.run_all_ht()

        self.cc.add_chromatogram(ticdat, color="gold", label="TIC*_HT")

        # tic2 = np.sum(self.eng.fullhstack, axis=(1, 2))
        # ticdat2 = np.transpose(np.vstack((self.eng.fulltime, tic2)))
        # self.cc.add_chromatogram(ticdat2, color="red", label="TIC*")

        self.showht = True
        self.plot_chromatograms()

    def run_all_mass_transform(self, e=None):
        """
        Button to trigger running mass transform on all data. Runs transform_stacks. Creates mass TICs.
        :param e: Unused event.
        :return: None
        """
        self.eng.transform_stacks()
        self.cc.add_chromatogram(self.eng.mass_tic, color="grey", label="Mass TIC")
        if self.eng.fullmstack_ht is not None:
            self.cc.add_chromatogram(self.eng.mass_tic_ht, color="goldenrod", label="Mass TIC HT")
        self.plot_chromatograms()

    def select_ht_range(self, range=None):
        """
        Select a time range and create a 2D m/z vs z sum from that time range post-HT.
        :param range: Time range to select. Default None, which should be all times
        :return: None
        """
        self.eng.select_ht_range(range=range)
        self.makeplot1()
        self.makeplot2()
        self.makeplot3()
        self.makeplot4()

    def make_charge_time_2dplot(self, e=None):
        """
        Creates a 2D plot of charge vs time
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.make_2d_plot(face=1)

    def make_mz_time_2dplot(self, e=None):
        """
        Creates a 2D plot of m/z vs time
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.make_2d_plot(face=2)

    def make_mass_time_2dplot(self, e=None):
        """
        Creates a 2D plot of mass vs time
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.make_2d_plot(face=3)

    def make_2d_plot(self, e=None, face=1):
        """
        Creates 2D plots of charge, m/z, or mass vs time. General 2D plot function called by others
        :param e: Unused event
        :param face: Which face to plot against time. 1=charge, 2=m/z, 3=mass
        :return: None
        """
        starttime = time.perf_counter()
        print("Starting 2D Plots...", )
        if self.eng.fullhstack is None:
            print("Need to create fullhstack")
            self.eng.process_data_scans()

        if self.eng.fullmstack is None and face == 3:
            print("Need to create fullmstack")
            self.run_all_mass_transform()

        if self.eng.config.HTxaxis == "Scans":
            xdat = self.eng.fullscans
            xlab = "Scan Number"
        else:
            xdat = self.eng.fulltime
            xlab = "Time"

        if face == 1:
            grid = np.sum(self.eng.fullhstack, axis=2)
            y = self.eng.zaxis[1:]
            ylab = "Charge"
            discrete = True
        elif face == 2:
            grid = np.sum(self.eng.fullhstack, axis=1)
            y = self.eng.mzaxis[1:]
            ylab = "m/z (Th)"
            discrete = True
        elif face == 3:
            grid = self.eng.fullmstack
            print(np.shape(grid))
            y = self.eng.massaxis
            ylab = "Mass (Da)"
            discrete = True
        else:
            return

        self.view.plot7.contourplot(xvals=xdat, yvals=y, zgrid=grid,
                                    xlab=xlab, ylab=ylab, config=self.eng.config, discrete=discrete)

        if self.eng.fullhstack_ht is None:
            return
        else:
            fullhstack_ht = np.clip(self.eng.fullhstack_ht, 0, np.amax(self.eng.fullhstack_ht))
        if face == 1:
            grid = np.sum(fullhstack_ht, axis=2)
        elif face == 2:
            grid = np.sum(fullhstack_ht, axis=1)
        elif face == 3:
            if self.eng.fullmstack_ht is None:
                return
            else:
                grid = np.clip(self.eng.fullmstack_ht, 0, np.amax(self.eng.fullmstack_ht))
        else:
            return
        self.view.plot8.contourplot(xvals=xdat, yvals=y, zgrid=grid,
                                    xlab=xlab, ylab=ylab, config=self.eng.config, discrete=discrete)

        print("Finished 2D Plot", (time.perf_counter() - starttime), " s")

    def make_mass_cube_plot(self, e=None):
        """
        Creates a cube plot of mass vs charge vs time
        :param e: Unused event
        :return: None
        """
        self.make_cube_plot(e, mass=True)

    def make_cube_plot(self, e=None, mass=False):
        """
        Creates a cube plot of mass or m/z vs charge vs time. Trys both regular and HT if possible.
        :param e: Unused event
        :param mass: Flag to plot mass vs charge vs time. Default False, which plots m/z vs charge vs time.
        :return: None
        """
        self.export_config(self.eng.config.confname)

        if mass and self.eng.fullmstack is None:
            self.run_all_mass_transform()

        if self.eng.config.HTxaxis == "Scans":
            xdat = self.eng.fullscans
            xlab = "Scan Number"
        else:
            xdat = self.eng.fulltime
            xlab = "Time"

        try:
            print(np.shape(self.eng.fullhstack))
            if self.eng.fullhstack is None:
                self.eng.process_data_scans()
            starttime = time.perf_counter()
            face1 = np.sum(self.eng.fullhstack, axis=2).transpose()
            if mass:
                face3 = np.reshape(self.eng.data.massgrid, (len(self.eng.massaxis), len(self.eng.ztab)))
                ydat = self.eng.massaxis
                ylab = "Mass (Da)"
                face2 = self.eng.fullmstack.transpose()
            else:
                # mz by z
                ydat = self.eng.mzaxis[1:]
                ylab = "m/z (Th)"
                face3 = np.sum(self.eng.fullhstack, axis=0).transpose()
                face2 = np.sum(self.eng.fullhstack, axis=1).transpose()
            self.view.plot9.cubeplot(ydat, self.eng.zaxis[1:], xdat,
                                     face3, face2, face1,
                                     xlab=ylab, ylab="Charge", zlab=xlab,
                                     cmap=self.eng.config.cmap)
            endtime = time.perf_counter()
            print("Finished Raw Cube in: ", (endtime - starttime), " s")
        except Exception as ex:
            print("Failed Raw cube", ex)
            return
            pass

        try:
            if self.eng.fullhstack_ht is None:
                return
            if mass and self.eng.fullmstack_ht is None:
                return
            fullhstack_ht = np.clip(self.eng.fullhstack_ht, 0, np.amax(self.eng.fullhstack_ht))
            starttime = time.perf_counter()
            face1 = np.sum(fullhstack_ht, axis=2).transpose()

            if mass:
                fullmstack_ht = np.clip(self.eng.fullmstack_ht, 0, np.amax(self.eng.fullmstack_ht))
                face3 = np.reshape(self.eng.data.massgrid, (len(self.eng.massaxis), len(self.eng.ztab)))
                face2 = fullmstack_ht.transpose()
            else:
                face3 = np.sum(self.eng.fullhstack_ht, axis=0).transpose()
                face2 = np.sum(self.eng.fullhstack_ht, axis=1).transpose()

            self.view.plot10.cubeplot(ydat, self.eng.zaxis[1:], xdat,
                                      face3, face2, face1,
                                      xlab=ylab, ylab="Charge", zlab=xlab,
                                      cmap=self.eng.config.cmap)
            endtime = time.perf_counter()
            print("Finished HT Cube in: ", (endtime - starttime), " s")
        except Exception as ex:
            print("Failed HT cube", ex)
            pass

    def on_auto_set_ct(self, e=None):
        """
        Automatically set the cycle index based on the autocorrelation.
        :param e: Unused event
        :return: None
        """
        self.eng.get_cycle_time()
        self.eng.config.HTcycleindex = self.eng.cycleindex
        self.import_config()

    def on_plot_kernel(self, e=None):
        """
        Plot the HT kernel
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        if not self.showht:
            self.eng.setup_ht()
        data = self.eng.htkernel

        if self.eng.config.HTxaxis == "Scans":
            xdat = self.eng.fullscans
            xlab = "Scan Number"
        else:
            xdat = self.eng.fulltime
            xlab = "Time"

        label = "HT Kernel"
        color = "orange"
        self.indexrange = [self.eng.padindex - self.eng.shiftindex, len(xdat) - self.eng.shiftindex]
        xdat = xdat[self.indexrange[0]:self.indexrange[1]]
        data *= np.amax(self.eng.fulltic) / np.amax(data)
        if not self.view.plottic.flag:
            self.view.plottic.plotrefreshtop(xdat, data, config=self.eng.config,
                                             zoomout=True, label=label, xlabel=xlab,
                                             color=color, nopaint=False)
        else:
            self.view.plottic.plotadd(xdat, data, colval=color, nopaint=False,
                                      newlabel=label)
        self.view.plottic.add_legend()

    def on_autocorr2(self, index):
        """
        Launch an autocorrelation window for a given chromatogram.
        :param index: Index of chromatogram in the list
        :return: None
        """
        data = self.cc.chromatograms[index].chromdat
        data[:, 0] = np.arange(0, len(data))
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, data, window=10)
        dlg.ShowModal()

    def on_export_arrays(self, e=None):
        """
        Export the chromatograms and HT kernel to text files.
        :param e: Unused event
        :return: None
        """
        self.export_config(self.eng.config.confname)
        self.export_arrays()

    def export_arrays(self):
        """
        Export the chromatograms and HT kernel to text files.
        :return: None
        """
        # Export Kernel Array
        np.savetxt(self.eng.config.outfname + "_htkernel.txt", self.eng.htkernel)
        # Export Chromatograms
        for c in self.cc.chromatograms:
            if c.ignore:
                continue
            newlabel = c.label.replace("/", "")
            newlabel = newlabel.replace(" ", "_")
            newlabel = newlabel.replace(":", "_")
            newlabel = os.path.join(os.path.split(self.eng.config.outfname)[0], newlabel)
            np.savetxt(newlabel + "_chrom.txt", c.chromdat)
        print("Saved Files")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = UniChromCDApp()
    app.start()

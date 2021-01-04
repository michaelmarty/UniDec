import time
import os
import subprocess
import atexit
import wx
import wx.html
import numpy as np
from metaunidec import mudview
from metaunidec import mudeng, metafft
import time
#
from pubsub import pub

import unidec_modules.unidectools as ud
from unidec_modules import Extract2D, masstools, miscwindows, \
    MassDefects, PlotAnimations, IM_functions, fft_window, AutocorrWindow
from unidec_modules.isolated_packages import FileDialogs
import datacollector
import multiprocessing
from unidec_modules.unidec_presbase import UniDecPres
from metaunidec import ultrameta
from metaunidec.mudhelp import *
from metaunidec.meta_import_wizard.meta_import_wizard import ImportWizard
from unidec_modules.plot_waterfall import WaterfallFrame

# import FileDialog  # Needed for pyinstaller

__author__ = 'Michael.Marty'


class MetaUniDecBase(UniDecPres):
    """
    Class for MetaUniDec functions that might want to be shared with UniChrom
    """

    def __init__(self, *args, **kwargs):
        UniDecPres.__init__(self, *args, **kwargs)
        atexit.register(self.repack_hdf5)

    def makeplot1(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        spectra = self.eng.data.get_spectra()
        if self.eng.config.datanorm == 0:
            sep = self.eng.config.separation * self.eng.data.get_max_data2()
        else:
            sep = self.eng.config.separation
        try:
            if len(spectra) > int(self.eng.config.crossover):
                mult = int(len(spectra) / self.eng.config.numtot)
                self.view.SetStatusText("Displaying subset of data", number=2)
            else:
                mult = 1
        except:
            mult = 1
        for i, s in enumerate(spectra[::mult]):
            if i == 0:
                self.view.plot1.plotrefreshtop(s.data2[:, 0], s.data2[:, 1], title="Processed Data", xlabel="m/z (Th)",
                                               ylabel="Intensity", label=s.name, config=self.eng.config, color=s.color,
                                               nopaint=True)
            else:
                self.view.plot1.plotadd(s.data2[:, 0], s.data2[:, 1] - i * sep, colval=s.color,
                                        newlabel=s.name)
        self.view.plot1.repaint()
        try:
            self.view.SetStatusText("Data Length: " + str(len(self.eng.data.data2)), number=2)
        except:
            pass

    def makeplot2_mud(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        if self.eng.config.datanorm == 0:
            sep = self.eng.config.separation * self.eng.data.get_max_data2()
        else:
            sep = self.eng.config.separation

        spectra = self.eng.data.get_spectra()
        if len(spectra) > self.eng.config.crossover:
            mult = int(len(spectra) / self.eng.config.numtot)
            self.view.SetStatusText("Displaying subset of data", number=2)
        else:
            mult = 1
        for i, s in enumerate(spectra[::mult]):
            if not ud.isempty(s.massdat):
                if i == 0:
                    self.view.plot2.plotrefreshtop(s.massdat[:, 0], s.massdat[:, 1], title="Zero-Charge Mass Spectrum",
                                                   xlabel="Mass (Da)",
                                                   ylabel="Intensity", label=s.name, config=self.eng.config,
                                                   color=s.color,
                                                   test_kda=True,
                                                   nopaint=True)
                else:
                    self.view.plot2.plotadd(s.massdat[:, 0], s.massdat[:, 1] - i * sep,
                                            colval=s.color, newlabel=s.name)
        self.view.plot2.repaint()
        try:
            self.makeplot9()
        except:
            pass

    def makeplot7(self, fitgrid=None):
        """
        Tested
        :return:
        """
        try:
            if fitgrid is None:
                fitgrid = self.eng.data.exgrid
            self.view.plot7.clear_plot()
            if not ud.isempty(self.eng.data.exgrid):
                ignore = self.eng.data.get_bool()
                var1 = np.array(self.eng.data.var1)[ignore]
                try:
                    var1 = var1.astype(float)
                    xlabel = self.eng.data.v1name
                except:
                    var1 = np.arange(0, len(self.eng.data.var1))[ignore]
                    xlabel = "Index"

                ylabel = self.view.extractlabels[self.eng.config.exchoice]
                self.view.plot7.clear_plot()
                self.view.plot7._axes = [0.2, 0.12, 0.75, 0.8]
                for i, p in enumerate(self.eng.pks.peaks):
                    if p.ignore == 0:
                        color = p.color

                        if not self.view.plot7.flag:
                            self.view.plot7.plotrefreshtop(var1, fitgrid[i][ignore],
                                                           title="Extracted Data", xlabel=xlabel
                                                           , ylabel=ylabel, color=color, test_kda=False)
                            self.view.plot7.plotadddot(var1, self.eng.data.exgrid[i][ignore], color, p.marker)
                        else:
                            self.view.plot7.plotadd(var1, fitgrid[i][ignore], color)
                            self.view.plot7.plotadddot(var1, self.eng.data.exgrid[i][ignore], color, p.marker)
                if self.eng.config.exnorm == 1:
                    self.view.plot7.subplot1.set_ylim([0, 1])
                self.view.plot7.repaint()
        except Exception as e:
            print(e)

    def make_waterfall_plots(self, e=None):
        print("Making Waterfall Plots")
        wt = WaterfallFrame(None)
        wt.make_plot(self.eng.data)
        # wt.draw()

    def on_plot_composite(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.export_config()
        self.eng.data.import_grids_and_peaks()
        if not ud.isempty(self.eng.data.massdat):
            self.view.plot2.clear_plot()

            self.view.plot2.plotrefreshtop(self.eng.data.massdat[:, 0], self.eng.data.massdat[:, 1],
                                           title="Zero-Charge Mass Spectrum",
                                           xlabel="Mass (Da)",
                                           ylabel="Intensity", label="Sum", config=self.eng.config,
                                           color="black",
                                           test_kda=True,
                                           nopaint=True)

            self.view.plot2.repaint()
        pass

    def plot_sums(self):
        """
        Tested
        :return:
        """
        try:
            maxsums = np.amax(self.eng.data.massdat[:, 1])
            maxpeaks = np.amax([p.height for p in self.eng.pks.peaks])
            self.eng.pks.norm = maxsums / maxpeaks
        except:
            self.eng.pks.norm = 1

        if self.eng.config.datanorm == 0:
            sep = self.eng.config.separation * self.eng.data.get_max_data2()
        else:
            sep = self.eng.config.separation

        if not ud.isempty(self.eng.data.massdat):
            self.view.plot2.plotadd(self.eng.data.massdat[:, 0],
                                    self.eng.data.massdat[:, 1] + sep,
                                    colval="black", newlabel="Sum")
            if self.eng.pks.plen > 0:
                for p in self.eng.pks.peaks:
                    if p.ignore == 0:
                        self.view.plot2.plotadddot(p.mass, p.height * self.eng.pks.norm + sep,
                                                   p.color, p.marker)
            self.view.plot2.repaint()

    def make2dplots(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.view.SetStatusText("Making 2D Plots...", number=5)
        self.export_config()
        self.eng.data.import_grids_and_peaks()
        self.makeplot3()
        self.makeplot5()
        self.view.SetStatusText("2D Plots Complete", number=5)

    def makeplot3(self):
        """
        Tested
        :return:
        """
        self.view.plot3.clear_plot()
        if not ud.isempty(self.eng.data.mzgrid):
            tstart = time.perf_counter()
            ignore = self.eng.data.get_bool()
            self.view.plot3.contourplot(
                xvals=self.eng.data.mzdat[:, 0], yvals=np.array(self.eng.data.var1)[ignore],
                zgrid=self.eng.data.mzgrid[ignore, :, 1].transpose(),
                config=self.eng.config, title="m/z vs. " + self.eng.data.v1name, test_kda=False, xlab="m/z (Th)",
                ylab=self.eng.data.v1name)
            self.view.plot3.repaint()
            tend = time.perf_counter()
            print("Plot 3: %.2gs" % (tend - tstart))
            pass

    def makeplot5(self):
        """
        Tested
        :return:
        """
        self.view.plot5.clear_plot()
        if not ud.isempty(self.eng.data.massgrid):
            tstart = time.perf_counter()
            ignore = self.eng.data.get_bool()
            self.view.plot5.contourplot(
                xvals=self.eng.data.massdat[:, 0], yvals=np.array(self.eng.data.var1)[ignore],
                zgrid=self.eng.data.massgrid[ignore, :, 1].transpose(),
                config=self.eng.config, title="Mass vs. " + self.eng.data.v1name, test_kda=True, xlab="Mass (Da)",
                ylab=self.eng.data.v1name)
            tend = time.perf_counter()
            print("Plot 5: %.2gs" % (tend - tstart))
            pass

    def recolor_spectra(self):
        self.view.ypanel.list.recolor()
        self.on_replot()

    def on_full(self, e=None):
        maxmz = np.amax(self.eng.data.spectra[0].rawdata[:, 0])
        minmz = np.amin(self.eng.data.spectra[0].rawdata[:, 0])
        self.view.controls.ctlminmz.SetValue(str(minmz))
        self.view.controls.ctlmaxmz.SetValue(str(maxmz))
        self.on_dataprep_button()

    def make_top(self, index=0):
        """
        Tested
        :param index:
        :return:
        """
        print("Top index is now:", index)
        self.eng.data.data2 = self.eng.data.spectra[index].data2

    def on_ignore(self, indexes):
        """
        Partly tested - Passed
        :param indexes:
        :return:
        """
        print("Ignoring:", indexes)
        spectra = self.eng.data.get_spectra()
        for i in indexes:
            spectra[i].ignore = 1
        self.on_replot()

    def on_isolate(self, indexes):
        """
        Partly tested - Passed
        :param indexes:
        :return:
        """
        print("Isolating:", indexes)
        spectra = self.eng.data.get_spectra()
        for i, s in enumerate(spectra):
            if np.any(np.array(indexes) == i):
                s.ignore = 0
            else:
                s.ignore = 1
        try:
            self.make_top(indexes[0])
        except:
            print("Failed to make top")
        self.on_replot(plotsums=False)

    def on_repopulate(self):
        """
        Partly tested - Passed
        :return:
        """
        for s in self.eng.data.spectra:
            s.ignore = 0
        self.on_replot()

    def on_color_change(self, item, color):
        """
        Manual Test - Passed
        :param item:
        :param color:
        :return:
        """
        spectra = self.eng.data.get_spectra()
        spectra[item].color = color
        self.on_replot()

    def rename_var1(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 1 Name", "Enter Variable 1 name:", defaultvalue=str(self.eng.data.v1name))
        dlg.ShowModal()
        self.eng.data.v1name = dlg.value
        self.view.ypanel.list.rename_column(1, self.eng.data.v1name)
        self.on_replot()

    def rename_var2(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 2 Name", "Enter Variable 2 name:", defaultvalue=str(self.eng.data.v2name))
        dlg.ShowModal()
        self.eng.data.v2name = dlg.value
        self.view.ypanel.list.rename_column(2, self.eng.data.v2name)
        self.on_replot()

    def export_vars(self):
        """
        Manual Test - Passed
        :return:
        """
        listdat = self.view.ypanel.list.get_list()
        for i, l in enumerate(listdat):
            self.eng.data.var1[i] = l[1]
            self.eng.data.var2[i] = l[2]
            self.eng.data.spectra[i].var1 = l[1]
            self.eng.data.spectra[i].var2 = l[2]
        self.eng.data.export_vars()

    def on_dataprep_button(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.view.SetStatusText("Preparing Data..", number=5)
        self.view.clear_plots()
        self.export_config()
        self.check_badness()
        self.eng.process_data()
        self.on_auto_peak_width(set=False)
        self.makeplot1()
        self.view.SetStatusText("Data Prep Complete", number=5)
        pass

    def on_unidec_button(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.view.SetStatusText("Running UniDec...", number=5)
        self.view.clear_plots()
        tstart = time.perf_counter()
        self.export_config()
        self.check_badness()
        self.eng.run_unidec()
        tend = time.perf_counter()
        self.eng.config.runtime = (tend - tstart)
        self.makeplot1()
        self.makeplot2_mud()
        print("Run Time:", self.eng.config.runtime)
        self.view.SetStatusText("UniDec Done %.2gs" % self.eng.config.runtime, number=5)
        pass

    def on_pick_peaks(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.view.SetStatusText("Picking Peaks...", number=5)
        self.export_config()
        self.eng.pick_peaks()
        self.view.peakpanel.add_data(self.eng.pks, show="dscore")
        self.view.peakpanel.meta = True
        self.peak_plots()
        self.view.SetStatusText("Peak Detection and Extraction Complete", number=5)
        pass

    def on_pick_scanpeaks(self, e=None):
        self.view.SetStatusText("Picking Peaks by Scan...", number=5)
        self.export_config()
        self.eng.pick_scanpeaks()
        self.view.peakpanel.add_data(self.eng.pks, show="dscore")
        self.view.peakpanel.meta = True
        self.peak_plots()
        self.view.SetStatusText("ScanPeak Detection and Extraction Complete", number=5)

    def on_filter_peaks_MUD(self, e=None):
        self.on_filter_peaks(e)
        self.view.SetStatusText("UniScore: " + str(round(self.eng.pks.uniscore * 100, 2)), number=3)
        self.peak_plots()

    def on_auto(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        self.export_config()
        self.on_dataprep_button()
        self.on_unidec_button()
        self.on_pick_peaks()
        pass

    def import_vars(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.import_config()
        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 1 Metadata Name", "Enter Variable 1 Metadata Name:",
                                 defaultvalue=str(self.eng.data.v1name))
        dlg.ShowModal()
        self.eng.data.v1name = dlg.value

        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 2 Metadata Name", "Enter Variable 2 Metadata Name:",
                                 defaultvalue=str(self.eng.data.v2name))
        dlg.ShowModal()
        self.eng.data.v2name = dlg.value

        self.eng.data.import_vars(get_vnames=False)
        self.view.ypanel.list.populate(self.eng.data)
        self.on_replot()

    def export_vars_dialog(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 1 Metadata Name", "Enter Variable 1 Metadata Name:",
                                 defaultvalue=str(self.eng.data.v1name))
        dlg.ShowModal()
        self.eng.data.v1name = dlg.value

        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface("Variable 2 Metadata Name", "Enter Variable 2 Metadata Name:",
                                 defaultvalue=str(self.eng.data.v2name))
        dlg.ShowModal()
        self.eng.data.v2name = dlg.value
        self.export_vars()

    def on_charge_states_mud(self, e=None):
        self.on_charge_states(plot=self.view.plot1)

    def on_mass_tools(self, e=None, show=True):
        """
        Opens masstools window. Manual Test - Passed

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
        if not ud.isempty(self.eng.config.matchlist) and result == 0:
            if len(self.eng.config.matchlist[3]) == self.eng.pks.plen:
                self.view.SetStatusText("Matching", number=5)
                np.savetxt(self.eng.config.matchfile, np.transpose(self.eng.config.matchlist), fmt='%s', delimiter=",")
                self.view.peakpanel.add_data(self.eng.pks)
                try:
                    self.makeplot6()
                    self.makeplot8()
                except:
                    pass
            else:
                self.eng.config.matchlist = []

        self.export_config()
        self.view.SetStatusText("Match Done", number=5)
        pass

    def on_autocorr2(self, index):
        """
        Manual Test - Passed
        :param index:
        :return:
        """
        spectra = self.eng.data.get_spectra()
        data = spectra[index].massdat
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, data)
        dlg.ShowModal()

    def on_autocorr_window_mud(self, e=None):
        if ud.isempty(self.eng.data.massdat):
            self.eng.data.import_grids_and_peaks()
        self.on_autocorr_window()

    def on_fft_window2(self, index):
        """
        Manual Test - Passed
        :param index:
        :return:
        """
        try:
            spectra = self.eng.data.get_spectra()
        except:
            spectra = self.eng.data.spectra
        data = spectra[index].data2
        fft_window.FFTWindow(self.view, data, self.eng.config)

    def on_kendrick(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.eng.data.import_grids_and_peaks()
        MassDefects.MassDefectWindow(self.view, self.eng.data.massgrid, self.eng.config, yvals=self.eng.data.var1,
                                     directory=self.eng.config.udir,
                                     value=self.eng.config.molig)
        pass

    def on_2d_grid(self, e):
        """
        Manual test - Passed
        :param e:
        :return:
        """
        self.eng.data.import_grids_and_peaks()
        exwindow = Extract2D.Extract2DPlot(self.view, self.eng.data.massgrid, self.eng.config, yvals=self.eng.data.var1,
                                           params=self.eng.config.gridparams,
                                           directory=os.path.split(self.eng.config.outfname)[0])
        self.eng.config.gridparams = exwindow.params

    def on_fft_window(self, e):
        """
        Manual test - Passed
        :param e:
        :return:
        """
        rawdatalist = [s.rawdata for s in self.eng.data.spectra]
        metafft.FFTWindow(self.view, rawdatalist, self.eng.data.var1, self.eng.config)
        pass

    def on_animate_mass(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.eng.sum_masses()
        PlotAnimations.AnimationWindow(self.view, self.eng.data.massgrid, self.eng.config, yvals=self.eng.data.var1,
                                       pksmode="mass")

    def on_animate_annotated_mass(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.eng.sum_masses()
        PlotAnimations.AnimationWindow(self.view, self.eng.data.massgrid, self.eng.config, pks=self.eng.pks,
                                       pksmode="mass", yvals=self.eng.data.var1)

    def on_animate_mz(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        # self.eng.sum_masses()
        newgrid = []
        for s in self.eng.data.spectra:
            newgrid.append(s.data2)

        PlotAnimations.AnimationWindow(self.view, newgrid, self.eng.config, yvals=self.eng.data.var1)

    def on_animate_annotated_mz(self, e=None):
        """
        :return:
        """
        if self.eng.config.rawflag > 1:
            self.warn("Unable to Animate: Need to turn off Fast Profile or Fast Centroid")
            return

        newgrid = []
        for s in self.eng.data.spectra:
            newgrid.append(s.data2)
        self.eng.peaks_heights()
        PlotAnimations.AnimationWindow(self.view, newgrid, self.eng.config, yvals=self.eng.data.var1,
                                       pks=self.eng.pks)

    def on_animate_2d(self, e=None, type="mass"):
        """
        Manual Test - Passed
        :param e:
        :param type:
        :return:
        """
        if self.eng.config.rawflag > 1:
            self.warn("Unable to Animate: Need to turn off Fast Profile or Fast Centroid")
            return

        self.eng.sum_masses()
        dlg = miscwindows.SingleInputDialog(self.view)
        dlg.initialize_interface(title="Set Compression", message="Number of x values to compress:", defaultvalue="10")
        dlg.ShowModal()
        try:
            compress = int(dlg.value)
            if compress > 1:
                print("Compressing Data by:", compress)
        except (ValueError, TypeError, AttributeError):
            print("Unrecognized compression value")
            compress = 0

        print("Loading 2D Data...")
        data2 = []
        for i, s in enumerate(self.eng.data.spectra):
            if type == "mz":
                igrid = s.mzgrid
                mdat = s.data2[:, 0]
            else:
                igrid = s.massgrid
                mdat = s.massdat[:, 0]
            igrid /= np.amax(igrid)
            mgrid, zgrid = np.meshgrid(mdat, s.ztab, indexing='ij')
            if compress > 1:
                igrid = np.reshape(igrid, mgrid.shape)
                m, z, d = IM_functions.compress_2d(mgrid, zgrid, igrid, compress)
                dat = np.transpose([np.ravel(m), np.ravel(z), np.ravel(d)])
            else:
                dat = np.transpose([np.ravel(mgrid), np.ravel(zgrid), np.ravel(igrid)])
            data2.append(dat)
            print(i, end=' ')
        data2 = np.array(data2)
        print("Loaded 2D Data", data2.shape)
        PlotAnimations.AnimationWindow(self.view, data2, self.eng.config, mode="2D", yvals=self.eng.data.var1,
                                       pksmode=type)

    def on_animate_2d_mass(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.on_animate_2d(type="mass")

    def on_animate_2d_mz(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        self.on_animate_2d(type="mz")

    def on_delete_spectrum(self, indexes=None):
        """
        Tested
        :param indexes:
        :return:
        """
        if indexes is not None:
            self.eng.data.remove_data(indexes)
            self.eng.data.import_vars()
            self.view.clear_plots()
            self.view.ypanel.list.populate(self.eng.data)
            try:
                self.eng.pick_peaks()
            except:
                pass
            self.on_replot()

    def on_export_params(self, e=None):
        """
        Runs self.eng.export_params(), which gets critical peak parameters and writes them to a file. Tested.
        :param e: event or arguments passed to self.eng.export_params()
        :return: None
        """
        self.eng.export_params(e)

    def repack_hdf5(self, e=None):
        if self.eng.config.hdf_file != 'default.hdf5':
            new_path = self.eng.config.hdf_file.replace(".hdf5", "temp.hdf5")
            if 0 == subprocess.call(
                    "\"" + self.eng.config.h5repackfile + "\" \"" + self.eng.config.hdf_file + "\" \"" + new_path + "\"") and os.path.isfile(
                new_path):
                os.remove(self.eng.config.hdf_file)
                os.rename(new_path, self.eng.config.hdf_file)

    def recursive_repack_hdf5(self, e=None):
        repack_dir = FileDialogs.open_dir_dialog(message="Select directory to repack")
        if repack_dir:
            for root, dirs, files in os.walk(repack_dir):
                for name in files:
                    if name.endswith('.hdf5'):
                        name = os.path.join(root, name)
                        new_path = name.replace(".hdf5", "temp.hdf5")
                        if 0 == subprocess.call(
                                "\"" + self.eng.config.h5repackfile + "\" \"" + name + "\" \"" + new_path + "\"") and os.path.isfile(
                            new_path):
                            os.remove(name)
                            os.rename(new_path, name)

    def on_undo(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        self.export_config()
        self.eng.undo()
        self.import_config()
        # print "Undo"

    def on_redo(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        self.export_config()
        self.eng.redo()
        self.import_config()
        # print("Redo")


class UniDecApp(MetaUniDecBase):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize App, Tested
        :param args:
        :param kwargs:
        :return: UniDecApp object
        """
        MetaUniDecBase.__init__(self, *args, **kwargs)
        self.init(*args, **kwargs)

        # self.on_open(0)
        try:
            if False:
                # testdir = "C:\Data\Others\\UniDec test data set"
                # testfile = "test.hdf5"
                testdir = "C:\\Data\\New"
                testfile = "20170209_P0B_dPOPC_POPC_ND_D1T0m_pos_ISTRAP_RAMP_0_275_25_1.hdf5"
                testpath = os.path.join(testdir, testfile)

                self.open_file(testpath)
                # self.on_pick_peaks()

                # self.on_match()
                # self.on_autoformat()

        except:
            pass
        # self.on_animate_annotated_mass()

    def init(self, *args, **kwargs):
        """
        Initialize Engine and View. Load defaults. Tested
        :param args:
        :param kwargs:
        :return:
        """
        pub.subscribe(self.on_get_mzlimits, 'mzlimits')
        pub.subscribe(self.on_left_click, 'left_click')

        self.eng = mudeng.MetaUniDec()

        self.view = mudview.Mainwindow(self, "MetaUniDec", self.eng.config)

        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)
        self.view.menu.update_recent()

    def on_open(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        print("Opening")
        dlg = wx.FileDialog(self.view, "Choose a data file in HDF5 format", '', "", "*.hdf5*")
        if dlg.ShowModal() == wx.ID_OK:
            self.view.SetStatusText("Opening", number=5)
            filename = dlg.GetFilename()
            print("Openening: ", filename)
            if os.path.splitext(filename)[1] != ".hdf5":
                print("Need HDF5 file")
                return
            dirname = dlg.GetDirectory()
            self.eng.config.dirname = dirname
            self.eng.config.hdf_file = os.path.join(dirname, filename)
            self.open_file(self.eng.config.hdf_file)
        dlg.Destroy()
        pass

    def open_file(self, path=None):
        """
        Tested
        :param path:
        :return:
        """
        tstart = time.perf_counter()
        self.view.clear_plots()
        if path is None:
            path = self.eng.config.hdf_file
        print("Opening:", path)
        self.eng.open(path)
        # print("1: %.2gs" % (time.perf_counter() - tstart))
        self.import_config()
        # print("2: %.2gs" % (time.perf_counter() - tstart))
        self.view.ypanel.list.populate(self.eng.data)
        # print("3: %.2gs" % (time.perf_counter() - tstart))
        # self.on_auto_peak_width(set=False)
        self.eng.get_auto_peak_width(set=False)
        # print("4: %.2gs" % (time.perf_counter() - tstart))
        self.makeplot1()
        self.makeplot2_mud()
        print("Load Time: %.2gs" % (time.perf_counter() - tstart))
        self.view.SetStatusText("File: " + self.eng.config.hdf_file, number=1)

        self.write_to_recent(self.eng.config.hdf_file)
        self.view.menu.update_recent()

    def makeplot9(self, e=None):
        """
        Tested
        :param e:
        :return:
        """
        spectra = self.eng.data.get_spectra()

        if self.eng.config.datanorm == 0:
            sep = self.eng.config.separation * self.eng.data.get_max_data2()
        else:
            sep = self.eng.config.separation

        if len(spectra) > self.eng.config.crossover:
            mult = int(len(spectra) / self.eng.config.numtot)
        else:
            mult = 1
        for i, s in enumerate(spectra[::mult]):
            if not ud.isempty(s.zdata):
                if i == 0:
                    self.view.plot9.plotrefreshtop(s.zdata[:, 0], s.zdata[:, 1], title="Total Charge Spectrum",
                                                   xlabel="Charge",
                                                   ylabel="Intensity", label=s.name, config=self.eng.config,
                                                   color=s.color,
                                                   test_kda=False,
                                                   nopaint=True)
                else:
                    self.view.plot9.plotadd(s.zdata[:, 0], s.zdata[:, 1] - i * sep,
                                            colval=s.color, newlabel=s.name)
        self.view.plot9.repaint()

    def makeplot6(self, e=None, show="height"):
        """
        Plots bar chart of peak heights or areas. Tested
        :param e: unused event
        :param show: What parameter to plot
        "height" will plot p.height for p in self.eng.pks.peaks
        "integral" will plot p.integral
        :return: None
        """
        if self.eng.pks.plen > 0:
            num = 0
            ints = []
            cols = []
            labs = []

            spectra = self.eng.data.get_spectra()
            if len(spectra) > self.eng.config.crossover:
                mult = int(len(spectra) / self.eng.config.numtot)
            else:
                mult = 1

            for i, s in enumerate(self.eng.data.spectra[::mult]):
                if s.ignore == 0:
                    num += 1
                    ints.append(0.0000000001)
                    cols.append([1, 1, 1, 1])
                    labs.append(s.index)

                    for p in self.eng.pks.peaks:
                        if p.ignore == 0:
                            num += 1
                            ints.append(p.extracts[i])
                            cols.append(p.color)
                            labs.append(p.label)

            num += 1
            ints.append(0.0000000001)
            cols.append([1, 1, 1, 1])
            labs.append("")
            self.view.plot6.barplottop(list(range(0, num)), ints, labs, cols, "Species", "Intensity",
                                       "Peak Intensities", repaint=False)
            num = 0
            for i, s in enumerate(self.eng.data.spectra[::mult]):
                if s.ignore == 0:
                    num += 1
                    for p in self.eng.pks.peaks:
                        if p.ignore == 0:
                            e = p.extracts[i]
                            num += 1
                            self.view.plot6.plotadddot(num - 1, e, p.color, p.marker)
            self.view.plot6.repaint()

    def makeplot8(self):
        """
        Tested
        :return:
        """
        try:
            self.view.plot8.clear_plot()
            if self.eng.pks.plen > 0:
                ignore = self.eng.data.get_bool()
                ignore2 = self.eng.pks.get_bool()
                zdat = self.eng.data.exgrid[ignore2, :]
                zdat = zdat[:, ignore]
                xvals = []
                for p in self.eng.pks.peaks:
                    if p.ignore == 0:
                        xvals.append(p.label)
                self.view.plot8._axes = [0.12, 0.12, 0.75, 0.8]
                var1 = np.array(self.eng.data.var1)[ignore]
                try:
                    var1 = var1.astype(float)
                    ylabel = self.eng.data.v1name
                except:
                    var1 = np.arange(0, len(self.eng.data.var1))[ignore]
                    ylabel = "Index"

                self.view.plot8.contourplot(xvals=np.arange(0, len(xvals)), yvals=var1,
                                            zgrid=zdat, normflag=0,
                                            normrange=[0, np.amax(zdat)], config=self.eng.config,
                                            xlab="Peaks", ylab=ylabel, discrete=1,
                                            ticloc=np.arange(0, len(xvals)),
                                            ticlab=xvals)
        except Exception as e:
            print(e)

    def peak_plots(self, e=None):
        """
        Called when peaks are picked or changed
        """
        self.makeplot2_mud()
        self.plot_sums()
        self.makeplot6()
        self.makeplot7()
        self.makeplot8()

    def on_replot(self, e=None, plotsums=True):
        """
        Tested
        :param e:
        :return:
        """
        self.view.SetStatusText("Replotting", number=5)
        self.export_config()
        self.makeplot1()
        self.makeplot2_mud()
        if plotsums:
            self.plot_sums()
        self.makeplot6()
        self.makeplot7()
        self.makeplot8()
        pass

    def on_delete(self, e=None):
        """
        Manual test - Passed
        :param e:
        :return:
        """
        # self.eng.pick_peaks()
        self.makeplot2_mud()
        self.makeplot6()
        self.makeplot7()
        self.makeplot8()
        self.plot_sums()
        pass



    def on_left_click(self, xpos, ypos):
        """
        Triggered by pubsub from plot windows.
        Gets a m/z peak near the click, stores it, and waits for another click.
        When two clicks has been performed, it tries to calculate the mass from their m/z value.
        Manual Test - Passed
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
                print(outstring)

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

    def on_import_mzml(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose ramp data files mzml or Thermo Raw format",
                                                       file_type="All Files|*.*| Thermo RAW files (*.RAW)|*.RAW|mzML files (*.mzML)|*.mzML")
        print("Paths:", paths)
        if paths is not None:
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Timestep", "Enter ramp timestep to compress in minutes:", defaultvalue=str(1.0))
            dlg.ShowModal()
            timestep = dlg.value
            self.eng.import_mzml(paths, timestep=timestep)

            print("Opening Just Created File:", self.eng.outpath)
            self.open_file(self.eng.outpath)

    def on_import_mzml_scans(self, e):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose ramp data files mzml or Thermo Raw format",
                                                       file_type="All Files|*.*| Thermo RAW files (*.RAW)|*.RAW|mzML files (*.mzML)|*.mzML")
        if paths is not None:
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Scan Step", "Enter number of scans to compress:", defaultvalue=str(1.0))
            dlg.ShowModal()
            scanstep = dlg.value
            self.eng.import_mzml(paths, scanstep=scanstep)

            print("Opening Just Created File:", self.eng.outpath)
            self.open_file(self.eng.outpath)

    def on_import_multiple_times(self, e):
        """
        Manual Test - Passed 2 RAW's, 2 mzml
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose ramp data files mzml or Thermo Raw format",
                                                       file_type="All Files|*.*| Thermo RAW files (*.RAW)|*.RAW|mzML files (*.mzML)|*.mzML")
        if paths is not None:
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Time point", "Enter start time point:", defaultvalue=str(0.0))
            dlg.ShowModal()
            starttp = dlg.value
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Time point", "Enter end time point:", defaultvalue=str(1.0))
            dlg.ShowModal()
            endtp = dlg.value

            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Name", "Enter name of final file (.hdf5 not required):", defaultvalue="")
            dlg.ShowModal()
            name = dlg.value
            # Dummy checks
            if float(endtp) < float(starttp) or float(endtp) < 0 or float(starttp) < 0:
                print("Bad values inputted")
            else:
                timestep = float(endtp) - float(starttp)
                self.eng.import_mzml(paths, timestep=timestep, name=name, starttp=float(starttp),
                                     endtp=float(endtp))

            print("Opening Just Created File:", self.eng.outpath)
            self.open_file(self.eng.outpath)

    def on_import_multiple_scans(self, e):
        paths = FileDialogs.open_multiple_files_dialog(message="Choose ramp data files mzml or Thermo Raw format",
                                                       file_type="All Files|*.*| Thermo RAW files (*.RAW)|*.RAW|mzML files (*.mzML)|*.mzML")
        if paths is not None:
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Scan Start", "Enter starting scan desired:", defaultvalue=str(1))
            dlg.ShowModal()
            startscan = dlg.value
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Scan End", "Enter ending scan desired:", defaultvalue=str(1))
            dlg.ShowModal()
            endscan = dlg.value
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg = miscwindows.SingleInputDialog(self.view)
            dlg.initialize_interface("Name", "Enter name of final file (.hdf5 not required):", defaultvalue="")
            dlg.ShowModal()
            name = dlg.value
            # + 1 including the endscan
            self.eng.import_mzml(paths, startscan=int(startscan), endscan=(int(endscan) + 1), name=name)

            print("Opening Just Created File:", self.eng.outpath)
            self.open_file(self.eng.outpath)

    def on_wizard(self, e=None):
        print("Launching Waters Converter Wizard")
        app = wx.App(False)
        frame = ImportWizard(None)
        frame.Show()
        app.MainLoop()

    def on_new_file(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        path = FileDialogs.save_file_dialog(message="Name New HDF5 Save File", file_types="*.hdf5",
                                            default_file="default.hdf5")
        if path is not None:
            self.new_file(path)

    def new_file(self, path):
        """
        Tested
        :param path:
        :return:
        """
        # self.eng = mudeng.MetaUniDec()
        # self.import_config()
        self.view.clear_plots()
        self.eng.data.new_file(path)
        self.open_file(path)

    def on_paste_spectrum(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        if self.eng.data.filename is None:
            self.on_new_file()

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
                self.eng.data.add_data(data)
                self.view.ypanel.list.populate(self.eng.data)
                self.makeplot1()
            else:
                print("Paste failed, got: ", data)
        except Exception as e:
            print(e)
            wx.MessageBox("Unable to open the clipboard", "Error")

    def add_files(self, paths):
        """
        Tested
        :param paths:
        :return:
        """
        self.view.clear_plots()
        for p in paths:
            dirname = os.path.dirname(p)
            filename = os.path.basename(p)
            self.eng.data.add_file(filename, dirname)
        self.view.ypanel.list.populate(self.eng.data)
        self.makeplot1()

    def on_add_file(self, e=None):
        """
        Manual Test - Passed
        :param e:
        :return:
        """
        if self.eng.data.filename is None:
            self.on_new_file()

        paths = FileDialogs.open_multiple_files_dialog(message="Choose data files in txt or mzml format")
        if paths is not None:
            print("Openening: ", paths)
            self.add_files(paths)

    def on_additional_parameters(self, e=None):
        """
        Open additional data parameter window, which hads some of the experimental and obscure variables.
        Window directly modifies self.eng.config.
        Saves config to file after window is closed.
        Manual Test - Unused
        :param e: unused event
        :return: None
        """
        dlg = miscwindows.AdditionalParameters(self.view)
        dlg.initialize_interface(self.eng.config)
        dlg.ShowModal()
        self.export_config(self.eng.config.confname)

    def on_data_collector(self, e=None):
        """
        Spawns separate DataCollector window.
        Manual Test - Failed
        :param e: unused event
        :return: None
        """
        dc = datacollector.DataCollector(None, "Data Collector", config=self.eng.config, pks=self.eng.pks,
                                         directory=self.eng.config.dirname, hdf_file=self.eng.config.hdf_file)

    def on_exp_fit(self, e=None):
        self.on_fit(fit="exp")

    def on_lin_fit(self, e=None):
        self.on_fit(fit="lin")

    def on_sig_fit(self, e=None):
        self.on_fit(fit="sig")

    def on_fit(self, fit="sig"):
        self.eng.fit_data(fit)
        self.makeplot7(fitgrid=self.eng.data.fitgrid)

    def on_ultra_meta(self, e=None):
        """
        Spawns separate DataCollector window.
        Manual Test - Failed
        :param e: unused event
        :return: None
        """
        dc = ultrameta.DataCollector(None, "Ultra Meta Data Collector", config=self.eng.config,
                                     directory=self.eng.config.dirname)

    def on_batch_config(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose HDF5 files to apply current config to",
                                                       file_type="*.hdf5")
        if paths is not None:
            self.view.SetStatusText("Applying Batch Configs...", number=5)
            self.eng.batch_set_config(paths)
            self.view.SetStatusText("Batch Configs Applied", number=5)

    def on_batch_run(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose HDF5 files to batch process.",
                                                       file_type="*.hdf5")
        if paths is not None:
            self.view.SetStatusText("Batch Setting Config...", number=5)
            self.eng.batch_set_config(paths)
            self.view.SetStatusText("Batch Running...", number=5)
            self.eng.batch_run_unidec(paths)
            self.view.SetStatusText("Batch Run Done", number=5)

    def on_batch_extract(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose HDF5 files to batch extract.",
                                                       file_type="*.hdf5")
        if paths is not None:
            self.view.SetStatusText("Batch Extracting...", number=5)
            self.eng.batch_extract(paths)
            self.view.SetStatusText("Batch Extract Done", number=5)

    def on_batch_cre(self, e=None):
        """
        Manual Test
        :param e:
        :return:
        """
        paths = FileDialogs.open_multiple_files_dialog(message="Choose HDF5 files to apply current config to",
                                                       file_type="*.hdf5")
        if paths is not None:
            self.view.SetStatusText("Batch Setting Config...", number=5)
            self.eng.batch_set_config(paths)
            self.view.SetStatusText("Batch Running and Extracting...", number=5)
            self.eng.batch_run_unidec(paths)
            self.view.SetStatusText("Batch Extract Done", number=5)

    def remake_mainwindow(self, tabbed=None):
        """
        Tested
        :param tabbed:
        :return:
        """
        iconfile = self.view.icon_path
        wx.GetApp().Yield()
        self.view.on_exit()
        self.view = []
        self.view = mudview.Mainwindow(self, "MetaUniDec", self.eng.config, iconfile=iconfile, tabbed=tabbed)
        self.view.Show()
        self.view.import_config_to_gui()

    def on_flip_tabbed(self, e=None):
        """
        Flips between tabbed plots and a single window of plots - Tested
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
        except Exception as exc:
            print("Failed to replot when making window:", exc)
        if self.view.tabbed == 1:
            print("Tabbed Mode")
        elif self.view.tabbed == 0:
            print("Single Plot Window Mode")

    def on_getting_started(self, e=None):
        helpDlg = HelpDlg(1)
        helpDlg.Show()

    def on_import_window(self, e=None):
        helpDlg = HelpDlg(2)
        helpDlg.Show()

    def on_plot_window(self, e=None):
        helpDlg = HelpDlg(3)
        helpDlg.Show()

    def on_peak_window(self, e=None):
        helpDlg = HelpDlg(4)
        helpDlg.Show()

    def on_data_processing(self, e=None):
        helpDlg = HelpDlg(5)
        helpDlg.Show()

    def on_unidec_parameters(self, e=None):
        helpDlg = HelpDlg(6)
        helpDlg.Show()

    def on_additional_filters(self, e=None):
        helpDlg = HelpDlg(7)
        helpDlg.Show()

    def on_peak_selection(self, e=None):
        helpDlg = HelpDlg(8)
        helpDlg.Show()

    def on_additional_plotting(self, e=None):
        helpDlg = HelpDlg(9)
        helpDlg.Show()

    def on_auto_import_help(self, e=None):
        helpDlg = HelpDlg(10)
        helpDlg.Show()

    def on_presets_help(self, e=None):
        helpDlg = HelpDlg(11)
        helpDlg.Show()

    def on_batch_help(self, e=None):
        helpDlg = HelpDlg(12)
        helpDlg.Show()

    def on_peak_width_tool_help(self, e=None):
        helpDlg = HelpDlg(13)
        helpDlg.Show()

    def on_oligomer_help(self, e=None):
        helpDlg = HelpDlg(14)
        helpDlg.Show()

    def on_auto_match_help(self, e=None):
        helpDlg = HelpDlg(15)
        helpDlg.Show()

    def on_animate_help(self, e=None):
        helpDlg = HelpDlg(16)
        helpDlg.Show()

    def on_autocorr_help(self, e=None):
        helpDlg = HelpDlg(17)
        helpDlg.Show()

    def on_fft_help(self, e=None):
        helpDlg = HelpDlg(18)
        helpDlg.Show()

    def on_baseline_help(self, e=None):
        helpDlg = HelpDlg(19)
        helpDlg.Show()

    def on_autoformat(self, e=None):
        self.on_match()
        self.eng.pks.auto_format()
        self.on_delete()
        self.view.peakpanel.add_data(self.eng.pks)



# Critical
# TODO: Thorough testing
# TODO: Better tuning and control of autobaseline

# Serious work
# TODO: Weighted average of charge states to calculate mass error
# TODO: Error as FWHM of peak
# TODO: Average charge state as extraction parameter

# Next
# TODO: Make Launcher Fancier
# TODO: Add DTIMS and TWIMS to Launcher
# TODO: Documentation

# Not critical
# TODO: Better MUD with data collector
# TODO: Light grey background for some list ctrls in MetaUniDec
# TODO: A notice that pops up when trying to open RAW files without multiplierz or MSFileReader


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = UniDecApp()
    app.start()

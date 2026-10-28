"""Standalone ion mobility-mass spectrometry GUI application."""

import multiprocessing
import os
import time

import numpy as np
import wx
from pubsub import pub

from unidec.GUniDec import UniDecApp as UniDecMSApp
from unidec.modules import IMEng
import unidec.modules.IM_functions as IM_func
import unidec.modules.IM_windows as IM_wind
from unidec.modules import IM_mainwindow
import unidec.tools as ud


class UniDecIMApp(UniDecMSApp):
    """Presenter for the dedicated UniDec IM-MS workflow."""

    def init(self, *args, **kwargs):
        self.eng = IMEng.UniDecIM(ignore_args=True)
        self.view = IM_mainwindow.IMMainwindow(
            self, "UniDec IM: Ion Mobility-Mass Spectrometry", self.eng.config
        )

        pub.subscribe(self.on_integrate, "integrate")
        pub.subscribe(self.on_smash, "smash")
        pub.subscribe(self.on_get_mzlimits, "mzlimits")
        pub.subscribe(self.on_left_click, "left_click")

        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)
        self.view.menu.update_recent()
        self.on_load_default(0)
        self.eng.config.imflag = 1
        if self.eng.config.mzbins == 0:
            self.eng.config.mzbins = 1
        self.import_config()

        path = kwargs.get("path", self.infile)
        if path is not None:
            directory, filename = os.path.split(path)
            self.on_open_file(filename, directory)

    def import_config(self, file_name=None):
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.read_hdf5(file_name)
            else:
                self.eng.config.config_import(file_name)
        self.eng.config.imflag = 1
        self.view.import_config_to_gui()
        if self.eng.config.filetype == 1:
            self.eng.config.write_hdf5()
        self.eng.update_history()

    def export_config(self, file_name=None):
        self.view.export_gui_to_config()
        self.eng.config.imflag = 1
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.write_hdf5(file_name)
            else:
                self.eng.config.config_export(file_name)
        self.eng.update_history()

    def on_open_file(self, filename, directory=None, skipengine=False, refresh=False, **kwargs):
        super().on_open_file(
            filename,
            directory=directory,
            skipengine=skipengine,
            refresh=refresh,
            **kwargs,
        )
        if self.eng.config.batchflag != 1 and len(self.eng.data.data3) > 0:
            self.view.controls.ctlmindt.SetValue(str(np.amin(self.eng.data.data3[:, 1])))
            self.view.controls.ctlmaxdt.SetValue(str(np.amax(self.eng.data.data3[:, 1])))

    def on_raw_open(self, e=None, dirname=None):
        self.export_config(self.eng.config.confname)
        if dirname is None:
            from unidec.modules.isolated_packages import FileDialogs

            dirname = FileDialogs.open_single_dir_dialog("Choose a raw file", "")
        if dirname is None:
            return

        if self.eng.config.compressflag == 1:
            binsize = str(self.eng.config.mzbins)
            print("Converting at resolution of: " + binsize)
        else:
            binsize = "0"
            print("Converting using full resolution")
        self.view.SetStatusText("Converting", number=5)
        dirname = os.path.abspath(dirname)
        filename, output_dir = self.eng.raw_process(dirname, True, binsize=binsize)
        if filename is not None:
            self.on_open_file(filename, output_dir)

    def after_unidec_run(self):
        self.view.SetStatusText("UniDec IM Plot", number=5)
        self.make_im_plots()
        self.view.SetStatusText("R\u00B2: " + str(self.eng.config.error), number=3)
        self.view.plot4.clear_plot()
        self.view.plot6.clear_plot()
        self.view.peakpanel.clear_list()

    def makeplot1(self, e=None, intthresh=False, imfit=True):
        self.eng.makeplot1(plot=self.view.plot1, intthresh=intthresh, imfit=imfit)
        if self.eng.config.batchflag == 0:
            self.eng.makeplot1im(
                plot1im=self.view.plot1im,
                plot1fit=self.view.plot1fit,
                imfit=imfit,
            )

    def make_im_plots(self):
        """Render the IM-MS data, CCS projections, and charge map."""
        if self.eng.config.batchflag != 0:
            return
        self.makeplot1(1)
        self.makeplot2(1)
        self.makeplot5(1)
        self.makeplot3(1)
        self.view.plot2ccs.plotrefreshtop(
            self.eng.data.ccsdata[:, 0],
            self.eng.data.ccsdata[:, 1],
            title="CCS Distribution",
            xlabel="CCS (${\\AA}$$^2$)",
            ylabel="Intensity",
            label="CCS Summation",
            config=self.eng.config,
            nopaint=False,
        )
        self.view.plot5mccs.contourplot(
            xvals=self.eng.data.massdat[:, 0],
            yvals=self.eng.data.ccsdata[:, 0],
            zgrid=np.ravel(self.eng.data.massccs),
            config=self.eng.config,
            ylab="CCS (${\\AA}$$^2$)",
            title="Mass vs. CCS",
            test_kda=True,
        )
        ccsgrid, zgrid = np.meshgrid(
            self.eng.data.ztab, self.eng.data.ccsdata[:, 0], sparse=False, indexing="ij"
        )
        self.view.plot5ccsz.contourplot(
            np.transpose([np.ravel(ccsgrid), np.ravel(zgrid), np.ravel(self.eng.data.ccsz)]),
            self.eng.config,
            xlab="Charge",
            ylab="CCS (${\\AA}$$^2$)",
            title="CCS vs. Charge",
        )
        try:
            self.view.plot3color.make_color_plot(
                self.eng.data.mztgrid,
                np.unique(self.eng.data.data3[:, 0]),
                np.unique(self.eng.data.data3[:, 1]),
                self.eng.data.ztab,
            )
        except Exception as exc:
            print("Color Plot Error", exc)

    def on_plot_nativeccs(self, e=None):
        if not ud.isempty(self.eng.data.massdat):
            ccses = [
                IM_func.calc_native_ccs(mass, self.eng.config.gasmass)
                for mass in self.eng.data.massdat[:, 0]
            ]
            self.view.plot5mccs.subplot1.plot(
                self.eng.data.massdat[:, 0] / self.view.plot5mccs.kdnorm,
                ccses,
                color="r",
            )
            self.view.plot5mccs.repaint()

    def on_replot(self, e=None):
        self.export_config(self.eng.config.confname)
        self.make_im_plots()
        self.makeplot4()
        self.makeplot6()
        if self.view.plot9.flag and self.view.plot10.flag:
            self.make_cube_plot()

    def make_cube_plot(self, event=None):
        self.export_config(self.eng.config.confname)
        try:
            start = time.perf_counter()
            self.view.plot9.cubeplot(
                np.unique(self.eng.data.data3[:, 0]),
                np.unique(self.eng.data.data3[:, 1]),
                self.eng.data.ztab,
                np.sum(self.eng.data.mztgrid, axis=2),
                np.sum(self.eng.data.mztgrid, axis=1),
                np.sum(self.eng.data.mztgrid, axis=0),
                xlab="m/z (Th)",
                ylab="Arrival Time (ms)",
                zlab="Charge",
                cmap=self.eng.config.cmap,
            )
            print("Finished m/z Cube in: ", time.perf_counter() - start, " s")
        except Exception as exc:
            print("Failed m/z cube", exc)
        try:
            start = time.perf_counter()
            self.view.plot10.cubeplot(
                self.eng.data.massdat[:, 0],
                self.eng.data.ccsdata[:, 0],
                self.eng.data.ztab,
                self.eng.data.massccs,
                self.eng.data.massgrid.reshape((len(self.eng.data.massdat), len(self.eng.data.ztab))),
                self.eng.data.ccsz.transpose(),
                xlab="Mass (Da)",
                ylab="CCS (${\\AA}$$^2$)",
                zlab="Charge",
                cmap=self.eng.config.cmap,
            )
            print("Finished Final Cube in: ", time.perf_counter() - start, " s")
        except Exception as exc:
            print("Failed final cube", exc)

    def on_im_tools(self, e=None):
        self.export_config()
        if ud.isempty(self.eng.data.data3):
            print("Load Data First")
            return
        dlg = IM_wind.IMTools(self.view)
        dlg.initialize_interface(self.eng.data.data3, self.eng.config)
        if dlg.ShowModal() == 0:
            self.import_config(None)
        else:
            self.export_config()

    def on_im_extract(self, e=None):
        if ud.isempty(self.eng.data.ccsdata):
            return
        print("Running UniDec to Generate Outputs")
        self.eng.config.zout = -1
        try:
            self.export_config(self.eng.config.confname)
            ud.unidec_call(self.eng.config)
            dlg = IM_wind.IMToolExtract(self.view)
            dlg.initialize_interface(
                self.eng.data.massdat,
                self.eng.data.ccsdata,
                self.eng.data.massccs,
                self.eng.config,
                self.eng.pks,
            )
            dlg.ShowModal()
        finally:
            self.eng.config.zout = 0

    def on_score(self, e=0):
        pass

    on_score2 = on_score
    on_score_window = on_score
    on_score_label = on_score
    on_score_FDR = on_score

    def on_pdf_report(self, e=None):
        self.view.on_save_figure_pdf(e)
        print("PDF Figures written.")

    def on_flip_twave(self, e):
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

    def remake_mainwindow(self, tabbed=None):
        iconfile = self.view.icon_path
        wx.GetApp().Yield()
        self.view.on_exit()
        self.view = IM_mainwindow.IMMainwindow(
            self,
            "UniDec IM: Ion Mobility-Mass Spectrometry",
            self.eng.config,
            iconfile=iconfile,
            tabbed=tabbed,
        )
        self.view.Show()
        self.view.import_config_to_gui()

    def on_gen_html_report(self, e=None):
        plots = [
            [self.view.plot4, self.view.plot2],
            [self.view.plot3, self.view.plot6],
            [self.view.plot5, self.view.plot1],
            [self.view.plot1im, self.view.plot1fit],
            [self.view.plot2ccs, self.view.plot5mccs],
            [self.view.plot5ccsz, self.view.plot3color],
            [self.view.plot9, self.view.plot10],
        ]
        self.eng.gen_html_report(plots=plots)


# Match the established presenter name while also exposing an explicit one.
UniDecApp = UniDecIMApp


def main():
    multiprocessing.freeze_support()
    app = UniDecIMApp()
    app.start()


if __name__ == "__main__":
    main()

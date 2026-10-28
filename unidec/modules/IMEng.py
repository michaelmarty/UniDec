"""Ion mobility-mass spectrometry engine.

The one-dimensional :mod:`unidec.engine` implementation owns the shared UniDec
workflow.  This subclass replaces only the data loading, processing, result
import, state restoration, and plotting steps whose data shapes differ for
IM-MS.
"""

from copy import deepcopy
import io
import fnmatch
import os
import shutil
import time
import zipfile

import numpy as np

from unidec.engine import UniDec as UniDecMS
import unidec.modules.IM_functions as IM_func
from unidec.modules import peakstructure
from unidec.modules.plotting import plot2d
import unidec.tools as ud


class UniDecIM(UniDecMS):
    """UniDec engine specialized for two-dimensional IM-MS data."""

    mode_imflag = 1

    def reset_config(self):
        super().reset_config()
        self.config.imflag = self.mode_imflag

    def load_config(self, f_name):
        super().load_config(f_name)
        self.config.imflag = self.mode_imflag

    def _load_raw_data(self, importer, time_range=None):
        self.data.rawdata = importer.get_imms_avg_scan(
            mzbins=self.config.mzbins, time_range=time_range
        )
        self.config.discreteplot = 1
        self.config.poolflag = 1
        self.data.rawdata3, self.data.rawdata = ud.unsparse(self.data.rawdata)
        print("Data Shape:", self.data.rawdata3.shape, self.data.rawdata.shape)
        self.data.data3 = self.data.rawdata3
        if getattr(self, "simple_output", False):
            return self.config.output_file("rawdata.npz"), self.data.rawdata3
        return self.config.outfname + "_imraw.txt", ud.sparse(self.data.rawdata3)

    def _load_managed_raw_cache(self, path):
        try:
            with np.load(path, allow_pickle=False) as cache:
                if int(cache["unidec_rawdata_version"]) != 1 or str(cache["mode"]) != "imms":
                    return None
                data = np.asarray(cache["data"])
                original_file = str(cache["original_file"])
        except (OSError, ValueError, KeyError, TypeError):
            return None
        if data.ndim != 2 or data.shape[1] < 3 or not np.isfinite(data[:, :3]).all():
            return None
        self.data.rawdata3, self.data.rawdata = ud.unsparse(data[:, :3])
        self.data.data3 = self.data.rawdata3
        return original_file

    def _restore_processed_data(self, refresh=False):
        # IM processing is inexpensive to identify but shape-sensitive; always
        # begin with the newly loaded 2D data and process it explicitly.
        self.data.data2 = self.data.rawdata
        self.config.procflag = 0

    def raw_process(self, dirname, inflag=False, binsize=1):
        """Convert a Waters RAW directory to sparse IM-MS text data."""
        self.config.dirname = dirname
        self.config.filename = os.path.split(dirname)[1]
        print("Opening: ", self.config.filename)

        extension = os.path.splitext(self.config.filename)[1]
        if extension == ".zip":
            print("Can't open zip, try Load State.")
            return None, None
        if extension.lower() == ".d" and self.config.system == "Windows":
            self.config.dirname = os.path.split(dirname)[0]
            return self.config.filename, self.config.dirname
        if extension.lower() != ".raw" or self.config.system != "Windows":
            print("Error in conversion or file type:", self.config.filename, dirname)
            return None, None

        basename = os.path.splitext(self.config.filename)[0]
        newfilename = basename + "_imraw.txt"
        output_dir = dirname if inflag else os.path.dirname(dirname)
        newfilepath = os.path.join(output_dir, newfilename)
        if os.path.isfile(newfilepath):
            print("Data already converted:", newfilename)
        else:
            start = time.perf_counter()
            call = [
                self.config.cdcreaderpath,
                "-r", dirname,
                "-m", newfilepath[:-10] + "_msraw.txt",
                "-i", newfilepath,
                "--ms_bin", binsize,
                "--ms_smooth_window", "0",
                "--ms_number_smooth", "0",
                "--im_bin", binsize,
                "--sparse", "1",
            ]
            result = ud.exe_call(call)
            print("Time: %.2gs" % (time.perf_counter() - start))
            if result != 0 or not os.path.isfile(newfilepath):
                print("Failed conversion to txt file. ", result, newfilepath)
                return None, None
            print("Converted IM data from raw to txt")

        self.config.filename = newfilename
        self.config.dirname = output_dir
        return self.config.filename, self.config.dirname

    def process_data(self, **kwargs):
        """Process the m/z and arrival-time dimensions of IM-MS data."""
        start = time.perf_counter()
        self.export_config()

        try:
            float(self.config.minmz)
        except ValueError:
            self.config.minmz = np.amin(self.data.rawdata[:, 0])
        try:
            float(self.config.maxmz)
        except ValueError:
            self.config.maxmz = np.amax(self.data.rawdata[:, 0])
        try:
            float(self.config.mindt)
        except ValueError:
            self.config.mindt = np.amin(self.data.rawdata3[:, 1])
        try:
            float(self.config.maxdt)
        except ValueError:
            self.config.maxdt = np.amax(self.data.rawdata3[:, 1])

        if self.check_badness() == 1:
            print("Badness found, aborting data prep")
            return 1

        mz, dt, intensity = IM_func.process_data_2d(
            self.data.rawdata3[:, 0],
            self.data.rawdata3[:, 1],
            self.data.rawdata3[:, 2],
            self.config,
        )
        self.data.data3 = np.transpose(
            [np.ravel(mz), np.ravel(dt), np.ravel(intensity)]
        )
        ud.dataexportbin(self.data.data3, self.config.infname)
        self.config.procflag = 1
        if not kwargs.get("silent"):
            print("Data Prep Time: %.2gs" % (time.perf_counter() - start))

    def unidec_imports(self, efficiency=False, everything=False):
        """Import the two-dimensional outputs produced by the UniDec core."""
        if everything:
            self.data.data3 = np.loadtxt(self.config.infname)
            intensity = self.data.data3[:, 2].reshape(
                (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1])))
            )
            self.data.data2 = np.transpose(
                [np.unique(self.data.data3[:, 0]), np.sum(intensity, axis=1)]
            )

        self.pks = peakstructure.Peaks()
        self.data.massdat = np.loadtxt(self.config.massdatfile)
        self.data.ztab = np.arange(self.config.startz, self.config.endz + 1)
        try:
            self.config.massdatnormtop = np.amax(self.data.massdat[:, 1])
        except Exception:
            self.data.massdat = np.array([self.data.massdat])
            self.config.massdatnormtop = np.amax(self.data.massdat[:, 1])

        if not efficiency:
            try:
                self.data.massgrid = np.fromfile(self.config.massgridfile, dtype=self.config.dtype)
            except Exception:
                pass
            self.data.fitdat = np.fromfile(self.config.fitdatfile, dtype=self.config.dtype)
            self.data.fitdat2d = deepcopy(self.data.data3)
            self.data.fitdat2d[:, 2] = self.data.fitdat
            self.data.fitdat = np.sum(
                self.data.fitdat.reshape(
                    (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1])))
                ),
                axis=1,
            )
            try:
                if self.config.aggressiveflag != 0:
                    self.data.baseline = np.fromfile(
                        self.config.output_file("baseline.bin"), dtype=self.config.dtype
                    )
                else:
                    self.data.baseline = np.array([])
            except Exception:
                self.data.baseline = np.array([])

        runstats = np.genfromtxt(self.config.errorfile, dtype="str")
        self.config.error = float(runstats[1])
        self.data.ccsdata = np.loadtxt(self.config.output_file("ccs.txt"))
        if not efficiency:
            masslen = len(self.data.massdat)
            ccslen = len(self.data.ccsdata)
            zlen = len(self.data.ztab)
            self.data.massccs = np.fromfile(
                self.config.output_file("massccs.bin"), dtype=self.config.dtype
            ).reshape((masslen, ccslen))
            self.data.ccsz = np.fromfile(
                self.config.output_file("ccsz.bin"), dtype=self.config.dtype
            ).reshape((zlen, ccslen))
            self.data.mztgrid = np.fromfile(
                self.config.output_file("mzgrid.bin"), dtype=self.config.dtype
            )
            self.data.mztgrid = np.clip(self.data.mztgrid, 0.0, np.amax(self.data.mztgrid))
            self.data.mztgrid = self.data.mztgrid.reshape(
                (
                    len(np.unique(self.data.data3[:, 0])),
                    len(np.unique(self.data.data3[:, 1])),
                    zlen,
                )
            )
            self.data.mzgrid = np.sum(self.data.mztgrid, axis=1)
            xv, yv = np.meshgrid(self.data.ztab, np.unique(self.data.data3[:, 0]))
            self.data.mzgrid = np.c_[np.c_[np.ravel(yv), np.ravel(xv)], np.ravel(self.data.mzgrid)]

    def convolve_peaks(self):
        return np.array(ud.makeconvspecies(self.data.data2, self.pks, self.config))

    def load_state(self, load_path):
        """Load an IM-MS state archive containing an ``_imraw`` data file."""
        print("Loading Zip File:", load_path)
        extension = "_imraw."
        with zipfile.ZipFile(load_path) as zipf:
            short_cache = next(
                (name for name in zipf.namelist() if os.path.basename(name) == "rawdata.npz"),
                None,
            )
            imfile = next(
                (name for name in zipf.namelist() if fnmatch.fnmatch(name, "*" + extension + "*")),
                None,
            )
            cached_data = None
            if short_cache is not None:
                with np.load(io.BytesIO(zipf.read(short_cache)), allow_pickle=False) as cache:
                    if str(cache["mode"]) != "imms":
                        print("Broken Save File. Raw cache is not IM-MS data")
                        return False
                    cached_data = np.asarray(cache["data"])
                    header = os.path.splitext(os.path.basename(str(cache["original_file"])))[0]
            elif imfile is None:
                print("Broken Save File. Unable to find _imraw")
                return False
            else:
                header = imfile[:-8]
            self.config.dirname = os.path.split(load_path)[0]
            basename = header.rsplit(sep="_", maxsplit=1)[0]
            print("Header:", basename, "Directory:", self.config.dirname)
            output_dir = os.path.join(self.config.dirname, basename + "_unidecfiles")
            os.makedirs(output_dir, exist_ok=True)
            zipf.extractall(output_dir)

        source = os.path.join(output_dir, basename + extension + "txt")
        filename = basename + ".txt"
        destination = os.path.join(self.config.dirname, filename)
        print("Data file:", source, destination)
        if cached_data is not None:
            np.savetxt(destination, cached_data)
        else:
            shutil.copy(source, destination)
        self.open_file(filename, self.config.dirname)

        if os.path.isfile(self.config.infname):
            self.data.data3 = np.loadtxt(self.config.infname)
            intensity = self.data.data3[:, 2].reshape(
                (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1])))
            )
            self.data.data2 = np.transpose(
                [np.unique(self.data.data3[:, 0]), np.sum(intensity, axis=1)]
            )
            self.config.procflag = 1
        else:
            self.config.procflag = 0
        if os.path.isfile(self.config.errorfile):
            self.unidec_imports()
        if os.path.isfile(self.config.peaksfile):
            self.pick_peaks()
        return True

    def makeplot1im(self, plot1im=None, plot1fit=None, imfit=False):
        """Plot IM-MS data and, optionally, its two-dimensional fit."""
        if plot1im is None:
            plot1im = plot2d.Plot2dBase()
        if plot1fit is None:
            plot1fit = plot2d.Plot2dBase()
        try:
            plot1im.contourplot(
                self.data.data3,
                self.config,
                xlab="m/z (Th)",
                ylab="Arrival Time (ms)",
                title="IM-MS Data",
            )
        except Exception:
            pass
        if imfit:
            try:
                plot1fit.contourplot(
                    self.data.fitdat2d,
                    self.config,
                    xlab="m/z (Th)",
                    ylab="Arrival Time (ms)",
                    title="IM-MS Fit",
                )
            except Exception:
                pass


# Conventional module-style and short aliases used by scripts.
UniDec = UniDecIM
IMEng = UniDecIM

# noinspection PyUnresolvedReferences
import os
import time
import shutil
import sys
import getopt
from copy import deepcopy
import zipfile
import fnmatch
import numpy as np

from unidec.modules import unidecstructure, peakstructure, MassFitter
import unidec.tools as ud
import unidec.modules.IM_functions as IM_func
import unidec.modules.MassSpecBuilder as MSBuild
from unidec.modules.unidec_enginebase import UniDecEngine
from unidec.modules.plotting import plot1d, plot2d
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

# import modules.DoubleDec as dd

__author__ = 'Michael.Marty'


def score_minimum(height, minimum):
    x2 = height
    x1 = height / 2.
    if minimum > x1:
        y = 1 - (minimum - x1) / (x2 - x1)
    else:
        y = 1
    return y


class UniDec(UniDecEngine):
    def __init__(self, *args, **kwargs):
        """
        UniDec Engine

        Consists of three main subclasses: Config, DataContiner, Peaks

        :return: None
        """
        UniDecEngine.__init__(self)
        self.autopeaks = None
        self.peakparams = None
        self.massfit = None
        self.massfitdat = None
        self.errorgrid = None
        self.infile = None
        self.outfile = None
        opts = None
        if "ignore_args" in kwargs:
            ignore_args = kwargs["ignore_args"]
        else:
            ignore_args = False
        if not ignore_args:
            try:
                opts, args = getopt.getopt(sys.argv[1:], "f:o:", ["file=", "out="])
            except getopt.GetoptError as e:
                print("Error in Argv. Likely unknown option: ", sys.argv, e)
                print("Known options: -f, -o")

            # print("ARGS:", args)
            # print("KWARGS:", kwargs)
            # print("OPTS:", opts)
            if opts is not None:
                for opt, arg in opts:
                    if opt in ("-f", "--file"):
                        self.infile = arg
                        print("Opening File:", self.infile)
                    if opt in ("-o", "--out"):
                        self.outfile = arg
                        print("Output File:", self.outfile)
            if ud.isempty(opts) or self.infile is None:
                if len(args) > 0:
                    maybe_file = args[0]
                    ext = os.path.splitext(maybe_file)[1]
                    if ext.lower()[1:] in ["mzml", "mzxml", "raw", "hdf5", "txt", "gz", "d"]:
                        self.infile = maybe_file
                        print("Opening File:", self.infile)
            if self.infile is not None:
                self.open_file(self.infile)
                self.autorun()
            pass

    def open_file(self, file_name, file_directory=None, time_range=None, refresh=False, load_results=False,
                 isodeceng=None, *args, **kwargs):
        """
        Open text or mzML file. Will create _unidecfiles directory if it does not exist.

        If it finds a _conf.dat file in _unidecfiles, it will import the old configuration.
        Otherwise, it keeps the existing configuration but resets file names.

        If silent=True is passed in kwargs, prints are suppressed.
        :param: file_name: Name of file to open.
        May be in  x y or x y z text format or in mzML format. May be tab or space delimited
        :param file_directory: Directory in which filename is located.
        Default is current directory. If None, will get from file_name.
        :param time_range: Array or tuple of (start, end) times to load from data file. Default is to load all.
        :param refresh: If True, will refresh the data file even if it already exists.
        :param load_results: If True, will load the results from the _unidecfiles directory.
        :return: None
        """
        tstart = time.perf_counter()
        self.pks = peakstructure.Peaks()
        self.data = unidecstructure.DataContainer()
        # Get the directory and file names
        if file_directory is None:
            file_directory = os.path.dirname(file_name)
            file_name = os.path.basename(file_name)
        if '.d' not in file_directory and '.raw' not in file_directory:
            file_path = os.path.join(file_directory, file_name)
        else:
            try:
                file_path = file_directory
            except:
                raise IOError("File path not found")
        # Handle Paths

        self.config.filename = file_path
        self.config.dirname = file_directory
        if "silent" not in kwargs or not kwargs["silent"]:
            print("Opening File: ", self.config.filename)

        if self.outfile is None:
            dirnew = os.path.splitext(file_path)[0] + "_unidecfiles"
            basename = os.path.split(os.path.splitext(file_name)[0])[1]
        else:
            splits = os.path.split(self.outfile)
            if splits[0] != '':
                dirnew = os.path.abspath(splits[0])
            else:
                dirnew = os.getcwd()
            basename = splits[1]

        if "clean" in kwargs and kwargs["clean"] and os.path.isdir(dirnew):
            shutil.rmtree(dirnew)

        print(self.config.udir)
        if not os.path.isdir(dirnew):
            os.mkdir(dirnew)
        self.config.udir = dirnew
        if "silent" not in kwargs or not kwargs["silent"]:
            print("Output Directory:", self.config.udir)
        self.config.outfname = os.path.join(self.config.udir, basename)
        self.config.extension = os.path.splitext(self.config.filename)[1]
        self.config.default_file_names()

        #All the magic happens here
        curr_importer = ImporterFactory.create_importer(self.config.filename)
        if isodeceng is not None:
            isodeceng.reader = curr_importer

        if self.config.imflag ==0:
            self.data.rawdata = curr_importer.get_avg_scan(time_range=time_range)

            newname = self.config.outfname + "_rawdata.txt"
            outputdata = self.data.rawdata

        else:
            self.data.rawdata = curr_importer.get_imms_avg_scan(mzbins=self.config.mzbins, time_range=time_range)
            self.config.discreteplot = 1
            self.config.poolflag = 1
            self.data.rawdata3, self.data.rawdata = ud.unsparse(self.data.rawdata)
            print("Data Shape:", self.data.rawdata3.shape, self.data.rawdata.shape)
            self.data.data3 = self.data.rawdata3

            newname = self.config.outfname + "_imraw.txt"
            outputdata = ud.sparse(self.data.rawdata3)

        if ud.isempty(self.data.rawdata):
            print("Error: Data Array is Empty")
            print("Likely an error with data conversion")
            raise ImportError

        if not os.path.isfile(newname):
            try:
                # shutil.copy(file_directory, newname)
                np.savetxt(newname, outputdata)
            except Exception as e:
                pass

        if os.path.isfile(self.config.infname) and not refresh and self.config.imflag == 0:
            try:
                self.data.data2 = np.loadtxt(self.config.infname)
                self.config.procflag = 1
            except:
                self.data.data2 = self.data.rawdata
                self.config.procflag = 0
        else:
            self.data.data2 = self.data.rawdata
            self.config.procflag = 0

        # Initialize Config
        if os.path.isfile(self.config.confname) and not refresh:
            self.load_config(self.config.confname)
        else:
            self.export_config()

        self.auto_polarity(importer = curr_importer)

        if load_results:
            self.unidec_imports(everything=True)
            self.pick_peaks()

        tend = time.perf_counter()
        if "silent" not in kwargs or not kwargs["silent"]:
            print("Loading Time: %.2gs" % (tend - tstart))

    def raw_process(self, dirname, inflag=False, binsize=1):
        """
        Processes Water's Raw files into .txt using external calls to:
        self.config.cdcreaderpath for IM-MS

        MS is processed by modules.waters_importer.Importer

        Default files are created with the header of the .raw file plus:
        _rawdata.txt for MS
        _imraw.txt for IM-MS
        :param dirname: .raw directory name
        :param inflag: If True, it will put the output .txt file inside the existing .raw directory. If False, it will
        put the file in the same directory that contains the .raw directory
        :param binsize: Parameter for IM-MS to specify the m/z bin size for conversion. If binsize=0, the conversion
        will be at full resolution (which is huge), so the default is every 1 m/z.
        :return: self.config.filename, self.config.dirname (name and location of created file)
        """
        self.config.dirname = dirname
        self.config.filename = os.path.split(self.config.dirname)[1]
        print("Opening: ", self.config.filename)
        if os.path.splitext(self.config.filename)[1] == ".zip":
            print("Can't open zip, try Load State.")
            return None, None

        elif os.path.splitext(self.config.filename)[1].lower() == ".d" and self.config.system == "Windows":
            self.config.dirname = os.path.split(self.config.dirname)[0]
            return self.config.filename, self.config.dirname

        elif os.path.splitext(self.config.filename)[1] == ".raw" and self.config.system == "Windows":
            basename = os.path.splitext(self.config.filename)[0]

            if self.config.imflag == 1:
                newfilename = basename + "_imraw.txt"
            else:
                newfilename = basename + "_rawdata.txt"

            if inflag:
                newfilepath = os.path.join(self.config.dirname, newfilename)
            else:
                newfilepath = os.path.join(os.path.dirname(self.config.dirname), newfilename)

            if os.path.isfile(newfilepath):
                self.config.filename = newfilename
                print("Data already converted:", newfilename)
            else:
                if self.config.system == "Windows":
                    if self.config.imflag == 1:
                        stime = time.perf_counter()
                        call = [self.config.cdcreaderpath, '-r', self.config.dirname, '-m',
                                newfilepath[:-10] + "_msraw.txt", '-i', newfilepath, '--ms_bin', binsize,
                                "--ms_smooth_window", "0", "--ms_number_smooth", "0", "--im_bin", binsize, "--sparse",
                                "1"]
                        result = ud.exe_call(call)
                        print("Time: %.2gs" % (stime - time.perf_counter()))

                        self.config.filename = newfilename
                        if result == 0 and os.path.isfile(newfilepath):
                            print("Converted IM data from raw to txt")
                        else:
                            print("Failed conversion to txt file. ", result, newfilepath)
                            return None, None
                else:
                    print("Sorry. Waters Raw converter only works on windows. Convert to txt file first.")
                    return None, None

            return self.config.filename, self.config.dirname

        else:
            print("Error in conversion or file type:", self.config.filename, self.config.dirname)
            return None, None


    def process_data(self, **kwargs):
        """
        Process data according to parameters in config.

        Checks certain parameters to make sure the limits make sense.
        Will accept silent=True kwarg to suppress printing.
        :return: None
        """
        tstart = time.perf_counter()
        self.export_config()

        try:
            float(self.config.minmz)
        except ValueError:
            self.config.minmz = np.amin(self.data.rawdata[:, 0])

        try:
            float(self.config.maxmz)
        except ValueError:
            self.config.maxmz = np.amax(self.data.rawdata[:, 0])
        # print("Min MZ: ", self.config.minmz, "Max MZ: ", self.config.maxmz)
        if self.config.imflag == 1:
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

        if self.config.imflag == 0:
            self.data.data2 = ud.dataprep(self.data.rawdata, self.config)
            if "scramble" in kwargs:
                if kwargs["scramble"]:
                    # np.random.shuffle(self.data.data2[:, 1])
                    self.data.data2[:, 1] = np.abs(
                        np.random.normal(0, 100 * np.amax(self.data.data2[:, 1]), len(self.data.data2)))
                    self.data.data2[:, 1] /= np.amax(self.data.data2[:, 1])
                    print("Added noise to data")
            ud.dataexport(self.data.data2, self.config.infname)
        else:
            tstart2 = time.perf_counter()
            mz, dt, i3 = IM_func.process_data_2d(self.data.rawdata3[:, 0], self.data.rawdata3[:, 1],
                                                 self.data.rawdata3[:, 2],
                                                 self.config)
            tend = time.perf_counter()
            if "silent" not in kwargs or not kwargs["silent"]:
                print("Time: %.2gs" % (tend - tstart2))
            self.data.data3 = np.transpose([np.ravel(mz), np.ravel(dt), np.ravel(i3)])
            self.data.data2 = np.transpose([np.unique(mz), np.sum(i3, axis=1)])
            ud.dataexportbin(self.data.data3, self.config.infname)
            pass

        self.config.procflag = 1
        tend = time.perf_counter()
        if "silent" not in kwargs or not kwargs["silent"]:
            print("Data Prep Time: %.2gs" % (tend - tstart))
        # self.get_spectrum_peaks()
        pass

    def run_unidec(self, silent=False, efficiency=False):
        """
        Runs unidec.

        Checks that everything is set to go and then places external call to:
            self.config.UniDecPath for MS
            self.config.UniDecIMPath for IM-MS

        If successful, calls self.unidec_imports()
        If not, prints the error code.
        :param silent: If True, it will suppress printing the output from unidec
        :param efficiency: Passed to self.unidec_imports()
        :return: out (stdout from external unidec call)
        """
        # Check to make sure everything is in order
        if self.config.procflag == 0:
            print("Need to process data first. Processing...")
            self.process_data()
        if self.check_badness() == 1:
            print("Badness found, aborting unidec run")
            return 1
        if self.config.doubledec:
            kpath = self.config.kernel
            try:
                with open(kpath, "r") as f:
                    pass
            except (IOError, FileNotFoundError) as err:
                print("Could not open kernel file.\nPlease select a valid kernel file to use DoubleDec")
                return 0
        # Export Config and Call
        self.export_config()
        tstart = time.perf_counter()

        out = ud.unidec_call(self.config, silent=silent)

        tend = time.perf_counter()
        self.config.runtime = (tend - tstart)
        if not silent:
            print("unidec run %.2gs" % self.config.runtime)
        # Import Results if Successful
        if out == 0:
            self.unidec_imports(efficiency)
            if not silent:
                print("File Name: ", self.config.filename, "R Squared: ", self.config.error)
            return out
        else:
            print("unidec Run Error:", out)
            return out

    def unidec_imports(self, efficiency=False, everything=False):
        """
        Imports files output from the unidec core executable into self.data.
        :param efficiency: If True, it will ignore the larger files to speed up the run.
        :param everything: If True, it will import all files, even the large ones.
        :return: None
        """
        if everything:
            self.data.data2 = np.loadtxt(self.config.infname)

        # Import Results
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
            try:
                if self.config.aggressiveflag != 0:
                    self.data.baseline = np.fromfile(self.config.outfname + "_baseline.bin", dtype=self.config.dtype)
                else:
                    self.data.baseline = np.array([])
            except Exception as e:
                self.data.baseline = np.array([])
                pass

            if self.config.imflag == 1:
                self.data.fitdat2d = deepcopy(self.data.data3)
                self.data.fitdat2d[:, 2] = self.data.fitdat
                self.data.fitdat = np.sum(self.data.fitdat.reshape(
                    (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1])))), axis=1)

        try:
            runstats = np.genfromtxt(self.config.errorfile, dtype='str')
        except Exception:
            runstats = []
            pass
        if self.config.imflag == 0:
            # Calculate Error
            try:
                sse = float(runstats[0, 2])
            except Exception as e:
                print("Error in error import", e)
                print(runstats)
                sse = 1

            mean = np.mean(self.data.data2[:, 1])
            self.config.error = 1 - sse / np.sum((self.data.data2[:, 1] - mean) ** 2)
            if not efficiency:
                # Import Grid
                try:
                    self.data.mzgrid = np.fromfile(self.config.mzgridfile, dtype=self.config.dtype)
                    xv, yv = np.meshgrid(self.data.ztab, self.data.data2[:, 0])
                    xv = np.c_[np.ravel(yv), np.ravel(xv)]
                    self.data.mzgrid = np.c_[xv, self.data.mzgrid]
                except Exception as e:
                    print("Error: Mismatched dimensions between processed and deconvolved data. ", e)
                    self.data.mzgrid = []

            for r in runstats:
                if r[0] == "avgscore":
                    self.config.avgscore = float(r[2])

        else:
            # Calculate Error
            self.config.error = float(runstats[1])

            self.data.ccsdata = np.loadtxt(self.config.outfname + "_ccs.txt")
            if not efficiency:
                # Import Grids and Reshape
                masslen = len(self.data.massdat)
                ccslen = len(self.data.ccsdata)
                zlen = len(self.data.ztab)

                self.data.massccs = np.fromfile(self.config.outfname + "_massccs.bin", dtype=self.config.dtype)
                self.data.ccsz = np.fromfile(self.config.outfname + "_ccsz.bin", dtype=self.config.dtype)
                self.data.mztgrid = np.fromfile(self.config.outfname + "_mzgrid.bin", dtype=self.config.dtype)

                self.data.massccs = self.data.massccs.reshape((masslen, ccslen))
                self.data.ccsz = self.data.ccsz.reshape((zlen, ccslen))
                self.data.mztgrid = np.clip(self.data.mztgrid, 0.0, np.amax(self.data.mztgrid))
                self.data.mztgrid = self.data.mztgrid.reshape(
                    (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1])), zlen))
                self.data.mzgrid = np.sum(self.data.mztgrid, axis=1)
                xv, yv = np.meshgrid(self.data.ztab, np.unique(self.data.data3[:, 0]))
                xv = np.c_[np.ravel(yv), np.ravel(xv)]
                self.data.mzgrid = np.c_[xv, np.ravel(self.data.mzgrid)]

    def pick_peaks(self, calc_dscore=True):
        """
        Detect, Normalize, and Output Peaks
        :return: None
        """
        self.export_config()
        # Detect Peaks and Normalize
        peaks = ud.peakdetect(self.data.massdat, self.config)
        if len(peaks) > 0:
            self.setup_peaks(peaks)
            if calc_dscore:
                try:
                    self.dscore()
                except:
                    pass
        else:
            print("No peaks detected", peaks, self.config.peakwindow, self.config.peakthresh)
            print("Mass Data:", self.data.massdat)
        return peaks

    def setup_peaks(self, peaks):
        """
        Setup Peaks from a list of peaks, mostly by normalizing everything but also generating mztab and error values.
        :param peaks: peaks to setup
        :return: normalization parameter
        """
        norm = 1
        if self.config.peaknorm == 1:
            norm = np.amax(peaks[:, 1]) / 100.
        elif self.config.peaknorm == 2:
            norm = np.sum(peaks[:, 1]) / 100.
        else:
            if self.config.massdatnormtop == 0:
                self.config.massdatnormtop = np.amax(self.data.massdat[:, 1])
            norm = np.amax(peaks[:, 1]) / self.config.massdatnormtop

        if norm == 0:
            norm = 1
            print("Warning: Normalization is 0. Setting to 1.")

        peaks[:, 1] = peaks[:, 1] / norm
        self.data.massdat[:, 1] = self.data.massdat[:, 1] / norm
        self.config.massdatnorm = self.config.massdatnormtop / np.amax(peaks[:, 1])

        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(peaks, massbins=self.config.massbins)
        self.pks.default_params(cmap=self.config.peakcmap)
        ud.dataexport(peaks, self.config.peaksfile)
        # Generate Intensities of Each Charge State for Each Peak
        try:
            mztab = ud.make_peaks_mztab(self.data.mzgrid, self.pks, self.config.adductmass)
        except Exception:
            mztab = []
            pass
        # Calculate errors for peaks with FWHM
        try:
            ud.peaks_error_FWHM(self.pks, self.data.massdat)
        except Exception as e:
            print("Error in FWHM calculation:", e)
            print(self.data.massdat)
        try:
            ud.peaks_error_mean(self.pks, self.data.massgrid, self.data.ztab, self.data.massdat, self.config)
        except Exception as e:
            print("Error in error calculations:", e)
        if self.config.batchflag == 0:
            try:
                ud.make_peaks_mztab_spectrum(self.data.mzgrid, self.pks, self.data.data2, mztab)
                self.export_config()
            except Exception:
                pass

        return norm

    def convolve_peaks(self):
        """
        Convolve Peaks with Peak Shape
        :return: None
        """
        if self.config.imflag == 1 or self.config.cdmsflag == 1:
            convdata = ud.makeconvspecies(self.data.data2, self.pks, self.config)
        else:
            # TODO: There's no reason this shouldn't work for CD-MS data, but we'd need to include a write to _grid.bin
            ud.unidec_call(self.config, conv=True)

            convdata = np.fromfile(self.config.outfname + "_conv.bin", dtype=self.config.dtype)
            xlen = len(self.data.data2)
            convdata = convdata.reshape((len(self.pks.peaks), xlen))

            self.pks.composite = np.zeros(xlen)
            for i in range(0, self.pks.plen):
                sd = convdata[i]
                self.pks.peaks[i].stickdat = sd
                self.pks.composite += np.array(sd)
            self.pks.convolved = True
        return np.array(convdata)

    def autorun(self, auto_peak_width=True, silent=False):
        self.process_data(silent=silent)
        if auto_peak_width:
            self.get_auto_peak_width()
        self.run_unidec(silent=silent)
        self.pick_peaks()
        self.autointegrate()
        self.export_params(silent=silent)

    def estimate_FDR(self):
        print("Starting FDR Estimate")
        dscores = []
        maxit = 50
        maxdscores = 300
        i = 0
        while i < maxit and len(dscores) < maxdscores:
            self.process_data(scramble="True", silent=False)
            self.run_unidec(silent=False)
            self.pick_peaks()
            self.dscore()
            for p in self.pks.peaks:
                dscores.append(p.dscore)
            i += 1
            print(i)

        self.autorun()

        print("Number of False Scores:", len(dscores))
        sd = np.sort(dscores)
        length = len(sd)
        fdrs = [0, 0.01, 0.05]
        cutoffs = []

        for fdr in fdrs:
            f = int(length * (1 - fdr)) + 1
            if f >= length:
                f = length - 1
            cutoffs.append(sd[f])
            print("FDR (%):", fdr * 100, "DScore Cut Off (%):", sd[f] * 100)

        self.fdrs = np.transpose([fdrs, cutoffs])

    def autocorrelation(self, massdat=None):
        """
        Performs autocorrelation on mass data. Result is stored as self.data.autocorr.
        Picks peaks greater than 0 using peak detection parameters in config file.
        Peaks are stored as a peak structure at self.autopeaks

        :param: massdat Data on which to run autocorrelation. Default is None, in which case self.data.massdat is used.
        :return: float. First peak in autocorrelation.
        """
        if massdat is None:
            massdat = self.data.massdat
        # corr=np.correlate(self.data.massdat[:,1],self.data.massdat[:,1],mode="same")
        self.data.autocorr, cpeaks = ud.autocorr(massdat, self.config)
        self.autopeaks = peakstructure.Peaks()
        self.autopeaks.add_peaks(cpeaks, massbins=self.config.massbins)
        self.autopeaks.default_params()
        print("Autocorrelation:", [p.mass for p in self.autopeaks.peaks])
        return self.autopeaks.peaks[0].mass

    def kendrick_peaks(self, kmass=None, centermode=1):
        """
        Run Kendrick analysis on peaks (self.pks object)
        :param kmass: Kendrick mass. Default is prior kendrick mass if it exists and is >0.
        Otherwise, default is oligomer mass (self.config.molig)
        :param centermode: Set range for normalization 1=(0,1),0=(-0.5,0.5)
        :return: Array of [mass,defect] for each peak in self.pks.
        """
        if kmass is not None:
            self.config.kendrickmass = kmass
        if not self.config.kendrickmass > 0:
            self.config.kendrickmass = self.config.molig
        if self.config.kendrickmass > 0:
            self.pks.get_mass_defects(self.config.kendrickmass, mode=centermode)
            return np.array([[p.mass, p.kendrickdefect] for p in self.pks.peaks])
        else:
            print("Need non-zero Kendrick mass")
            return None

    def kendrick_continuous(self, ref_mass=None, centermode=0, nbins=50, transformmode=0, xaxistype=1):
        """
        Runs continuous Kendrick analysis on self.data.massdat
        :param ref_mass: Kendrick mass. Default is self.config.kendrickmass if it is already set and >0.
        Otherwise, default is oligomer mass (self.config.molig)
        :param centermode: Set range for normalization 0=(0,1),1=(-0.5,0.5). Default is 0.
        :param nbins: Set mass defect axis density. Default is 50 bins.
        :param transformmode: Set type of transformation. 0=Interpolation. 1=Integration. Default is 0.
        :param xaxistype: Set x-axis dimensions. 0=Kendrick Mass Number, 1=Mass Number * Kendrick Mass. Default is 1.
        :return: mass grid, mass defect grid, intensity grid. All with shape (len(self.data.massdat),nbins)
        """
        if ref_mass is not None:
            self.config.kendrickmass = ref_mass
        if not self.config.kendrickmass > 0:
            self.config.kendrickmass = self.config.molig
        if self.config.kendrickmass > 0:

            data1, data2, m1grid, m2grid, igrid = ud.kendrick_analysis(self.data.massdat, self.config.kendrickmass,
                                                                       centermode=centermode, nbins=nbins,
                                                                       transformmode=transformmode,
                                                                       xaxistype=xaxistype)
            # Write outputs
            outfile2 = os.path.join(self.config.outfname + "_2D_Mass_Defects.txt")
            outfile1 = os.path.join(self.config.outfname + "_1D_Mass_Defects.txt")
            np.savetxt(outfile2, data2)
            np.savetxt(outfile1, data1)
            print("Saved Mass Defects:", outfile2, outfile1)
            return m1grid, m2grid, igrid
        else:
            print("Need non-zero Kendrick mass")
            return None, None, None

    def mass_grid_to_f_grid(self):
        """
        Convert the mass vs charge grid to a mass vs charge offset grid.

        Calculates the charge offset for each (mass,charge) point, creates a new axis of regularly spaced charge
        offsets (oaxis), and the interpolates a new grid of (mass, offset) from oaxis, which is output as outgrid.
        :return: oxais, outgrid: offset axis (N) and offset grid (M x N)
        """
        mgrid, zgrid = np.meshgrid(self.data.massdat[:, 0], np.array(self.data.ztab), indexing="ij")
        ogrid = ud.get_z_offset(mgrid, zgrid)
        oaxis = np.arange(np.amin(ogrid), np.amax(ogrid), 0.5)
        mgrid2, ogrid2 = np.meshgrid(self.data.massdat[:, 0], oaxis, indexing="ij")
        massgrid = self.data.massgrid.reshape((len(self.data.massdat[:, 0]), len(self.data.ztab)))
        outgrid = ud.mergedata2d(mgrid2, ogrid2, mgrid, ogrid, massgrid)
        outgrid -= np.amin(outgrid)
        outgrid /= np.amax(outgrid)
        return oaxis, outgrid

    def autointegrate(self, ztab=None):
        """
        Perform automatic integration of peaks.

        If self.config.integrateup is empty, the upperbound becomes self.config.peakwindow.
        If self.config.integratelb is empty, the lowerbound becomes -self.config.peakwindow.

        Integral range for each peak is set to peak.integralrange.
        Integral value is set to peak.integral.

        If ztab parameter is set to a list of charge states, it will integrate the mass vs charge grid at each
        individual charge state. Otherwise, this is ignored.
        :param ztab: List of charge states (default = None)
        :return: zarea: P x Z array where P is the number of peaks and Z is the number of charge states.
        Each value of the array is the integral of peak P at charge state Z.
        """
        if self.config.integrateub == "":
            ub = self.config.peakwindow
        else:
            ub = self.config.integrateub
        if self.config.integratelb == "":
            lb = -self.config.peakwindow
        else:
            lb = self.config.integratelb
        zarea = []
        for p in self.pks.peaks:
            p.integralrange = [p.mass + lb, p.mass + ub]
            p.integral = self.integrate(p.integralrange)[0]
            zlist = []
            if ztab is not None:
                for i in range(0, len(ztab)):
                    integral = self.integrate(p.integralrange,
                                              data=np.reshape(self.data.massgrid, (len(self.data.massdat), len(ztab)))[
                                                   :, i])[0]
                    zlist.append(integral)
            zarea.append(zlist)

        self.normalize_peaks()
        return np.array(zarea)

    def export_params(self, e=None, silent=False):
        """
        Export a number of different parameters about the peaks into different text files.
        :param e: if e is "PostFit", it will output mass fit parameters as well
        :param silent: if True, it will not print the output file names
        :return: None
        """
        if self.pks.plen > 0:
            # Export Peaks Height by Charge Grid
            mztab = np.array([p.mztab for p in self.pks.peaks])
            ud.dataexport(mztab[:, :, 1], self.config.outfname + "_chargedata.dat")
            if not silent:
                print("Exported data to " + self.config.outfname + "_chargedata.dat")

            ud.dataexport(mztab[:, :, 0], self.config.outfname + "_mzpeakdata.dat")
            if not silent:
                print("Exported data to " + self.config.outfname + "_mzpeakdata.dat")

            # Export Peaks Integral by Charge Grid
            if self.config.batchflag == 0:
                try:
                    chargeareas = self.autointegrate(ztab=self.data.ztab)
                    ud.dataexport(chargeareas, self.config.outfname + "_chargedata_areas.dat")
                except (IndexError, ValueError, AttributeError, ZeroDivisionError):
                    print("Unable to autointegrate")

            # Get Params
            peaks = np.array([[p.mass, p.height] for p in self.pks.peaks])
            try:
                self.autointegrate()
                areas = [p.integral for p in self.pks.peaks]
            except (IndexError, ValueError, AttributeError, ZeroDivisionError):
                areas = peaks[:, 1]
                print("Failed to integrate. Substituting heights for areas.")

            peakparams = []
            for i in range(0, len(peaks)):
                try:
                    avg = np.average(self.data.ztab, weights=mztab[i, :, 1])
                    std = np.sqrt(np.average((np.array(self.data.ztab) - avg) ** 2, weights=mztab[i, :, 1]))
                except:
                    avg = -1
                    std = -1
                p = self.pks.peaks[i]
                if e == "PostFit":
                    peakparams.append(
                        [peaks[i, 0], self.config.mzsig * avg, avg, std, peaks[i, 1] / np.sum(peaks[:, 1]),
                         self.massfit[i, 1], self.massfit[i, 2] / np.sum(self.massfit[:, 2]),
                         p.centroid, p.errorFWHM, p.errorreplicate, p.dscore])
                else:
                    peakparams.append([peaks[i, 0], self.config.mzsig * avg, avg, std, peaks[i, 1], areas[i],
                                       p.centroid, p.errorFWHM, p.errormean, p.dscore])
            self.peakparams = np.array(peakparams)

            header = "Mass MassStdGuess AvgCharge StdDevCharge Height Area MassCentroid " \
                     "MassFWHM MassErrorBetweenZ DScore"
            outfile = self.config.outfname + "_peakparam.dat"
            ud.dataexport(self.peakparams, outfile, header)
            if not silent:
                print(header)
                np.set_printoptions(precision=2, formatter={'float': '{: 0.2f}'.format})
                print(self.peakparams)
                np.set_printoptions()
                print("Peak Parameters (Saved To", outfile, ")")
        else:
            print("Pick Peaks First")

            # TODO: Streamline to remove multiple integration steps
            # TODO: Rework params into peakstructure
            # TODO: Better docstring

    def process_mass_data(self):
        """
        Apply the same parameters used to process the data to process the mass distribution. Linearization parameters
        are ignored, but smoothing, baseline subtraction, normalization, and intensity threshold all apply.
        :return: None
        """
        self.pks = peakstructure.Peaks()
        if self.config.smooth > 0:
            self.data.massdat = ud.gsmooth(self.data.massdat, self.config.smooth)
        # Baseline Subtraction
        buff = abs(self.config.subbuff)
        subtype = self.config.subtype
        if subtype == 1 and buff != 0:
            self.data.massdat = ud.datasimpsub(self.data.massdat, buff)
        elif subtype == 2 and buff != 0:
            self.data.massdat = ud.datacompsub(self.data.massdat, buff)
        elif subtype == 0 and buff != 0:
            self.data.massdat[:, 1] = self.data.massdat[:, 1] - np.amin(self.data.massdat[:, 1])
        # Normalization
        self.data.massdat = ud.normalize(self.data.massdat)
        # Intensity Threshold
        self.data.massdat = ud.intensitythresh(self.data.massdat, self.config.intthresh)  # thresh

    def center_of_mass(self, data=None, limits=None):
        """
        Return the center of mass and weighted standard deviation for data within some limits. If data is None,
        self.data.massdat is used. If limits is None, the whole range is used.
        :param data: mass data to determine center of mass
        :param limits: limits to restrict the calculation
        :return: com, std (center of mass, weighted standard deviation)
        """
        if data is None:
            data = self.data.massdat
        if limits is None:
            com = np.average(data[:, 0], weights=data[:, 1])
            std = ud.weighted_std(data[:, 0], data[:, 1])
        else:
            com, std = ud.center_of_mass(data, limits[0], limits[1])
        return com, std

    def fit_all_masses(self):
        """
        Fit all masses to a series of peaks, with initial guesses defined by the peak parameters.
        :return: self.massfitdat, self.massfit (fit to data, fit parameters)
        """
        self.massfitdat, self.massfit = MassFitter.MassFitter(self.data.massdat, self.peakparams,
                                                              self.config.psfun).perform_fit("nonorm", "sort")
        return self.massfitdat, self.massfit

    def get_charge_peaks(self):
        """
        Determines total charge distribution. Imports each charge state as a peak in self.pks.
        Will overwrite mass peaks.
        :return: cpeaks (Z x 2 array of (charge state, intensity))
        """
        if not ud.isempty(self.data.mzgrid):
            dat = self.data.mzgrid
            c = dat[:, 2]
            xlen = len(np.unique(dat[:, 0]))
            ylen = len(np.unique(dat[:, 1]))
            newgrid = np.reshape(c, (xlen, ylen))

            cint = np.sum(newgrid, axis=0)
            if self.config.peaknorm == 1:
                cint = cint / np.amax(cint) * 100.
            elif self.config.peaknorm == 2:
                cint = cint / np.sum(cint) * 100.

            cpeaks = np.transpose([self.data.ztab, cint])
            np.savetxt(self.config.outfname + "_chargepeaks.txt", cpeaks)
            # com, std = self.center_of_mass(data=cpeaks)
            self.pks = peakstructure.Peaks()
            self.pks.add_peaks(cpeaks, massbins=1)
            self.pks.default_params(self.config.peakcmap)
            for i, p in enumerate(self.pks.peaks):
                p.stickdat = newgrid[:, i]
                p.label = str(int(self.data.ztab[i]))
            return cpeaks
        else:
            print("Error: no m/z grid.")
            return None

    def save_state(self, file_name):
        ud.zip_folder(file_name, directory=self.config.udir)

    def load_state(self, load_path):
        """
        Load unidec state from a zip save file.

        Note: save_state is located under unidectools (ud.savestate)
        :param load_path: .zip file to load
        :return: True is successful, False if failed
        """
        print("Loading Zip File:", load_path)
        # Set up extensions
        extension = "_rawdata."
        extension2 = "_imraw."
        # In zip file, search for correct files
        zipf = zipfile.ZipFile(load_path)
        imfile = None
        msfile = None
        for file_path in zipf.namelist():
            if fnmatch.fnmatch(file_path, '*' + extension + "*") or fnmatch.fnmatch(file_path, '*' + extension2 + "*"):
                if fnmatch.fnmatch(file_path, '*' + extension + "*"):
                    msfile = file_path
                elif fnmatch.fnmatch(file_path, '*' + extension2 + "*"):
                    imfile = file_path

        # Set file and extension
        header = None
        if imfile is not None:
            if msfile is None or self.config.imflag == 1:
                header = imfile[:-8]
                extension = "_imraw."
        elif msfile is not None:
            if imfile is None or self.config.imflag == 0:
                header = msfile[:-10]
                extension = "_rawdata."
        else:
            print("Broken Save File. Unable to find _rawdata or _imraw")
            return False

        # Get directory, filename, and header
        self.config.dirname = os.path.split(load_path)[0]
        basename = header.rsplit(sep="_", maxsplit=1)[0]
        print("Header:", basename, "Directory:", self.config.dirname)

        # Setup default file names, unidecfile directory, and extract there
        dirnew = os.path.join(self.config.dirname, basename + "_unidecfiles")
        flag = os.path.isdir(dirnew)
        if not flag:
            os.mkdir(dirnew)
        zipf.extractall(dirnew)

        # Copy data file from unidecfiles to directory above it
        file_name = os.path.join(dirnew, basename + extension + "txt")
        # if not os.path.isfile(file_name):
        #    file_name = self.config.outfname + extension + "dat"
        filename2 = basename + ".txt"
        load_path = os.path.join(self.config.dirname, filename2)
        print("Data file:", file_name, load_path)
        shutil.copy(file_name, load_path)

        # Open File
        self.open_file(filename2, self.config.dirname)

        # Import Processed Data
        if os.path.isfile(self.config.infname):
            if self.config.imflag == 0:
                self.data.data2 = np.loadtxt(self.config.infname)
            else:
                self.data.data3 = np.loadtxt(self.config.infname)
                i3 = self.data.data3[:, 2].reshape(
                    (len(np.unique(self.data.data3[:, 0])), len(np.unique(self.data.data3[:, 1]))))
                self.data.data2 = np.transpose([np.unique(self.data.data3[:, 0]), np.sum(i3, axis=1)])
            self.config.procflag = 1
        else:
            self.config.procflag = 0

        # Import unidec Results
        if os.path.isfile(self.config.errorfile):
            self.unidec_imports()

        # Import Peaks
        if os.path.isfile(self.config.peaksfile):
            self.pick_peaks()

        return True
        # TODO: Import Matches, others things in state?

    def normalize_peaks(self):
        """
        Noamlize everything in the peaks accoring to the type set in self.config.peaknorm
        0 = No normalization
        1 = Normalize the max value to 1
        2 = Normalize the sum to 1
        :return: None
        """
        if len(self.pks.peaks) == 0:
            print("No Peaks Detected")
            return

        integrals = np.array([p.integral for p in self.pks.peaks])
        heights = np.array([p.height for p in self.pks.peaks])
        # corrints = np.array([p.corrint for p in self.pks.peaks])
        # fitareas = np.array([p.fitarea for p in self.pks.peaks])
        if self.config.peaknorm == 1:
            inorm = np.amax(integrals) / 100.
            hnorm = np.amax(heights) / 100.
            # cnorm = np.amax(corrints) / 100.
            # fnorm = np.amax(fitareas) / 100.
        elif self.config.peaknorm == 2:
            inorm = np.sum(integrals) / 100.
            hnorm = np.sum(heights) / 100.
            # cnorm = np.sum(corrints) / 100.
            # fnorm = np.sum(fitareas) / 100.
        else:
            inorm = 1.
            hnorm = 1.
            # cnorm = 1.
            # fnorm = 1.

        if inorm != 0:
            for p in self.pks.peaks:
                p.integral /= inorm

        if hnorm != 0:
            for p in self.pks.peaks:
                p.height /= hnorm

        # if cnorm != 0:
        #    for p in self.pks.peaks:
        #        p.corrint /= cnorm
        #        p.correrr /= cnorm

        # if fnorm != 0:
        #    for p in self.pks.peaks:
        #        p.fitarea /= fnorm
        #        p.fitareaerr /= fnorm

    def get_zstack(self, xfwhm=1):
        zarr = np.reshape(self.data.massgrid, (len(self.data.massdat), len(self.data.ztab)))
        # zarr = zarr / np.amax(np.sum(zarr, axis=1)) * np.amax(self.data.massdat[:, 1])

        for i, p in enumerate(self.pks.peaks):
            # fwhm = p.errorFWHM
            # if fwhm == 0:
            #    fwhm = self.config.massbins
            # width = fwhm / 2. * xfwhm
            # interval = [p.mass - width, p.mass + width]

            interval = np.abs(np.array(p.intervalFWHM) - p.mass) * xfwhm
            if interval[0] == 0:
                interval[0] = self.config.massbins * xfwhm
            if interval[1] == 0:
                interval[1] = self.config.massbins * xfwhm
            interval = np.array([p.mass - interval[0], p.mass + interval[1]])

            boo1 = self.data.massdat[:, 0] < interval[1]
            boo2 = self.data.massdat[:, 0] > interval[0]
            boo3 = np.all([boo1, boo2], axis=0)
            top = self.data.massdat[boo3]

            stack = [top[:, 0]]
            for j, z in enumerate(self.data.ztab):
                stack.append(zarr[boo3, j])
            p.zstack = np.array(stack)

    def pks_mscore(self, xfwhm=2, power=2):
        self.get_zstack(xfwhm=xfwhm)
        for p in self.pks.peaks:
            m = p.zstack[0]
            ints = p.zstack[1:]  # mzgrid
            msum = np.sum(ints, axis=0)
            zsum = np.sum(ints, axis=1)
            msum /= np.amax(msum)
            zsum /= np.amax(zsum)
            p.mdist = np.transpose([m, msum])
            # mmax = np.argmax(msum)
            # sumz = np.sum(ints, axis=1)
            # zsum = deepcopy(ints[:, mmax])
            p.zdist = np.transpose([self.data.ztab, zsum])

            rats = []
            weights = []
            for i, z in enumerate(self.data.ztab):
                Y = ints[i]  # mzgrid
                sY = np.sum(Y)
                X = msum * sY / np.sum(msum)  # summed decon
                sX = np.sum(X)
                sae = np.sum(np.abs(Y - X))
                M = 1 - ud.safedivide1(sae, sX)
                rats.append(M)
                weights.append(sY)
            rats = np.array(rats)
            weights = np.array(weights)
            avg = ud.weighted_avg(rats, weights ** power)
            # print("Peak Mass:", p.mass, "Peak Shape Score", avg, p.mscore)
            p.mscore = avg

    def pks_csscore(self, xfwhm=2):
        try:
            if len(self.pks.peaks[0].zstack) < 1:
                self.get_zstack(xfwhm=xfwhm)
        except Exception as e:
            print("Error in z score: Make sure you get peaks first", e)

        for p in self.pks.peaks:
            ints = p.zstack[1:]
            msum = np.sum(ints, axis=0)
            # mmax = np.argmax(msum)
            sumz = np.sum(ints, axis=1)
            # sumz = deepcopy(ints[:, mmax]) # Switch to peak center only
            sumz /= np.amax(sumz)
            # sumz = deepcopy(p.pca_zdist[:, 1])
            # sumz -= np.amin(sumz)
            zm = np.argmax(sumz)

            zs = np.sum(sumz)
            badarea = 0

            index = zm
            low = sumz[zm]
            while index < len(sumz) - 1:
                index += 1
                val = sumz[index]
                if val < low:
                    low = val
                else:
                    badarea += val - low

            index = zm
            low = sumz[zm]
            while index > 0:
                index -= 1
                val = sumz[index]
                if val < low:
                    low = val
                else:
                    badarea += val - low

            p.cs_score = 1 - ud.safedivide1(badarea, zs)
            # print(badarea, zs, p.cs_score)

    def get_mzstack(self, xfwhm=2):
        zarr = np.reshape(self.data.mzgrid[:, 2], (len(self.data.data2), len(self.data.ztab)))
        # zarr = zarr / np.amax(np.sum(zarr, axis=1)) * np.amax(self.data.data2[:, 1])
        for i, p in enumerate(self.pks.peaks):
            # fwhm = p.errorFWHM
            # if fwhm == 0:
            #    fwhm = self.config.massbins
            # width = fwhm / 2. * xfwhm
            # interval = np.array([p.mass - width, p.mass + width])

            interval = np.abs(np.array(p.intervalFWHM) - p.mass) * xfwhm
            if interval[0] == 0:
                interval[0] = self.config.massbins * xfwhm
            if interval[1] == 0:
                interval[1] = self.config.massbins * xfwhm
            interval = np.array([p.mass - interval[0], p.mass + interval[1]])

            stack = []
            bvals = []
            for j, z in enumerate(self.data.ztab):
                interval2 = (interval + z * self.config.adductmass) / z
                boo1 = self.data.data2[:, 0] <= interval2[1]
                boo2 = self.data.data2[:, 0] >= interval2[0]
                boo3 = np.all([boo1, boo2], axis=0)
                top = self.data.data2[boo3]
                stack.append(np.transpose([top[:, 0], top[:, 1], zarr[boo3, j]]))
                bvals.append(boo3)
            btot = np.any(bvals, axis=0)
            fit = self.data.fitdat[btot]
            top = self.data.data2[btot]
            sse = np.sum((fit - top[:, 1]) ** 2)
            denom = np.sum((top[:, 1] - np.mean(top[:, 1])) ** 2)
            p.rsquared = 1 - ud.safedivide1(sse, denom)
            p.mzstack = np.array(stack, dtype='object')

    def pks_uscore(self, xfwhm=2, power=1):
        self.get_mzstack(xfwhm=xfwhm)
        for p in self.pks.peaks:
            rats = []
            weights = []
            for i, zval in enumerate(self.data.ztab):
                v = p.mzstack[i]
                X = v[:, 1]  # spectrum
                Y = v[:, 2]  # mzgrid
                if self.config.orbimode == 1:
                    Y = Y * zval
                sx = np.sum(X)
                sy = np.sum(Y)
                # r = 1 - ud.safedivide1(np.abs(sy - sz), sy)
                sae = np.sum(np.abs(X - Y))  # ** 2)
                U = 1 - ud.safedivide1(sae, sx)
                rats.append(U)
                weights.append(sy)
                # try:
                #    print(sse, r, sz, sy, z[0], z[-1])
                # except:
                # //    print(z)
            rats = np.array(rats)

            weights = np.array(weights)
            avg = ud.weighted_avg(rats, weights ** power)
            # print("Peak Mass:", p.mass, "Uniqueness Score", avg)
            p.uscore = avg

    def pks_fscore(self):
        for p in self.pks.peaks:
            # Tests if the FWHM interval is highly asymetric
            index = ud.nearest(self.data.massdat[:, 0], p.mass)
            height = self.data.massdat[index, 1]

            interval = np.array(p.intervalFWHM)
            diff = np.abs(np.array(p.intervalFWHM) - p.mass)

            p.fscore = 1
            if p.badFWHM:
                if diff[0] > diff[1]:
                    dcut1 = ud.datachop(self.data.massdat, interval[0] - self.config.massbins, p.mass)
                    lmin = np.amin(dcut1[:, 1])
                    fscore1 = score_minimum(height, lmin)
                    # print("Fscore2", fscore1, lmin, height)
                    p.fscore *= fscore1
                else:
                    dcut2 = ud.datachop(self.data.massdat, p.mass, interval[1] + self.config.massbins)
                    umin = np.amin(dcut2[:, 1])
                    fscore2 = score_minimum(height, umin)
                    # print("Fscore3", fscore2, umin, height)
                    p.fscore *= fscore2

                # diff = np.abs(interval - p.mass)
                # r = ud.safedivide1(diff, fwhm)
                # if np.any(r > thresh):
                #    p.fscore = 1 - (np.amax(r) - thresh) / (1 - thresh)

            # Test if the FWHM dip isn't low enough
            masses = self.pks.masses
            lower = np.all([interval[0] < masses, masses < p.mass], axis=0)
            upper = np.all([p.mass < masses, masses < interval[1]], axis=0)

            if np.any(lower):
                badmasses = masses[lower]
                for m in badmasses:
                    dcut1 = ud.datachop(self.data.massdat, m, p.mass)
                    lmin = np.amin(dcut1[:, 1])
                    fscore1 = score_minimum(height, lmin)
                    # print("Fscore4", fscore1, lmin, height)
                    p.fscore *= fscore1

            if np.any(upper):
                badmasses = masses[upper]
                for m in badmasses:
                    dcut2 = ud.datachop(self.data.massdat, p.mass, m)
                    umin = np.amin(dcut2[:, 1])
                    fscore2 = score_minimum(height, umin)
                    # print("Fscore5", fscore2, umin, height)
                    p.fscore *= fscore2

    '''
    def tscore(self):
        try:
            tscore = 1 - np.sum(np.abs(self.data.fitdat - self.data.data2[:, 1])) / np.sum(self.data.data2[:, 1])
        except Exception as e:
            print("Error in Tscore: ", e)
            tscore = 0
        self.data.tscore = tscore
        return tscore'''

    def dscore(self, xfwhm=2, power=2):
        """
        Combined score should include:

        1) Fit to data -> Uniqueness Score
        2) Background -> Uniqueness Score
        3) Unique evidence for each peak -> Uniqueness Score
        4) Noise -> M Score
        5) Consistent peak shape -> M Score
        6) Charge state distribution -> Z_Score
        7) FWHM -> Penalties for poor FWHM

        :return:
        """
        if len(self.pks.peaks) > 0:
            self.pks_mscore(xfwhm=xfwhm, power=power)
            self.pks_uscore(power=power, xfwhm=xfwhm)
            self.pks_csscore(xfwhm=xfwhm)
            self.pks_fscore()
            # self.tscore()
            dscores = []
            ints = []
            for p in self.pks.peaks:
                scores = np.array([p.mscore, p.uscore, p.fscore, p.cs_score])  # p.rsquared,
                p.dscore = np.prod(scores)
                p.lscore = np.prod(scores[:-1])
                dscores.append(p.dscore)
                ints.append(p.height)
                '''
                print("Mass:", p.mass,
                      "Peak Shape:", round(p.mscore * 100, 2),
                      "Uniqueness:", round(p.uscore * 100, 2),
                      # "Fitting R^2", round(p.rsquared * 100, 2),
                      "Charge:", round(p.cs_score * 100, 2),
                      "FWHM:", round(p.fscore * 100, 2),
                      "Combined:", round(p.dscore * 100, 2))
                # print(p.intervalFWHM)'''
            ints = np.array(ints)
            self.pks.uniscore = ud.weighted_avg(dscores, ints ** power) * self.config.error
            print("R Squared:", self.config.error)
            # print("TScore:", self.data.tscore)
            print("Average Peaks Score (UniScore):", self.pks.uniscore)
        else:
            print("No peaks present")

    def filter_peaks(self, minscore=0.4):
        self.pick_peaks()
        self.dscore()

        newpeaks = []
        for p in self.pks.peaks:
            if p.dscore > minscore:
                newpeaks.append([p.mass, p.height])
        newpeaks = np.array(newpeaks)
        if len(newpeaks) > 0:
            self.setup_peaks(newpeaks)
            self.dscore()
        else:
            print("No Peaks Found with Scores Above", minscore)
            self.pks = peakstructure.Peaks()

    def open_test_spectrum(self, masslist=None, n=1, **kwargs):
        data, ztab = MSBuild.simple_spectrum(masslist, **kwargs)
        testdir = os.path.join(self.config.UniDecDir, "../TestSpectra")
        if not os.path.isdir(testdir):
            os.mkdir(testdir)

        fname = "test_" + str(n) + ".txt"
        np.savetxt(os.path.join(testdir, fname), data)

        self.open_file(fname, testdir, **kwargs)
        self.config.maxmz, self.config.minmz = "", ""
        self.data.ztab = ztab
        pass

    def pass_data_in(self, data, n=1, dirname=None, fname=None, **kwargs):
        # Specify the directory for the new file
        # If one is not specified, then use the TestSpectra location
        if dirname is None:
            # Find parent of parent of UniDecDir
            parent = os.path.split(self.config.UniDecDir)[0]
            parent = os.path.split(parent)[0]
            testdir = os.path.join(parent, "TestSpectra")
        else:
            testdir = dirname
        # Create the directory if it doesn't exist
        if not os.path.isdir(testdir):
            os.mkdir(testdir)

        # If the fname is not passed in, create it.
        if fname is None:
            fname = "test_" + str(n) + ".txt"

        # Create the file in the new path
        path = os.path.join(testdir, fname)
        if "silent" in kwargs and not kwargs["silent"]:
            print("Creating file:", path)
        np.savetxt(path, data)

        # Open it
        self.open_file(fname, testdir, **kwargs)
        self.config.maxmz, self.config.minmz = "", ""
        # self.config.default_isotopic_res()

    def get_spectrum_peaks(self, threshold=0.05, window=None):
        if window is None:
            window = self.config.mzsig

        peaks = ud.peakdetect(self.data.data2, None, window, threshold)
        self.data.data2 = np.array(peaks)
        ud.dataexport(self.data.data2, self.config.infname)

        print(self.config.dirname)
        return peaks

    def estimate_areas(self):
        gauss_coeff = np.sqrt(np.pi / np.log(2)) / 2
        adjusted_coeff = ((0.5 * gauss_coeff) + (np.pi / 4))
        for p in self.pks.peaks:
            if self.config.psfun == 0:  # Gaussian
                p.estimatedarea = p.height * p.errorFWHM * gauss_coeff
            elif self.config.psfun == 1:  # Lorentzian
                p.estimatedarea = p.height * p.errorFWHM * np.pi / 2
            elif self.config.psfun == 2:  # Split G/L
                p.estimatedarea = p.height * p.errorFWHM * adjusted_coeff

    def make_plot(self, massrange=None, autorange=True, show=True):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax1 = plt.subplot(121)
        try:
            plt.plot(self.data.data2[:, 0], self.data.data2[:, 1], color="k")
        except (IndexError, AttributeError):
            plt.plot(self.data.rawdata[:, 0], self.data.rawdata[:, 1], color="k")
        ax2 = plt.subplot(122)
        try:
            plt.plot(self.data.massdat[:, 0], self.data.massdat[:, 1], color="k")
        except (IndexError, AttributeError):
            pass

        if massrange is not None:
            plt.xlim(massrange)
        elif autorange:
            plt.xlim(np.amin(self.pks.masses) - 2000, np.amax(self.pks.masses) + 2000)
        if show:
            plt.show()
        return fig, ax1, ax2

    def label_plot_correct(self, ax, matches):
        for i, match in enumerate(matches):
            p = self.pks.peaks[i]
            m = p.mass
            h = p.height

            if "Correct" in match:
                color = "g"
                ax.text(m, h + 5, match)
            elif match == "":
                color = "y"
            else:
                color = "r"
                ax.text(m, h + 5, match)

            ax.plot(m, h, color=color, marker="o")

    def makeplot1im(self, plot1im=None, plot1fit=None, imfit=False):
        if plot1im is None:
            plot1im = plot2d.Plot2dBase()
        if plot1fit is None:
            plot1fit = plot2d.Plot2dBase()

        if self.config.imflag == 1:
            try:
                plot1im.contourplot(self.data.data3, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)",
                                    title="IM-MS Data")
            except Exception:
                pass
            if imfit:
                try:
                    plot1fit.contourplot(self.data.fitdat2d, self.config, xlab="m/z (Th)", ylab="Arrival Time (ms)",
                                         title="IM-MS Fit")
                except Exception:
                    pass

    def makeplot1(self, plot=None, intthresh=False, imfit=True, config=None):
        """
        Plot data and fit in self.view.plot1 and optionally in plot1fit
        :param plot: plot to use
        :param intthresh: if True, plot intensity threshold
        :param imfit: if True, plot IM fit
        :param config: config to use, default is self.config
        :return: plot object
        """
        if config is None:
            config = self.config
        if plot is None:
            plot = plot1d.Plot1dBase()

        if self.config.batchflag == 0:
            tstart = time.perf_counter()
            leg = False
            # Honestly, I don't remember what this was supposed to do. Probably something with peak mode.
            try:
                test = float(self.config.reductionpercent)
            except Exception as e:
                self.config.reductionpercent = 0

            if self.config.reductionpercent < 0:
                print("Making Dot Plot")
                data2 = ud.dataprep(self.data.rawdata, self.config, peaks=False, intthresh=False)
                plot.plotrefreshtop(data2[:, 0], data2[:, 1], "Data and UniDec Fit",
                                    "m/z (Th)", "Normalized Intensity", "Data", self.config, nopaint=True)
                plot.plotadddot(self.data.data2[:, 0], self.data.data2[:, 1], 'blue', "o", "Peaks")

                # Add red line if there is a fitdata
                try:
                    if len(self.data.fitdat) > 0 and imfit:
                        plot.plotadddot(self.data.data2[:, 0], self.data.fitdat, 'red', "s", "Fit Data")
                        leg = True
                    pass
                except Exception:
                    pass

            else:
                # Add Processed Data Plot
                plot.plotrefreshtop(self.data.data2[:, 0], self.data.data2[:, 1],
                                    "Data and UniDec Fit",
                                    "m/z (Th)", "Normalized Intensity", "Data", self.config,
                                    nopaint=True)

                # Add red line if there is a threshold
                if self.config.intthresh != 0 and self.config.imflag == 0 and intthresh:
                    plot.plotadd(self.data.data2[:, 0], np.zeros_like(self.data.data2[:, 1]) + self.config.intthresh,
                                 "red", "Noise Threshold")
                    leg = True

                # Add red line for fit
                try:
                    if len(self.data.fitdat) > 0 and imfit:
                        plot.plotadd(self.data.data2[:, 0], self.data.fitdat, 'red', "Fit Data")
                        leg = True
                    pass
                except Exception as e:
                    pass

            # Add baseline if the flag is set
            if self.config.aggressiveflag != 0 and len(self.data.baseline) == len(self.data.data2):
                plot.plotadd(self.data.data2[:, 0], self.data.baseline, 'blue', "Baseline")

            # Add Legend
            if leg:
                plot.add_legend()
            # Repaint
            plot.repaint()
            print("Plot 1: %.2gs" % (time.perf_counter() - tstart))
        return plot



# Optional Run
if __name__ == "__main__":
    print("Running unidec Command Line")
    eng = UniDec()

    filename = "BSA.txt"
    path = "bin/Example Data"
    # files=["150611_ND_AmtB_05_100.txt"]
    # files = ["0.txt"]  # ,"250313_AQPZ_POPC_100_imraw.txt"]
    # files=["Traptavidin_20mer+B4F_new_40Colli.txt"]
    # eng.read_hdf5()
    # eng.config.mzbins = 1
    # eng.config.psfun=2
    # eng.config.peakwindow=2000.
    # eng.process_data()
    # eng.run_unidec(silent=False)

    test = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw"
    test1 = "C:\\Users\\MartyLabsOfficePC\\OneDrive - University of Arizona\\Desktop\\20230816_Myoglobin 0666 ugmL 01.wiff"
    dat = eng.open_file(test1)
    for i in eng.data.rawdata:
        print(i)
    exit()
    eng.config.imflag = 1
    eng.config.mzbins = 1
    eng.open_file(test1)
    # eng.match()
    eng.run_unidec()
    # plot = eng.makeplot2()
    # os.chdir(path)
    # plot.save_figure("test.png")
    # eng.gen_html_report()
    # eng.pick_peaks()

    # eng.pick_peaks()
    # eng.write_hdf5()

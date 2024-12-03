import numpy as np
import os
import scipy
import unidec.tools as ud
from unidec.modules.unidec_enginebase import UniDecEngine
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from copy import deepcopy
import matplotlib.pyplot as plt
import scipy.fft as fft
from unidec.modules import unidecstructure, peakstructure, IM_functions, fitting, i2ms_importer
import time
from unidec import engine
from scipy.optimize import curve_fit

xp = np


def Gmax2(darray):
    ints = darray[:, 1]
    boo1 = ints < np.median(ints) * 2

    noisedata = darray[boo1]

    hist, edges = np.histogram(noisedata[:, 1], bins=100)
    edges = edges[1:] - (edges[1] - edges[0]) / 2.

    smoothhist = scipy.ndimage.gaussian_filter1d(hist, 1)

    grad = np.gradient(smoothhist)
    grad2 = np.gradient(grad)

    # gradmin = edges[np.argmin(grad)]

    maxindex = np.argmax(hist)
    grad2maxindex = np.argmax(grad2[maxindex:])
    gradmax2 = edges[maxindex + grad2maxindex]
    # fit = np.polyfit(noisedata[:, 0], noisedata[:, 1], 1)
    # print(fit)
    print("Gradmax2:", gradmax2)
    return gradmax2


def fft_fun(a):
    return fft.rfft2(a)


def ifft_fun(A, shape):
    return fft.irfft2(A, shape)


def safedivide(a, b):
    return ud.safedivide(a, b)


def ndis(x, y, s):
    return np.exp(-(x - y) * (x - y) / (2.0 * s * s))


def cconv2D_preB(a, B):
    A = fft_fun(a)
    C = A * B
    c = ifft_fun(C, a.shape)
    return xp.abs(c)


def softmax(I, beta):
    numz = len(I)
    E = xp.exp(beta * I)
    Min2 = xp.amin(E, axis=0)
    Sum1 = xp.sum(I, axis=0)
    Sum2 = xp.sum(E, axis=0)
    factor = safedivide(Sum1, Sum2 - Min2 * numz)
    I = (E - Min2) * factor
    return I


def filter_centroid_const(d, mz_threshold=1):
    """
    """
    try:
        boo1 = np.diff(d[:, 0]) > mz_threshold
        boo2 = np.append(boo1[0], boo1)
        boo3 = np.append(boo1, boo1[-1])
        boo1 = np.logical_and(boo2, boo3)
    except:
        print("No Centroid Overlaps Found")
        return d
    return d[boo1]


def slopefunc(x, slope):
    return slope * x


def linearmodel(x, slope, intercept):
    return slope * x + intercept


def quadraticmodel(x, a, b, c):
    return a * x ** 2 + b * x + c


class UniDecCD(engine.UniDec):
    def __init__(self):
        """
        UniDecCD Engine

        Consists of three main subclasses: Config, DataContiner, Peaks

        :return: None
        """
        UniDecEngine.__init__(self)
        # Creates the initial array objects
        self.darray = None
        self.farray = None
        self.zarray = None
        self.thermodata = False
        self.harray = []
        # Sets a few default config parameters
        self.config.mzbins = 1
        self.config.rawflag = 1
        self.config.poolflag = 1
        self.exemode = True
        self.massaxis = None
        self.invinjtime = None
        pass

    def exe_mode(self, exemode=True):
        self.exemode = exemode

    def open_file(self, path, refresh=False):
        """
        Passthrough function for opening CD-MS files. Calls self.open_cdms_file.
        :param path: File path
        :param refresh: Whether to refresh file
        :return: None
        """
        self.open_cdms_file(path, refresh=refresh)

    def open_cdms_file(self, path, refresh=False):
        """
        Main function for opening CD-MS files. Supports Thermo Raw and DMT, mzML, text, binary, and numpy compressed.

        Calls self.before_open() first to setup UniDec files and directories. Searches first for previously opened
        numpy compressed data. Switches path to that if found. Otherwise it opens the specified path.

        Thermo Data via .Raw or .mzML needs to be corrected first by multiplying by the injection times in units of s.

        Text, binary, and numpy files are opened assuming the first column is m/z, the second is intensity, the third is scans.
        If no third column is provided, it assigns a scan number of 1 to each.
        The m/z, intensity, and scan data is saved in self.darray (the raw data array).

        Next, it filters out any data with intensities of zero.
        The noise is calculated as the maximum of the second derivative of the intensity distribution, also known as the Gmax2.
        The self.farray (filtered array) is created from self.darray.

        An Import Error is raised if the darray is empty

        The min and max m/z values are set from the darray.

        Then, saves a numpy compressed file if it doesn't exist for faster imports in the future.
        Finally, it opens the config file if found in the _unidecfiles folder.

        :param path: File path to open.
        :return: None
        """
        # Setup
        starttime = time.perf_counter()
        self.path = path
        # Set up unidec paths
        self.before_open(refresh=refresh)

        if os.path.isdir(self.path):
            self.path = self.convert_stori(self.path)

        # Get the extension
        extension = os.path.splitext(self.path)[1]
        self.invinjtime = None

        if extension.lower() == ".raw":
            # Import Thermo Raw file using ThermoDataImporter
            self.TDI = ImporterFactory.create_importer(self.path)
            # Get the data
            data = self.TDI.grab_data(threshold=self.config.CDprethresh)
            # Get the scans
            self.scans = self.TDI.scans
            # Get the injection times for each scan
            self.it = self.TDI.get_inj_time_array()
            # Get the overall resolution
            self.res = self.TDI.msrun.resolution
            # Set flag for correcting injection times later
            self.thermodata = True
            self.invinjtime = 1. / self.it

        elif extension.lower() == ".i2ms" or extension.lower() == ".dmt":
            # Import Thermo Raw file using ThermoDataImporter
            self.I2MSI = i2ms_importer.I2MSImporter(self.path)
            # Get the data
            data = self.I2MSI.grab_data(threshold=self.config.CDprethresh)
            mz = data[:, 0]
            intensity = data[:, 1]
            # Get the scans
            scans = self.I2MSI.scans
            # Set flag for correcting injection times later
            self.thermodata = False
            # Set the inverse injection time array for DMT data
            self.invinjtime = self.I2MSI.invinjtime

        elif extension.lower() == ".mzml" or extension.lower() == ".gz":
            # Import mzML data, scans, and injection time
            # Note, resolution info is not transferred to mzML to my knowledge
            self.MLI = ImporterFactory.create_importer(self.path)
            data = self.MLI.grab_data(threshold=self.config.CDprethresh)
            self.scans = self.MLI.scans
            self.it = self.MLI.get_inj_time_array()
            self.invinjtime = 1. / self.it
            self.thermodata = True

        elif extension.lower() == ".txt":
            # Simple Text File
            try:
                # Load text file
                data = np.loadtxt(self.path)

                # Filter pre-threshold
                b1 = data[:, 1] > self.config.CDprethresh
                data = data[b1]

                # Assume m/z is in column 1 and intensity column 2
                mz = data[:, 0]
                intensity = data[:, 1]

                # Look for scans in column 3
                try:
                    scans = data[:, 2]
                except Exception:
                    # If not found, assume all are scan 1
                    scans = np.ones_like(mz)
                try:
                    self.invinjtime = data[:, 3]
                except Exception:
                    self.invinjtime = None
            except:
                # More complex text file from Jarrold Lab
                data = np.genfromtxt(path, dtype=np.str, comments=None)
                mzkey = "m/z_Lin_Reg"
                intkey = "Charge"
                mzcol = np.argmax(data[0] == mzkey)
                intcol = np.argmax(data[0] == intkey)
                mz = data[1:, mzcol].astype(float)
                intensity = data[1:, intcol].astype(float)
                scans = np.arange(len(mz))

                # Filter pre-threshold
                b1 = intensity > self.config.CDprethresh
                mz = mz[b1]
                intensity = intensity[b1]
                scans = scans[b1]

                self.invinjtime = None
            # Don't do post-processing for thermo data
            self.thermodata = False

        elif extension.lower() == ".csv":
            # Load text file
            data = np.genfromtxt(self.path, delimiter=",")
            # Filter pre-threshold
            b1 = data[:, 1] > self.config.CDprethresh
            data = data[b1]

            # Assume m/z is in column 1 and intensity column 2
            mz = data[:, 0]
            intensity = data[:, 1]
            # Look for scans in column 3
            try:
                scans = data[:, 2]
            except Exception:
                # If not found, assume all are scan 1
                scans = np.ones_like(mz)
            try:
                self.invinjtime = data[:, 3]
            except Exception:
                self.invinjtime = None
            # Don't do post-processing for thermo data
            self.thermodata = False

        elif extension.lower() == ".bin":
            # Load numpy compressed file
            data = np.fromfile(self.path)
            try:
                data = data.reshape((int(len(data) / 3), 3))
            except Exception:
                data = data.reshape((int(len(data) / 2)), 2)

            # Filter pre-threshold
            b1 = data[:, 1] > self.config.CDprethresh
            data = data[b1]
            # Assume m/z is in column 1 and intensity column 2
            mz = data[:, 0]
            intensity = data[:, 1]
            # Look for scans in column 3
            try:
                scans = data[:, 2]
            except:
                # If not found, assume all are scan 1
                scans = np.ones_like(mz)
            try:
                self.invinjtime = data[:, 3]
            except Exception:
                self.invinjtime = None
            # Don't do post-processing for thermo data
            self.thermodata = False

        elif extension.lower() == ".npz":
            # Load numpy compressed file
            data = np.load(self.path, allow_pickle=True)['data']

            # Filter pre-threshold
            b1 = data[:, 1] > self.config.CDprethresh
            data = data[b1]

            # Assume m/z is in column 1 and intensity column 2
            mz = data[:, 0]
            intensity = data[:, 1]
            # Look for scans in column 3
            try:
                scans = data[:, 2]
            except Exception:
                # If not found, assume all are scan 1
                scans = np.ones_like(mz)
            try:
                self.invinjtime = data[:, 3]
            except Exception:
                self.invinjtime = None
            # Don't do post-processing for thermo data
            self.thermodata = False

        else:
            print("Unrecognized file type:", self.path)
            return 0

        print("File Read. Length: ", len(data))
        # Post processing if data is from raw or mzML
        # Ignored if text file input
        if self.thermodata:
            scans = np.concatenate([s * np.ones(len(data[i])) for i, s in enumerate(self.scans)])

            mz = np.concatenate([d[:, 0] for d in data])
            try:
                intensity = np.concatenate([d[:, 1] * self.it[i] / 1000. for i, d in enumerate(data)])
            except Exception as e:
                print(e, self.it)
                intensity = np.concatenate([d[:, 1] for i, d in enumerate(data)])

            try:
                it = np.concatenate([it * np.ones(len(data[i])) for i, it in enumerate(self.invinjtime)])
                self.invinjtime = it
            except Exception as e:
                print("Error with injection time correction:", e)
                print(mz.shape, intensity.shape, scans.shape, self.invinjtime.shape)
                self.invinjtime = None

        if self.invinjtime is None:
            self.invinjtime = np.ones_like(scans)

        # Create data array
        self.darray = np.transpose([mz, intensity, scans, self.invinjtime])
        # Filter out only the data with positive intensities
        boo1 = self.darray[:, 1] > 0
        self.darray = self.darray[boo1]
        print("Filtered 0 intensity values. Length: ", len(self.darray))
        # Define the noise level as the median of the data intensities
        self.noise = Gmax2(self.darray)
        # Create the filtered array (farray) object
        self.farray = deepcopy(self.darray)
        self.data.rawdata = self.darray[:, :2]

        # Check for empty data
        if ud.isempty(self.darray):
            print("Error: Data Array is Empty")
            print("Likely an error with data conversion")
            raise ImportError

        # Get the min and max m/z value
        self.config.minmz = np.amin(mz)
        self.config.maxmz = np.amax(mz)

        # Create the npz file of the extracted values if it doens't exist
        if not os.path.isfile(self.config.cdrawextracts) or refresh:
            try:
                np.savez_compressed(self.config.cdrawextracts, data=self.darray)
            except Exception as e:
                pass

        self.config.cdmsflag = 1
        # Load the config if you can find it
        if os.path.isfile(self.config.confname):
            self.load_config(self.config.confname)
        else:
            self.export_config()

        print("Read File Length: ", len(self.farray), "Noise: ", self.noise, "Time:", time.perf_counter() - starttime)

    def before_open(self, refresh=False):
        """
        Creates the initial blank self.pks and self.data objects. Creates the _unidecfiles folder if needed.
        Sets the default file names. If it finds raw data already in the numpy compressed format, it switches the path to that.
        :return: None
        """
        # Create blank peaks and data objects
        self.pks = peakstructure.Peaks()
        self.data = unidecstructure.DataContainer()

        # Get the directory and file names
        file_directory = os.path.dirname(self.path)
        file_name = os.path.basename(self.path)
        # Handle Paths
        self.config.filename = self.path
        self.config.dirname = file_directory

        # Change paths to unidecfiles folder
        dirnew = os.path.splitext(self.path)[0] + "_unidecfiles"
        if not os.path.isdir(dirnew):
            os.mkdir(dirnew)
        self.config.udir = dirnew
        # Default file names
        basename = os.path.split(os.path.splitext(file_name)[0])[1]
        self.config.outfname = os.path.join(self.config.udir, basename)
        self.config.extension = os.path.splitext(self.config.filename)[1]
        self.config.default_file_names()

        # Look for already processed data in the form of an npz file and load it if it exists for speed.
        if os.path.isfile(self.path) and os.path.isfile(self.config.cdrawextracts) and not refresh:
            print("Raw data found:", self.config.cdrawextracts)
            self.path = self.config.cdrawextracts

    def open_stori(self, path):
        starttime = time.perf_counter()
        files = os.listdir(path)
        alldata = []
        for i, f in enumerate(files):
            if os.path.splitext(f)[1] == ".csv":
                p = os.path.join(path, f)
                data = np.genfromtxt(p, delimiter="\t", skip_header=1)
                alldata.append(data)

                # if ud.isempty(alldata):
                #    alldata = data
                # else:
                #    alldata = np.append(alldata, data, axis=0)
        alldata = np.concatenate(alldata, axis=0)
        print(alldata.shape)
        print("Opening Time: ", time.perf_counter() - starttime)
        # Scan	Ionnumber	Segnumber	M/Z	Frequency	Slope	Slope R Squared	Time of Birth	Time of Death
        mz = alldata[:, 3]
        intensity = alldata[:, 5]
        scan = alldata[:, 0]
        return mz, intensity, scan

    def convert_stori(self, path):
        mz, intensity, scans = self.open_stori(path)
        darray = np.transpose([mz, intensity, scans])
        # Filter out only the data with positive intensities
        boo1 = darray[:, 1] > 0
        darray = darray[boo1]
        outfile = path + "_combined.npz"
        np.savez_compressed(outfile, data=darray)
        return outfile

    def process_data(self, transform=True):
        """
        Main function for processing CDMS data, includes filtering, converting intensity to charge, histogramming,
        processing the histograms, and transforming (if specified). Transforming is sometimes unnecessary, which is why
        it can be turned off for speed.

        :param transform: Sets whether to transform the data from m/z to mass. Default True.
        :return: None
        """
        starttime = time.perf_counter()
        # Copy filtered array and delete peaks.
        self.farray = deepcopy(self.darray)
        self.pks = peakstructure.Peaks()

        # Filter Scans
        try:
            if int(self.config.CDScanStart) > 0:
                self.farray = self.farray[self.farray[:, 2] > self.config.CDScanStart]
        except:
            pass
        try:
            if int(self.config.CDScanEnd) > 0:
                self.farray = self.farray[self.farray[:, 2] < self.config.CDScanEnd]
        except:
            pass

        # Compress Scans
        if self.config.CDScanCompress > 1:
            self.farray[:, 2] = np.floor(self.farray[:, 2] / self.config.CDScanCompress)

        # Filter m/z
        print("Filtering m/z range:", self.config.minmz, self.config.maxmz, "Start Length:", len(self.farray))
        self.filter_mz(mzrange=np.abs([self.config.minmz, self.config.maxmz]))

        # Filter Centroids
        print("Filtering centroids:", self.config.CDres, "Start Length:", len(self.farray))
        self.filter_centroid_all(self.config.CDres)

        # Filter Charge States
        print("Filtering Charge range:", self.config.startz, self.config.endz, "Start Length:", len(self.farray))
        self.filter_z(zrange=np.abs([self.config.startz, self.config.endz + 1]))

        # Convert intensity to charge
        print("Converting From Intensity to Charge. Slope:", self.config.CDslope, "Start Length:", len(self.farray))
        self.convert(slope=self.config.CDslope)

        # Create Histogram
        print("Creating Histogram Bins:", self.config.mzbins, self.config.CDzbins, "Start Length:", len(self.farray))
        self.histogram(mzbins=self.config.mzbins, zbins=self.config.CDzbins)

        if len(self.harray) > 0:
            self.harray = self.hist_data_prep()
            self.harray_process(transform=transform)

        else:
            print("ERROR: Empty histogram array on process")
            return 0
        print("Process Time:", time.perf_counter() - starttime)

    def harray_process(self, transform=True):
        self.data.data2 = np.transpose([self.mz, np.sum(self.harray, axis=0)])
        np.savetxt(self.config.infname, self.data.data2)

        print("Transforming m/z to mass:", self.config.massbins, "Start Length:", len(self.farray))
        if transform:
            self.transform()
            np.savetxt(self.config.massdatfile, self.data.massdat)
            self.unprocessed = deepcopy(self.data.massdat)

    def filter_int(self, int_range):
        """
        Filter the self.farray to include only intensities within the int_range.

        :param int_range: list or array with int_range[0] as the minimum intensity and int_range[1] as the maximum
        :return: None
        """
        if int_range is not None:
            boo1 = self.farray[:, 1] > int_range[0]
            boo2 = self.farray[:, 1] < int_range[1]
            boo3 = np.logical_and(boo1, boo2)
            self.farray = self.farray[boo3]
        return self.farray

    def filter_z(self, zrange, slope=None):
        """
        Function to filter intensities based on their eventual charge state. Calculates int_range from zrange based on the slope.
        Calls self.filter_int from the calculated int_range.

        :param zrange: list or array with zrange[0] as the minimum charge state and zrange[1] as the maximum.
        :param slope: Optional parameter to specify the slope. If None (default), will use self.config.CDslope
        :return: None
        """
        if slope is not None:
            self.config.CDslope = slope
        else:
            slope = self.config.CDslope

        # If slope is negative (old dev settings) or subtype is 0, multiply the slope by the noise ratio.
        if slope < 0 or self.config.subtype == 0:
            slope = np.abs(slope) * self.noise
        # The int_range is the zrange * the slope
        int_range = np.array(zrange) * slope
        # Run filter_int with the int_range calculated from zrange
        self.filter_int(int_range=int_range)

    def filter_mz(self, mzrange=None):
        """
        Filter the self.farray to include only m/z value within the mzrange.

        :param mzrange: list or array with mzrange[0] as the minimum m/z and mzrange[1] as the maximum
        :return: self.farray, the filtered data array
        """
        if mzrange is not None:
            boo1 = self.farray[:, 0] > mzrange[0]
            boo2 = self.farray[:, 0] < mzrange[1]
            boo3 = np.logical_and(boo1, boo2)
            self.farray = self.farray[boo3]
        return self.farray

    def filter_centroid_all(self, mz_threshold=1):
        """
        Filter the self.farray to remove ions with nearby ions. See Worner et Al., Nat. Methods, 2020.

        Removed all the resolution settings of Worner et Al in favor of a simple single mz_threshold.

        :param mz_threshold: Distance cutoff. Two points in the same scan closer than this will be removed.
        :return: None
        """
        if mz_threshold > 0:
            # farray2 = []
            # for i, s in enumerate(self.scans):
            #    sbool = self.farray[:, 2] == s
            #    farray2.append(filter_centroid_scan(self.farray[sbool], mz_threshold=mz_threshold))
            # self.farray = np.concatenate(farray2)
            self.farray = filter_centroid_const(self.farray, mz_threshold=mz_threshold)
        return self.farray

    def simp_convert(self, y):
        if self.config.CDslope > 0 and self.config.subtype == 1:
            y = y / self.config.CDslope
        elif self.config.CDslope < 0 or self.config.subtype == 0:
            y = y / (np.abs(self.config.CDslope) * self.noise)
        return y

    def convert(self, slope=None):
        if slope is not None:
            self.config.CDslope = slope
        else:
            slope = self.config.CDslope

        if slope > 0 and self.config.subtype == 1:
            self.zarray = self.farray[:, 1] / slope
        elif slope < 0 or self.config.subtype == 0:
            self.zarray = self.farray[:, 1] / (np.abs(slope) * self.noise)

    def histogram(self, mzbins=1, zbins=1, x=None, y=None, mzrange=None, zrange=None):
        if x is None:
            x = self.farray[:, 0]

        if y is None:
            y = self.zarray

        if mzbins < 0.001:
            print("Error, mzbins too small. Changing to 1", mzbins)
            mzbins = 1
            self.config.mzbins = 1
        self.config.mzbins = mzbins
        self.config.CDzbins = zbins
        if len(x) == 0:
            print("ERROR: Empty Filtered Array, check settings")
            self.harray = []
            return 0

        if mzrange is None:
            mzrange = np.array([self.config.minmz, self.config.maxmz])
            if np.any(mzrange < 0):
                mzrange = np.abs(mzrange)
            else:
                mzrange = [np.floor(np.amin(x)), np.amax(x)]
        if zrange is None:
            zrange = np.array([self.config.startz, self.config.endz + 1])
            if np.any(zrange < 0):
                zrange = np.abs(zrange)
            else:
                zrange = np.array([np.floor(np.amin(y)), np.ceil(np.amax(y))])

        mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + mzbins / 2, mzbins)
        # Weird fix to make this axis even is necessary for CuPy fft for some reason...
        if len(mzaxis) % 2 == 1:
            mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + 3 * mzbins / 2, mzbins)
        zaxis = np.arange(zrange[0] - zbins / 2., zrange[1] + zbins / 2, zbins)

        if self.config.CDiitflag and self.invinjtime is not None:
            weights = self.farray[:, 3]
            print("Using weighted hist", np.mean(weights), np.amin(weights), np.amax(weights))
        else:
            weights = None

        self.harray, self.mz, self.ztab = np.histogram2d(x, y, [mzaxis, zaxis], weights=weights)

        self.mz = self.mz[1:] - mzbins / 2.
        self.ztab = self.ztab[1:] - zbins / 2.
        self.data.ztab = self.ztab
        self.harray = np.transpose(self.harray)

        if self.config.datanorm == 1:
            maxval = np.amax(self.harray)
            if maxval > 0:
                self.harray /= maxval

        self.X, self.Y = np.meshgrid(self.mz, self.ztab, indexing='xy')

        self.mass = (self.X - self.config.adductmass) * self.Y

        self.data.data3 = np.transpose([np.ravel(self.X, order="F"), np.ravel(self.Y, order="F"),
                                        np.ravel(self.harray, order="F")])

    def hist_data_prep(self, harray=None):
        if harray is None:
            harray = self.harray
        if self.config.smooth > 0 or self.config.smoothdt > 0:
            print("Histogram Smoothing:", self.config.smooth, self.config.smoothdt)
            harray = IM_functions.smooth_2d(harray, self.config.smoothdt, self.config.smooth)

        if self.config.intthresh > 0:
            print("Histogram Intensity Threshold:", self.config.intthresh)
            harray = self.hist_int_threshold(harray, self.config.intthresh)
        if self.config.reductionpercent > 0:
            print("Histogram Data Reduction:", self.config.reductionpercent)
            harray = self.hist_datareduction(harray, self.config.reductionpercent)

        if self.config.subbuff > 0 or self.config.subbufdt > 0:
            print("Histogram Background Subtraction:", self.config.subbuff, self.config.subbufdt)
            harray = IM_functions.subtract_complex_2d(harray.transpose(), self.config).transpose()

        harray = self.hist_filter_smash(harray)

        return harray

    def hist_int_threshold(self, harray, int_threshold):
        """
        Filter the harray to include only intensities above the int_threshold.

        :param harray: 2D numpy array to be filtered
        :param int_threshold: Intensity threshold
        :return: harray, the filtered 2D numpy array
        """
        boo1 = harray > int_threshold
        harray *= boo1
        return harray

    def hist_datareduction(self, harray, red_per):
        sdat = np.sort(np.ravel(harray))
        index = round(len(sdat) * red_per / 100.)
        cutoff = sdat[index]
        harray = self.hist_int_threshold(harray, cutoff)
        return harray

    def hist_mass_filter(self, harray, massrange=None):
        # Get values from config if not supplied
        if massrange is None:
            massrange = [self.config.masslb, self.config.massub]

        # Filter values
        boo1 = self.mass < massrange[0]
        boo2 = self.mass > massrange[1]
        boo3 = np.logical_or(boo1, boo2)
        # Set values outside range to 0
        harray[boo3] = 0
        return harray

    def hist_filter_smash(self, harray, smashrange=None):
        """
        Smashes the range region to 0. Used to eliminate unwanted peaks.
        :return: None
        """
        if smashrange is None:
            smashrange = self.config.smashrange
        else:
            self.config.smashrange = smashrange

        if not ud.isempty(smashrange):
            print("Smashing: ", smashrange)
            # Filter values
            boo1 = self.X > smashrange[0]
            boo2 = self.X < smashrange[1]
            boo3 = np.logical_and(boo1, boo2)
            print(boo3.shape)
            # Set values outside range to 0
            harray[boo3] = 0

        if not ud.isempty(self.config.smashlist) and self.config.smashflag == 1:
            print("Smashing: ", self.config.smashlist)
            for i in range(len(self.config.smashlist)):
                smashrange = self.config.smashlist[i]
                # Filter values
                boo1 = self.X > smashrange[0] - smashrange[1]
                boo2 = self.X < smashrange[0] + smashrange[1]
                boo3 = np.logical_and(boo1, boo2)

                boo1 = self.Y > smashrange[2] - smashrange[3]
                boo2 = self.Y < smashrange[2] + smashrange[3]
                boo4 = np.logical_and(boo1, boo2)

                boo3 = np.logical_and(boo3, boo4)

                # Set values outside range to 0
                harray[boo3] = 0
        return harray

    def hist_nativeZ_filter(self, nativeZrange=None):
        # Get values from config if not supplied
        if nativeZrange is None:
            nativeZrange = [self.config.nativezlb, self.config.nativezub]

        if nativeZrange != [-1000, 1000]:
            # Find predicted native charge state
            nativeZ = ud.predict_charge(self.mass)
            # Find distance from predicted charge state
            offset = self.Y - nativeZ

            # Filter values
            boo1 = offset < nativeZrange[0]
            boo2 = offset > nativeZrange[1]
            boo3 = np.logical_or(boo1, boo2)
            # Set values outside range to 0
            self.harray[boo3] = 0

    def create_mass_axis(self, harray=None, mass=None):
        if harray is None:
            harray = self.harray
        if mass is None:
            mass = self.mass

        # filter out zeros
        harray = np.array(harray)
        boo1 = harray > 0
        mass = mass[boo1]

        # Test for if array is empty and error if so
        if ud.isempty(mass):
            print("ERROR: Empty histogram array on transform")
            return 0
        # Find the max and min and create the new linear mass axis
        minval = np.amax([np.amin(mass) - self.config.massbins * 3, self.config.masslb])
        maxval = np.amin([np.amax(mass) + self.config.massbins * 3, self.config.massub])
        minval = round(minval / self.config.massbins) * self.config.massbins  # To prevent weird decimals
        massaxis = np.arange(minval, maxval, self.config.massbins)
        return massaxis

    def transform(self, harray=None, dataobj=None, ztab=None, mass=None, mz=None):
        if harray is None:
            harray = self.harray
        if dataobj is None:
            dataobj = self.data
            flag = True
        else:
            flag = False
        if ztab is None:
            ztab = self.ztab
        if mass is None:
            mass = self.mass
        if mz is None:
            mz = self.mz
        # Test for if array is empty and error if so
        if len(harray) == 0:
            print("ERROR: Empty histogram array on transform")
            return 0

        massaxis = self.create_mass_axis(harray, mass=mass)

        # Create the mass grid
        dataobj.massgrid = []
        for i in range(len(ztab)):
            d = harray[i]
            if self.config.poolflag == 1:
                newdata = np.atleast_2d(np.transpose([mass[i], d]))
                massdata = ud.linterpolate(newdata, massaxis)
            else:
                boo1 = d > 0
                newdata = np.transpose([mass[i][boo1], d[boo1]])
                if len(newdata) < 2:
                    massdata = np.transpose([massaxis, massaxis * 0])
                else:
                    massdata = ud.lintegrate(newdata, massaxis)

            dataobj.massgrid.append(massdata[:, 1])
        dataobj.massgrid = np.transpose(dataobj.massgrid)
        # Create the linearized mass data by integrating everything into the new linear axis
        dataobj.massdat = np.transpose([massaxis, np.sum(dataobj.massgrid, axis=1)])
        if flag:
            self.config.massdatnormtop = np.amax(dataobj.massdat[:, 1])
        # Ravel the massgrid to make the format match unidec
        dataobj.massgrid = np.ravel(dataobj.massgrid)
        # Create the data2 and mzgrid objects from the histogram array for compatiblity with unidec functions
        dataobj.data2 = np.transpose([mz, np.sum(harray, axis=0)])
        dataobj.zdat = np.transpose([ztab, np.sum(harray, axis=1)])
        if self.X.shape != harray.shape:
            X, Y = np.meshgrid(mz, ztab, indexing='xy')
        else:
            X = self.X
            Y = self.Y

        dataobj.mzgrid = np.transpose(
            [np.ravel(X.transpose()), np.ravel(Y.transpose()), np.ravel(harray.transpose())])
        if flag:
            self.massaxis = massaxis
        return dataobj

    def transform_mzmass(self):
        # Test for if array is empty and error if so
        if len(self.harray) == 0:
            print("ERROR: Empty histogram array on transform")
            return 0
        massaxis = self.data.massdat[:, 0]
        print("m/z Length:", len(self.mz))
        print("Mass Length:", len(massaxis))
        # Create the mass grid
        self.data.mzmassgrid = []
        for i in range(len(self.mz)):
            d = self.harray[:, i]
            boo1 = d >= 0
            newdata = np.transpose([self.mass[:, i][boo1], d[boo1]])
            # if :
            #    self.data.mzmassgrid.append(massaxis * 0)
            # else:
            if self.config.poolflag == 1 and len(newdata) >= 2:
                massdata = ud.linterpolate(newdata, massaxis)
            else:
                massdata = ud.lintegrate(newdata, massaxis)
            self.data.mzmassgrid.append(massdata[:, 1])
        self.data.mzmassgrid = np.transpose(self.data.mzmassgrid)
        print("Created m/z vs. Mass: ", self.data.mzmassgrid.shape)
        pass

    def make_kernel(self, mzsig=None, zsig=None):
        """
        Sets up the deconvolution kernels based on the provided mzsig and zsig values, which specify the peak FWHM in m/z and z.
        :param mzsig: FWHM for m/z in units of Th.
        :param zsig: FWHM for charge. This is also the uncertainty for the charge assignment.
        :return: None
        """
        if mzsig is not None:
            self.config.mzsig = mzsig
        else:
            mzsig = self.config.mzsig
        if zsig is not None:
            self.config.csig = zsig
        else:
            zsig = self.config.csig

        X, Y = np.meshgrid(self.mz, self.ztab, indexing='xy')
        self.kernel = np.zeros_like(X)
        self.kernel[0, 0] = 1
        if mzsig > 0:
            self.mkernel = np.zeros_like(X)
            self.mkernel[0] = fitting.psfit(X[0], mzsig, np.amin(self.mz), psfun=self.config.psfun) \
                              + fitting.psfit(X[0], mzsig, np.amax(self.mz) + self.config.mzbins,
                                              psfun=self.config.psfun)
            if zsig > 0:
                self.kernel = fitting.psfit(X, mzsig, np.amin(self.mz), psfun=self.config.psfun) \
                              + fitting.psfit(X, mzsig, np.amax(self.mz) + self.config.mzbins, psfun=self.config.psfun)
                self.kernel *= fitting.psfit(-Y, zsig, -np.amin(self.ztab), psfun=self.config.psfunz) \
                               + fitting.psfit(-Y, zsig, -(np.amax(self.ztab) + self.config.CDzbins),
                                               psfun=self.config.psfunz)
            else:
                self.kernel = deepcopy(self.mkernel)
        elif zsig > 0 >= mzsig:
            self.mkernel = None
            self.kernel[:, 0] = fitting.psfit(-Y[:, 0], zsig, -np.amin(self.ztab), psfun=self.config.psfunz) \
                                + fitting.psfit(-Y[:, 0], zsig, -(np.amax(self.ztab) + self.config.CDzbins),
                                                psfun=self.config.psfunz)
        else:
            self.mkernel = None

        if self.config.psig > 0:
            self.mkernel2 = np.zeros_like(X)
            self.mkernel2[0] = fitting.psfit(X[0], self.config.psig, np.amin(self.mz), psfun=self.config.psfun) \
                               + fitting.psfit(X[0], self.config.psig, np.amax(self.mz) + self.config.mzbins,
                                               psfun=self.config.psfun)
            self.mkernel2 /= np.sum(self.mkernel2)

        self.kernel /= np.sum(self.kernel)
        # Create the flipped kernel for the correlation. Roll so that the 0,0 point is still the highest
        self.ckernel = np.roll(np.flip(self.kernel), (1, 1), axis=(0, 1))

        if mzsig > 0:
            self.mkernel /= np.sum(self.mkernel)

    def setup_zsmooth(self):
        # Setup Grids for  Z+1, and Z-1
        X, uY = np.meshgrid(self.mz, self.ztab + 1, indexing='xy')
        X, lY = np.meshgrid(self.mz, self.ztab - 1, indexing='xy')

        # Calculate m/z values for Z+1 and Z-1
        uppermz = (self.mass + uY * self.config.adductmass) / uY
        lowermz = ud.safedivide((self.mass + lY * self.config.adductmass), lY)  # In case the starting charge state is 1

        # Calculate the indexes for where to find the Z+1 and Z-1 m/z values
        m1 = self.mz[0]
        m2 = self.mz[-1]
        lm = len(self.mz)

        # Calculate normal indexes
        indexes = np.arange(0, lm)
        normindexes = np.array([indexes for i in self.ztab]).astype(int)

        # Calculate upperindexes for Z+1
        upperindex = np.round(((uppermz - m1) / (m2 - m1)) * (lm - 1))
        upperindex[-1] = normindexes[-1]
        upperindex[upperindex < 0] = normindexes[upperindex < 0]
        upperindex[upperindex >= lm] = normindexes[upperindex >= lm]

        # Calculate lowerindexes for Z-1
        lowerindex = np.round(((lowermz - m1) / (m2 - m1)) * (lm - 1))
        lowerindex[0] = normindexes[0]
        lowerindex[lowerindex < 0] = normindexes[lowerindex < 0]
        lowerindex[lowerindex >= lm] = normindexes[lowerindex >= lm]
        # For both, use the normal index if it is out of range

        self.upperindex = np.array(upperindex, dtype=int)
        self.lowerindex = np.array(lowerindex, dtype=int)

    def filter_zdist(self, I, setup=True):
        if setup:
            self.setup_zsmooth()

        # Get the intensities of the Z+1 values
        upperints = xp.zeros_like(I)
        for i, row in enumerate(self.upperindex):
            if i + 1 != len(I):
                upperints[i] = I[i + 1, row]
            else:
                upperints[i] = I[i, row]

        # Get the intensities of the Z-1 values
        lowerints = xp.zeros_like(I)
        for i, row in enumerate(self.lowerindex):
            if i != 0:
                lowerints[i] = I[i - 1, row]
            else:
                lowerints[i] = I[i, row]

        floor = self.config.zzsig
        if floor > 0:
            I = xp.clip(xp.exp(
                xp.mean(xp.asarray([xp.log(upperints + floor), xp.log(lowerints + floor), xp.log(I + floor)]), axis=0))
                        - floor, 0, None)
        else:
            ratio = xp.abs(floor)
            I = (upperints * ratio + I + lowerints * ratio) / 3.
        return I

    def setup_msmooth(self):
        mzoffsets = self.config.molig / self.ztab
        indexoffsets = np.round(mzoffsets / self.config.mzbins)
        lmz = len(self.mz)
        indexes = np.arange(0, lmz)
        upperindexes = np.array([indexes + i for i in indexoffsets]).astype(int)
        lowerindexes = np.array([indexes - i for i in indexoffsets]).astype(int)
        normindexes = np.array([indexes for i in indexoffsets]).astype(int)

        upperindexes[upperindexes < 0] = normindexes[upperindexes < 0]
        lowerindexes[lowerindexes < 0] = normindexes[lowerindexes < 0]
        upperindexes[upperindexes >= lmz] = normindexes[upperindexes >= lmz]
        lowerindexes[lowerindexes >= lmz] = normindexes[lowerindexes >= lmz]

        self.mupperindexes = upperindexes
        self.mlowerindexes = lowerindexes

    def filter_mdist(self, I, setup=True):
        if setup:
            self.setup_msmooth()

        upperints = xp.array([d[self.mupperindexes[i]] for i, d in enumerate(I)])
        lowerints = xp.array([d[self.mlowerindexes[i]] for i, d in enumerate(I)])

        floor = self.config.msig
        if floor > 0:
            I = xp.clip(xp.exp(
                xp.mean(xp.array([xp.log(upperints + floor), xp.log(lowerints + floor), xp.log(I + floor)]),
                        axis=0)) - floor,
                        0, None)
        else:
            ratio = xp.abs(floor)
            I = (upperints * ratio + I + lowerints * ratio) / 3.

        return I

    def decon_core(self):
        """
        Core unidec deconvolution function. Needs kernels already in place.
        :return: self.harray, the intensity array.
        """
        # Create a working array of intensity values
        I = deepcopy(self.harray)
        D = deepcopy(self.harray)

        # Perform the FFTs for convolution and correlation kernels
        ftk = fft_fun(self.kernel)
        ftck = fft_fun(self.ckernel)
        if self.config.psig != 0:
            ftmk = fft_fun(self.mkernel2)

        # Set up the while loop
        i = 0
        diff = 1
        # Continue until i = max number of iterations or the relative difference is less than 0.0001
        while i < self.config.numit and diff > 0.0001:
            # If beta is not 0, add in Softmax with beta value
            if self.config.beta > 0:
                I = softmax(I, self.config.beta)

            # Point Smoothing
            if self.config.psig > 0:
                I = cconv2D_preB(I, ftmk)

            # Run the smooth charge state filter and smooth mass filter. Set up the dist if needed.
            if i == 0:
                setup = True
            else:
                setup = False

            if self.config.zzsig != 0:
                if self.config.CDzbins == 1:  # Zdist smoothing currently only defined for unit charge bins
                    I = self.filter_zdist(I, setup)
                else:
                    print("Error: Charge Bins Size must be 1 for Charge State Smoothing")
            if self.config.msig != 0:
                I = self.filter_mdist(I, setup)

            # Classic Richardson-Lucy Algorithm here. Most the magic happens in this one line...
            if self.config.mzsig != 0 or self.config.csig != 0:
                newI = I * cconv2D_preB(safedivide(D, cconv2D_preB(I, ftk)), ftck)
                # newI = I * safedivide(D, cconv2D_preB(I, ftk))

                if i > 10:
                    # Calculate the difference and increment the counter to halt the while loop if needed
                    diff = np.sum((I - newI) ** 2) / np.sum(I)
                I = newI
            i += 1
            # print(i)
        print("Deconvolution iterations: ", i)
        if self.config.datanorm == 1:
            # Normalize the final array
            I /= xp.amax(I)

        # Get the reconvolved data
        recon = cconv2D_preB(I, ftk)

        # Get the fit data in 1D for the DScore calc
        self.data.fitdat = np.sum(recon, axis=0)
        if self.config.datanorm == 1:
            self.data.fitdat /= np.amax(self.data.data2[:, 1]) / np.amax(self.data.fitdat)

        if self.config.mzsig > 0 and self.config.rawflag == 0:
            # Reconvolved/Profile: Reconvolves with the peak shape in the mass dimension only
            ftmk = fft_fun(self.mkernel)
            recon2 = cconv2D_preB(I, ftmk)
            self.harray = recon2
        else:
            # Raw/Centroid: Takes the deconvolved data straight
            self.harray = I

        return self.harray

    def run_deconvolution(self, process_data=True):
        """
        Function for running the full deconvolution sequence, including setup, pre-, and post-processing.

        :return: None
        """
        if process_data:
            # Process data but don't transform
            self.process_data(transform=False)
        # Filter histogram to remove masses that are not allowed
        self.harray = self.hist_mass_filter(self.harray)
        # Filter histogram to remove charge states that aren't allowed based on the native charge state filter
        self.hist_nativeZ_filter()

        print("Running Deconvolution", self.config.mzsig, self.config.csig)
        starttime = time.perf_counter()
        if self.exemode:
            # Run the deconvolution core by calling on C external
            self.decon_external_call()
        else:
            # Make kernels for convolutions based on peak shapes
            self.make_kernel(self.config.mzsig, self.config.csig)
            # Run deconvolution
            self.decon_core()
        print("Deconvolution Time:", time.perf_counter() - starttime)
        # Transform m/z to mass
        self.transform()
        np.savetxt(self.config.massdatfile, self.data.massdat)

    def decon_external_call(self):
        self.export_config()
        # Check for this
        if self.config.CDzbins != 1 and self.config.zzsig != 0:
            print("ERROR: Charge smoothing is only define for when charges are binned to unit charge")
            self.harray = [[]]
            return
        # Output input data
        self.harray = np.array(self.harray)
        X, Y = np.meshgrid(self.mz, self.ztab, indexing='ij')
        outarray = self.harray.transpose()
        startdims = np.shape(outarray)
        outdat = np.transpose([np.ravel(X), np.ravel(Y), np.ravel(outarray)])
        np.savetxt(self.config.infname, outdat)
        print("Saved Input File:", self.config.infname)

        # Make the call
        ud.unidec_call(self.config)

        # Load in deconvolved data
        self.harray = np.loadtxt(self.config.deconfile)
        print("Loaded Output File:", self.config.deconfile)
        self.harray = self.harray.reshape(startdims).transpose()

        # Load in fit data, needed for scoring
        self.data.fitdat = np.fromfile(self.config.fitdatfile, dtype=self.config.dtype)
        self.data.fitdat = self.data.fitdat.reshape(startdims).transpose()
        self.data.fitdat = np.sum(self.data.fitdat, axis=0)
        print("Loaded Output File:", self.config.deconfile)

    def extract_intensities(self, mass, minz, maxz, window=25, sdmult=2, noise_mult=0):
        ztab = np.arange(minz, maxz + 1)
        mztab = (mass + ztab * self.config.adductmass) / ztab
        extracts = []
        windows = [-window, window]
        for i, z in enumerate(ztab):
            self.farray = deepcopy(self.darray)
            self.filter_mz(mzrange=windows + mztab[i])
            # ext = self.filter_centroid_all(1)
            int_range = [self.noise * noise_mult, 100000000000000]
            ext = self.filter_int(int_range)
            if not ud.isempty(self.farray):
                median = np.median(self.farray[:, 1])
                stddev = np.std(self.farray[:, 1]) * sdmult
                # stddev = median * medrange
                int_range = np.array([-stddev, stddev]) + median
                ext = self.filter_int(int_range)

                ext = ext[:, 1]
                ext = np.transpose([np.ones_like(ext) * z, ext])
                if i == 0:
                    extracts = ext
                else:
                    try:
                        extracts = np.concatenate((extracts, ext), axis=0)
                    except:
                        print("Error in extraction. Likely no data extracted.", z, len(ext))
                pass
                print("Filtering mz: ", mztab[i], len(ext))
        snextracts = deepcopy(extracts)
        try:
            snextracts[:, 1] = snextracts[:, 1] / self.noise
        except:
            print("No Intensity Found", snextracts)
        return extracts, snextracts

    def get_fit(self, extracts):
        fits = curve_fit(slopefunc, extracts[:, 0], extracts[:, 1])[0]
        x = np.unique(extracts[:, 0])
        fitdat = x * fits[0]
        return fits, np.transpose([x, fitdat])

    def get_fit_linear(self, extracts):
        fits = curve_fit(linearmodel, extracts[:, 0], extracts[:, 1])[0]
        x = np.unique(extracts[:, 0])
        fitdat = linearmodel(x, *fits)
        return fits, np.transpose([x, fitdat])

    def get_fit_quadratic(self, extracts):
        fits = curve_fit(quadraticmodel, extracts[:, 0], extracts[:, 1])[0]
        x = np.unique(extracts[:, 0])
        fitdat = quadraticmodel(x, *fits)
        return fits, np.transpose([x, fitdat])

    def plot_add(self):
        plt.subplot(132)
        plt.plot(self.mz, np.sum(self.harray, axis=0) / np.amax(np.sum(self.harray, axis=0)))
        plt.subplot(133)
        plt.plot(self.ztab, np.sum(self.harray, axis=1) / np.amax(np.sum(self.harray, axis=1)))

    def plot_hist(self):
        plt.subplot(131)
        plt.contourf(self.mz, self.ztab, self.harray, 100)
        plt.subplot(132)
        plt.plot(self.mz, np.sum(self.harray, axis=0) / np.amax(np.sum(self.harray, axis=0)))
        plt.subplot(133)
        plt.plot(self.ztab, np.sum(self.harray, axis=1) / np.amax(np.sum(self.harray, axis=1)))
        plt.show()

    def plot_mzmass_hist(self):
        plt.subplot(131)
        plt.contourf(self.mz, self.data.massdat[:, 0], self.data.mzmassgrid, 100)
        plt.subplot(132)
        plt.plot(self.mz, np.sum(self.data.mzmassgrid, axis=0) / np.amax(np.sum(self.data.mzmassgrid, axis=0)))
        plt.subplot(133)
        plt.plot(self.data.massdat[:, 0],
                 np.sum(self.data.mzmassgrid, axis=1) / np.amax(np.sum(self.data.mzmassgrid, axis=1)))
        plt.show()

    def sim_dist(self):
        self.config.mzsig = 1000
        self.config.csig = 0.2
        self.config.rawflag = 1
        mzmean = self.mz[round(len(self.mz) / 2)]
        zmean = self.ztab[round(len(self.ztab) / 2)]
        mzsig = self.config.mzsig
        zsig = self.config.csig
        X, Y = np.meshgrid(self.mz, self.ztab, indexing='ij')
        self.harray = fitting.psfit(X, mzsig, mzmean, psfun=self.config.psfun)
        self.harray *= fitting.psfit(Y, zsig, zmean, psfun=self.config.psfunz)
        self.harray = np.transpose(self.harray)
        self.harray /= np.amax(self.harray)


if __name__ == '__main__':
    eng = UniDecCD()
    # path = "C:\Data\CDMS\AqpZ_STORI\AqpZ_STORI\\072621AquaZ_low-high_noDi_IST10_processed"
    # eng.open_file(path)

    # exit()
    path = "C:\\Data\\CDMS\\spike trimer CDMS data.csv"
    path = "/unidec/bin\\Example Data\\CDMS\\GroEL_CDMS_1.RAW"
    eng.open_file(path)
    eng.process_data()
    eng.run_deconvolution()
    # eng.sim_dist()
    # eng.plot_add()
    # maxtup = np.unravel_index(np.argmax(eng.harray, axis=None), eng.harray.shape)
    # print(maxtup)
    # eng.make_kernel(eng.config.mzsig, eng.config.csig)
    # eng.harray = np.roll(eng.ckernel, maxtup, axis=(0, 1))
    # eng.plot_hist()
    # exit()
    # eng.decon_core()
    # print(np.unravel_index(np.argmax(eng.harray, axis=None), eng.harray.shape))
    # eng.plot_mzmass_hist()
    # eng.plot_hist()
    exit()

    eng = UniDecCD()
    eng.mz = np.array([5000, 6000, 7000, 8000, 9000])
    eng.ztab = np.array([12, 13, 14])
    eng.make_kernel(10, 2)
    print(eng.kernel)

    exit()
    path = "C:\\Data\\CDMS\\20210309_MK_ADH_pos_CDMS_512ms_5min_50ms_pressure01.RAW"
    eng = UniDecCD()
    eng.open_file(path)

    extracts, snextracts = eng.extract_intensities(147720, 21, 27)

    fits, fitdat = eng.get_fit(extracts)
    print(fits)
    snfits, snfitdat = eng.get_fit(snextracts)
    print(snfits)

    plt.subplot(121)
    plt.hist2d(extracts[:, 0], extracts[:, 1], cmap="nipy_spectral")
    plt.plot(fitdat[:, 0], fitdat[:, 1], color="red")
    plt.subplot(122)
    plt.hist2d(snextracts[:, 0], snextracts[:, 1], cmap="nipy_spectral")
    plt.plot(snfitdat[:, 0], snfitdat[:, 1], color="red")
    plt.show()

    # eng.run_deconvolution()

    exit()
    eng = UniDecCD()
    eng.mz = np.array([5000, 6000, 7000, 8000, 9000])
    eng.ztab = np.array([12, 13, 14])

    eng.config.mzbins = 1000
    eng.setup_zsmooth()
    eng.harray = eng.mass
    # eng.filter_zdist()
    exit()

import os
import numpy as np
import unidec_modules.unidectools as ud
from unidec_modules.mzMLimporter import mzMLimporter
from metaunidec.mudeng import MetaUniDec
from unidec import UniDec
from copy import deepcopy

chrom_file_exts = [".mzML", ".raw", ".Raw", ".RAW", ".d"]


class ChromEngine(MetaUniDec):
    """
    Imports mzML data files.
    """

    def __init__(self):
        MetaUniDec.__init__(self)
        self.filename = None
        self.dirname = None
        self.path = None
        self.chromdat = None
        self.tic = None
        self.ticdat = None
        self.spectra = None
        self.massdat = None
        self.mzdata = None
        self.config.default_high_res()
        self.unidec_eng = UniDec()

    def open_chrom(self, path, load_hdf5=True):
        if os.path.splitext(path)[1] == ".hdf5":
            self.open_hdf5_file(path)
            return True
        else:
            self.dirname, self.filename = os.path.split(path)

            print("Opening: ", self.filename)
            self.path = path
            load_hdf5 = self.load_mzml(self.path, load_hdf5)
            return load_hdf5

    def open_hdf5_file(self, path):
        print("Opening HDF5 File: ", path)
        header = os.path.splitext(path)[0]
        for ext in chrom_file_exts:
            testpath = header + ext
            if os.path.isfile(testpath) or os.path.isdir(testpath):
                try:
                    self.open_chrom(testpath)
                    # TODO: Make it so that you could match HDF5 file and data with different name
                    return True, testpath
                except Exception as e:
                    print("Error opening: ", testpath)
                    print(e)
        print("Attempted to open: ", path)
        print("Unable to find chromatograph file matching the header: ", header)
        print("Please rename the files such that the raw data matches the HDF5 file name")
        return False, header

    def load_mzml(self, path, load_hdf5=True, *args, **kwargs):
        self.path = path
        name = os.path.splitext(path)[0]
        self.outpath = name + ".hdf5"
        self.setup_filenames(self.outpath)
        self.data.filename = self.outpath
        hdf5 = False
        self.clear()
        if os.path.isfile(self.outpath) and load_hdf5:
            print('Opening HDF5 File:', self.outpath)
            try:
                self.open(self.outpath)
                hdf5 = True
            except Exception as e:
                print("Error opening prior hdf5 file:", e)
        if not os.path.isfile(self.outpath):
            self.data.new_file(self.outpath)
            self.open(self.outpath)

        self.update_history()

        self.chromdat = ud.get_importer(path)
        self.tic = self.chromdat.get_tic()
        self.ticdat = np.array(self.tic)

        return hdf5

    def get_scans(self, scan_range=None):
        self.mzdata = self.chromdat.get_data(scan_range)
        return self.mzdata

    def reset_vars(self, export=True):
        self.data.update_var_array()
        if export:
            self.data.export_hdf5(delete=True)
        self.data.v1name = "timestart"
        self.data.v2name = "scanstart"
        self.data.import_vars(get_vnames=False)

    def get_chrom_peaks(self, window=None):
        # Cleanup TIC Data
        ticdat = deepcopy(self.ticdat)
        ticdat = ud.gsmooth(ticdat, 2)
        ticdat[:, 1] -= np.amin(ticdat[:, 1])
        #ticdat = ud.gaussian_backgroud_subtract(ticdat, 100)
        maxval = np.amax(ticdat[:, 1])
        ticdat[:, 1] /= maxval
        maxt = np.amax(ticdat[:, 0])
        mint = np.amin(ticdat[:, 0])

        # Set Window
        if window is None:
            window = self.config.chrom_peak_width

        # Set Threshold
        noise = ud.noise_level2(ticdat, percent=0.50)
        print("Noise Level:", noise, "Window:", window)

        # Detect Peaks
        peaks = ud.peakdetect_nonlinear(ticdat, window=window, threshold=noise)

        # Filter Peaks
        goodpeaks = []
        tranges=[]
        diffs = np.diff(ticdat[:,0])
        for p in peaks:
            fwhm, range = ud.calc_FWHM(p[0], ticdat)
            index = ud.nearest(ticdat[:,0], p[0])
            if index >= len(diffs):
                index=len(diffs)-1
            localdiff = diffs[index]
            if p[0] - fwhm / 2. < mint or p[0] + fwhm / 2. > maxt or fwhm > 4 * window or fwhm < localdiff * 2 or range[
                0] == p[0] or range[1] == p[0]:
                print("Bad Peak", p, fwhm, range)
                pass
            else:
                print(p[0], fwhm)
                goodpeaks.append(p)
                tranges.append(range)
        self.chrompeaks = goodpeaks
        self.chrompeaks_tranges=tranges
        return goodpeaks, tranges
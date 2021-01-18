import os
import numpy as np
import unidec_modules.unidectools as ud
from metaunidec.mudeng import MetaUniDec
from unidec import UniDec
from copy import deepcopy

chrom_file_exts = [".raw", ".Raw", ".RAW", ".d", ".mzML.gz", ".mzML" ]


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
        self.procdata = None
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
        raise FileNotFoundError
        # return False, header

    def load_mzml(self, path, load_hdf5=True, *args, **kwargs):
        self.path = path
        name = os.path.splitext(path)[0]
        if name[-5:].lower() == ".mzml":
            name = name[:-5]
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

    def get_data_from_scans(self, scan_range=None):
        self.mzdata = self.chromdat.get_data(scan_range)
        self.procdata = None
        return self.mzdata

    def get_data_from_times(self, min, max):
        minscan = ud.nearest(self.ticdat[:, 0], min)
        if self.ticdat[minscan, 0] < min:
            minscan += 1
        maxscan = ud.nearest(self.ticdat[:, 0], max)
        if self.ticdat[maxscan, 0] > max:
            maxscan -= 1
        if maxscan <= minscan:
            maxscan = minscan + 1
        self.scans = [minscan, maxscan, min, max]

        attrs = {"timestart": min, "timeend": max,
                 "timemid": (min + max) / 2.,
                 "scanstart": minscan, "scanend": maxscan,
                 "scanmid": (minscan + maxscan) / 2.}
        self.attrs = attrs

        self.get_data_from_scans([minscan, maxscan])
        return self.mzdata

    def get_times(self):
        starts = []
        ends = []
        for s in self.data.spectra:
            starts.append(s.attrs["timestart"])
            ends.append(s.attrs["timeend"])
        return np.array(starts), np.array(ends)

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
        # ticdat = ud.gaussian_backgroud_subtract(ticdat, 100)
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
        tranges = []
        diffs = np.diff(ticdat[:, 0])
        for p in peaks:
            fwhm, range = ud.calc_FWHM(p[0], ticdat)
            index = ud.nearest(ticdat[:, 0], p[0])
            if index >= len(diffs):
                index = len(diffs) - 1
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
        self.chrompeaks_tranges = tranges
        return goodpeaks, tranges

    def add_manual_selection(self):
        self.data.add_data(self.mzdata, name=str(self.scans[2]), attrs=self.attrs, export=False)

    def add_regular_times(self):
        times = np.arange(0, np.amax(self.ticdat[:, 0]), self.config.time_window)
        self.data.clear()
        for i, t in enumerate(times):
            data = self.get_data_from_times(t, t + self.config.time_window)
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)

    def add_chrom_peaks(self):
        self.get_chrom_peaks()
        times = self.chrompeaks_tranges
        self.data.clear()
        for i, t in enumerate(times):
            data = self.get_data_from_times(t[0], t[1])
            self.data.add_data(data, name=str(t[0]), attrs=self.attrs, export=False)

    def add_sliding_window(self):
        if self.config.sw_scan_offset < 1:
            self.config.sw_scan_offset = 1
        tindex = np.arange(0, len(self.ticdat), int(self.config.sw_scan_offset))
        self.data.clear()
        for i in tindex:
            t = self.ticdat[i, 0]
            data = self.get_data_from_times(t, t + self.config.sw_time_window)
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)
        pass

    def add_list_times(self, starts, ends):
        self.data.clear()
        for i, t in enumerate(starts):
            data = self.get_data_from_times(t, ends[i])
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)

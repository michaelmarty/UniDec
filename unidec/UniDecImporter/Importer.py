import time
import os
import numpy as np
import unidec.tools as ud
from copy import deepcopy
from unidec.UniDecImporter.ImportTools import merge_spectra, merge_im_spectra, IndexedFile, IndexedScan
import subprocess
class Importer:
    def __init__(self, file_path, **kwargs):
        self._file_path = file_path
        self._params = kwargs
        self.ext = os.path.splitext(file_path)[1]
        self.filesize = os.stat(file_path).st_size
        self.scans = None
        self.times = None
        self.levels = None
        self.centroided = False
        self.polarity = "Positive"
        self.scan_range = None # Note, we will use inclusive scan range such that the last scan is self.scan_range[1]
        self.data = None
        self.immsdata = None
        self.cdms_support = False
        self.imms_support = False
        self.chrom_support = False
        self.centroid_threshold = 0.8
        self.thermo_support = False

        self.indexed_file = None


    def get_polarity(self, scan=None):
        return self.polarity

    def get_ms_order(self, scan=1):
        if self.levels is not None:
            index = self.get_scan_index(scan)
            return self.levels[index]
        else:
            return 1

    # list of all scans
    def get_all_scans(self):
        pass

    # averaged scans
    def get_avg_scan(self, scan_range=None, time_range=None):
        pass

    # Single scan data
    def get_single_scan(self, scan):
        pass

    def get_max_time(self):
        return self.times[-1]

    def get_max_scan(self):
        return self.scans[-1]

    def get_scans_from_times(self, time_range):
        mins = self.get_time_scan(time_range[0])
        if time_range[1] - time_range[0] > 0:
            maxs = self.get_time_scan(time_range[1])
            return [mins, maxs]
        else:
            return [mins, mins]

    def get_times_from_scans(self, scan_range):
        mint = self.get_scan_time(scan_range[0])
        if scan_range[1] - scan_range[0] > 1:
            maxt = self.get_scan_time(scan_range[1])
            avgt = (mint + maxt) / 2.
            return [mint, avgt, maxt]
        elif scan_range[1] - scan_range[0] == 0 or scan_range[1] - scan_range[0] == 1:
            maxt = self.get_scan_time(scan_range[1]+1)
            avgt = (mint + maxt) / 2.
            return [mint, avgt, maxt]
        else:
            # Assume it's a single scan
            return [mint, mint, mint]

    def get_tic(self):
        if not self.chrom_support:
            print("TIC data not supported for this file type:", self._file_path)
            raise Exception
        pass

    def index_scans(self, min_mz, bin_width):
        if not self.chrom_support:
            print("Indexing not supported for this file type:", self._file_path)
            raise Exception
        else:
            self.indexed_file = IndexedFile()
            scans = self.get_all_scans()
            for scan_number in self.scans:
                if self.get_ms_order(scan_number) == 1:
                    scan = scans[self.get_scan_index(scan_number)]
                    rt = self.get_scan_time(scan_number)
                    ind_scan = IndexedScan(scan, rt, scan_number, min_mz, bin_width)
                    self.indexed_file.indexed_scans.append(ind_scan)
        return

    def get_eic(self, mass, mz_tol, rt_range=None):
        if not self.chrom_support:
            print("EIC data not supported for this file type:", self._file_path)
            raise Exception
        else:
            if self.indexed_file is None:
                self.index_scans(0, 1)
            return self.indexed_file.extract_xic(mass, mz_tol, rt_range)

    def get_scan_time(self, scan):
        # Find index in self.scans
        index = self.get_scan_index(scan)
        t = self.times[index]
        return t

    def get_scan_index(self, scan):
        if scan in self.scans:
            index = np.where(np.array(self.scans) == scan)[0][0]
        elif scan < self.scans[0]:
            index = 0
        elif scan > self.scans[-1]:
            index = len(self.scans) - 1
        else:
            index = np.searchsorted(self.scans, scan)
        return index

    def get_time_scan(self, time):
        index = np.argmin(np.abs(self.times - time))
        scan = self.scans[index]
        return scan

    def check_centroided(self):
        scan_data = self.get_single_scan(self.scan_range[0])
        ratio = ud.get_autocorr_ratio(scan_data)
        print(f"Autocorr ratio {round(ratio, 3)}")
        if ratio < self.centroid_threshold:
            self.centroided = True
        else:
            self.centroided = False
        return self.centroided

    def scan_range_from_inputs(self, scan_range=None, time_range=None):
        if scan_range is None and time_range is None:
            scan_range = self.scan_range
        elif time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        scan_range = np.array(scan_range)
        if scan_range[0] < np.amin(self.scans):
            scan_range[0] = np.amin(self.scans)
        if scan_range[0] > np.amax(self.scans):
            scan_range[0] = np.amax(self.scans)
        if scan_range[1] > np.amax(self.scans) or scan_range[1] < scan_range[0]:
            scan_range[1] = np.amax(self.scans)

        print("Scan Range:", scan_range)
        return scan_range

    def avg_fast(self, scan_range=None, time_range=None):
        if self.data is None:
            self.get_all_scans()

        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        data = deepcopy(self.data)
        if scan_range is not None:
            startindex = self.get_scan_index(scan_range[0])
            endindex = self.get_scan_index(scan_range[1])
            data = data[startindex:endindex + 1]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans))

        if data is None or ud.isempty(data):
            print("Error: Empty Data Object")
            return None
        elif len(data) > 1:
            try:
                data = merge_spectra(data)
            except Exception as e:
                print("Merge Spectra Error 2:", e)
                print(data)
        else:
            data = data[0]

        return data

    def get_cdms_data(self, scan_range=None):
        if not self.cdms_support:
            print("CDMS data not supported for this file type:", self._file_path)
            raise Exception

        raw_dat = self.get_all_scans()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])
        scans = np.ones(len(mz))  # Update this
        it = np.ones(len(mz))

        data_array = np.transpose([mz, intensity, scans, it])

        return data_array

    def get_imms_avg_scan(self, scan_range=None, time_range=None, mzbins=1):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception
        start_time = time.perf_counter()
        if self.immsdata is None:
            self.get_all_imms_scans()

        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        data = deepcopy(self.immsdata)
        if scan_range is not None:
            startindex = self.get_scan_index(scan_range[0])
            endindex = self.get_scan_index(scan_range[1])
            data = data[startindex:endindex + 1]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans))

        if data is None or ud.isempty(data):
            print("Error: Empty Data Object")
            return None

        # Need to merge to get 2D from sparse array
        data = merge_im_spectra(data, mzbins=mzbins)
        # try:
        #     data = merge_im_spectra(data, mzbins=mzbins)
        # except Exception as e:
        #     concat = np.concatenate(data)
        #     sort = concat[concat[:, 0].argsort()]
        #     data = ud.removeduplicates(sort)
        #     print("Mark2", e)

        print("Import Time:", time.perf_counter() - start_time, "s")
        return data

    def get_imms_scan(self, s):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception
        start_time = time.perf_counter()
        if self.immsdata is None:
            self.get_all_imms_scans()

        index = self.get_scan_index(s)
        data = self.immsdata[index]
        print("Import Time:", time.perf_counter() - start_time, "s")
        return data

    def get_all_imms_scans(self):
        if not self.imms_support:
            print("IMMS data not supported for this file type:", self._file_path)
            raise Exception

    def close(self):
        pass



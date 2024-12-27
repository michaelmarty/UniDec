import numpy as np
import unidec.tools as ud
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.ImportTools import merge_spectra
from unidec.UniDecImporter.Waters import MassLynxRawScanReader as MLRSR, MassLynxRawInfoReader as MLRIR, \
    MassLynxRawChromatogramReader as MLCR

# Note, waters scans are 0 indexed, so we need to subtract 1 from the scan number

class WatersDataImporter(Importer):
    """
    Imports Waters data files.

    Note: Waters scans are 0 indexed, so the first scan is scan 0, not scan 1.
    """
    def __init__(self, path, function=0, *args, **kwargs):
        super().__init__(path, **kwargs)
        print("Reading Waters Data:", path)
        self.function = int(function)
        self.reader = MLRIR.MassLynxRawInfoReader(path)
        self.readerMS = MLRSR.MassLynxRawScanReader(path)
        self.init_scans()

        self.imms_support = True
        self.cdms_support = False
        self.chrom_support = True


    def init_scans(self):
        self.nfunc = self.reader.GetNumberofFunctions()
        try:
            fn = self.reader.GetFunctionType(self.function)
        except:
            try:
                print("Function number not found in Raw file.", self.function, "Defaulting to:", 0)
                self.function = 0
                fn = self.reader.GetFunctionType(self.function)
            except:
                try:
                    print("Function number not found in Raw file.", self.function, "Defaulting to:", 1)
                    self.function = 1
                    fn = self.reader.GetFunctionType(self.function)
                except:
                    print("Function number not found in Raw file:", self.function)
                    raise IOError

        print("Waters Data Type:", self.reader.GetFunctionTypeString(fn), "Function:", self.function)
        self.maxscans = self.reader.GetScansInFunction(self.function)
        self.scan_range = [1, self.maxscans]
        print("Waters Scan Range", self.scan_range)
        self.scans = np.arange(self.scan_range[0], self.scan_range[1]+1)
        self.times = []
        for s in self.scans:
            self.times.append(self.reader.GetRetentionTime(self.function, s-1))
        self.times = np.array(self.times)

    def get_stats(self):
        self.reader = MLRIR.MassLynxRawInfoReader(self._file_path)
        self.stat_nums = self.reader.GetItemsInFunction(self.function, 0)
        self.stat_vals = self.reader.GetScanItems(self.function, 0, self.stat_nums)
        self.stat_names = self.reader.GetScanItemString(self.stat_nums)

    def get_stat_code(self, stat_code):
        self.get_stats()
        if stat_code in self.stat_nums:
            index = ud.nearestunsorted(self.stat_nums, stat_code)
            value = self.stat_vals[index]
            print(self.stat_names[index] + ": ", value)
        else:
            print("Stat Code not found in data")
            value = 0
        return value

    def get_stat_name(self, stat_name):
        self.get_stats()
        for i, sn in enumerate(self.stat_names):
            if sn.lower() == stat_name.lower():
                print(sn, ":", self.stat_vals[i])
                return self.stat_vals[i]
        print("Could not find stat in stat names")
        print("Names: ", self.stat_names)
        return 0


    def get_avg_scan(self, scan_range=None, time_range=None, mzbins=None):
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        mzs, ivals = self.readerMS.combineScan(self.function, np.arange(scan_range[0]-1, scan_range[1]))
        data = np.transpose([mzs, ivals])
        if mzbins is None or float(mzbins) == 0:
            return data
        else:
            data = merge_spectra([data], mzbins=mzbins, type="Integrate")
            return data

    def get_all_scans(self):
        self.data = []
        for s in self.scans:
            impdat = np.transpose(self.readerMS.ReadScan(self.function, s-1))
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
        return self.data

    def get_single_scan(self, scan):
        #This should return a format which is similar to the numpy 2 column format
        dat = self.readerMS.ReadScan(self.function, scan-1)
        mz, intensity = dat
        res = np.column_stack((mz, intensity))
        return res

    def get_tic(self):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self._file_path)
        tic = np.transpose(self.readerLC.ReadTIC(self.function))
        # self.readerLC.__del__()
        return tic

    def get_bpi(self):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self._file_path)
        tic = np.transpose(self.readerLC.ReadBPI(self.function))
        return tic

    def get_eic(self, mass=811, tolerance=0.10):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self._file_path)
        tic = np.transpose(self.readerLC.ReadMassChromatogram(self.function, mass, tolerance, False))
        return tic

    def get_all_imms_scans(self):
        # self.mindrift = self.get_stat_name('Minimum Drift Time Channel')
        # self.maxdrift = self.get_stat_name('Maximum Drift Time Channel')
        self.trf = self.get_stat_name("Transport RF")
        self.pusher = np.floor((1. / float(self.trf)) * 1000 * 1000)
        print("Pusher Freq:", self.pusher)
        self.immsdata = []
        for s in self.scans:
            self.immsdata.append(self.get_imms_scan(s))
        return self.immsdata

    def get_imms_scan(self, s):
        scandat = []
        for i in range(0, 200):
            o = self.readerMS.ReadDriftScan(self.function, s - 1, i)
            if len(o) == 0:
                continue
            dts = np.ones(len(o[0])) * i
            dtdat = np.transpose([o[0], dts, o[1]])
            scandat.extend(dtdat)
        return np.array(scandat)

    def get_polarity(self, scan=None):
        line = self.reader.GetIonModeString(self.function)
        if "+" in line:
            print("Polarity: Positive")
            return "Positive"
        if "-" in line:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")
        return None

    def get_ms_order(self, scan=None):
        ms_level =  self.reader.GetFunctionTypeString(self.function)
        count = ms_level.count("MS")
        return count

    def close(self):
        del self.reader
        del self.readerMS


if __name__ == "__main__":
    test = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw"
    d = WatersDataImporter(test)
    d.close()






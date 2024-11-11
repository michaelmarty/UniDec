import numpy as np

from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.MZML.mzML import merge_spectra
from unidec.UniDecImporter.Waters import MassLynxRawScanReader as MLRSR, MassLynxRawInfoReader as MLRIR, \
    MassLynxRawChromatogramReader as MLCR

class WatersDataImporter(Importer):
    def __init__(self, path, do_import = False, function = 0, *args, **kwargs):
        super().__init__(path, **kwargs)
        print("Reading Data:", path)
        self.path = path

        self.reader = MLRIR.MassLynxRawInfoReader(path)
        self.readerMS = MLRSR.MassLynxRawScanReader(path)

        self.nfunc = self.reader.GetNumberofFunctions()
        self.function = int(function)
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
        self.scanrange = [0, self.maxscans]
        self.scans = np.arange(self.scanrange[0], self.scanrange[1])
        self.times = []
        for s in self.scans:
            self.times.append(self.reader.GetRetentionTime(self.function, s))
        self.times = np.array(self.times)
        if do_import:
            self.slow_get_data()
        else:
            self.data = None
        # self.readerMS.__del__()
        # self.reader.__del__()

    def slow_get_data(self):
        self.data = []
        for s in self.scans:
            impdat = np.transpose(self.readerMS.ReadScan(self.function, s))
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
        self.data = np.array(self.data)

    def slow_get_data_scans(self, scan_range, mzbins=None):
        if self.data is None:
            self.slow_get_data()
        data = deepcopy(self.data)
        if scan_range is not None:
            data = data[scan_range[0]:scan_range[1]]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans))

        if len(data) > 1:
            try:
                data = merge_spectra(data, mzbins=mzbins)
            except Exception as e:
                concat = np.concatenate(data)
                sort = concat[concat[:, 0].argsort()]
                data = ud.removeduplicates(sort)
                print(e)
        elif len(data) == 1:
            data = data[0]
        else:
            data = data
        return data

    def fast_get_data_scans(self, scan_range, mzbins=None):
        scan_range = np.array(scan_range)
        # if scan_range[0]<1:
        #    scan_range[0]=1
        if scan_range[1] > self.maxscans:
            scan_range[1] = self.maxscans

        mzs, ivals = self.readerMS.CombineScan(self.function, np.arange(scan_range[0], scan_range[1]))
        data = np.transpose([mzs, ivals])
        if mzbins is None or float(mzbins) == 0:
            return data
        else:
            data = merge_spectra([data], mzbins=mzbins, type="Integrate")
            return data

    def get_data(self, scan_range=None, time_range=None, mzbins=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """

        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        if scan_range is None:
            scan_range = self.scanrange

        if True:
            data = self.fast_get_data_scans(scan_range, mzbins)
        '''
        try:
            data=self.fast_get_data_scans(scan_range, mzbins)
        except Exception as e:
            print("ERROR with fast Waters combineScans, using slow method: ", e)
            data=self.slow_get_data_scans(scan_range, mzbins)'''

        return data

    def get_tic(self):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self.path)
        tic = np.transpose(self.readerLC.ReadTIC(self.function))
        # self.readerLC.__del__()
        return tic

    def get_bpi(self):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self.path)
        tic = np.transpose(self.readerLC.ReadBPI(self.function))
        # self.readerLC.__del__()
        return tic

    def get_eic(self, mass=811, tolerance=0.10):
        self.readerLC = MLCR.MassLynxRawChromatogramReader(self.path)
        tic = np.transpose(self.readerLC.ReadMassChromatogram(self.function, mass, tolerance, False))
        return tic

    def get_max_time(self):
        return self.times[len(self.times) - 1]

    def get_max_scans(self):
        return self.scanrange[1] - 1

    def get_scans_from_times(self, time_range):
        boo1 = self.times >= time_range[0]
        boo2 = self.times < time_range[1]
        try:
            min = np.amin(self.scans[boo1])
            max = np.amax(self.scans[boo2])
        except:
            min = -1
            max = -1
        return [min, max]

    def get_times_from_scans(self, scan_range):
        boo1 = self.scans >= scan_range[0]
        boo2 = self.scans < scan_range[1]
        boo3 = np.logical_and(boo1, boo2)
        min = np.amin(self.times[boo1])
        max = np.amax(self.times[boo2])
        try:
            avg = np.mean(self.times[boo3])
        except:
            avg = min
        return [min, avg, max]

    def get_stats(self):
        self.reader = MLRIR.MassLynxRawInfoReader(self.path)
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

    def get_IMMS_data(self):
        self.mindrift = self.get_stat_name('Minimum Drift Time Channel')
        self.maxdrift = self.get_stat_name('Maximum Drift Time Channel')
        self.trf = self.get_stat_name("Transport RF")
        self.pusher = np.floor((1. / float(self.trf)) * 1000 * 1000)
        print(self.pusher)
        self.immsdata = []
        for s in self.scans:
            scandat = []
            for i in range(0, 200):
                o = np.array(self.readerMS.ReadDriftScan(self.function, s, i))
                scandat.append(o)
                '''
                mz = o[0]
                try:
                    print(np.amax(mz))
                except:
                    pass'''
            self.immsdata.append(scandat)
        return self.immsdata

    def get_polarity(self):
        line = self.reader.GetIonModeString(self.function)
        print(line)
        if "+" in line:
            print("Polarity: Positive")
            return "Positive"
        if "-" in line:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")
        return None

if __name__ == "__main__":
    test = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw"
    importer = WatersDataImporter(test)
    print(type(importer))
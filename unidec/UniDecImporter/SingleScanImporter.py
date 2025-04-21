import os
import numpy as np
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.ImportTools import header_test


class SingleScanImporter(Importer):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.scans = [1]
        self.times = [0]
        self.scan_range = [1, 1]
        self.centroided = False
        self.polarity = "Positive"
        self.scan_number = 1
        self.injection_time = 1
        self.load_data()
        self.cdms_support = True
        self.imms_support = True
        self.chrom_support = False

        # Check if file exists
        if not os.path.isfile(filename):
            print("File not found:", filename)
            raise FileNotFoundError

    def __len__(self):
        return len(self.data)

    def load_data(self):
        if self.data is not None:
            return self.data
        if self.ext == ".txt" or self.ext == ".dat":
            self.data = np.loadtxt(self._file_path, skiprows=header_test(self._file_path))
        elif self.ext == ".csv":
            self.data = np.loadtxt(self._file_path, delimiter=",", skiprows=header_test(self._file_path),
                                   usecols=(0, 1))
        elif self.ext == '.npz':
            self.data = np.load(self._file_path, allow_pickle=True)['data']
        elif self.ext == '.bin':
            self.data = np.fromfile(self._file_path, dtype=np.float64)
        else:
            print("Unsupported file extension in SingleScanImporter", )
            return None
        return self.data


    def get_all_scans(self):
        return [self.load_data()]

    def get_single_scan(self, scan=None):
        return self.load_data()

    def get_avg_scan(self, scan_range=None, time_range=None):
        return self.load_data()

    def get_cdms_data(self):
        ext = os.path.splitext(self._file_path)[1]
        raw_dat = self.get_all_scans()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])

        try:
            scans = np.concatenate([d[:, 2] for d in raw_dat])
        except Exception as e:
            print("No scan data in CSV file, populating with 1's")
            scans = np.ones_like(intensity)
        try:
            it = np.concatenate([d[:, 3] for d in raw_dat])
        except Exception as e:
            print("No injection time data in CSV file, populating with 1's")
            it = np.ones_like(mz)
        # elif ext.lower() == '.bin':
        #     try:
        #         raw_dat = raw_dat.reshape((int(len(raw_dat) / 3), 3))
        #     except Exception as e:
        #         raw_dat = raw_dat.reshape((int(len(raw_dat) / 2)), 2)
        #     mz = raw_dat[:, 0] > 0
        #     intensity = raw_dat[:, 1] > 0
        #     it = np.ones_like(mz)
        #     scans = np.ones_like(intensity)
        # else:
        #     it = np.ones_like(mz)
        #     scans = np.ones_like(intensity)

        return np.transpose([mz, intensity, scans, it])

    def get_imms_avg_scan(self, scan_range=None, time_range=None, mzbins=None):
        self.immsdata = np.loadtxt(self._file_path, skiprows=header_test(self._file_path))
        return self.immsdata

    def get_all_imms_scans(self):
        return [self.get_imms_avg_scan()]

    def get_imms_scan(self, s):
        return self.get_imms_avg_scan()


if __name__ == "__main__":
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\CDMS\\test_csv_cdms.csv"
    path = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Skippy\\HT SEC paper data\\Bgal GroEL\\20240412 Bgal GroEL bit5 zp10 3_2024-04-16-11-38-01_unidecfiles\\20240412 Bgal GroEL bit5 zp10 3_2024-04-16-11-38-01_rawdata.npz"
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\SingleScan\\test_csv.csv"
    importer = SingleScanImporter(path)
    print(importer.get_scan_time(12))








import os
import numpy as np
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.ImportTools import header_test


class SingleScanImporter(Importer):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.data = None
        self.load_data()
        self.scans = [1]
        self.times = [0]
        self.centroided = True
        self.polarity = "Positive"

    def __len__(self):
        return len(self.data)
    # This needs to actually aggregrate the data if it is not already aggregated
    # Need way to determine if data is already merged or not
    def load_data(self):
        if self.ext == ".txt" or self.ext == ".dat":
            self.data = np.loadtxt(self._file_path, skiprows=header_test(self._file_path))
        elif self.ext == ".csv":
            self.data = np.loadtxt(self._file_path, delimiter=",", skiprows=header_test(self._file_path),
                                   usecols=(0, 1))
        elif self.ext == '.npz':
            self.data = np.load(self._file_path, allow_pickle=True)['data']
        else:
            print("Unsupported file extension in SingleScanImporter", )
            return None

    def get_data(self):
        return self.data

    def grab_data(self):
        return self.data

    def grab_scan_data(self, scan=None):
        return self.data

    def get_scan_time(self, scan=None):
        return self.times[0]

    def get_ms_order(self, scan=None):
        return 1

    def get_polarity(self, scan=None):
        return self.polarity


if __name__ == "__main__":
    path = "C:\\Python\\UniDec3\\TestSpectra\\test_data.csv"
    importer = SingleScanImporter(path)
    print(len(importer))
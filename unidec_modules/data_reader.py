from multiplierz.mzAPI import mzFile
from mzMLimporter import merge_spectra
from copy import deepcopy
import unidectools as ud
import numpy as np


class DataImporter:
    """
    Imports mzML data files.
    """

    def __init__(self, path, *args, **kwargs):
        """
        Imports mzML file, adds the chromatogram into a single spectrum.
        :param path: .mzML file path
        :param args: arguments (unused)
        :param kwargs: keywords (unused)
        :return: mzMLimporter object
        """
        print "Reading Data:", path
        self.msrun = mzFile(path)
        self.scanrange = self.msrun.scan_range()
        self.scans = np.arange(self.scanrange[0], self.scanrange[1] + 1)
        self.times = []
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.scan(s))
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
            self.times.append(self.msrun.scan_time_from_scan_name(s))
        self.times = np.array(self.times)
        self.data = np.array(self.data)

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        data = deepcopy(self.data)
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print "Getting times:", time_range

        if scan_range is not None:
            data = data[scan_range[0]:scan_range[1]]
            print "Getting scans:", scan_range
        else:
            print "Getting all scans, length:", len(self.scans)

        if len(data) > 1:
            try:
                data = merge_spectra(data)
            except Exception, e:
                concat = np.concatenate(data)
                sort = concat[concat[:, 0].argsort()]
                data = ud.removeduplicates(sort)
                print e
        elif len(data) == 1:
            data = data[0]
        else:
            data = data

        return data

    def get_tic(self):
        return self.msrun.xic()

    def get_max_time(self):
        times = self.msrun.time_range()
        return times[1]

    def get_max_scans(self):
        scans = self.msrun.scan_range()
        return scans[1]

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
        try:
            min = np.amin(self.times[boo1])
            max = np.amax(self.times[boo2])
            avg = np.mean(self.times[boo3])
        except:
            min = -1
            max = -1
        return [min, avg, max]


if __name__ == "__main__":
    test = "C:\\Data\\test.RAW"
    d = DataImporter(test)

    print d.get_times_from_scans([15, 30])

    exit()

    d = DataImporter(test).get_data(time_range=(0, 1))
    import matplotlib.pyplot as plt

    plt.plot(d[:, 0], d[:, 1])
    plt.show()

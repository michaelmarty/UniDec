import numpy as np
#import unidec_modules.thermo_reader.RawFileReader
import time

class ThermoDataImporter:
    """
    Imports Thermo data files.
    """

    def __init__(self, path, *args, **kwargs):
        """
        Imports Thermo file.
        :param path: Thermo .Raw file path
        :param args: arguments (unused)
        :param kwargs: keywords (unused)
        :return: mzMLimporter object
        """
        from unidec_modules.thermo_reader.RawFileReader import RawFileReader as rr
        print("Reading Thermo Data:", path)
        self.msrun = rr(path)
        self.scanrange = self.msrun.scan_range()
        # print(self.scanrange)
        self.scans = np.arange(self.scanrange[0], self.scanrange[1])
        self.times = []
        self.data = None
        for s in self.scans:
            self.times.append(self.msrun.scan_time_from_scan_name(s))
        self.times = np.array(self.times)
        # print(len(self.data), len(self.times), len(self.scans))

    def grab_data(self):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetSpectrum(s))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
        self.data = np.array(self.data)
        return self.data

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if scan_range is not None:
            scan_range = np.array(scan_range, dtype=np.int)
            scan_range = scan_range + 1
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)
        if scan_range is None:
            scan_range = [np.amin(self.scans), np.amax(self.scans)]
        scan_range = np.array(scan_range, dtype=np.int)
        print("Scan Range:", scan_range)

        if scan_range[0] < np.amin(self.scans):
            scan_range[0] = np.amin(self.scans)
        if scan_range[1] > np.amax(self.scans):
            scan_range[1] = np.amin(self.scans)

        if scan_range[1] - scan_range[0] > 1:
            data = np.array(list(self.msrun.GetAverageSpectrum(scan_range)))
        else:
            impdat = np.array(self.msrun.GetSpectrum(scan_range[0]))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            data = impdat

        return data

    def get_tic(self):
        return self.msrun.GetChromatogram()

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
        if scan_range[1] - scan_range[0] > 1:
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
        else:
            t = self.times[scan_range[0]]
            return [t, t, t]


if __name__ == "__main__":
    test = u"C:\Python\\UniDec3\TestSpectra\\test.RAW"
    tstart = time.perf_counter()
    d = ThermoDataImporter(test).get_data()
    # d.get_data(time_range=(0, 10))
    print("ImportData: %.2gs" % (time.perf_counter() - tstart))
    # import matplotlib.pyplot as plt
    # plt.plot(d[:, 0], d[:, 1])
    # plt.show()

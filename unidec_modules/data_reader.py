from multiplierz.mzAPI import mzFile
from unidec_modules import unidectools as ud
from unidec_modules.mzMLimporter import merge_spectra
from copy import deepcopy
import numpy as np
import time

# import sys
# import wx

def register():
    print("Trying to Register Interfaces...")
    from multiplierz.mzAPI.management import registerInterfaces
    try:
        print("Registering...")
        registerInterfaces()
    except Exception as e:
        print("Failed Interface Registration:", e)
        print("NOTE: TRY RUNNING AS ADMINISTRATOR")
        pass


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
        # TODO make this work
        """
        del sys.modules[mzFile]
        print "Breaking"
        if mzFile not in sys.modules or merge_spectra not in sys.module:
            dlg = dlg = wx.MessageDialog(parent=None, message='Please install multiplierz and MSFileReader',
                               caption='Error', style=wx.OK)
            dlg.ShowModal()
            return
            """
        print("Reading Data:", path)
        try:
            self.msrun = mzFile(path)
        except:
            register()
            self.msrun = mzFile(path)
        self.scanrange = self.msrun.scan_range()
        # print(self.scanrange)
        self.scans = np.arange(self.scanrange[0], self.scanrange[1])
        self.times = []
        self.data = None
        for s in self.scans:
            s=s-1
            try:
                self.times.append(self.msrun.scan_time_from_scan_name(s))
            except Exception as e:
                try:
                    t = self.msrun.info[s][0]
                    self.times.append(t)
                except Exception as e2:
                    try:
                        t = self.msrun.scan_info()[s][0]
                        self.times.append(t)
                    except Exception as e3:
                        print("Error getting scan times:", e, e2, e3)
                        print("Using Scan rather than Time")
                        self.times.append(s)
        self.times = np.array(self.times)
        # print(self.times)
        # print(len(self.data), len(self.times), len(self.scans))

    def grab_data(self):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.scan(s))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
        self.data = np.array(self.data, dtype=object)
        return self.data

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        try:
            if scan_range is not None:
                scan_range = np.array(scan_range, dtype=np.int)
                scan_range = scan_range+1
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

            if scan_range[1]-scan_range[0] > 1:
                data = np.array(list(self.msrun.average_scan(int(scan_range[0]), int(scan_range[1]), filter="Full")))
            else:
                impdat = np.array(self.msrun.scan(scan_range[0]))  # May want to test this.
                impdat = impdat[impdat[:, 0] > 10]
                data = impdat

        except Exception as e:
            print("Failed native raw averaging. Using Python averaging.")
            if self.data is None:
                self.grab_data()
            data = deepcopy(self.data)
            if time_range is not None:
                scan_range = self.get_scans_from_times(time_range)
                print("Getting times:", time_range)

            if scan_range is not None:
                data = data[scan_range[0]:scan_range[1]]
                print("Getting scans:", scan_range)
            else:
                print("Getting all scans, length:", len(self.scans))

            if len(data) > 1:
                try:
                    data = merge_spectra(data)
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

    def get_tic(self):
        return self.msrun.xic(filter="Full")

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
    test = "Z:\Group Share\Group\\Archive\\Scott\\test.RAW"
    test = "C:\Data\Others\Agilent\\2019_05_15_bsa_ccs_02.d"
    test = u"C:\Python\\UniDec3\TestSpectra\\test.RAW"
    #tstart = time.perf_counter()

    #d = DataImporter(test)
    #print(d.get_times_from_scans([0, 10]))
    #d.get_data(time_range=(0, 10))
    #print("ImportData: %.2gs" % (time.perf_counter() - tstart))
    #import matplotlib.pyplot as plt

    #plt.plot(d[:, 0], d[:, 1])
    #plt.show()

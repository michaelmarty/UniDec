import os.path

from unidec.UniDecImporter.Agilent.MZFILE import MZFile as mzFile, MZFile
from unidec.UniDecImporter.Importer import *
from unidec import tools as ud
from unidec.UniDecImporter.MZML.mzML import merge_spectra
from copy import deepcopy
import numpy as np

import time


# When is the msrun scaninfo even being populated ever


class AgilentImporter(Importer):
    def __init__(self, path, *args, **kwargs):
        print(f"Path exists: {os.path.exists(path)}")
        print(f"Is file: {os.path.isfile(path)}")
        print("Reading Data:", path)
        self.msrun = mzFile(path)
        self.scanrange = self.msrun.scan_range()
        self.scans = np.arange(self.scanrange[0], self.scanrange[1])
        self.times = []
        # self.times = [self.get_times_from_scans(self.scanrange)]
        self.data = None

        self.datascans = []
        curr_info = self.msrun.scan_info()
        for i in curr_info:
            self.datascans.append(self.msrun.scan(i[2]))
            self.times.append(i[0])
        self.times = np.array(self.times)

    def grab_data(self):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.scan(int(s)))
            try:
                impdat = impdat[impdat[:, 0] > 10]
            except:
                print("Error Importing Data, Scan:", s, "Data:", impdat)
            self.data.append(impdat)
        self.data = np.array(self.data, dtype=object)
        return self.data


    def grab_scan_data(self, scan):
        impdat = np.array(self.msrun.scan(int(scan)))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def get_data(self, scan_range=None, time_range=None, mzbins=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        try:
            if scan_range is not None:
                scan_range = np.array(scan_range, dtype=int)
                scan_range = scan_range
            if time_range is not None:
                scan_range = self.get_scans_from_times(time_range)
                print("Getting times:", time_range)
            if scan_range is None:
                scan_range = [np.amin(self.scans), np.amax(self.scans)]
            print("Scan Range:", scan_range)

            if scan_range[0] < np.amin(self.scans):
                scan_range[0] = np.amin(self.scans)
            if scan_range[1] > np.amax(self.scans):
                scan_range[1] = np.amin(self.scans)

            if scan_range[1] - scan_range[0] > 1:
                # noinspection PyUnresolvedReferences
                # print("I got here")
                # data = np.array(list(self.msrun.average_scan(int(scan_range[0]), int(scan_range[1]))))
                data = np.array(list(self.msrun.average_scan(self.scanrange, self.datascans)))
                # print("made it past before dying")
            else:
                impdat = np.array(self.msrun.scan(scan_range[0]))
                impdat = impdat[impdat[:, 0] > 10]
                data = impdat


        except Exception as e:
            print(e)
            print("Failed native raw averaging. Using Python averaging.")
            if self.data is None:
                self.grab_data()
            data = deepcopy(self.data)
            if time_range is not None:
                scan_range = self.get_scans_from_times(time_range)
                print("Getting times:", time_range)
            if scan_range is not None and len(scan_range) == 2:
                data = data[scan_range[0]:scan_range[1]]
            elif scan_range is not None:
                print("Getting scan:", scan_range[0])
                data = data[scan_range[0]]
                return data
            else:
                print("Getting all scans, length:", len(self.scans))

            if len(data) > 1:
                try:
                    data = merge_spectra(data, mzbins=mzbins, type="Integrate")
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

    # noinspection PyArgumentList
    def get_tic(self):
        try:
            xic = self.msrun.xic(filter="Full")
        except:
            xic = self.msrun.xic()
        return xic

    def get_max_time(self):
        times = self.msrun.time_range()
        return times[1]

    def get_max_scans(self):
        # noinspection PyUnresolvedReferences
        scans = self.msrun.scan_range()
        return scans[1]

    def get_scans_from_times(self, time_range):
        boo1 = self.times >= time_range[0]
        boo2 = self.times < time_range[1]
        try:
            minimum = np.amin(self.scans[boo1])
            maximum = np.amax(self.scans[boo2])
        except:
            minimum = -1
            maximum = -1
        return [minimum, maximum]

    def get_times_from_scans(self, scan_range):
        if scan_range[1] - scan_range[0] > 1:
            boo1 = self.scans >= scan_range[0]
            boo2 = self.scans < scan_range[1]
            boo3 = np.logical_and(boo1, boo2)
            minimum = np.amin(self.times[boo1])
            maximum = np.amax(self.times[boo2])
            try:
                avg = np.mean(self.times[boo3])
            except:
                avg = minimum
            return [minimum, avg, maximum]
        else:
            t = self.times[scan_range[0]]
            return [t, t, t]

    def get_polarity(self, scan=0):

        self.scan_info = self.msrun.scan_info()
        line = self.scan_info[scan]
        # previous checking for +
        if "Positive" in line:
            print("Polarity: Positive")
            return "Positive"
        # previous checking for -
        if "Negative" in line:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")
        return None
        pass

# if __name__ == '__main__':
#     path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d"
#     test = AgilentImporter(path)
#     print(type(test))
#     res = test.grab_data()
#     for i in res:
#         print(i)
#     print("Finished")

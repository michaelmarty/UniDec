from copy import deepcopy
import numpy as np
from unidec.UniDecImporter.Agilent.MZFILE import MZFile as mzFile
from unidec.UniDecImporter.Importer import *
from unidec import tools as ud
from unidec.UniDecImporter.ImportTools import merge_spectra

# When is the msrun scaninfo even being populated ever


class AgilentImporter(Importer):
    """
    Imports Agilent data files.

    Note: Agilent scans are 0 indexed, so the first scan is scan 0, not scan 1.
    """
    def __init__(self, path, **kwargs):
        super().__init__(path, **kwargs)
        print("Reading Agilent Data:", path)
        self.msrun = mzFile(path)
        self.init_scans()

        self.cdms_support = False
        self.imms_support = False
        self.chrom_support = True

    def init_scans(self):
        self.times = []
        self.scans = []

        curr_info = self.msrun.scan_info()
        for i in curr_info:
            self.scans.append(i[2]+1)
            self.times.append(i[0])

        self.times = np.array(self.times)
        self.scans = np.array(self.scans)
        self.scan_range = [np.amin(self.scans), np.amax(self.scans)]

    def get_all_scans(self):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.scan(int(s)-1))
            try:
                impdat = impdat[impdat[:, 0] > 10]
            except Exception as e:
                print("Error Importing Data, Scan:", s, "Data:", impdat, "Error:", e)
            self.data.append(impdat)
        return self.data


    def get_single_scan(self, scan, threshold=-1):
        try:
            impdat = np.array(self.msrun.scan(int(scan)-1))
            impdat = np.transpose([impdat[:, 0], impdat[:, 1]])
            impdat = impdat[impdat[:, 0] > 10]
            if threshold >= 0:
                impdat = impdat[impdat[:, 1] > threshold]
            return impdat
        except Exception as e:
            print(f"Error in get_single_scan for scan {scan}: {e}")
            return None


    def get_avg_scan(self, scan_range=None, time_range=None, mzbins=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        self.datascans = []
        for s in self.scans:
            self.datascans.append(self.msrun.scan(int(s-1)))

        if scan_range[1] - scan_range[0] > 1:
            # noinspection PyUnresolvedReferences
            data = np.array(list(self.msrun.average_scan(list(np.array(scan_range)-1), self.datascans)))
        else:
            impdat = np.array(self.msrun.scan(int(scan_range[0]-1)))
            impdat = impdat[impdat[:, 0] > 10]
            data = impdat

        return data

    # noinspection PyArgumentList
    def get_tic(self):
        try:
            xic = self.msrun.xic(filter="Full")
        except:
            xic = self.msrun.xic()
        return xic


    def get_polarity(self, scan=1):
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

    def get_ms_order(self, scan=1):
        scan_info = self.msrun.scan_info()
        for line in scan_info:
            if scan == int(line[2]):
                if "MS1" in line:
                    return 1
                elif "MS2" in line:
                    return 2
        print("MS Order: Unknown")
        return 1

    def close(self):
        self.msrun.close()


if __name__ == '__main__':
    import time
    # path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d"
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\test_agilent.d"
    #path = "C:\\Data\\DataTypeCollection\\test_agilent.d"
    d = AgilentImporter(path)
    sing = d.get_single_scan(1)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use("WxAgg")
    plt.plot(sing[:, 0], sing[:, 1])
    plt.show()




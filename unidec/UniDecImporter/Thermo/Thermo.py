import numpy as np
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.Thermo.RawFileReader import RawFileReader as rr


class ThermoImporter(Importer):
    """
    Imports Thermo data files.

    Note: Thermo scans are 1 indexed, so the first scan is scan 1, not scan 0.
    """
    def __init__(self, path, silent=False, **kwargs):
        super().__init__(path, **kwargs)

        if not silent:
            print("Launching Thermo Importer. If it fails, try to unblock the zip before unzipping")
        print("Reading Thermo Data:", path)
        self.msrun = rr(path)
        self.init_scans()

        self.cdms_support = True
        self.imms_support = False
        self.chrom_support = True

    def init_scans(self):
        self.scan_range = self.msrun.scan_range
        self.scans = np.arange(self.scan_range[0], self.scan_range[1] + 1)
        self.times = []
        for s in self.scans:
            self.times.append(self.msrun.scan_time_from_scan_name(s))

        self.times = np.array(self.times)
        print("Number of Scans", len(self.scans))

    def get_all_scans(self, threshold=-1):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetSpectrum(s))
            impdat = impdat[impdat[:, 0] > 10]
            if threshold >= 0:
                impdat = impdat[impdat[:, 1] > threshold]
            self.data.append(impdat)
        return self.data

    def get_single_scan(self, s):
        impdat = np.array(self.msrun.GetSpectrum(s))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def grab_centroid_data(self, s):
        impdat = np.array(self.msrun.GetCentroidSpectrum(s))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def get_avg_scan(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        if scan_range[1] - scan_range[0] > 1:
            print("Getting Data from Scans:", scan_range)
            scan_range = [scan_range[0], scan_range[1]]
            data = np.array(list(self.msrun.GetAverageSpectrum(scan_range)))
        else:
            print("Getting Data from Scan:", scan_range[0])
            impdat = self.get_single_scan(scan_range[0])
            data = impdat

        return data

    def get_tic(self):
        return self.msrun.GetChromatogram()

    def get_eic(self, mass_range=None, scan_range=None):
        if mass_range is None:
            mass_range = [0, 1000000]
        return self.msrun.Get_EIC(massrange=mass_range, scanrange=scan_range)

    def get_inj_time_array(self):
        its = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            try:
                it = float(it)
            except:
                print("Error in scan header:", i, s, it)
                it = 1
            its.append(it)
        return np.array(its)

    def get_analog_voltage1(self):
        vs = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            vs.append(an1)
        return np.array(vs)

    def get_analog_voltage2(self):
        vs = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            vs.append(an2)
        return np.array(vs)

    def get_polarity(self, scan=1):
        # print(dir(self.msrun.source))
        scan_mode = self.msrun.source.GetScanEventStringForScanNumber(scan)
        if "+" in scan_mode:
            print("Polarity: Positive")
            return "Positive"
        if "-" in scan_mode[:10]:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")
        return None

    def get_ms_order(self, scan=1):
        order = self.msrun.GetMSOrder(scan)
        return order

    def get_isolation_mz_width(self, s):
        scanFilter = self.msrun.GetScanFilter(s)
        reaction = scanFilter.GetReaction(0)
        mz = reaction.PrecursorMass
        width = reaction.IsolationWidth
        return mz, width

    def get_cdms_data(self, scan_range=None):
        raw_dat = self.get_all_scans(threshold=0)
        scans = self.scans

        it = 1. / self.get_inj_time_array()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        scans = np.concatenate([s * np.ones(len(raw_dat[i])) for i, s in enumerate(self.scans)])
        try:
            intensity = np.concatenate([d[:, 1] * it[i] / 1000. for i, d in enumerate(raw_dat)])
        except Exception as e:
            print("Mark1:", e, it)
            intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])
        try:
            it = np.concatenate([it * np.ones(len(raw_dat[i])) for i, it in enumerate(it)])
        except Exception as e:
            print("Error with injection time correction:", e)

        data_array = np.transpose([mz, intensity, scans, it])
        return data_array

    def close(self):
        self.msrun.Close()
        return


if __name__ == "__main__":
    # import matplotlib.pyplot as plt
    test = "C:\\Python\\UniDec3\\TestSpectra\\test.raw"
    d = ThermoImporter(test, silent=False)
    # plt.plot(data[:, 0], data[:, 1])
    # plt.show()
    exit()
    cdms_dat = importer.get_cdms_data()
    for i in cdms_dat:
        print(i)

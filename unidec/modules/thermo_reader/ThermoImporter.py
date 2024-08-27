import numpy as np
# import modules.thermo_reader.RawFileReader
import time


class ThermoDataImporter:
    """
    Imports Thermo data files.
    """

    def __init__(self, path, silent=False, *args, **kwargs):
        """
        Imports Thermo file.
        :param path: Thermo .Raw file path
        :param args: arguments (unused)
        :param kwargs: keywords (unused)
        :return: mzMLimporter object
        """
        if not silent:
            print("Launching Thermo Importer. If it fails after this step, try this:")
            print("Delete your whole UniDec folder but keep the zip file.")
            print("Right click on the zip file and open properties. You should see a box to Unblock it. Check that.")
            print("Click ok. Unzip it again. Try it once more.")
        from unidec.modules.thermo_reader.RawFileReader import RawFileReader as rr
        print("Reading Thermo Data:", path)
        self.msrun = rr(path)
        self.scanrange = self.msrun.scan_range()
        # print(self.scanrange)
        self.scans = np.arange(self.scanrange[0], self.scanrange[1]+1)
        self.times = []
        self.data = None
        for s in self.scans:
            self.times.append(self.msrun.scan_time_from_scan_name(s))
        self.times = np.array(self.times)
        print("Number of Scans", len(self.scans))
        # print(self.times)

    def grab_data(self, threshold=-1):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetSpectrum(s))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            if threshold >= 0:
                impdat = impdat[impdat[:, 1] > threshold]
            self.data.append(impdat)
        self.data = np.array(self.data, dtype=object)
        return self.data

    def grab_centroid_data(self):
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetCentroidArray(s))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
        self.data = np.array(self.data, dtype=object)
        return self.data

    def grab_scan_data(self, s):
        impdat = np.array(self.msrun.GetSpectrum(s))  # May want to test this.
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def grab_centroid_data(self, s):
        impdat = np.array(self.msrun.GetCentroidSpectrum(s))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def get_ms_order(self, s):
        order = self.msrun.GetMSOrder(s)
        return order

    def get_scan_time(self, s):
        return self.msrun.scan_time_from_scan_name(s)

    def get_isolation_mz_width(self, s):
        scanFilter =  self.msrun.GetScanFilter(s)
        reaction = scanFilter.GetReaction(0)
        mz = reaction.PrecursorMass
        width = reaction.IsolationWidth
        return mz, width

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if scan_range is not None:
            scan_range = np.array(scan_range, dtype=int)
            scan_range = scan_range + 1
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)
        if scan_range is None:
            try:
                scan_range = [np.amin(self.scans), np.amax(self.scans)]
            except:
                scan_range = [1, 2]
        scan_range = np.array(scan_range, dtype=int)
        print("Thermo Scan Range:", scan_range)

        try:
            if scan_range[0] < np.amin(self.scans) or scan_range[0] == -1:
                scan_range[0] = np.amin(self.scans)
            if scan_range[1] > np.amax(self.scans) or scan_range[1] == -1:
                scan_range[1] = np.amax(self.scans)
        except:
            scan_range = [1, 2]

        if scan_range[1] - scan_range[0] > 1:
            data = np.array(list(self.msrun.GetAverageSpectrum(scan_range)))
        else:
            impdat = np.array(self.msrun.GetSpectrum(scan_range[0]))  # May want to test this.
            impdat = impdat[impdat[:, 0] > 10]
            data = impdat

        return data

    def get_tic(self):
        return self.msrun.GetChromatogram()

    def get_eic(self, mass_range=None, scan_range=None):
        if mass_range is None:
            mass_range = [0, 1000000]
        return self.msrun.Get_EIC(massrange=mass_range, scanrange=scan_range)

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

    def get_sid_voltage(self, scan=1):
        try:
            scan_mode = self.msrun.source.GetScanEventStringForScanNumber(scan)
            #Find where sid= is in the string
            sid_index = scan_mode.find("sid=")
            #Find the next space after sid=
            space_index = scan_mode.find(" ", sid_index)
            #Get the sid value
            sid_value = scan_mode[sid_index+4:space_index]
            #Convert to float
            sid_value = float(sid_value)
        except:
            sid_value = 0

        return sid_value



    def get_polarity(self, scan=1):
        #print(dir(self.msrun.source))
        """
        im = self.msrun.source.GetInstrumentMethod(0)
        print(im)
        for line in im.split("\n"):
            if "Polarity" in line:
                if "Positive" in line:
                    print("Polarity: Positive")
                    return "Positive"
                if "Negative" in line:
                    print("Polarity: Negative")
                    return "Negative"
        print("Polarity: Unknown")
        return None
        # exit()"""
        scan_mode = self.msrun.source.GetScanEventStringForScanNumber(scan)
        if "+" in scan_mode:
            print("Polarity: Positive")
            return "Positive"
        if "-" in scan_mode[:10]:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")

        return None

if __name__ == "__main__":
    test = u"C:\\Python\\UniDec3\\TestSpectra\\test.RAW"
    #test = "Z:\\Group Share\\Levi\\MS DATA\\vt_ESI data\\DMPG LL37 ramps\\18to1\\20210707_LB_DMPG3_LL37_18to1_RAMP_16_37_3.RAW"
    #test = "Z:\Group Share\Group\Archive\JamesKeener Keener\AqpZ mix lipid ND\\20190226_JEK_AQPZ_E3T0_PGPC_GC_NEG.RAW"

    tstart = time.perf_counter()
    d = ThermoDataImporter(test)
    dat = d.grab_centroid_data(1)
    print(len(dat))

    import matplotlib.pyplot as plt
    plt.plot(dat[:,0], dat[:,1])
    plt.show()

    exit()

    d.get_polarity()
    exit()
    vdata = d.get_analog_voltage1()
    times = d.get_tic()[1:, 0]
    data = d.get_data()
    # vdata = (-34.48*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata)))+(263.91*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata)))-(811.83*(vdata-np.amin(vdata))*(vdata-np.amin(vdata))*(vdata-np.amin(vdata)))+(1258.4*(vdata-np.amin(vdata))*(vdata-np.amin(vdata)))-(1032.3*(vdata-np.amin(vdata)))+409.12

    vdata = (-44.115 * vdata * vdata * vdata) + (201.67 * vdata * vdata) + (-347.15 * vdata) + 242.19

    import matplotlib.pyplot as plt

    plt.plot(times, vdata)
    plt.show()

    # d.get_data(time_range=(0, 10))
    print("ImportData: %.2gs" % (time.perf_counter() - tstart))
    # import matplotlib.pyplot as plt
    # plt.plot(d[:, 0], d[:, 1])
    # plt.show()

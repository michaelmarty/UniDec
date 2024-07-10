from unidec.modules.mzMLimporter import *
from pyteomics import mzxml

__author__ = 'Michael.Marty'

def get_data_from_spectrum(spectrum, threshold=-1):
    impdat = np.transpose([spectrum['m/z array'], spectrum['intensity array']])
    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 1] > threshold]
    return impdat



class mzXMLimporter:
    """
    Imports mzXML data files.
    """

    def __init__(self, path, *args, **kwargs):
        """
        Imports mzXML file, adds the chromatogram into a single spectrum.
        :param path: .mzXML file path
        :param args: arguments (unused)
        :param kwargs: keywords (unused)
        :return: mzXMLimporter object
        """
        print("Reading mzXML:", path)
        self.path = path
        self.filesize = os.stat(path).st_size
        self.msrun = mzxml.read(path)
        self.data = None
        self.polarity = None
        # self.scans = []
        self.times = []
        self.ids = []
        for i, spectrum in enumerate(self.msrun):
            try:
                if spectrum["msLevel"] is None:
                    continue
            except:
                pass
            try:
                t = spectrum["retentionTime"]
                id = spectrum["num"]
                self.times.append(float(t))
            except Exception as e:
                self.times.append(-1)
                id = -1
                # print("1", spectrum, e)
            # self.scans.append(i)
            self.ids.append(id)
            # print(i, end=" ")
        self.times = np.array(self.times)
        self.ids = np.array(self.ids)
        self.scans = np.arange(0, len(self.ids))
        # print("Reading Complete")

    def grab_scan_data(self, scan):
        data = get_data_from_spectrum(self.msrun[self.ids[scan]])
        return data

    def get_data_memory_safe(self, scan_range=None, time_range=None):
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)
        if scan_range is None:
            scan_range = [np.amin(self.scans), np.amax(self.scans)]
        print("Scan Range:", scan_range)
        data = get_data_from_spectrum(self.msrun[self.ids[0]])

        resolution = get_resolution(data)
        axis = ud.nonlinear_axis(np.amin(data[:, 0]), np.amax(data[:, 0]), resolution)
        template = np.transpose([axis, np.zeros_like(axis)])
        # print("Length merge axis:", len(template))
        newdat = ud.mergedata(template, data)
        template[:, 1] += newdat[:, 1]

        for i in range(int(scan_range[0]) + 1, scan_range[1] + 1):
            try:
                data = get_data_from_spectrum(self.msrun[self.ids[i]])
                newdat = ud.mergedata(template, data)
                template[:, 1] += newdat[:, 1]
            except Exception as e:
                print("Error", e, "With scan number:", i)
        return template

    def grab_data(self, threshold=-1):
        newtimes = []
        # newscans = []
        newids = []
        self.data = []
        for i, s in enumerate(self.ids):
            try:
                impdat = get_data_from_spectrum(self.msrun[s], threshold=threshold)
                self.data.append(impdat)
                newtimes.append(self.times[i])
                # newscans.append(self.scans[i])
                newids.append(s)
            except Exception as e:
                print("mzML import error")
                print(e)
        # self.scans = np.array(newscans)
        self.times = np.array(newtimes)
        self.ids = np.array(newids)
        self.scans = np.arange(0, len(self.ids))
        return self.data

    def get_data_fast_memory_heavy(self, scan_range=None, time_range=None):
        if self.data is None:
            self.grab_data()

        data = deepcopy(self.data)
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        if scan_range is not None:
            data = data[int(scan_range[0]):int(scan_range[1] + 1)]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans), len(data))

        if data is None or ud.isempty(data):
            print("Error: Empty Data Object")
            return None

        if len(data) > 1:
            try:
                data = merge_spectra(data)
            except Exception as e:
                concat = np.concatenate(data)
                sort = concat[concat[:, 0].argsort()]
                data = ud.removeduplicates(sort)
                print("2", e)
        elif len(data) == 1:
            data = data[0]
        else:
            data = data
        # plt.figure()
        # plt.plot(data)
        # plt.show()
        return data

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if self.filesize > 1000000000 and self.data is None:
            try:
                data = self.get_data_memory_safe(scan_range, time_range)
            except Exception as e:
                print("Error in Memory Safe mzML, trying memory heavy method")
                data = self.get_data_fast_memory_heavy(scan_range, time_range)
        else:
            data = self.get_data_fast_memory_heavy(scan_range, time_range)
        return data

    def get_tic(self):
        tic = []
        print("Imported Data. Constructing TIC")
        for i, s in enumerate(self.ids):
            spectrum = self.msrun[s]
            try:
                tot = float(spectrum['totIonCurrent'])
            except:
                tot = 0
            tic.append(tot)
        t = self.times
        ticdat = np.transpose([t, tic])
        return ticdat

    def get_scans_from_times(self, time_range):
        boo1 = self.times >= time_range[0]
        boo2 = self.times < time_range[1]
        try:
            min = np.amin(self.scans[boo1])
            max = np.amax(self.scans[boo2])
        except Exception as e:
            min = -1
            max = -1
            print("3", e)
        return [min, max]

    def get_times_from_scans(self, scan_range):
        boo1 = self.scans >= scan_range[0]
        boo2 = self.scans < scan_range[1]
        boo3 = np.logical_and(boo1, boo2)
        min = np.amin(self.times[boo1])
        max = np.amax(self.times[boo2])
        try:
            avg = np.mean(self.times[boo3])
        except Exception as e:
            avg = min
            print("4", e)
        return [min, avg, max]

    def get_max_time(self):
        return np.amax(self.times)

    def get_max_scans(self):
        return np.amax(self.scans)

    def get_polarity(self):
        return self.polarity


if __name__ == "__main__":
    test = u"C:\\Data\\Wilson_Genentech\\Raw Files for Michael Marty\\B4.mzXML"
    test = "D:\Data\TopDown\denatured_discovery_input_topfd.mzxml"
    import time

    tstart = time.perf_counter()

    d = mzXMLimporter(test)
    d.grab_scan_data(10)
    exit()

    tic = d.get_tic()
    #print(len(tic))
    import matplotlib.pyplot as plt

    plt.plot(tic[:, 0], tic[:, 1])
    plt.show()

    exit()
    print(len(d.scans))

    #data = d.get_data_memory_safe()
    data = d.get_data()
    tend = time.perf_counter()
    # print(call, out)
    print("Execution Time:", (tend - tstart))

    exit()
    import matplotlib.pyplot as plt

    plt.plot(data[:, 0], data[:, 1])
    plt.show()

    print(len(data))
    exit()
    # get_data_from_spectrum(d.msrun[239])
    # exit()
    data = d.get_data()

    print(data)


    # print d.get_times_from_scans([15, 30])

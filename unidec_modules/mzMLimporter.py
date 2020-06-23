import numpy as np
import pymzml
from unidec_modules import unidectools as ud
import os
from copy import deepcopy

__author__ = 'Michael.Marty'


def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = ud.safedivide(diffs[:, 0], diffs[:, 1])
    # popt, pcov = scipy.optimize.curve_fit(fit_line, diffs[:, 0], resolutions, maxfev=1000000)
    # fitLine = fit_line(diffs[:, 0], *popt)
    # fitMax = np.max(fitLine)
    # fitMin = np.min(fitLine)
    # diffs_new = diffs[(1.2 * fitMin < resolutions) & (resolutions < 1.2 * fitMax)]
    # resolutions_new = resolutions[(1.2 * fitMin < resolutions) & (resolutions < 1.2 * fitMax)]
    # popt2, pcov2 = scipy.optimize.curve_fit(fit_line, diffs_new[:, 0], resolutions_new, maxfev=1000000)
    # plt.figure()
    # plt.plot(diffs[:,0], resolutions)
    # plt.plot(diffs[:, 0], fit_line(diffs[:, 0], *popt2), 'r-')
    # plt.show()
    # Currently use A * m ^1.5 (0.5?)
    # Maybe use a*M^b
    # return popt2
    return np.median(resolutions)


def fit_line(x, a, b):
    return a * x ** b


def merge_spectra(datalist, mzbins=None):
    """
    Merge together a list of data.
    Interpolates each data set in the lit to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    # print "Median Resolution:", resolution
    # axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    if mzbins is None or mzbins == 0:
        resolution = get_resolution(datalist[0])
        axis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        axis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))
    template = np.transpose([axis, np.zeros_like(axis)])
    print("Length merge axis:", len(template))
    for d in datalist:
        if len(d) > 1:
            newdat = ud.mergedata(template, d)
            # newdat=ud.lintegrate(d,axis)
            template[:, 1] += newdat[:, 1]
    return template


def nonlinear_axis(start, end, res):
    """
    Creates a nonlinear axis with the m/z values spaced with a defined and constant resolution.
    :param start: Minimum m/z value
    :param end: Maximum m/z value
    :param res: Resolution of the axis ( m / delta m)
    :return: One dimensional array of the nonlinear axis.
    """
    axis = []
    i = start
    axis.append(i)
    i += i / fit_line(i, res[0], res[1])
    while i < end:
        axis.append(i)
        i += i / fit_line(i, res[0], res[1])
    return np.array(axis)


def get_data_from_spectrum(spectrum):
    impdat = np.transpose([spectrum.mz, spectrum.i])
    impdat = impdat[impdat[:, 0] > 10]
    return impdat


class mzMLimporter:
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
        print("Reading mzML:", path)
        self.filesize = os.stat(path).st_size
        self.msrun = pymzml.run.Reader(path)
        self.data=None
        self.scans = []
        self.times = []
        for i, spectrum in enumerate(self.msrun):
            if '_scan_time' in list(spectrum.__dict__.keys()):
                try:
                    t = spectrum.scan_time_in_minutes()
                    self.times.append(float(t))
                except Exception as e:
                    self.times.append(-1)
                    print("1", spectrum, e)
                self.scans.append(i)
                # print(i, end=" ")
            else:
                print("Scan time not found", i)
        self.times = np.array(self.times)
        self.scans = np.array(self.scans)
        # print("Reading Complete")

    def get_data_memory_safe(self, scan_range=None, time_range=None):
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)
        if scan_range is None:
            scan_range = [np.amin(self.scans), np.amax(self.scans)]
        print("Scan Range:", scan_range)
        data = get_data_from_spectrum(self.msrun[int(scan_range[0] + 1)])

        resolution = get_resolution(data)
        axis = ud.nonlinear_axis(np.amin(data[:, 0]), np.amax(data[:, 0]), resolution)
        template = np.transpose([axis, np.zeros_like(axis)])
        # print("Length merge axis:", len(template))
        newdat = ud.mergedata(template, data)
        template[:, 1] += newdat[:, 1]

        for i in range(int(scan_range[0]) + 1, scan_range[1] + 1):
            print("Importing Scan:", i)
            data = get_data_from_spectrum(self.msrun[i + 1])
            newdat = ud.mergedata(template, data)
            template[:, 1] += newdat[:, 1]
        return template

    def grab_data(self):
        self.data = []
        for i, spectrum in enumerate(self.msrun):
            if '_scan_time' in list(spectrum.__dict__.keys()):
                impdat = np.transpose([spectrum.mz, spectrum.i])
                impdat = impdat[impdat[:, 0] > 10]
                self.data.append(impdat)
        self.data = np.array(self.data)

    def get_data_fast_memory_heavy(self, scan_range=None, time_range=None):
        if self.data is None:
            self.grab_data()

        data = deepcopy(self.data)
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        if scan_range is not None:
            data = data[int(scan_range[0]):int(scan_range[1])]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans), data.shape)

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
        if self.filesize>1000000000:
            try:
                data=self.get_data_memory_safe(scan_range, time_range)
            except Exception as e:
                print("Error in Memory Safe mzML, trying memory heavy metho")
                data=self.get_data_fast_memory_heavy(scan_range, time_range)
        else:
            data=self.get_data_fast_memory_heavy(scan_range, time_range)
        return data

    def get_tic(self):
        return self.msrun["TIC"]

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


if __name__ == "__main__":
    test = u"C:\Python\\UniDec\TestSpectra\JAW.mzML"
    d = mzMLimporter(test).get_data()
    # print d.get_times_from_scans([15, 30])

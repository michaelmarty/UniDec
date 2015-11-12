import os
import time

import numpy as np
import pymzml
import unidectools as ud

__author__ = 'Michael.Marty'


def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = diffs[:, 0] / diffs[:, 1]
    # plt.figure()
    # plt.plot(diffs[:,0],resolutions)
    # plt.show()
    return np.median(resolutions)


def merge_spectra(datalist):
    """
    Merge together a list of data.
    Interpolates each data set in the lit to a new nonlinear axis with the median resolution of the first element.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    resolution = get_resolution(datalist[0])
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    print "Median Resolution:", resolution
    axis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    template = np.transpose([axis, np.zeros_like(axis)])
    print "Length merge axis:", len(template)

    for d in datalist:
        if len(d) > 1:
            newdat = ud.mergedata(template, d)
            # newdat=ud.lintegrate(d,axis)
            template[:, 1] += newdat[:, 1]
    return template


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
        print "Importing mzML"
        self.msrun = pymzml.run.Reader(path)
        i = 0
        self.data = []
        for spectrum in self.msrun:
            impdat = np.transpose([spectrum.mz, spectrum.i])
            impdat = impdat[impdat[:, 0] > 10]
            self.data.append(impdat)
            i += 1
        print "Number of scans:", i
        self.data = np.array(self.data)
        if len(self.data) > 1:
            try:
                self.data = merge_spectra(self.data)
            except Exception, e:
                concat = np.concatenate(self.data)
                sort = concat[concat[:, 0].argsort()]
                self.data = ud.removeduplicates(sort)
                print e

    def get_data(self):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        return self.data


if __name__ == "__main__":
    dir = "C:\\cprog\\Joe"
    file = "20150304_myo_70r_5u_full-ALL.mzML"
    dir = "C:\\Users\\michael.marty\\Hugh\Data\\290415\\mzML"
    file = "20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_2_35000.mzML"
    file = "20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_1_35000.mzML"
    # file="20150429_ND_POPG_PG02_RF50_ETHMR_0_HCD_24_35000.mzML"
    tstart = time.time()
    data = mzMLimporter(os.path.join(dir, file)).get_data()
    tend = time.time()
    print data
    print "Done %.2gs" % float(tend - tstart)

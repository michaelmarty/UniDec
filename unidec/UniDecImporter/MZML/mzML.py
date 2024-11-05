import gc

import numpy as np
import pymzml
from unidec import tools as ud
import os
from copy import deepcopy
from pymzml.utils.utils import index_gzip
import pymzml.obo
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
__author__ = 'Michael.Marty'

from unidec.UniDecImporter.Importer import Importer


def gzip_files(mzml_path, out_path):
    """
    Create and indexed gzip mzML file from a plain mzML.
    """
    with open(mzml_path) as fin:
        fin.seek(0, 2)
        max_offset_len = fin.tell()
        max_spec_no = pymzml.run.Reader(mzml_path).get_spectrum_count() + 32

    index_gzip(
        mzml_path, out_path, max_idx=max_spec_no, idx_len=len(str(max_offset_len))
    )


def auto_gzip(mzml_path):
    out_path = mzml_path + ".gz"
    gzip_files(mzml_path, out_path)
    return out_path

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


def get_longest_index(datalist):
    lengths = [len(x) for x in datalist]
    return np.argmax(lengths)


def merge_spectra(datalist, mzbins=None, type="Interpolate"):
    """
    Merge together a list of data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Filter out junk spectra that are empty
    datalist = [x for x in datalist if len(x) > 0]

    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    # print "Median Resolution:", resolution
    # axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    #for d in datalist:
    #    print(d)
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.
    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution(datalist[maxlenpos])
        if resolution < 0:
            print("ERROR with auto resolution:", resolution, maxlenpos, datalist[maxlenpos])
            print("Using ABS")
            resolution = np.abs(resolution)
        elif resolution == 0:
            print("ERROR, resolution is 0, using 20000.", maxlenpos, datalist[maxlenpos])
            resolution = 20000

        axis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        axis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))
    template = np.transpose([axis, np.zeros_like(axis)])

    print("Length merge axis:", len(template))

    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = ud.mergedata(template, d)
            elif type == "Integrate":
                newdat = ud.lintegrate(d, axis)
            else:
                print("ERROR: unrecognized trdtrmerge spectra type:", type)

            template[:, 1] += newdat[:, 1]

    # Trying to catch if the data is backwards
    try:
        if template[1, 0] < template[0, 0]:
            template = template[::-1]
    except:
        pass
    return template


def merge_im_spectra(datalist, mzbins=None, type="Integrate"):
    """
    Merge together a list of ion mobility data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    if len(datalist) > 1:
        concat = np.concatenate(datalist)
    else:
        concat = np.array(datalist[0])
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.
    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution(datalist[maxlenpos])
        mzaxis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        mzaxis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))

    # For drift time, use just unique drift time values. May need to make this fancier.
    dtaxis = np.sort(np.unique(concat[:, 1]))

    # Create the mesh grid from the new axes
    X, Y = np.meshgrid(mzaxis, dtaxis, indexing="ij")

    template = np.transpose([np.ravel(X), np.ravel(Y), np.ravel(np.zeros_like(X))])
    print("Shape merge axis:", X.shape)
    xbins = deepcopy(mzaxis)
    xbins[1:] -= np.diff(xbins)
    xbins = np.append(xbins, xbins[-1] + np.diff(xbins)[-1])
    ybins = deepcopy(dtaxis)
    ybins[1:] -= np.diff(ybins) / 2.
    ybins = np.append(ybins, ybins[-1] + np.diff(ybins)[-1])
    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = ud.mergedata2d(template[:, 0], template[:, 1], d[:, 0], d[:, 1], d[:, 2])
            elif type == "Integrate":
                newdat, xedges, yedges = np.histogram2d(d[:, 0], d[:, 1], bins=[xbins, ybins], weights=d[:, 2])
            else:
                print("ERROR: unrecognized merge spectra type:", type)
            template[:, 2] += np.ravel(newdat)
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


def get_data_from_spectrum(spectrum, threshold=-1):
    impdat = np.transpose([spectrum.mz, spectrum.i])
    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 1] > threshold]
    # print(len(impdat))
    # with np.printoptions(threshold=10000):
    #    print(np.sum(impdat[:,1]>0))
    # print(impdat)
    # exit()
    return impdat


def get_im_data_from_spectrum(spectrum, threshold=-1):
    array_params = spectrum._get_encoding_parameters("raw ion mobility array")
    dtarray = spectrum._decode(*array_params)

    impdat = np.transpose([spectrum.mz, dtarray, spectrum.i])

    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 2] > threshold]
    return impdat


def search_by_id(obo, id):
    key = "MS:{0}".format(id)
    return_value = ""
    for lookup in obo.lookups:
        if key in lookup:
            if obo.MS_tag_regex.match(key):
                for fn in FIELDNAMES:
                    if fn in lookup[key].keys():
                        return_value += "{0}\n".format(lookup[key][fn])
    return return_value


FIELDNAMES = ["id", "name", "def", "is_a"]


import os
import numpy as np
from pyteomics import mzml

import os
import numpy as np

class MZMLImporter(Importer):
    """
    Imports mzML data files.
    """
    def __init__(self, path, *args, **kwargs):
        self.filesize = os.stat(path).st_size
        self.times = []
        self.ids = []
        self.scans = []
        self.msrun = pymzml.run.Reader(path)
        self.data = None
        self.global_counter = 0
        self.process_scan()


    def process_scan(self):
        for i, spectrum in enumerate(self.msrun):
            if '_scan_time' in list(spectrum.__dict__.keys()):
                try:
                    if spectrum.ms_level is None:
                        continue
                    t = spectrum.scan_time_in_minutes()
                    id = spectrum.ID
                    self.times.append(float(t))
                    self.ids.append(id)
                except Exception as e:
                    continue
                    # self.times.append(-1)
                    # self.ids.append(-1)
            else:
                print("Scan time not found for spectrum ID:", spectrum.ID)
        self.times = np.array(self.times)
        self.ids = np.array(self.ids)
        self.scans = np.arange(len(self.ids))


    def grab_scan_data(self, scan):
        try:
            data = get_data_from_spectrum(self.msrun[self.ids[scan]])
        except Exception as e:
            print("Error in grab_scan_data:", e)
            data = None
        return data

    def get_data_memory_safe(self, scan_range=None, time_range=None):
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)
        if scan_range is None:
            scan_range = [int(np.amin(self.scans)), int(np.amax(self.scans))]
        print("Scan Range:", scan_range)
        data = get_data_from_spectrum(self.msrun[self.ids[scan_range[0]]])

        resolution = get_resolution(data)
        axis = ud.nonlinear_axis(np.amin(data[:, 0]), np.amax(data[:, 0]), resolution)
        template = np.transpose([axis, np.zeros_like(axis)])
        newdat = ud.mergedata(template, data)
        template[:, 1] += newdat[:, 1]
        if self.filesize < 1e6:
            index = 0
            while index <= scan_range[1] - scan_range[0]:

                try:
                    spec = self.msrun.next()
                except:
                    break
                if spec.ID in self.ids:
                    index += 1
                    if scan_range[0] <= index <= scan_range[1]:
                        # try:
                        data = get_data_from_spectrum(spec)
                        newdat = ud.mergedata(template, data)
                        template[:, 1] += newdat[:, 1]
                        # except Exception as e:
                        #     print("Error", e, "With scan number:", index)
            return template
        else:
            thread_count = 8
            index_range = list(range(scan_range[0], scan_range[1]))
            batches = self.generate_batches(index_range, thread_count)
            for batch in batches:
                print(batch)
            with ThreadPoolExecutor(max_workers=thread_count) as executor:
                results = list(executor.map(lambda batch: self.process_batch_and_merge(batch, template), batches))

                for local_template in results:
                    template[:, 1] += local_template[:, 1]

            return template

    def process_batch_and_merge(self, batch, template):
        local_template = template.copy()
        for index in batch:
                if index < len(self.ids):
                    try:

                        spec = self.msrun[self.ids[index]]
                        data = get_data_from_spectrum(spec)
                        newdat = ud.mergedata(local_template, data)
                        local_template[:, 1] += newdat[:, 1]
                    except Exception as e:
                        print("Error", e, "With scan number:", index)
        return local_template

    def generate_batches(self, scan_range, thread_count):
        total_length = scan_range[-1]-scan_range[0]+1
        batch_size = total_length // thread_count
        batches = []


        for i in range(thread_count):
            start_index = i * batch_size
            if i == thread_count - 1:  # Last batch gets any remaining items
                batches.append(list(range(start_index, scan_range[1])))
            else:
                end_index = start_index + batch_size
                batches.append(scan_range[start_index:end_index])

        return batches

    def grab_data(self):
        for spectrum in self.msrun:
            if 'scanList' in spectrum:
                try:
                    scan_info = spectrum['scanList']['scan'][0]
                    t = scan_info['scan start time']
                    id = spectrum['id']

                    raw_data = list(spectrum.get('m/z array', []))
                    intensity_data = list(spectrum.get('intensity array', []))

                    if raw_data and intensity_data:
                        self.data.append((raw_data, intensity_data))
                        print(
                            f"Data appended for ID {id}: {raw_data[:5]}, {intensity_data[:5]}")
                    else:
                        print(f"Missing data for spectrum ID {id}")

                except Exception as e:
                    continue
        #             print("Error processing spectrum:", e)
        #     else:
        #         print("Scan list not found for spectrum ID:", spectrum['id'])
        #
        # print("Final data length:", len(self.data))

    def get_data_fast_memory_heavy(self, scan_range=None, time_range=None):

        if self.data is None:
            self.grab_data()

        data = deepcopy(self.data)
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        if scan_range is not None:
            data = data[int(scan_range[0]):int(scan_range[1] + 1)]
            print("Getting scans1:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans), data.shape)

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
        return data

    def get_data(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if self.filesize > 1e9 and self.data is None:
            # try:

            data = self.get_data_memory_safe(scan_range, time_range)
            # except Exception as e:
            #     print("Error in Memory Safe mzML, trying memory heavy method")
            #     data = self.get_data_fast_memory_heavy(scan_range, time_range)
        else:
            data = self.get_data_memory_safe(scan_range, time_range)
        return data

    def get_tic(self):
        try:
            tic = self.msrun["TIC"]
            ticdat = np.transpose([tic.time, tic.i])
            if len(ticdat) != len(self.scans):
                print("TIC too long. Likely extra non-MS scans", len(ticdat), len(self.scans))
                raise Exception
        except:
            print("Error getting TIC in mzML; trying to make it...")

            tic = []
            self.grab_data()
            print("Imported Data. Constructing TIC")
            for d in self.data:
                try:
                    tot = np.sum(d[:, 1])
                except:
                    tot = 0
                tic.append(tot)
            t = self.times
            ticdat = np.transpose([t, tic])
            print("Done")
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

    def get_inj_time(self, spectrum):
        element = spectrum.element
        it = 1
        for child in element.iter():
            if 'name' in child.attrib:
                if child.attrib['name'] == 'ion injection time':
                    it = child.attrib['value']
                    try:
                        it = float(it)
                    except:
                        it = 1
        return it

    def get_property(self, s, name):
        element = self.msrun[s].element
        it = 1
        for child in element.iter():
            if 'name' in child.attrib:
                if child.attrib['name'] == name:
                    it = child.attrib['value']
                    try:
                        it = float(it)
                    except:
                        it = 1
        return it

    def get_dts(self):
        dts = []
        for i, s in enumerate(self.ids):
            try:
                dt = self.get_property(s, 'ion mobility drift time')
                dts.append(dt)
            except:
                dts.append(-1)
        return np.array(dts)

    def get_inj_time_array(self):
        its = []
        for i, s in enumerate(self.ids):
            it = self.get_inj_time(self.msrun[s])
            try:
                it = float(it)
            except:
                print("Error in scan header:", i, s, it)
                it = 1
            its.append(it)
        return np.array(its)

    def grab_im_data(self):
        newtimes = []
        newids = []
        self.data = []
        """
        # Old Slow Method
        for i, s in enumerate(self.ids):
            try:
                array = get_im_data_from_spectrum(self.msrun[s])
                self.data.append(np.array(array))
                newtimes.append(self.times[i])
                newids.append(s)
            except:
                pass"""

        # New Faster Method
        for n, spec in enumerate(self.msrun):
            if spec.ID in self.ids:
                try:
                    impdat = get_im_data_from_spectrum(spec)
                    self.data.append(impdat)
                    newtimes.append(self.times[n])
                    newids.append(spec.ID)
                except Exception as e:
                    print("mzML import error")
                    print(e)

        self.data = np.array(self.data, dtype='object')
        self.times = np.array(newtimes)
        self.ids = np.array(newids)
        self.scans = np.arange(0, len(self.ids))
        return self.data

    def get_im_data(self, scan_range=None, time_range=None, mzbins=None):
        start_time = time.perf_counter()
        if self.data is None:
            self.grab_im_data()

        data = deepcopy(self.data)
        if time_range is not None:
            scan_range = self.get_scans_from_times(time_range)
            print("Getting times:", time_range)

        if scan_range is not None:
            data = data[int(scan_range[0]):int(scan_range[1] + 1)]
            print("Getting scans:", scan_range)
        else:
            print("Getting all scans, length:", len(self.scans), data.shape)

        if data is None or ud.isempty(data):
            print("Error: Empty Data Object")
            return None

        # Need to merge to get 2D from sparse array
        try:
            data = merge_im_spectra(data, mzbins=mzbins)
        except Exception as e:
            concat = np.concatenate(data)
            sort = concat[concat[:, 0].argsort()]
            data = ud.removeduplicates(sort)
            print("2", e)

        # plt.figure()
        # plt.plot(data)
        # plt.show()
        print("Import Time:", time.perf_counter() - start_time)
        return data

    def get_polarity(self, scan=1):
        for s, spec in enumerate(self.msrun):
            if s == scan:
                # spec = self.msrun[scan]
                negative_polarity = spec["negative scan"]
                if negative_polarity == "" or negative_polarity:
                    negative_polarity = True
                    print("Polarity: Negative")
                    return "Negative"

                positive_polarity = spec["positive scan"]
                if positive_polarity == "" or positive_polarity:
                    print("Polarity: Positive")
                    return "Positive"

                print(positive_polarity, negative_polarity)
                print("Polarity: Unknown")
                return None

    def get_ms_order(self, s):
        order = self.msrun[s-1].ms_level
        #s is the scan number, this is 0 indexed, so we subtract 1 to access the correct scan.
        return order

    def get_scan_time(self, s):
        scantime = self.msrun[s-1].scan_time_in_minutes()
        return scantime

    def get_isolation_mz_width(self, s):
        scan = self.msrun[s-1]
        precursors = scan.selected_precursors
        if len(precursors) == 0:
            return None, None
        else:
            precursor = precursors[0]
            mz = precursor['mz']
            width = 3
            return mz, width



if __name__ == "__main__":
    test = u"C:\\Python\\UniDec3\\TestSpectra\\JAW.mzML"
    # test = "C:\Data\CytC_Intact_MMarty_Share\\221114_STD_Pro_CytC_2ug_r1.mzML.gz"
    # test = "C:\Data\CytC_Intact_MMarty_Share\\221114_STD_Pro_CytC_2ug_r1_2.mzML"
    # test = "C:\Data\IMS Example Data\imstest2.mzML"
    test = u"C:\\Data\\CytC_Intact_MMarty_Share\\20221215_MMarty_Share\\SHA_1598_9.mzML.gz"
    #test = "C:\\Users\\marty\\OneDrive - University of Arizona\\Attachments\\S203.mzML.gz"
    test = "C:\\Users\\marty\\Downloads\\BSA1_scan214_full.mzML"
    import time

    tstart = time.perf_counter()

    d = MZMLImporter(test)
    data = d.get_data()
    print(len(data))
    tend = time.perf_counter()
    # print(call, out)
    print("Execution Time:", (tend - tstart))
    exit()
    '''
    spectrum = d.msrun[10]
    # it = it.get("ion inject time")
    element = spectrum.element
    for child in element.iter():
        if 'name' in child.attrib:
            if child.attrib['name'] == 'ion injection time':
                it = child.attrib['value']
                print(it)
    exit()
    tic = d.get_tic()
    print(len(tic))
    print(len(d.scans))

    exit()'''
    # data = d.get_data_memory_safe()
    data = d.get_data(time_range=(3, 5))
    tend = time.perf_counter()
    # print(call, out)
    print("Execution Time:", (tend - tstart))

    print(len(data))
    # exit()
    # get_data_from_spectrum(d.msrun[239])
    # exit()
    # data = d.get_data()

    print(data)
    import matplotlib.pyplot as plt

    plt.plot(data[:, 0], data[:, 1])
    plt.show()

    # print d.get_times_from_scans([15, 30])

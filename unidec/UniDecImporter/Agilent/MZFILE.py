# Rip of multiplierz's mzFile
import os
import sys
import clr
from numpy import average
from itertools import chain

from unidec.UniDecImporter.Thermo.RawFileReader import pathtothisfile

'''

pathtothisfile = os.path.dirname(__file__)
# print(pathtothisfile)
dlls = ['ThermoFisher.CommonCore.Data', 'ThermoFisher.CommonCore.RawFileReader',
        'ThermoFisher.CommonCore.BackgroundSubtraction', 'ThermoFisher.CommonCore.MassPrecisionEstimator']

print(pathtothisfile)

for dll in dlls:
    testpath = os.path.join(pathtothisfile, dll) + ".dll"
    # print(testpath)
    if os.path.isfile(testpath):
        # print("1")
        clr.AddReference(testpath)
    else:
        try:
            # print("2")
            import sys

            sys.path.append(pathtothisfile)
            clr.AddReference(dll)
        except:
            # print("3")
            clr.AddReference(dll)'''

#dll_path = "C:\\Python\\UniDec3\\unidec\\UniDecImporter\\Agilent"
pathtothisfile = os.path.dirname(__file__)
print(pathtothisfile)

dlls = ['MassSpecDataReader', 'BaseCommon', 'BaseDataAccess']
for dll in dlls:
    testpath = os.path.join(pathtothisfile, dll) + ".dll"
    if os.path.isfile(testpath):
        clr.AddReference(testpath)
    else:
        try:
            sys.path.append(pathtothisfile)
            clr.addReference(dll)
        except:
            clr.AddReference(dll)
clr.AddReference("System.Collections")

import Agilent
from Agilent.MassSpectrometry.DataAnalysis import (MassSpecDataReader, BDAChromFilter, MsdrPeakFilter,
                                                   MsdrChargeStateAssignmentFilter, IMsdrDataReader, IBDAChromFilter,
                                                   IMsdrChargeStateAssignmentFilter,
                                                   DesiredMSStorageType, ChromType, MinMaxRange, MSLevel)

for item in dir(IMsdrDataReader):
    try:
        setattr(MassSpecDataReader, item, getattr(IMsdrDataReader, item))
    except:
        pass
for item in dir(IBDAChromFilter):
    try:
        setattr(BDAChromFilter, item, getattr(IBDAChromFilter, item))
    except:
        pass
for item in dir(IMsdrChargeStateAssignmentFilter):
    try:
        setattr(MsdrChargeStateAssignmentFilter, item, getattr(IMsdrChargeStateAssignmentFilter, item))
    except:
        pass


class MZFile:
    def __init__(self, datafile, **kwargs):
        self.file_type = '.d'
        self.data_file = datafile
        self._filters = None
        self.ticObj = None
        self.info = None
        self.scanRange = None
        self.noFilter = MsdrPeakFilter()
        self.source = MassSpecDataReader()
        s = self.source.OpenDataFile(self.source, datafile)
        if not s: raise IOError("Error opening %s" % datafile)

    def close(self):
        self.source.CloseDataFile(self.source)

    def time_range(self):
        if not self.ticObj:
            self.ticObj = self.source.GetTIC(self.source)
        assert self.ticObj.AcquiredTimeRange.Length == 1, "Multiple time ranges are not supported"
        return self.ranges.GetValue(0).Min, self.ranges.GetValue(0).Max

    def scan_range(self):
        return 0, self.source.FileInformation.MSScanFileInformation.TotalScansPresent

    def scan_info(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None):
        if self.info == None:
            self.info = []
            for index in range(self.source.FileInformation.MSScanFileInformation.TotalScansPresent):
                infoObj = self.source.GetScanRecord(self.source, index)
                rt = infoObj.RetentionTime
                mz = infoObj.MZOfInterest
                if not rt: break
                if start_time != None and rt <= start_time: continue
                if stop_time != None and rt >= stop_time: continue
                if start_mz != None and mz <= start_mz: continue
                if stop_mz != None and mz >= stop_mz: continue
                level = 'MS%d' % int(infoObj.MSLevel)
                polarity = str(infoObj.IonPolarity)
                # scanType = str(infoObj.MSScanType)
                self.info.append((rt, mz, index, level, polarity))
        return self.info

    def headers(self):
        return self.scan_info()

    def scan(self, index, mode=None):

        """
        Returns a spectrum from the specified scan index, where type is
        controlled by mode argument; defaults to prioiritize Peak (Centroid)
        and then try Profile. If the requested mode is not present, an empty
        spectrum  will be returned. Alternative modes are 'profile', 'centroid'
        or  'profileElsePeak'.
        """

        if mode != None: mode = mode.lower()
        if mode == None or mode == 'peakelseprofile':
            mode = DesiredMSStorageType.PeakElseProfile
        elif mode == 'profileelsepeak':
            mode = DesiredMSStorageType.ProfileElsePeak
        elif mode == 'profile':
            mode = DesiredMSStorageType.Profile
        elif mode == 'peak':
            mode = DesiredMSStorageType.Peak
        else:
            return []
        scanObj = self.source.GetSpectrum(self.source, index, self.noFilter, self.noFilter, mode)
        return list(zip(scanObj.XArray, scanObj.YArray))

    def cscan(self, index):
        """
               Calculates a centroided scan from profile-mode data. If a profile
               mode copy of the specified scan is not available, this raises an
               exception; in that case, you can use mzFile.scan() with mode =
               'centroid' to return the machine-centroided scan.
               """

        mode = DesiredMSStorageType.Profile
        scanObj = self.source.GetSpectrum(self.source, index, self.noFilter, self.noFilter, mode)
        mzs, ints = list(scanObj.XArray), list(scanObj.YArray)
        if not mzs:
            raise IOError("Profile data for scan index %s not available." % index)

        threshold = average(ints)
        peaks = []
        peak = []
        for pt in zip(mzs, ints):
            if pt[1] > threshold:
                peak.append(pt)
            elif peak:
                centroid = average(list(zip(*peak))[0], weights=list(zip(*peak))[1]), max(list(zip(*peak))[1])
                peaks.append(centroid)
                peak = []

        return peaks

    def xic(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None, filter=None, UV=False):
        if filter:
            assert filter.strip().lower() == 'full ms', 'Thermo-style XIC filters are not supported for Agilent files.'

        # A full XIC can be performed with an existing TIC object
        if self.ticObj and not any([start_time, stop_time, start_mz, stop_mz]):
            return list(zip(self.ticObj.XArray, self.ticObj.YArray))

        if start_time == None: start_time = 0
        if stop_time == None: stop_time = 999999
        if start_mz == None: start_mz = 0
        if stop_mz == None: stop_mz = 999999

        chromFilter = BDAChromFilter()

        chromFilter.set_MSLevelFilter(chromFilter, MSLevel.MS)  # Alt value is MSLevel.MSMS
        if not UV:
            chromFilter.set_ChromatogramType(chromFilter, ChromType.ExtractedIon)
        else:
            chromFilter.set_ChromatogramType(chromFilter, ChromType.ExtractedWavelength)
        chromFilter.set_SingleChromatogramForAllMasses(chromFilter, True)

        mzRange = MinMaxRange()
        mzRange.set_Min(start_mz)
        mzRange.set_Max(stop_mz)
        chromFilter.set_IncludeMassRanges(chromFilter, (mzRange,))

        rtRange = MinMaxRange()
        rtRange.set_Min(start_time)
        rtRange.set_Max(stop_time)
        chromFilter.set_ScanRange(chromFilter, rtRange)

        xic = self.source.GetChromatogram(self.source, chromFilter).Get(0)
        return list(zip(xic.XArray, xic.YArray))

    def uv_trace(self):

        # Cannot verify functionality, so leaving a couple potential methods here
        nonmsDevs = self.source.GetNonmsDevices()
        if not nonmsDevs.Length: raise IOError("No NonmsDevices were available")
        return self.source.GetTWC(nonmsDevs[0])

    def average_scan(self, scan_range):
        """
        Average the scan data across multiple scans.
        :param scan_range: Tuple or list defining the scan range.
        :return: Averaged scan data as a sorted list of tuples (average mz, average intensity)
        """
        # Assuming `scans` is a list of lists, where each sublist contains tuples (mz, intensity)
        dists = list(chain(
            *[[s[i + 1][0] - s[i][0] for i in range(len(s) - 1)] for s in self.scans[scan_range[0]:scan_range[1]]]))

        # Find the smallest distance, adjust by a tiny bit to avoid rounding issues
        max_width = min(dists) - 0.000001

        # Aggregate points based on the max_width found above
        aggregated_scan = self.aggregate_points(list(chain(*self.scans[scan_range[0]:scan_range[1]])),
                                                MAX_WIDTH=max_width)

        if aggregated_scan:
            avg_scan = []
            for agg_pts in aggregated_scan:
                avg_mz = self.average([x[0] for x in agg_pts])  # Average mz values
                avg_int = sum([x[1] for x in agg_pts]) / len(self.scans)  # Average intensity values across scans
                avg_scan.append((avg_mz, avg_int))

            return sorted(avg_scan)
        else:
            print("Error: Aggregated scan data is empty")
            return []


def average(xs, weights=None):
    if weights:
        assert len(weights) == len(xs), "Weights must be same length as input sequence."
        return (sum([x * w for x, w in zip(xs, weights)])
                /
                float(len(xs)) * reduce((lambda x, y: x * y), weights, 1))
    else:
        return sum(xs) / float(len(xs))


# arg 3 initial=_initial_missing
# if initial is _initial_missing:
#     try:
#         value = next(it)
#     except StopIteration:
#         raise TypeError(
#             "reduce() of empty iterable with no initial value") from None
# else:
#     value = initial





def aggregate_points(pointlist, distance_function=None, MAX_WIDTH=0.025):
    """
    Given a list of points (which may be floats or tuples with the
    position-relevant value in the first position), returns a list of
    grouped points such that all points within a group are within MAX_WIDTH
    of each other.
    """
    peaks = []
    agg = []

    if distance_function is not None:
        distance = distance_function
    elif isinstance(pointlist[0], tuple):
        def distance(x, y):
            return x[0] - y[0]

        pointlist.sort()
    else:
        def distance(x, y):
            return x - y

        pointlist.sort()
    return 0

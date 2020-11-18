import os
import numpy as np
# require pythonnet, pip install pythonnet
import clr

pathtothisfile = os.path.dirname(__file__)
# print(pathtothisfile)
dlls = ['ThermoFisher.CommonCore.Data', 'ThermoFisher.CommonCore.RawFileReader',
        'ThermoFisher.CommonCore.BackgroundSubtraction', 'ThermoFisher.CommonCore.MassPrecisionEstimator']

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
            clr.AddReference(dll)

clr.AddReference('System.Collections')
from System import *
from System.Collections.Generic import *

import ThermoFisher
from ThermoFisher.CommonCore.Data import ToleranceUnits
from ThermoFisher.CommonCore.Data import Extensions
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, ChromatogramTraceSettings, DataUnits, Device, \
    GenericDataTypes, SampleType, Scan, TraceType
from ThermoFisher.CommonCore.Data.FilterEnums import IonizationModeType, MSOrderType
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings, IScanEventBase, IScanFilter, \
    RawFileClassification
from ThermoFisher.CommonCore.MassPrecisionEstimator import PrecisionEstimate
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter, SequenceFileReader

""""""
'''
APIs are similar to pymsfilereader(https://github.com/frallain/pymsfilereader), but some APIs have not be implemented yet."
Borrowed some code from pyRawRileReader from pDeep3: https://github.com/pFindStudio/pDeep3/blob/master/pDeep/pyRawFileReader/RawFileReader.py
'''


def DotNetArrayToNPArray(arr, dtype=np.float):
    return np.array(list(arr), dtype=dtype)


class RawFileReader(object):
    def __init__(self, filename, **kwargs):
        # static class members
        self.sampleType = {0: 'Unknown',
                           1: 'Blank',
                           2: 'QC',
                           3: 'Standard Clear (None)',
                           4: 'Standard Update (None)',
                           5: 'Standard Bracket (Open)',
                           6: 'Standard Bracket Start (multiple brackets)',
                           7: 'Standard Bracket End (multiple brackets)'}

        self.controllerType = {-1: 'No device',
                               0: 'MS',
                               1: 'Analog',
                               2: 'A/D card',
                               3: 'PDA',
                               4: 'UV',
                               'No device': -1,
                               'MS': 0,
                               'Analog': 1,
                               'A/D card': 2,
                               'PDA': 3,
                               'UV': 4}

        self.massAnalyzerType = {'ITMS': 0,
                                 'TQMS': 1,
                                 'SQMS': 2,
                                 'TOFMS': 3,
                                 'FTMS': 4,
                                 'Sector': 5,
                                 0: 'ITMS',
                                 1: 'TQMS',
                                 2: 'SQMS',
                                 3: 'TOFMS',
                                 4: 'FTMS',
                                 5: 'Sector'}
        self.activationType = {'CID': 0,
                               'MPD': 1,
                               'ECD': 2,
                               'PQD': 3,
                               'ETD': 4,
                               'HCD': 5,
                               'Any activation type': 6,
                               'SA': 7,
                               'PTR': 8,
                               'NETD': 9,
                               'NPTR': 10,
                               'UVPD': 11,
                               'ETHCD': 12,  # not Thermo's build-in activation types
                               'ETCID': 13,  # not Thermo's build-in activation types
                               0: 'CID',
                               1: 'MPD',
                               2: 'ECD',
                               3: 'PQD',
                               4: 'ETD',
                               5: 'HCD',
                               6: 'Any activation type',
                               7: 'SA',
                               8: 'PTR',
                               9: 'NETD',
                               10: 'NPTR',
                               11: 'UVPD',
                               12: 'ETHCD',  # not Thermo's build-in activation types
                               13: 'ETCID',  # not Thermo's build-in activation types
                               }

        self.detectorType = {'Valid': 0,
                             'Any': 1,
                             'NotValid': 2,
                             0: 'Valid',
                             1: 'Any',
                             2: 'NotValid',
                             }

        self.scanDataType = {'Centroid': 0,
                             'Profile': 1,
                             'Any': 2,
                             0: 'Centroid',
                             1: 'Profile',
                             2: 'Any',
                             }

        self.scanType = {'Full': 0,
                         'Zoom': 1,
                         'SIM': 2,
                         'SRM': 3,
                         'CRM': 4,
                         'Any': 5,
                         'Q1MS': 6,
                         'Q3MS': 7,
                         0: 'Full',
                         1: 'SIM',
                         2: 'Zoom',
                         3: 'SRM',
                         4: 'CRM',
                         5: 'Any',
                         6: 'Q1MS',
                         7: 'Q3MS',
                         }

        self.filename = os.path.abspath(filename)
        self.filename = os.path.normpath(self.filename)

        # Check to see if the specified RAW file exists
        if not os.path.isfile(self.filename):
            print('The file doesn\'t exist in the specified location - {}'.format(filename))
            return

        self.source = ThermoFisher.CommonCore.RawFileReader.RawFileReaderAdapter.FileFactory(self.filename)

        # Check for any errors in the RAW file
        if self.source.IsError:
            print('Error opening ({}) - {}'.format(self.source.FileError, filename))

        if not self.source.IsOpen or self.source.IsError:
            raise IOError(
                "RAWfile {0} could not be opened, is the file accessible ?".format(
                    self.filename))

        # Check if the RAW file is being acquired
        if self.source.InAcquisition:
            print('RAW file still being acquired - {}'.format(filename))
            # May be able to do some cool stuff here...

        self.source.SelectInstrument(ThermoFisher.CommonCore.Data.Business.Device.MS, 1)

        self.StartTime = self.source.RunHeaderEx.StartTime
        self.EndTime = self.source.RunHeaderEx.EndTime
        self.FirstSpectrumNumber = self.source.RunHeaderEx.FirstSpectrum
        self.LastSpectrumNumber = self.source.RunHeaderEx.LastSpectrum
        self.LowMass = self.source.RunHeaderEx.LowMass
        self.HighMass = self.source.RunHeaderEx.HighMass
        self.MassResolution = self.source.RunHeaderEx.MassResolution
        self.NumSpectra = self.source.RunHeaderEx.SpectraCount
        self.scanrange = [self.FirstSpectrumNumber, self.LastSpectrumNumber]
        self.timerange = [self.StartTime, self.EndTime]

        try:
            self.header = {}
            extra_header_info = self.source.GetTrailerExtraHeaderInformation()
            extra_header_values = self.source.GetTrailerExtraValues(2, True)
            extra_header_values = DotNetArrayToNPArray(extra_header_values, dtype=str)
            for i in range(len(extra_header_info)):
                item = extra_header_info[i]
                self.header[item.Label[:-1]] = extra_header_values[i]

            try:
                self.injection_time = float(self.Get_Header_Item('Ion Injection Time (ms)'))
                self.resolution = float(self.Get_Header_Item('FT Resolution'))
            except:
                self.injection_time = None
                self.resolution = None
        except:
            self.header = None
            print("Error getting header")

    def Get_Header_Item(self, item):
        return self.header[item]

    def Close(self):
        '''Closes a raw file and frees the associated memory.'''
        self.source.Dispose()

    def scan_range(self):
        return self.scanrange

    def time_range(self):
        return self.timerange

    def scan_time_from_scan_name(self, scan):
        return self.source.RetentionTimeFromScanNumber(scan)

    def GetFilters(self):
        """Returns the list of unique scan filters for the raw file. This function is only supported for MS
        device controllers."""
        return list(self.source.GetFilters())

    def PrintInfo(self):
        # Get some information from the header portions of the RAW file and
        # display that information.  The information is general information
        # pertaining to the RAW file.
        print('General File Information:')
        print('   RAW file: {}'.format(self.source.FileName))
        print('   RAW file version: {}'.format(self.source.FileHeader.Revision))
        print('   Creation date: {}'.format(self.source.FileHeader.CreationDate))
        print('   Operator: {}'.format(self.source.FileHeader.WhoCreatedId))
        print('   Number of instruments: {}'.format(self.source.InstrumentCount))
        print('   Description: {}'.format(self.source.FileHeader.FileDescription))
        print('   Instrument model: {}'.format(self.source.GetInstrumentData().Model))
        print('   Instrument name: {}'.format(self.source.GetInstrumentData().Name))
        print('   Serial number: {}'.format(self.source.GetInstrumentData().SerialNumber))
        print('   Software version: {}'.format(self.source.GetInstrumentData().SoftwareVersion))
        print('   Firmware version: {}'.format(self.source.GetInstrumentData().HardwareVersion))
        print('   Units: {}'.format(Enum.GetName(DataUnits, self.source.GetInstrumentData().Units)))
        print('   Mass resolution: {:.3f}'.format(self.source.RunHeaderEx.MassResolution))
        print('   Number of scans: {}'.format(self.source.RunHeaderEx.SpectraCount))
        print('   Scan range: {} - {}'.format(self.FirstSpectrumNumber, self.LastSpectrumNumber))
        print('   Time range: {:.2f} - {:.2f}'.format(self.StartTime, self.EndTime))
        print(
            '   Mass range: {:.4f} - {:.4f}'.format(self.source.RunHeaderEx.LowMass, self.source.RunHeaderEx.HighMass))
        print()

        # Get information related to the sample that was processed
        print('Sample Information:')
        print('   Sample name: {}'.format(self.source.SampleInformation.SampleName))
        print('   Sample id: {}'.format(self.source.SampleInformation.SampleId))
        print('   Sample type: {}'.format(Enum.GetName(SampleType, self.source.SampleInformation.SampleType)))
        print('   Sample comment: {}'.format(self.source.SampleInformation.Comment))
        print('   Sample vial: {}'.format(self.source.SampleInformation.Vial))
        print('   Sample volume: {}'.format(self.source.SampleInformation.SampleVolume))
        print('   Sample injection volume: {}'.format(self.source.SampleInformation.InjectionVolume))
        print('   Sample row number: {}'.format(self.source.SampleInformation.RowNumber))
        print('   Sample dilution factor: {}'.format(self.source.SampleInformation.DilutionFactor))
        print()

        # Read the first instrument method (most likely for the MS portion of
        # the instrument).  NOTE: This method reads the instrument methods
        # from the RAW file but the underlying code uses some Microsoft code
        # that hasn't been ported to Linux or MacOS.  Therefore this method
        # won't work on those platforms therefore the check for Windows.
        if 'Windows' in str(Environment.OSVersion):
            deviceNames = self.source.GetAllInstrumentNamesFromInstrumentMethod()
            for device in deviceNames:
                print('Instrument method: {}'.format(device))
            print()

    def ListTrailerExtraFields(self, scanrange=None):
        '''Reads and reports the trailer extra data fields present in the RAW
        file.        '''

        if scanrange is None:
            scanrange = [self.FirstSpectrumNumber, self.LastSpectrumNumber]

        # Get the Trailer Extra data fields present in the RAW file
        trailerFields = self.source.GetTrailerExtraHeaderInformation()

        # Display each value
        i = 0
        print('Trailer Extra Data Information:')

        for field in trailerFields:
            print('   Field {} = {} storing data of type {}'.format(
                i, field.Label, Enum.GetName(GenericDataTypes, field.DataType)))
            i += 1
        print()
        # Get the number of filters present in the RAW file
        numberFilters = len(self.source.GetFilters())

        # Get the scan filter for the first and last spectrum in the RAW file
        firstFilter = IScanFilter(self.source.GetFilterForScanNumber(scanrange[0]))
        lastFilter = IScanFilter(self.source.GetFilterForScanNumber(scanrange[1]))

        print('Filter Information:')
        print('   Scan filter (first scan): {}'.format(firstFilter.ToString()))
        print('   Scan filter (last scan): {}'.format(lastFilter.ToString()))
        print('   Total number of filters: {}'.format(numberFilters))
        print()

    def GetChromatogram(self, scanrange=None, trace_number=0):
        '''Reads the base peak chromatogram for the RAW file.
        Args:
            scanrange = [ start scan for the chromatogram, end scan for the chromatogram.]
        '''
        if scanrange is None:
            scanrange = [self.FirstSpectrumNumber, self.LastSpectrumNumber]
        # Define the settings for getting the Base Peak chromatogram #TraceType.BasePeak
        # Define the settings for getting the TIC chromatogram #TraceType.TIC
        settings = ChromatogramTraceSettings(TraceType.TIC)
        # Get the chromatogram from the RAW file.
        data = self.source.GetChromatogramData([settings], int(scanrange[0]), int(scanrange[1]))
        # Split the data into the chromatograms
        trace = ChromatogramSignal.FromChromatogramData(data)
        # Convert to Numpy Array
        self.ticdat = np.transpose(
            [DotNetArrayToNPArray(trace[trace_number].Times), DotNetArrayToNPArray(trace[trace_number].Intensities)])
        return self.ticdat

    def ReadScanInformation(self, scanrange=None, outputData=True):
        '''Reads the general scan information for each scan in the RAW file
        using the scan filter object and also the trailer extra data
        section for that same scan.
        Args:
            scanrange = [ start scan for the chromatogram, end scan for the chromatogram.]
            outputData (bool): the output data flag.
        '''
        if scanrange is None:
            scanrange = [self.FirstSpectrumNumber, self.LastSpectrumNumber]
        # Read each scan in the RAW File
        for scan in range(scanrange[0], scanrange[1]):
            # Get the retention time for this scan number.  This is one of
            # two comparable functions that will convert between retention
            # time and scan number.
            time = self.source.RetentionTimeFromScanNumber(scan)

            # Get the scan filter for this scan number
            scanFilter = IScanFilter(self.source.GetFilterForScanNumber(scan))

            # Get the scan event for this scan number
            scanEvent = IScanEventBase(self.source.GetScanEventForScanNumber(scan))

            # Get the ionizationMode, MS2 precursor mass, collision
            # energy, and isolation width for each scan
            if scanFilter.MSOrder == MSOrderType.Ms2:
                # Get the reaction information for the first precursor
                reaction = scanEvent.GetReaction(0)

                precursorMass = reaction.PrecursorMass
                collisionEnergy = reaction.CollisionEnergy
                isolationWidth = reaction.IsolationWidth
                monoisotopicMass = 0.0
                masterScan = 0
                ionizationMode = scanFilter.IonizationMode
                order = scanFilter.MSOrder

                # Get the trailer extra data for this scan and then look
                # for the monoisotopic m/z value in the trailer extra data
                # list
                trailerData = self.source.GetTrailerExtraInformation(scan)

                for i in range(trailerData.Length):
                    if trailerData.Labels[i] == 'Monoisotopic M/Z:':
                        monoisotopicMass = float(trailerData.Values[i])
                    elif trailerData.Labels[i] in ('Master Scan Number:', 'Master Scan Number', 'Master Index:'):
                        masterScan = int(trailerData.Values[i])

                if outputData:
                    print(
                        '''Scan number {} @ time {:.2f} - Master scan = {}, Ionization mode={},\
                        MS Order={}, Precursor mass={:.4f}, Monoisotopic Mass = {:.4f},\
                        Collision energy={:.2f}, Isolation width={:.2f}'''.format(
                            scan, time, masterScan, Enum.GetName(IonizationModeType, ionizationMode),
                            Enum.GetName(MSOrderType, order), precursorMass, monoisotopicMass,
                            collisionEnergy, isolationWidth))

            elif scanFilter.MSOrder == MSOrderType.Ms:
                scanDependents = self.source.GetScanDependents(scan, 5)

                print(
                    'Scan number {} @ time {:.2f} - Instrument type={}, Number dependent scans={}'.format(
                        scan, time, Enum.GetName(RawFileClassification, scanDependents.RawFileInstrumentType),
                        scanDependents.ScanDependentDetailArray.Length))

    def GetSpectrum(self, scanNumber, scanFilter=None, outputData=False):
        '''Gets the spectrum from the RAW file.

        Args:
            scanNumber (int): the scan number being read.
            scanFilter (str): the scan filter for that scan.
            outputData (bool): the output data flag.
        '''
        if scanFilter is None:
            scanFilter = IScanFilter(self.source.GetFilterForScanNumber(scanNumber))
        # Check for a valid scan filter
        if not scanFilter:
            return

        # Get the scan statistics from the RAW file for this scan number
        scanStatistics = self.source.GetScanStatsForScanNumber(scanNumber)

        # Check to see if the scan has centroid data or profile data.  Depending upon the
        # type of data, different methods will be used to read the data.  While the ReadAllSpectra
        # method demonstrates reading the data using the Scan.FromFile method, generating the
        # Scan object takes more time and memory to do, so that method isn't optimum.
        if scanStatistics.IsCentroidScan:
            # Get the centroid (label) data from the RAW file for this
            # scan
            centroidStream = self.source.GetCentroidStream(scanNumber, False)

            # Print the spectral data (mass, intensity, charge values).
            # Not all of the information in the high resolution centroid
            # (label data) object is reported in this example.  Please
            # check the documentation for more information about what is
            # available in high resolution centroid (label) data.
            if outputData:
                print('Spectrum (centroid/label) {} - {} points'.format(scanNumber, centroidStream.Length))
                for i in range(centroidStream.Length):
                    print('  {} - {:.4f}, {:.0f}, {:.0f}'.format(
                        i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]))
                print()
            self.data = np.transpose(
                [DotNetArrayToNPArray(centroidStream.Masses), DotNetArrayToNPArray(centroidStream.Intensities)])
        else:
            # Get the segmented (low res and profile) scan data
            segmentedScan = self.source.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics)

            # Print the spectral data (mass, intensity values)
            if outputData:
                print('Spectrum (normal data) {} - {} points'.format(scanNumber, segmentedScan.Positions.Length))
                for i in range(segmentedScan.Positions.Length):
                    print('  {} - {:.4f}, {:.0f}'.format(
                        i, segmentedScan.Positions[i], segmentedScan.Intensities[i]))
                print()
            self.data = np.transpose(
                [DotNetArrayToNPArray(segmentedScan.Positions), DotNetArrayToNPArray(segmentedScan.Intensities)])
        return self.data

    def GetAverageSpectrum(self, scanrange=None, outputData=False, filter="FTMS"):
        '''Gets the average spectrum from the RAW file.

        Args:
            scanrange = [ start scan for the chromatogram, end scan for the chromatogram.]
            outputData (bool): the output data flag.
        '''

        if scanrange is None:
            scanrange = [self.FirstSpectrumNumber, self.LastSpectrumNumber]

        # Create the mass options object that will be used when averaging
        # the scans
        options = Extensions.DefaultMassOptions(self.source)
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = 5.0

        # Get the scan filter for the first scan.  This scan filter will be used to located
        # scans within the given scan range of the same type
        if filter is None:
            scanFilter = IScanFilter(self.source.GetFilterForScanNumber(scanrange[0]))
        else:
            filterhelper = Extensions.BuildFilterHelper(self.source, filter)
            scanFilter = filterhelper.Filter
        scanStatistics = self.source.GetScanStatsForScanNumber(scanrange[0])

        # Get the average mass spectrum for the provided scan range. In addition to getting the
        # average scan using a scan range, the library also provides a similar method that takes
        # a time range.

        averageScan = Extensions.AverageScansInScanRange(
            self.source, scanrange[0], scanrange[1], scanFilter, options)

        if averageScan is None:
            filterhelper = Extensions.BuildFilterHelper(self.source, "Full")
            scanFilter = filterhelper.Filter
            averageScan = Extensions.AverageScansInScanRange(
                self.source, scanrange[0], scanrange[1], scanFilter, options)
        # This example uses a different method to get the same average spectrum that was calculated in the
        # previous portion of this method.  Instead of passing the start and end scan, a list of scans will
        # be passed to the GetAveragedMassSpectrum function.
        # scans = List[int]([1, 6, 7, 9, 11, 12, 14])
        # averageScan = Extensions.AverageScans(self.source, scans, options)

        # Check to see if the scan has centroid data or profile data.  Depending upon the
        # type of data, different methods will be used to read the data.  While the ReadAllSpectra
        # method demonstrates reading the data using the Scan.FromFile method, generating the
        # Scan object takes more time and memory to do, so that method isn't optimum.
        if scanStatistics.IsCentroidScan:
            # Get the centroid (label) data from the RAW file for this
            # scan
            centroidStream = averageScan.CentroidScan

            if outputData:
                for i in range(centroidStream.Length):
                    print('  {} - {:.4f}, {:.0f}, {:.0f}'.format(
                        i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]))
                print()
            self.data = np.transpose(
                [DotNetArrayToNPArray(centroidStream.Masses), DotNetArrayToNPArray(centroidStream.Intensities)])
        else:
            # Get the segmented (low res and profile) scan data
            segmentedScan = averageScan.SegmentedScan
            if outputData:
                for i in range(segmentedScan.Positions.Length):
                    print('  {} - {:.4f}, {:.0f}'.format(
                        i, segmentedScan.Positions[i], segmentedScan.Intensities[i]))
                print()
            self.data = np.transpose(
                [DotNetArrayToNPArray(segmentedScan.Positions), DotNetArrayToNPArray(segmentedScan.Intensities)])
        return self.data

    def CalculateMassPrecision(self, scanNumber=1):
        '''Calculates the mass precision for a spectrum.

        Args:
            scanNumber (int): the scan to process.
        '''

        # Get the scan from the RAW file
        scan = Scan.FromFile(self.source, scanNumber)

        # Get the scan event and from the scan event get the analyzer type for this scan
        scanEvent = IScanEventBase(self.source.GetScanEventForScanNumber(scanNumber))

        # Get the trailer extra data to get the ion time for this file
        logEntry = self.source.GetTrailerExtraInformation(scanNumber)

        trailerHeadings = List[str]()
        trailerValues = List[str]()
        for i in range(logEntry.Length):
            trailerHeadings.Add(logEntry.Labels[i])
            trailerValues.Add(logEntry.Values[i])

        # Create the mass precision estimate object
        precisionEstimate = PrecisionEstimate()

        # Get the ion time from the trailer extra data values
        ionTime = precisionEstimate.GetIonTime(scanEvent.MassAnalyzer, scan, trailerHeadings, trailerValues)

        # Calculate the mass precision for the scan
        listResults = precisionEstimate.GetMassPrecisionEstimate(
            scan, scanEvent.MassAnalyzer, ionTime, self.source.RunHeader.MassResolution)

        # Output the mass precision results
        ''''''
        # if listResults.Count:
        #    print('Mass Precision Results:')

        for result in listResults:
            print('Mass {:.5f}, mmu = {:.3f}, ppm = {:.2f}'.format(
                result.Mass, result.MassAccuracyInMmu, result.MassAccuracyInPpm))

    def GetSeqInfo(self):
        """Returns the sequence row number for this sample in an acquired sequence. The numbering
        starts at 1.
        NOTE : XCALIBUR INTERFACE "View/Report/Sample Information" part
        """
        # result = c_long()
        seq = SequenceFileReader.OpenSequence(self.filename)
        print(seq.Samples)
        # print(str(self.source.GetCompoundNames()))
        # error = self.source.SequenceFileReader()
        # if error:
        #    raise IOError("GetSeqRowNumber error : ", error)
        # return error


if __name__ == "__main__":
    test = u"C:\Python\\UniDec3\TestSpectra\\test.RAW"
    rr = RawFileReader(test)
    rr.GetSeqInfo()
    rr.Close()

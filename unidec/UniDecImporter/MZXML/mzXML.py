from copy import deepcopy
import numpy as np
import unidec.tools as ud
import lxml.etree as ET
import time
from pyteomics import mzxml

from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.ImportTools import get_resolution, merge_spectra, IndexedFile, IndexedScan


def get_data_from_spectrum(spectrum, threshold=-1):
    impdat = np.transpose([spectrum['m/z array'], spectrum['intensity array']])
    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 1] > threshold]
    return impdat


class MZXMLImporter(Importer):
    """
    Imports mzXML data files.

    Note: mzXML are 1 indexed, so the first scan is scan 1, not scan 0.
    """

    def __init__(self, path, **kwargs):
        """
        Imports mzXML file, adds the chromatogram into a single spectrum.
        :param path: .mzXML file path
        :param args: arguments (unused)
        :param kwargs: keywords (unused)
        :return: mzXMLimporter object
        """
        print("Reading mzXML:", path)
        super().__init__(path, **kwargs)
        self.msrun = mzxml.read(path)
        self.init_scans()

        self.cdms_support = True
        self.imms_support = False
        self.chrom_support = True

    def init_scans(self):
        self.times = []
        self.scans = []
        self.levels = []
        for i, spectrum in enumerate(self.msrun):
            try:
                level = spectrum["msLevel"]
                if level is None:
                    continue
                self.levels.append(level)
                t = spectrum["retentionTime"]
                id = spectrum["num"]
                self.times.append(float(t))
                self.scans.append(int(id))
            except Exception as e:
                print("Error Importing Data, Scan:", i, "Error:", e)

        self.times = np.array(self.times)
        self.scans = np.array(self.scans)
        self.levels = np.array(self.levels)
        self.scan_range = [np.amin(self.scans), np.amax(self.scans)]

    def get_single_scan(self, scan):
        # Iterate until we hit the desired scan, then return it, and reset the iterator
        # self.msrun.reset()
        # s = 0
        # # Iterate up until one scan minus the desired scan
        # while s < scan - 1:
        #     spectrum = self.msrun.next()
        #     s = int(spectrum['num'])
        spectrum = self.msrun.get_by_id(str(scan))
        # if int(spectrum['num']) != scan:
        #     raise ValueError(f"Scan {scan} not found in mzXML file")
        dat = get_data_from_spectrum(spectrum)
        # reset the reader
        # self.msrun.reset()
        return dat

    def avg_safe(self, scan_range=None, time_range=None):
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        self.msrun.reset()
        spec = self.msrun.get_by_id(str(scan_range[0]))
        data = get_data_from_spectrum(spec)
        resolution = get_resolution(data)
        axis = ud.nonlinear_axis(np.amin(data[:, 0]), np.amax(data[:, 0]), resolution)
        template = np.transpose([axis, np.zeros_like(axis)])
        # print("Length merge axis:", len(template))
        newdat = ud.mergedata(template, data)
        template[:, 1] += newdat[:, 1]

        # Get other data points
        index = 0
        while index <= scan_range[1] - scan_range[0]:
            try:
                spec = self.msrun.next()
            except:
                break
            index += 1
            s = int(spec['num'])
            if s in self.scans:
                if scan_range[0] < s <= scan_range[1]:
                    data = get_data_from_spectrum(spec)
                    newdat = ud.mergedata(template, data)
                    template[:, 1] += newdat[:, 1]
                elif s <= scan_range[0]:
                    pass
                else:
                    break
        self.msrun.reset()
        return template

    def get_all_scans(self, threshold=-1):
        newtimes = []
        newids = []
        self.data = []
        self.msrun.reset()
        for i, s in enumerate(self.scans):
            try:
                impdat = get_data_from_spectrum(self.msrun.next(), threshold=threshold)
                self.data.append(impdat)
                newtimes.append(self.times[i])
                newids.append(s)
            except Exception as e:
                print("mzXML import error")
                print(e)
        self.times = np.array(newtimes)
        self.scans = np.array(newids)
        return self.data


    def get_avg_scan(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        if self.filesize > 1e9 and self.data is None:
            data = self.avg_safe(scan_range, time_range)
        else:
            data = self.avg_fast(scan_range, time_range)
        return data

    def get_tic(self):
        tic = []
        print("Constructing TIC")
        self.msrun.reset()
        for i, s in enumerate(self.scans):
            spectrum = self.msrun.next()
            try:
                tot = float(spectrum['totIonCurrent'])
            except:
                tot = np.sum(spectrum['intensity array'])
            tic.append(tot)
        t = self.times
        ticdat = np.transpose([t, tic])
        return ticdat

    def get_polarity(self, scan=None):
        tree = ET.parse(self._file_path)
        root = tree.getroot()
        polarity = None
        # Define the namespaces used in the XML
        namespaces = {'mz': 'http://sashimi.sourceforge.net/schema_revision/mzXML_3.2'}
        scans = root.xpath('//mz:scan', namespaces=namespaces)
        if scans is not None:
            # Extract the polarity attribute from the scan element
            polarity = scans[0].attrib.get('polarity', None)
        if polarity == '+':
            print("Polarity: Positive")
            return "Positive"
        if polarity == '-':
            print("Polarity: Negative")
            return "Negative"
        else:
            print("Polarity: Unknown")
            return None

    def get_ms_order(self, scan=None):
        index = self.get_scan_index(scan)
        ms_order = self.levels[index]
        return ms_order


    def close(self):
        self.msrun.close()



if __name__ == "__main__":
    import time

    test = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.mzXML"
    tstart = time.perf_counter()
    d = MZXMLImporter(test)

    exit()
    d.get_single_scan(1000)
    print("Time:", time.perf_counter() - tstart)


    scan_range = [1, 1000]

    starttime = time.perf_counter()
    print(len(d.avg_fast(scan_range=scan_range)))
    print("Time:", time.perf_counter() - starttime)

    # test it again
    starttime = time.perf_counter()
    print(len(d.avg_fast(scan_range=scan_range)))
    print("Time:", time.perf_counter() - starttime)

    d = MZXMLImporter(test)
    startime = time.perf_counter()
    print(len(d.avg_safe(scan_range=scan_range)))
    print("Time:", time.perf_counter() - startime)

    # test it again
    startime = time.perf_counter()
    print(len(d.avg_safe(scan_range=scan_range)))
    print("Time:", time.perf_counter() - startime)
    exit()


    d.get_avg_scan(scan_range=[1, 5])
    print("Time:", time.perf_counter() - tstart)
    exit()

    d.get_tic()

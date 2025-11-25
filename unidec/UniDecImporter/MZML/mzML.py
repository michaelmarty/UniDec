import time
import numpy as np
import os
from copy import deepcopy
from unidec import tools as ud
import pymzml
from pymzml.utils.utils import index_gzip
import pymzml.obo
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.ImportTools import get_resolution, IndexedFile, IndexedScan
import re

__author__ = 'Michael.Marty'


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


def get_data_from_spectrum(spectrum, threshold=-1, check_duplicates=False):
    if spectrum == None:
        return
    impdat = np.transpose([spectrum.mz, spectrum.i])
    if check_duplicates:
        mz = np.unique(impdat[:, 0])
        if len(mz) < len(impdat):
            # remove duplicates
            # Sort mz and intensities based on mz values
            sarray = impdat[np.argsort(impdat[:, 0])]
            # Sum intensities for duplicate mz values
            mz, starts = np.unique(sarray[:, 0], return_index=True)
            inten = [np.sum(sarray[starts[i]:starts[i + 1], 1]) for i in range(len(starts) - 1)]
            inten.append(np.sum(sarray[starts[-1]:, 1]))
            impdat = np.transpose([mz, inten])
            pass
    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 1] > threshold]
    return impdat


def get_im_data_from_spectrum(spectrum, threshold=-1):
    array_params = spectrum._get_encoding_parameters("raw ion mobility array")
    dtarray = spectrum._decode(*array_params)

    impdat = np.transpose([spectrum.mz, dtarray, spectrum.i])

    impdat = impdat[impdat[:, 0] > 10]
    if threshold >= 0:
        impdat = impdat[impdat[:, 2] > threshold]
    return impdat


def get_inj_time(spectrum):
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

def correct_mzML_ID(mzmlpath):
    """
    Some mzML files have weird scan IDs, which pymzML cannot handle.
    This function creates a temporary mzML file with corrected integer IDs.
    """
    temp_path = mzmlpath + "_corrected.mzML"
    with open(mzmlpath, 'r') as fin, open(temp_path, 'w') as fout:
        for line in fin:
            if "<spectrum " not in line:
                fout.write(line)
                continue
            match = re.search(r'id="([^"]+)"', line)
            if match:
                # Find "index=" in the ID and extract the number after it
                index = re.search(r'index="(\d+)"', line)
                if index:
                    indexval = index.group(1)
                else:
                    raise Exception("Could not find index= in mzML id")
                try:
                    int_id = int(indexval) + 1
                    print(match.group(1), "->", int_id)
                    line = line.replace(f'id="{match.group(1)}"', f'id="scan={int_id}"')
                    print(line)
                except ValueError:
                    pass
            fout.write(line)
    # Rename original file and move corrected file to original name, overwriting if needed
    os.replace(mzmlpath, mzmlpath[:-5] + "_original.mzML")
    os.replace(temp_path, mzmlpath)
    return mzmlpath

FIELDNAMES = ["id", "name", "def", "is_a"]


class MZMLImporter(Importer):
    """
    Imports mzML data files.

    Note: Some mzML files are 1 indexed, so the first scan is 1, not 0. Others are not, which is frustrating...
    """
    def __init__(self, path, *args, **kwargs):
        super().__init__(path, **kwargs)

        self.msrun = pymzml.run.Reader(path) #, build_index_from_scratch=True)
        self.levels = None
        self.init_scans()

        self.cdms_support = True
        self.imms_support = True
        self.chrom_support = True
        print("Reading mzML file:", path)

    def init_scans(self):
        get_time = lambda s: s.scan_time_in_minutes()

        times, scans, levels = [], [], []
        for spectrum in self.msrun:
            level = getattr(spectrum, 'ms_level', None)
            if level is None:
                continue
            try:
                times.append(float(get_time(spectrum)))
                scans.append(spectrum.ID)
                levels.append(level)
            except Exception:
                continue

        self.times = np.array(times, dtype=np.float32)
        self.scans = np.array(scans)
        # Check if self.scans is all the same, if so, call correct_mzML_ID
        if len(np.unique(self.scans)) < len(self.scans) and len(self.scans) > 1:
            self.msrun.close()
            print("Duplicate Scan IDs Detected. Attempting to correct mzML IDs.")
            corrected_path = correct_mzML_ID(self._file_path)

            self.msrun = pymzml.run.Reader(corrected_path)
            # Re-initialize scans
            times, scans, levels = [], [], []
            for spectrum in self.msrun:
                level = getattr(spectrum, 'ms_level', None)
                if level is None:
                    continue
                try:
                    times.append(float(get_time(spectrum)))
                    scans.append(spectrum.ID)
                    levels.append(level)
                except Exception:
                    continue

            self.times = np.array(times, dtype=np.float32)
            self.scans = np.array(scans)
        self.levels = np.array(levels)
        self.scan_range = [int(self.scans.min()), int(self.scans.max())]
        self.reset_reader()

    def get_single_scan(self, scan):
        try:
            data = get_data_from_spectrum(self.msrun[scan])
        except Exception as e:
            print("Error in get_single_scan:", e)
            data = None
        return data

    #This is the merging function
    def avg_safe(self, scan_range=None, time_range=None):
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        self.reset_reader()
        # Get first data point and interpolate it to match template axis
        data = get_data_from_spectrum(self.msrun[scan_range[0]])
        resolution = get_resolution(data)
        axis = ud.nonlinear_axis(np.amin(data[:, 0]), np.amax(data[:, 0]), resolution)
        template = np.transpose([axis, np.zeros_like(axis)])
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
            if spec.ID in self.scans:
                if scan_range[0] < spec.ID <= scan_range[1]:
                    data = get_data_from_spectrum(spec)
                    newdat = ud.mergedata(template, data)
                    template[:, 1] += newdat[:, 1]
                elif spec.ID <= scan_range[0]:
                    pass
                else:
                    break
        self.reset_reader()
        return template

    def reset_reader(self):
        self.msrun.close()
        self.msrun = pymzml.run.Reader(self._file_path)

    def get_all_scans(self, threshold=-1):
        self.reset_reader()
        newtimes = []
        newids = []
        self.data = []
        for n, spec in enumerate(self.msrun):
            # print("Processing scan:", spec.ID)
            if spec.ID in self.scans:
                try:
                    impdat = get_data_from_spectrum(spec, threshold=threshold)
                    self.data.append(impdat)
                    newtimes.append(self.times[n])
                    newids.append(spec.ID)
                except Exception as e:
                    print("Failed to get all scans:", e)
        self.times = np.array(newtimes)
        self.scans = np.array(newids)
        self.reset_reader()
        return self.data

    def get_avg_scan(self, scan_range=None, time_range=None):
        """
        Returns merged 1D MS data from mzML import
        :return: merged data
        """
        starttime = time.perf_counter()
        if self.filesize > 1e9 and self.data is None:
            data = self.avg_safe(scan_range, time_range)
        else:
            data = self.avg_fast(scan_range, time_range)
        print("Import Time:", time.perf_counter() - starttime)
        return data

    def get_tic(self):
        try:
            tic = self.msrun["TIC"]
            ticdat = np.transpose([tic.time, tic.i])
            if len(ticdat) != len(self.scans):
                print("TIC wrong size. Likely extra non-MS scans", len(ticdat), len(self.scans))
                raise Exception
            if not np.isclose(ticdat[-1,0],  self.times[-1]):
                print("TIC end time not equal to scan end time")
                raise Exception
        except:
            print("Error getting TIC in mzML; trying to make it...")
            tic = []
            self.get_all_scans()
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

    def get_inj_time_array(self):
        self.reset_reader()
        its = []
        for i, s in enumerate(self.scans):
            it = get_inj_time(self.msrun[s])
            try:
                it = float(it)
            except:
                print("Error in scan header:", i, s, it)
                it = 1
            its.append(it)
        self.reset_reader()
        return np.array(its)

    def get_all_imms_scans(self):
        self.reset_reader()
        newtimes = []
        newids = []
        self.immsdata = []
        for n, spec in enumerate(self.msrun):
            if spec.ID in self.scans:
                try:
                    impdat = get_im_data_from_spectrum(spec)
                    self.immsdata.append(impdat)
                    newtimes.append(self.times[n])
                    newids.append(spec.ID)
                except Exception as e:
                    print(f"Error processing spectrum ID {spec.ID}: {e}")

        self.times = np.array(newtimes)
        self.scans = np.array(newids)
        self.reset_reader()
        return self.immsdata

    # grab pol info from single scan. Similar to old impl but use direct idxing
    def get_polarity(self, scan=1):
        if scan == -1:
            scan = self.scans[0]
        # Directly access the scan at the provided index
        try:
            spec = self.msrun[scan]
        except Exception as e:
            try:
                print(e)
                print("Error getting scan", scan, "trying to get scan", self.scans[1])
                scan = int(self.scans[1])
                spec = self.msrun[scan]
            except Exception as e:
                print("Error getting scan", scan, e)
                print("Defaulting to positive polarity")
                return "Positive"  # Default to positive if error occurs
        # to string gets the raw byte xml format of the scan
        comp = ""
        xml = spec.to_string()
        # convert the byte xml to a string
        for i in xml:
            comp += chr(i)
        # if we want to look at the raw xml data
        # look for the string telling us what polarity is
        if "negative scan" in comp:
            print("Polarity: Negative")
            self.polarity = "Negative"
            return "Negative"
        if "positive scan" in comp:
            print("Polarity: Positive")
            self.polarity = "Positive"
            return "Positive"
        return None


    def get_isolation_mz_width(self, s):
        scan = self.msrun[s]
        precursors = scan.selected_precursors
        if len(precursors) == 0:
            return None, None
        else:
            precursor = precursors[0]
            mz = precursor['mz']
            width = 3
            return mz, width

    def get_cdms_data(self, scan_range=None):
        raw_dat = self.get_all_scans(threshold=0)

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
        self.msrun.close()

if __name__ == "__main__":
    test = "Z:\\Group Share\\JGP\\DataForJoe\\TF_centroided.mzML"
    # test = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\test_mzml.mzML"
    test = "Z:\\Group Share\\JGP\\DataForJoe\\TF_centroided.mzML"
    # test = "C:\\Data\\Sergei\\50.mzML"
    test = r"C:\Data\DataTypeCollection\mzML_collection\20250523_HBD_200ng_NR_1.mzML"
    test = r"C:\Data\DataTypeCollection\mzML_collection\20250901_MB_NOPFBS_IMS_10eV.mzML"
    importer = MZMLImporter(test)
    print(importer.scans)
    # print(importer.times)
    # print(importer.get_tic())
    # print(importer.msrun["cycle=10"])
    # s1 = importer.get_all_scans()[10]
    # s2 = importer.get_all_imms_scans()[10]
    # print(np.sum(s1[:,1]), np.sum(s2[:,2]))

    import matplotlib.pyplot as plt


    mzdata = importer.get_avg_scan()
    mzdata = ud.linearize(mzdata, 1, 0)
    plt.plot(mzdata[:,0], mzdata[:,1])
    imzdata = importer.get_imms_avg_scan()
    # Reshape imzdata into 2D array for plotting
    mz_values = np.unique(imzdata[:, 0])
    im_values = np.unique(imzdata[:, 1])
    intensity_matrix = imzdata[:, 2].reshape(len(mz_values), len(im_values))
    mzdata2 = np.transpose([mz_values, np.sum(intensity_matrix, axis=1)])
    plt.plot(mzdata2[:,0], mzdata2[:,1]/np.amax(mzdata2[:,1]) * np.amax(mzdata[:,1]), 'r--')

    print(np.sum(mzdata[:,1]), np.sum(mzdata2[:,1]), np.sum(mzdata[:,1])/np.sum(mzdata2[:,1]))
    plt.show()


    exit()
    #test = "C:\\Data\\RileyLab\\exportMGF_10spectra.mzML"
    #test = "C:\\Data\\DataTypeCollection\\IMMS\\test_agilentimms_mzml.mzML"
    # test = "C:\\Data\\DataTypeCollection\\IMMS\\test_watersimms_mzml.mzML"
    # importer = MZMLImporter(test)
    # # importer.msrun[1]
    # importer.get_single_scan(0)
    # importer.get_all_scans()

    msrun = pymzml.run.Reader(test)
    msrun[0]
    msrun[1]
    msrun[0]
    msrun[9]

    exit()
    test1 = "Z:\\Group Share\\JGP\\DataForJoe\\TF.mzML"
    cimp = MZMLImporter(test1)
    rdat = cimp.get_avg_scan()
    print("Length of normal", len(rdat))
    print("Length of centroided", len(cdat))
    exit()
    scan_range = [1, 5]
    d.get_imms_avg_scan(scan_range=scan_range)

    exit()



    starttime = time.perf_counter()
    print(len(d.avg_fast(scan_range=scan_range)))
    print("Time:", time.perf_counter() - starttime)

    # test it again
    starttime = time.perf_counter()
    print(len(d.avg_fast(scan_range=scan_range)))
    print("Time:", time.perf_counter() - starttime)

    d = MZMLImporter(test)
    startime = time.perf_counter()
    print(len(d.avg_safe(scan_range=scan_range)))
    print("Time:", time.perf_counter() - startime)

    # test it again
    startime = time.perf_counter()
    print(len(d.avg_safe(scan_range=scan_range)))
    print("Time:", time.perf_counter() - startime)

    # print(d.get_cdms_data())


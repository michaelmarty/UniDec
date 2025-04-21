from pyopenms import *
import os
import math
import time
from operator import itemgetter
from numba import njit
from collections import defaultdict

@njit
def compute_bin_indices(mzs, intensities, min_mz, bin_width):
    indices = ((mzs - min_mz) / bin_width).astype(np.int32)
    return indices, mzs, intensities

class XIC:
    """
    Takes a centroided MSExperiment and extracts the XIC of a given precursor.
    """
    def __init__(self, indexed_file, precursor, mz_width):
        self.indexed_file = indexed_file
        self.precursor = precursor
        self.mz_width = mz_width
        self.xic = self.extract_xic()

class IndexedScan:
    """
    A class to store an in individual scan with the peaks indexed by m/z.
    """
    def __init__(self, scan, rt, scan_num, indexed_file):
        self.scan = scan
        self.retention_time = rt
        self.scan_num = scan_num
        self.min_mz = indexed_file.min_mz
        self.max_mz = indexed_file.max_mz
        self.bin_width = indexed_file.bin_width
        self.indexed_peaks = self.index_peaks(indexed_file.min_mz, indexed_file.bin_width)


    def index_peaks(self, min_mz, bin_width):
        peaks = self.scan.get_peaks()
        mzs = np.array(peaks[0])
        intensities = np.array(peaks[1])
        indices, mzs, intensities = compute_bin_indices(mzs, intensities, min_mz, bin_width)

        indexed_peaks = defaultdict(list)
        for mz, intensity, idx in zip(mzs, intensities, indices):
            indexed_peaks[idx].append([mz, intensity])

        return indexed_peaks

    def get_intensity(self, mz, mz_tolerance):
        #this is easy if the mz+/- the tolerance is within a single bin..
        min_index = math.trunc((mz - mz_tolerance - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tolerance - self.min_mz) / self.bin_width)
        if min_index == max_index:
            sum_intensity = 0
            #sum the intensities of all peaks in the bin within the tolerance of the mz
            if self.indexed_peaks.get(min_index, 0) == 0:
                return 0
            else:
                bin_peaks = self.indexed_peaks.get(min_index)
                for i in range(len(bin_peaks)):
                    if abs(bin_peaks[i][0] - mz) <= mz_tolerance:
                        sum_intensity += bin_peaks[i][1]
            return sum_intensity
        else:
            #check peaks from both bins..
            sum_intensity = 0
            indices = range(min_index, max_index + 1, 1)
            for i in range(len(indices)):
                if self.indexed_peaks.get(indices[i], 0) == 0:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(indices[i])
                    for j in range(len(bin_peaks)):
                        if abs(bin_peaks[j][0] - mz) <= mz_tolerance:
                            sum_intensity += bin_peaks[j][1]
            return sum_intensity

    def get_exp_mz(self, mz, mz_tol):
        min_index = math.trunc((mz - mz_tol - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tol - self.min_mz) / self.bin_width)
        if min_index == max_index:
            peaks_sorted = sorted(self.indexed_peaks.get(min_index), key=itemgetter(1), reverse=True)
        else:
            indices = range(min_index, max_index+1, 1)
            peaks_within_tol = []
            for s in indices:
                if self.indexed_peaks.get(s, 0) == 0:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(s)
                    for p in bin_peaks:
                        if abs(p[0] - mz) <= mz_tol:
                            peaks_within_tol.append(p)
            peaks_sorted = sorted(peaks_within_tol, key=itemgetter(1), reverse=True)
        return peaks_sorted[0][0]

class IndexedFile:
    """
    A class to store a set of indexed scans originating from a single mzML file.
    """
    def __init__(self, exp):
        self.exp = exp
        self.min_mz = 10
        self.max_mz = 2500
        self.bin_width = 1
        self.indexed_scans = self.index_scans()


    def index_scans(self):
        starttime = time.perf_counter()
        indexed_scans = []
        for i, s in enumerate(self.exp):
            if s.getMSLevel() != 1:
                continue
            indexed_scan = IndexedScan(s, s.getRT(), i+1 , self)
            indexed_scans.append(indexed_scan)
        print("Indexing took", time.perf_counter() - starttime, "seconds.")
        return indexed_scans

    def extract_xic(self, mz, mz_width, t_min, t_max):
        rts = []
        scans = []
        intensities = []
        for i in range(len(self.indexed_scans)):
            if self.indexed_scans[i].retention_time/60 < t_min or self.indexed_scans[i].retention_time/60 > t_max:
                continue
            else:
                intensity = self.indexed_scans[i].get_intensity(mz, mz_width)
                scans.append(self.indexed_scans[i].scan_num)
                rts.append(self.indexed_scans[i].retention_time)
                intensities.append(intensity)
        return np.transpose([scans, rts, intensities])


    def get_indexed_spectrum_atRt(self, rt):
        index_guess = int(round(len(self.indexed_scans)/2))
        finding_spec = True
        min_error = self.indexed_scans[index_guess].retention_time - rt
        while finding_spec:
            current_rt = self.indexed_scans[index_guess].retention_time
            current_error = current_rt - rt

            if abs(current_error) <= abs(min_error):
                min_error = current_error

                if min_error < 0:
                    index_guess += 1
                if min_error > 0:
                    index_guess -= 1
                if min_error == 0:
                    finding_spec = False


            else:
                finding_spec = False


        return self.indexed_scans[index_guess]
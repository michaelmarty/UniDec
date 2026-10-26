import os
import numpy as np
import unidec.tools as ud
from UniDecImporter.ImporterFactory import ImporterFactory
from unidec.metaunidec.mudeng import MetaUniDec, metaunidec_call
from copy import deepcopy

chrom_file_exts = [".raw", ".Raw", ".RAW", ".d", ".mzML.gz", ".mzML"]


def chromatographic_peak_scores(xvalues, intensities):
    """Return FWHM symmetry and within-FWHM monotonicity scores."""
    xvalues = np.asarray(xvalues, dtype=float)
    intensities = np.asarray(intensities, dtype=float)
    if (xvalues.ndim != 1 or intensities.ndim != 1 or len(xvalues) != len(intensities)
            or len(xvalues) < 3 or not np.all(np.isfinite(xvalues))
            or not np.all(np.isfinite(intensities)) or np.any(np.diff(xvalues) <= 0)):
        return 0., 0.

    intensities = np.maximum(intensities, 0)
    height = np.max(intensities)
    if height <= 0:
        return 0., 0.

    apex_start = int(np.argmax(intensities))
    apex_end = apex_start
    while apex_start > 0 and intensities[apex_start - 1] == height:
        apex_start -= 1
    while apex_end < len(intensities) - 1 and intensities[apex_end + 1] == height:
        apex_end += 1
    apex_x = (xvalues[apex_start] + xvalues[apex_end]) / 2.
    half_height = height / 2.

    left_index = apex_start - 1
    while left_index >= 0 and intensities[left_index] > half_height:
        left_index -= 1
    right_index = apex_end + 1
    while right_index < len(intensities) and intensities[right_index] > half_height:
        right_index += 1
    if left_index < 0 or right_index >= len(intensities):
        return 0., 0.

    def crossing(index1, index2):
        y1, y2 = intensities[index1], intensities[index2]
        if y1 == half_height:
            return xvalues[index1]
        return xvalues[index1] + (half_height - y1) * (xvalues[index2] - xvalues[index1]) / (y2 - y1)

    left_x = crossing(left_index, left_index + 1)
    right_x = crossing(right_index, right_index - 1)
    left_width = apex_x - left_x
    right_width = right_x - apex_x
    if left_width <= 0 or right_width <= 0:
        return 0., 0.
    symmetry = min(left_width, right_width) / max(left_width, right_width)

    left_mask = np.logical_and(xvalues < apex_x, xvalues > left_x)
    left_trace_x = np.concatenate(([apex_x], xvalues[left_mask][::-1], [left_x]))
    left_trace_y = np.concatenate(([height], intensities[left_mask][::-1], [half_height]))
    right_mask = np.logical_and(xvalues > apex_x, xvalues < right_x)
    right_trace_x = np.concatenate(([apex_x], xvalues[right_mask], [right_x]))
    right_trace_y = np.concatenate(([height], intensities[right_mask], [half_height]))

    def side_areas(side_x, side_y):
        distance = np.abs(side_x - apex_x)
        minimum = side_y[0]
        violations = np.zeros(len(side_y))
        for i, value in enumerate(side_y[1:], 1):
            if value > minimum:
                violations[i] = value - minimum
            else:
                minimum = value
        widths = np.diff(distance) / 2.
        violation_area = np.sum((violations[:-1] + violations[1:]) * widths)
        trace_area = np.sum((side_y[:-1] + side_y[1:]) * widths)
        return violation_area, trace_area

    left_violation, left_area = side_areas(left_trace_x, left_trace_y)
    right_violation, right_area = side_areas(right_trace_x, right_trace_y)
    total_area = left_area + right_area
    if total_area <= 0 or not np.isfinite(total_area):
        return 0., 0.
    monotonicity = np.clip(1. - (left_violation + right_violation) / total_area, 0., 1.)
    return float(np.clip(symmetry, 0., 1.)), float(monotonicity)


class ChromEngine(MetaUniDec):
    """
    UniChrom Engine.
    """

    def __init__(self):
        MetaUniDec.__init__(self)
        self.polarity = None
        self.chrompeaks = None
        self.chrompeaks_tranges = None
        self.outpath = None
        self.attrs = None
        self.scans = None
        self.fullscans = None
        self.filename = None
        self.dirname = None
        self.path = None
        self.chromdat = None
        self.tic = None
        self.ticdat = None
        self.bpc = None
        self.bpcdat = None
        self.spectra = None
        self.massdat = None
        self.mzdata = None
        self.procdata = None
        self.config.default_high_res()
        # Zero preserves the original scan-wise MetaUniDec workflow. Positive
        # values enable coupled chromatographic deconvolution in the C engine.
        self.config.dtsig = 0.0
    def run_unidec(self):
        """Select scan-wise MetaUniDec or coupled UniChrom using ``dtsig``."""
        if not self.check_badness():
            self.pks.peaks = []
            self.config.write_hdf5()
            self.out = metaunidec_call(self.config)
            if self.out == 0:
                self.data.import_hdf5()
                self.update_history()
            return self.out
        return 1

    def pick_peaks(self, refresh=True):
        super().pick_peaks(refresh=refresh)
        if not self.config.calc_unichrom_scores:
            return
        current_peaks = np.asarray(self.data.peaks)
        base_scores = np.array([p.dscore for p in self.pks.peaks], dtype=float)
        saved_scores = np.asarray(self.data.unichrom_scores)
        if (not refresh and saved_scores.ndim == 2 and saved_scores.shape == (len(base_scores), 5)
                and current_peaks.ndim == 2 and current_peaks.shape[1] >= 3
                and np.allclose(saved_scores[:, 0], current_peaks[:, 0])
                and np.allclose(saved_scores[:, 4], current_peaks[:, 2])):
            base_scores = saved_scores[:, 1]

        score_rows = []
        dscores = []
        weights = []
        for i, peak in enumerate(self.pks.peaks):
            peak.base_dscore = float(np.clip(base_scores[i], 0., 1.)) if np.isfinite(base_scores[i]) else 0.
            peak.chrom_symmetry_score, peak.chrom_monotonicity_score = chromatographic_peak_scores(
                self.data.var1, peak.extracts)
            peak.dscore = (peak.base_dscore * peak.chrom_symmetry_score
                           * peak.chrom_monotonicity_score)
            if current_peaks.ndim == 2 and current_peaks.shape[1] >= 3:
                current_peaks[i, 2] = peak.dscore
            score_rows.append([peak.mass, peak.base_dscore, peak.chrom_symmetry_score,
                               peak.chrom_monotonicity_score, peak.dscore])
            if np.isfinite(peak.height) and peak.height > 0:
                dscores.append(peak.dscore)
                weights.append(peak.height ** 2)

        self.pks.uniscore = float(np.clip(np.average(dscores, weights=weights), 0., 1.)) if weights else 0.
        self.data.peaks = current_peaks
        self.data.export_unichrom_scores(np.asarray(score_rows).reshape((-1, 5)), self.pks.uniscore)

    def open_chrom(self, path, load_hdf5=True, refresh=False):
        if os.path.splitext(path)[1] == ".hdf5":
            self.open_hdf5_file(path)
            return True
        else:
            self.dirname, self.filename = os.path.split(path)

            print("Opening: ", self.filename)
            self.path = path
            load_hdf5 = self.load_mzml(self.path, load_hdf5, refresh=refresh)
            return load_hdf5

    def open_hdf5_file(self, path):
        print("Opening HDF5 File: ", path)
        header = os.path.splitext(path)[0]
        for ext in chrom_file_exts:
            testpath = header + ext
            if os.path.isfile(testpath) or os.path.isdir(testpath):
                try:
                    self.open_chrom(testpath)
                    # TODO: Make it so that you could match HDF5 file and data with different name
                    return True, testpath
                except Exception as e:
                    print("Error opening: ", testpath)
                    print(e)
        print("Attempted to open: ", path)
        print("Unable to find chromatograph file matching the header: ", header)
        print("Please rename the files such that the raw data matches the HDF5 file name")
        raise FileNotFoundError
        # return False, header

    def load_mzml(self, path, load_hdf5=True, refresh=False, *args, **kwargs):
        self.path = path
        name = os.path.splitext(path)[0]
        if name[-5:].lower() == ".mzml":
            name = name[:-5]
        self.outpath = name + ".hdf5"
        self.setup_filenames(self.outpath)
        self.data.filename = self.outpath
        hdf5 = False
        self.clear()
        if os.path.isfile(self.outpath) and load_hdf5 and not refresh:
            print('Opening HDF5 File:', self.outpath)
            try:
                self.open(self.outpath)
                hdf5 = True
            except Exception as e:
                print("Error opening prior hdf5 file:", e)
        if not os.path.isfile(self.outpath):
            self.data.new_file(self.outpath)
            self.open(self.outpath)

        self.update_history()

        self.chromdat = ImporterFactory.create_importer(path)
        self.polarity = self.chromdat.get_polarity()
        self.tic = self.chromdat.get_tic()
        self.ticdat = np.array(self.tic)
        self.bpc = self.chromdat.get_bpc()
        self.bpcdat = np.array(self.bpc)
        self.fullscans = self.chromdat.scans
        return hdf5

    def get_selected_chromatogram(self):
        """Return the cached chromatogram selected for plotting and parsing."""
        if self.config.chromatogram_type == 1:
            return self.bpcdat, "Base Peak Chromatogram"
        return self.ticdat, "Total Ion Chromatogram"

    def get_data_from_scans(self, scan_range=None):
        self.mzdata = self.chromdat.get_avg_scan(scan_range=scan_range)
        self.procdata = None

        minscan, maxscan = scan_range
        time_range = self.chromdat.get_times_from_scans(scan_range)
        minval, midval, maxval = time_range
        self.scans = [minscan, maxscan, minval, maxval]
        print("Scans:", self.scans)

        attrs = {"timestart": minval, "timeend": maxval,
                 "timemid": midval, "retention_time": midval,
                 "scanstart": minscan, "scanend": maxscan,
                 "scanmid": (minscan + maxscan) / 2.}
        self.attrs = attrs

        return self.mzdata

    def get_data_from_times(self, minval, maxval):
        scan_range = self.chromdat.get_scans_from_times([minval, maxval])
        self.get_data_from_scans(scan_range)
        return self.mzdata

    def get_minmax_times(self):
        chromatogram, _ = self.get_selected_chromatogram()
        return np.amin(chromatogram[:, 0]), np.amax(chromatogram[:, 0])

    def get_times(self):
        starts = []
        ends = []
        for s in self.data.spectra:
            starts.append(s.attrs["timestart"])
            ends.append(s.attrs["timeend"])
        return np.array(starts), np.array(ends)

    def reset_vars(self, export=True):
        self.data.update_var_array()
        if export:
            self.data.export_hdf5(delete=True)
        self.data.v1name = "timemid"
        self.data.v2name = "scanstart"
        self.data.import_vars(get_vnames=False)

    def get_chrom_peaks(self, window=None):
        # Clean up the selected TIC or BPC data.
        chromatogram, _ = self.get_selected_chromatogram()
        chromdat = deepcopy(chromatogram)
        chromdat = ud.gsmooth(chromdat, 2)
        chromdat[:, 1] -= np.amin(chromdat[:, 1])
        # chromdat = ud.gaussian_backgroud_subtract(chromdat, 100)
        maxval = np.amax(chromdat[:, 1])
        chromdat[:, 1] /= maxval
        maxt = np.amax(chromdat[:, 0])
        mint = np.amin(chromdat[:, 0])

        # Set Window
        if window is None:
            window = self.config.chrom_peak_width

        # Set Threshold
        noise = ud.noise_level2(chromdat, percent=0.50)
        print("Noise Level:", noise, "Window:", window)

        # Detect Peaks
        peaks = ud.peakdetect_nonlinear(chromdat, window=window, threshold=noise)

        # Filter Peaks
        goodpeaks = []
        tranges = []
        diffs = np.diff(chromdat[:, 0])
        for p in peaks:
            fwhm, fwhmrange = ud.calc_FWHM(p[0], chromdat)
            index = ud.nearest(chromdat[:, 0], p[0])
            if index >= len(diffs):
                index = len(diffs) - 1
            localdiff = diffs[index]
            if p[0] - fwhm / 2. < mint or p[0] + fwhm / 2. > maxt or fwhm > 4 * window or fwhm < localdiff * 2 \
                    or fwhmrange[0] == p[0] or fwhmrange[1] == p[0]:
                print("Bad Peak", p, fwhm, fwhmrange)
                pass
            else:
                print(p[0], fwhm)
                goodpeaks.append(p)
                tranges.append(fwhmrange)
        self.chrompeaks = goodpeaks
        self.chrompeaks_tranges = tranges
        return goodpeaks, tranges

    def add_manual_selection(self):
        self.data.add_data(self.mzdata, name=str(self.scans[2]), attrs=self.attrs, export=False)

    def add_regular_times(self):
        chromatogram, _ = self.get_selected_chromatogram()
        times = np.arange(np.amin(chromatogram[:, 0]), np.amax(chromatogram[:, 0]), self.config.time_window)

        if self.config.time_start is not None and self.config.time_end is not None:
            times = np.arange(self.config.time_start, np.amax(chromatogram[:, 0]), self.config.time_window)
            boo1 = times >= self.config.time_start
            boo2 = times <= self.config.time_end
            boo3 = np.logical_and(boo1, boo2)
            times = times[boo3]

        self.data.clear()
        for i, t in enumerate(times):
            print(t)
            data = self.get_data_from_times(t, t + self.config.time_window)
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)

    def add_regular_scans(self):
        if self.config.time_start is not None and self.config.time_end is not None:
            # get scan for self.config.time_start and self.config.time_end
            scan_range = self.chromdat.get_scans_from_times([self.config.time_start, self.config.time_end])
            scans = np.arange(scan_range[0], scan_range[1]+self.config.scan_window, self.config.scan_window)
        else:
            scans = np.arange(np.amin(self.chromdat.scans), len(self.chromdat.scans), self.config.scan_window)


        self.data.clear()
        for i, t in enumerate(scans):
            print(t)
            data = self.get_data_from_scans([t, t + self.config.scan_window])
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)

    def add_chrom_peaks(self):
        self.get_chrom_peaks()
        times = np.array(self.chrompeaks_tranges)

        if self.config.time_start is not None and self.config.time_end is not None:
            boo1 = times[:, 0] >= self.config.time_start
            boo2 = times[:, 1] <= self.config.time_end
            boo3 = np.logical_and(boo1, boo2)
            times = times[boo3]

        self.data.clear()
        for i, t in enumerate(times):
            data = self.get_data_from_times(t[0], t[1])
            self.data.add_data(data, name=str(t[0]), attrs=self.attrs, export=False)

    def add_sliding_window(self):
        if self.config.sw_scan_offset < 1:
            self.config.sw_scan_offset = 1
        chromatogram, _ = self.get_selected_chromatogram()
        tindex = np.arange(0, len(chromatogram), int(self.config.sw_scan_offset))
        self.data.clear()
        for i in tindex:
            t = chromatogram[i, 0]

            if self.config.time_start is not None:
                if t < self.config.time_start:
                    continue
            if self.config.time_end is not None:
                if t + self.config.sw_time_window > self.config.time_end:
                    continue

            data = self.get_data_from_times(t, t + self.config.sw_time_window)
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)
        pass

    def add_list_times(self, starts, ends):
        self.data.clear()
        for i, t in enumerate(starts):
            data = self.get_data_from_times(t, ends[i])
            self.data.add_data(data, name=str(t), attrs=self.attrs, export=False)


if __name__ == "__main__":
    inpath = "C:\\Python\\UniDec3\\unidec\\bin\\Example Data\\UniChrom\\SEC_Native_Herceptin.raw"
    eng = ChromEngine()
    eng.open_chrom(inpath)
    #eng.gen_html_report()
    pass

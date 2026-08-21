"""Standalone configuration for IsoDec deconvolution."""


class IsoDecConfig:
    """Parameters used by the native and Python IsoDec engines.

    The class intentionally remains mutable for compatibility with the
    configuration object historically supplied by UniDec.
    """

    def __init__(self):
        self.adductmass = 1.007276467
        self.verbose = 0
        self.filepath = ""
        self.batch_size = 32
        self.test_batch_size = 2048
        self.current_KD_round = 0
        self.activescan = -1
        self.activescanrt = -1
        self.activescanorder = -1
        self.meanpeakspacing_thresh = 0.01
        self.background_subtraction = 0

        self.mass_diff_c = 1.0033
        self.peakwindow = 80
        self.phaseres = 8
        self.matchtol = 5
        self.minpeaks = 3
        self.peakthresh = 0.0001
        self.css_thresh = 0.7
        self.maxshift = 3
        self.mzwindowlb = -1.05
        self.mzwindowub = 4.05
        self.plusoneintwindowlb = 0.1
        self.plusoneintwindowub = 0.6
        self.knockdown_rounds = 5
        self.min_score_diff = 0.1
        self.minareacovered = 0.20
        self.minusoneaszero = 1
        self.isotopethreshold = 0.01
        self.datathreshold = 0.05
        self.zscore_threshold = 0.95

        self.avgpeakmasses = 0
        self.report_multiple_monoisos = 1
        self.write_scans_without_precs = 1
        self.write_msalign = 0
        self.write_tsv = 1
        self.analyte_type = "Peptide"

    def set_scan_info(self, scan, reader=None):
        """Record scan metadata from a spectrum reader when available."""
        self.activescan = scan
        if reader is None:
            return
        self.activescanrt = reader.get_scan_time(scan)
        self.activescanorder = reader.get_ms_order(scan)
        self.adductmass = -1.007276467 if reader.polarity == "Negative" else 1.007276467


__all__ = ["IsoDecConfig"]

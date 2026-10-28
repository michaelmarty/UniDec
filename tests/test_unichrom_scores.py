import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import h5py
import numpy as np

from unidec.metaunidec.mudeng import MetaUniDec
from unidec.metaunidec.mudstruct import MetaDataSet
from unidec.modules.ChromEng import ChromEngine, chromatographic_peak_scores
from unidec.modules.peakstructure import Peaks


class TestChromatographicPeakScores(unittest.TestCase):
    def test_symmetric_monotonic_peak_scores_one(self):
        scores = chromatographic_peak_scores([0, 1, 2, 3, 4], [0, .5, 1, .5, 0])
        np.testing.assert_allclose(scores, [1, 1])

    def test_known_asymmetric_peak_uses_half_width_ratio(self):
        scores = chromatographic_peak_scores([0, 1, 2, 3, 4, 5], [0, .5, 1, .75, .5, 0])
        np.testing.assert_allclose(scores, [.5, 1])

    def test_shoulder_penalizes_monotonicity(self):
        scores = chromatographic_peak_scores([0, 1, 2, 3, 4, 5, 6], [0, .5, 1, .7, .8, .5, 0])
        np.testing.assert_allclose(scores, [1 / 3, 1 - .1 / 3])

    def test_uneven_axis_uses_interpolated_half_height_crossing(self):
        scores = chromatographic_peak_scores([0, 1, 2, 4, 7], [0, .25, 1, .5, 0])
        np.testing.assert_allclose(scores, [1 / 3, 1])

    def test_flat_top_uses_plateau_midpoint(self):
        scores = chromatographic_peak_scores([0, 1, 2, 3, 4, 5], [0, .5, 1, 1, .5, 0])
        np.testing.assert_allclose(scores, [1, 1])

    def test_scores_are_intensity_scale_invariant(self):
        xvalues = np.arange(7)
        intensities = np.array([0, .5, 1, .7, .8, .5, 0])
        np.testing.assert_allclose(
            chromatographic_peak_scores(xvalues, intensities),
            chromatographic_peak_scores(xvalues, intensities * 37),
        )

    def test_invalid_traces_score_zero(self):
        invalid = [
            ([], []),
            ([0], [1]),
            ([0, 1], [0, 1]),
            ([0, 1, 2], [0, 0, 0]),
            ([0, 1, 2], [0, np.nan, 0]),
            ([0, 1, 1], [0, 1, 0]),
            ([0, 1, 2], [0, 1]),
            ([0, 1, 2], [1, .5, 0]),
            ([0, 1, 2], [0, .5, 1]),
        ]
        for xvalues, intensities in invalid:
            with self.subTest(xvalues=xvalues, intensities=intensities):
                self.assertEqual(chromatographic_peak_scores(xvalues, intensities), (0., 0.))


class TestUniChromScoreIntegration(unittest.TestCase):
    @staticmethod
    def make_engine(path=None, masses=(100, 200), heights=(2, 1), dscores=(.8, .2)):
        engine = object.__new__(ChromEngine)
        engine.config = SimpleNamespace(dtsig=0, calc_unichrom_scores=True)
        engine.data = MetaDataSet(engine)
        engine.data.filename = path
        engine.data.var1 = np.arange(5, dtype=float)
        engine.data.peaks = np.column_stack((masses, heights, dscores)).astype(float)
        engine.pks = Peaks()
        engine.pks.add_peaks(engine.data.peaks, scores_included=True)
        for peak in engine.pks.peaks:
            peak.extracts = np.array([0, .5, 1, .5, 0])
        return engine

    def test_engine_combines_scores_and_calculates_uniscore_for_both_modes(self):
        for dtsig in (0, 2.5):
            with self.subTest(dtsig=dtsig):
                engine = self.make_engine()
                engine.config.dtsig = dtsig
                with patch.object(MetaUniDec, "pick_peaks"):
                    engine.pick_peaks()

                np.testing.assert_allclose(
                    [[p.base_dscore, p.chrom_symmetry_score, p.chrom_monotonicity_score, p.dscore]
                     for p in engine.pks.peaks],
                    [[.8, 1, 1, .8], [.2, 1, 1, .2]],
                )
                self.assertAlmostEqual(engine.pks.uniscore, .68)
                np.testing.assert_allclose(engine.data.peaks[:, 2], [.8, .2])

    def test_zero_height_peak_does_not_make_uniscore_nonfinite(self):
        engine = self.make_engine(masses=(100,), heights=(0,), dscores=(.8,))
        with patch.object(MetaUniDec, "pick_peaks"):
            engine.pick_peaks()
        self.assertEqual(engine.pks.uniscore, 0)

    def test_scores_are_disabled_by_default(self):
        engine = self.make_engine(masses=(100,), heights=(1,), dscores=(.8,))
        engine.config.calc_unichrom_scores = False
        engine.data.export_unichrom_scores = Mock()
        with patch.object(MetaUniDec, "pick_peaks"):
            engine.pick_peaks()
        self.assertEqual(engine.pks.peaks[0].dscore, .8)
        self.assertEqual(engine.data.unichrom_scores, [])
        engine.data.export_unichrom_scores.assert_not_called()

    def test_hdf5_scores_round_trip_and_old_files_remain_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "scores.hdf5")
            with h5py.File(path, "w") as hdf:
                peaks = hdf.require_group("peaks")
                peaks.create_dataset("peakdata", data=[[100, 2, .8]])
                peaks.create_dataset("extracts", data=[[0], [.5], [1], [.5], [0]])

            data = MetaDataSet(SimpleNamespace())
            data.filename = str(path)
            data.import_peaks()
            self.assertEqual(data.uniscore, 0)
            self.assertEqual(data.unichrom_scores.size, 0)

            data.export_unichrom_scores([[100, .8, 1, 1, .8]], .8)
            restored = MetaDataSet(SimpleNamespace())
            restored.filename = str(path)
            restored.import_peaks()
            np.testing.assert_allclose(restored.unichrom_scores, [[100, .8, 1, 1, .8]])
            self.assertEqual(restored.uniscore, .8)

    def test_refresh_false_is_idempotent(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory, "scores.hdf5"))
            engine = self.make_engine(path=path, masses=(100,), heights=(1,), dscores=(.8,))
            with patch.object(MetaUniDec, "pick_peaks"):
                engine.pick_peaks(refresh=False)
                first_score = engine.pks.peaks[0].dscore
                engine.pick_peaks(refresh=False)
            self.assertEqual(engine.pks.peaks[0].base_dscore, .8)
            self.assertEqual(engine.pks.peaks[0].dscore, first_score)

    def test_stale_component_rows_are_ignored(self):
        engine = self.make_engine(masses=(100,), heights=(1,), dscores=(.8,))
        engine.data.unichrom_scores = np.array([[999, .2, 1, 1, .2]])
        with patch.object(MetaUniDec, "pick_peaks"):
            engine.pick_peaks(refresh=False)
        self.assertEqual(engine.pks.peaks[0].base_dscore, .8)


if __name__ == "__main__":
    unittest.main()

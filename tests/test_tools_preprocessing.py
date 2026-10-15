import unittest
from types import SimpleNamespace

import numpy as np

from unidec import tools


class TestNormalization(unittest.TestCase):
    def test_normalize_preserves_axis_and_mutates_float_input(self):
        data = np.array([[10.0, 2.0], [11.0, 4.0], [12.0, 1.0]])
        result = tools.normalize(data)

        self.assertIs(result, data)
        np.testing.assert_allclose(result[:, 0], [10.0, 11.0, 12.0])
        np.testing.assert_allclose(result[:, 1], [0.5, 1.0, 0.25])

    def test_normalize_leaves_zero_and_singleton_inputs_finite(self):
        cases = (
            np.array([[10.0, 0.0], [11.0, 0.0]]),
            np.array([[10.0, 0.0]]),
        )
        for data in cases:
            with self.subTest(shape=data.shape):
                result = tools.normalize(data)
                np.testing.assert_array_equal(result[:, 1], np.zeros(len(data)))
                self.assertTrue(np.all(np.isfinite(result)))


class TestLinearizationAndDuplicates(unittest.TestCase):
    def test_linearize_all_documented_modes(self):
        data = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0], [5.0, 50.0]])
        cases = {
            0: np.array([[1.0, 0.0], [2.0, 20.0], [3.0, 30.0], [4.0, 0.0]]),
            1: np.array([[1.0, 0.0], [2.0, 35.0], [4.0, 15.0]]),
            3: np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]]),
            4: np.array([[1.0, 10.0], [2.0, 20.0], [4.0, 40.0]]),
        }
        for mode, expected in cases.items():
            with self.subTest(mode=mode):
                np.testing.assert_allclose(tools.linearize(data, 1.0, mode), expected)

    def test_linear_interpolation_reproduces_linear_signal_with_uneven_input(self):
        data = np.array([[0.0, 1.0], [0.7, 2.4], [2.2, 5.4], [4.0, 9.0]])
        result = tools.linearize(data, 1.0, 3)
        np.testing.assert_allclose(result[:, 1], 2.0 * result[:, 0] + 1.0)

    def test_linearize_empty_and_singleton_inputs(self):
        for data in (np.empty((0, 2)), np.array([[10.0, 3.0]])):
            for mode in (0, 1, 3, 4):
                with self.subTest(shape=data.shape, mode=mode):
                    result = tools.linearize(data, 1.0, mode)
                    np.testing.assert_array_equal(result, data)
                    self.assertIsNot(result, data)

    def test_remove_duplicates_sums_every_group_including_last(self):
        data = np.array([[1.0, 2.0], [1.0, 3.0], [2.0, 5.0], [3.0, 7.0], [3.0, 11.0]])
        expected = np.array([[1.0, 5.0], [2.0, 5.0], [3.0, 18.0]])
        np.testing.assert_allclose(tools.removeduplicates(data), expected)


class TestThresholdCropAndPipeline(unittest.TestCase):
    def test_crop_includes_both_boundaries(self):
        data = np.column_stack((np.arange(5.0), np.arange(5.0)))
        np.testing.assert_array_equal(tools.datachop(data, 1.0, 3.0), data[1:4])

    def test_threshold_keeps_value_equal_to_cutoff(self):
        data = np.array([[0.0, 0.9], [1.0, 1.0], [2.0, 1.1]])
        result = tools.intensitythresh(data, 1.0)
        np.testing.assert_allclose(result[:, 1], [0.0, 1.0, 1.1])

    @staticmethod
    def _config(**overrides):
        values = dict(
            minmz="", maxmz="", subbuff=0, smooth=0, mzbins=0,
            intthresh=0, subtype=0, detectoreffva=0, linflag=0,
            reductionpercent=0, smashflag=0, smashlist=[], datanorm=0,
        )
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_dataprep_crop_minimum_subtraction_normalization_and_threshold_order(self):
        raw = np.array([[0.0, 100.0], [1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 100.0]])
        original = raw.copy()
        config = self._config(minmz=1.0, maxmz=3.0, subbuff=1, subtype=0, datanorm=1, intthresh=0.5)

        result = tools.dataprep(raw, config, silent=True)

        np.testing.assert_array_equal(raw, original)
        np.testing.assert_allclose(result, [[1.0, 0.0], [2.0, 0.5], [3.0, 1.0]])

    def test_dataprep_without_crop_does_not_mutate_when_no_in_place_step_is_enabled(self):
        raw = np.array([[1.0, 2.0], [2.0, 4.0]])
        original = raw.copy()
        result = tools.dataprep(raw, self._config(), silent=True)
        np.testing.assert_array_equal(result, original)
        np.testing.assert_array_equal(raw, original)

    def test_dataprep_without_crop_preserves_raw_input_during_normalization(self):
        raw = np.array([[1.0, 2.0], [2.0, 4.0]])
        original = raw.copy()
        result = tools.dataprep(raw, self._config(datanorm=1), silent=True)
        np.testing.assert_allclose(result[:, 1], [0.5, 1.0])
        np.testing.assert_array_equal(raw, original)

    def test_detector_correction_matches_documented_efficiency_formula(self):
        data = np.array([[1000.0, 2.0], [2000.0, 4.0], [4000.0, 8.0]])
        original_x = data[:, 0].copy()
        voltage = 2.0
        efficiency = 1.0 - np.exp(-1620.0 * (voltage / original_x) ** 1.75)

        result = tools.detectoreff(data, voltage)

        np.testing.assert_array_equal(result[:, 0], original_x)
        np.testing.assert_allclose(result[:, 1], np.array([2.0, 4.0, 8.0]) / efficiency)

    def test_gaussian_smoothing_preserves_a_constant_signal(self):
        data = np.column_stack((np.arange(9.0), np.full(9, 3.5)))
        result = tools.gsmooth(data, 1.5)
        np.testing.assert_allclose(result[:, 1], np.full(9, 3.5), atol=1e-12)

    def test_smash_zeros_only_points_strictly_inside_window(self):
        data = np.column_stack((np.arange(5.0), np.ones(5)))
        result = tools.smash(data, midpoint=2.0, window=1.0)
        np.testing.assert_array_equal(result[:, 1], [1.0, 1.0, 0.0, 1.0, 1.0])

    def test_positive_percent_reduction_zeros_values_below_percentile_cutoff(self):
        data = np.column_stack((np.arange(5.0), np.arange(1.0, 6.0)))
        result = tools.remove_noise(data, percent=40)
        np.testing.assert_array_equal(result[:, 1], [0.0, 0.0, 3.0, 4.0, 5.0])

    def test_negative_percent_reduction_returns_detected_peaks_in_dataprep(self):
        raw = np.column_stack((np.arange(9.0), [0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0]))
        result = tools.dataprep(
            raw,
            self._config(reductionpercent=-1),
            silent=True,
        )
        np.testing.assert_array_equal(result, [[1.0, 1.0], [4.0, 2.0], [7.0, 1.0]])


class TestBackgroundSubtraction(unittest.TestCase):
    def setUp(self):
        x = np.arange(101.0)
        self.constant = np.column_stack((x, np.full_like(x, 5.0)))
        self.with_peak = self.constant.copy()
        self.with_peak[50, 1] = 10.0

    def test_constant_baseline_and_peak_for_each_algorithm(self):
        algorithms = {
            "minimum": lambda data: self._minimum_subtract(data),
            "line": lambda data: tools.datasimpsub(data, 5),
            "curved": lambda data: tools.datacompsub(data, 5),
            "relative_constant": lambda data: tools.intensitythresh_sub(data, 0.5),
            "polynomial": lambda data: tools.polynomial_background_subtract(data, 1),
            "savgol": lambda data: tools.savgol_background_subtract(data, 5),
        }
        for name, subtract in algorithms.items():
            with self.subTest(name=name, signal="constant"):
                result = subtract(self.constant.copy())
                expected = 2.5 if name == "relative_constant" else 0.0
                np.testing.assert_allclose(result[:, 1], expected, atol=1e-12)
                self.assertTrue(np.all(np.isfinite(result)))
            with self.subTest(name=name, signal="peak"):
                result = subtract(self.with_peak.copy())
                self.assertAlmostEqual(result[0, 1], 0.0, places=10)
                self.assertAlmostEqual(result[-1, 1], 0.0, places=10)
                if name == "savgol":
                    self.assertEqual(np.argmax(result[:, 1]), 50)
                    self.assertGreater(result[50, 1], 0.0)
                    self.assertLess(result[50, 1], 5.0)
                else:
                    self.assertAlmostEqual(result[50, 1], 5.0, places=10)

    @staticmethod
    def _minimum_subtract(data):
        data[:, 1] -= np.amin(data[:, 1])
        return data


if __name__ == "__main__":
    unittest.main()

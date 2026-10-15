import unittest

import numpy as np

from unidec import tools


class TestDataExtraction(unittest.TestCase):
    def setUp(self):
        self.data = np.array([[10.0, 0.0], [11.0, 1.0], [12.0, 4.0], [13.0, 1.0], [14.0, 0.0]])

    def test_hand_computable_windowed_modes(self):
        expected = {
            0: 4.0,
            1: 4.0,
            2: 5.0,
            3: 12.0,
            4: 12.0,
            5: 12.0,
            6: 12.0,
            7: 12.0,
            8: 12.0,
            9: 12.0,
            10: 12.0,
            12: 5.0,
        }
        for mode, value in expected.items():
            with self.subTest(mode=mode):
                self.assertAlmostEqual(tools.data_extract(self.data, 12.0, mode, window=1.5), value)

    def test_thresholded_center_modes_without_window_use_filtered_length(self):
        for mode in (3, 4, 5, 6, 7, 8, 9, 10):
            with self.subTest(mode=mode):
                self.assertEqual(tools.data_extract(self.data, 12.0, mode), 12.0)

    def test_unknown_mode_returns_zero(self):
        self.assertEqual(tools.data_extract(self.data, 12.0, 999, window=1.0), 0)

    def test_automatic_width_area_matches_analytic_unit_fwhm_gaussian(self):
        x = np.linspace(-5.0, 5.0, 1001)
        intensity = np.exp(-4.0 * np.log(2.0) * x ** 2)
        data = np.column_stack((x, intensity))
        expected_area = np.sqrt(np.pi / np.log(2.0)) / 2.0
        self.assertAlmostEqual(
            tools.data_extract(data, 0.0, extract_method=11, window=4.0),
            expected_area,
            places=10,
        )


class TestExtractNormalization(unittest.TestCase):
    def test_all_normalization_modes_and_mutation_contract(self):
        source = np.array([[1.0, 0.0], [3.0, 1.0]])
        cases = {
            0: np.array([[1.0, 0.0], [3.0, 1.0]]),
            1: np.array([[1.0 / 3.0, 0.0], [1.0, 1.0]]),
            2: np.array([[0.25, 0.0], [0.75, 1.0]]),
            3: np.array([[1.0, 0.0], [1.0, 1.0 / 3.0]]),
            4: np.array([[1.0, 0.0], [0.75, 0.25]]),
        }
        for mode, expected in cases.items():
            with self.subTest(mode=mode):
                grid = source.copy()
                result = tools.normalize_extracts(grid, mode)
                self.assertIs(result, grid)
                np.testing.assert_allclose(result, expected)

    def test_zero_rows_and_columns_remain_finite(self):
        for mode in (1, 2, 3, 4):
            with self.subTest(mode=mode):
                grid = np.zeros((2, 2))
                np.testing.assert_array_equal(tools.normalize_extracts(grid, mode), np.zeros((2, 2)))
                self.assertTrue(np.all(np.isfinite(grid)))


if __name__ == "__main__":
    unittest.main()

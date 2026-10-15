import unittest

import numpy as np

from unidec import tools


class TestPeakDetection(unittest.TestCase):
    def test_isolated_peaks_and_edges(self):
        data = np.array([[0.0, 3.0], [1.0, 0.0], [2.0, 2.0], [3.0, 0.0], [4.0, 4.0]])
        expected = np.array([[0.0, 3.0], [2.0, 2.0], [4.0, 4.0]])
        np.testing.assert_array_equal(tools.peakdetect(data, window=1, norm=False), expected)

    def test_equal_height_neighbors_select_first(self):
        data = np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 2.0], [3.0, 0.0]])
        np.testing.assert_array_equal(tools.peakdetect(data, window=1, norm=False), [[1.0, 2.0]])

    def test_threshold_is_strict(self):
        data = np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 0.0], [3.0, 4.0], [4.0, 0.0]])
        np.testing.assert_array_equal(
            tools.peakdetect(data, window=1, threshold=0.5, norm=True),
            [[3.0, 4.0]],
        )

    def test_window_suppresses_smaller_neighbor(self):
        data = np.array([[0.0, 0.0], [1.0, 3.0], [2.0, 0.0], [3.0, 4.0], [4.0, 0.0]])
        np.testing.assert_array_equal(tools.peakdetect(data, window=1, norm=False), [[1.0, 3.0], [3.0, 4.0]])
        np.testing.assert_array_equal(tools.peakdetect(data, window=2, norm=False), [[3.0, 4.0]])

    def test_empty_input(self):
        self.assertEqual(tools.peakdetect(np.empty((0, 2))).shape, (0,))

    def test_singleton_and_all_zero_inputs(self):
        np.testing.assert_array_equal(tools.peakdetect(np.array([[10.0, 2.0]])), [[10.0, 2.0]])
        self.assertEqual(tools.peakdetect(np.array([[10.0, 0.0], [11.0, 0.0]])).shape, (0,))


if __name__ == "__main__":
    unittest.main()

import unittest

import numpy as np
from sqlalchemy.sql.sqltypes import NULLTYPE

from unidec.IsoDec.datatools import fastnearest, fastwithin_abstol, fastpeakdetect, fastcalc_FWHM


class FastNearest(unittest.TestCase):
    def test_accuracy_exact(self):
        numbers = np.array([i for i in range(1, 20)], dtype=np.float32)
        target = 10.0
        val = fastnearest(numbers, target)
        self.assertEqual(val, 9)

    def test_nearest_index(self):
        list_one = [1, 2, 3, 4, 5, 6, 7,12,13,14,15,16]
        target = 9.0
        # Should return index at value 7(6), before index 8(12)
        numbers = np.array(list_one, dtype=np.float32)
        val = fastnearest(numbers, target)
        self.assertEqual(val, 6)

    def test_higher_or_lower(self):
        list_one = [1, 2, 3, 4, 5, 6, 7,9,10,11,12,13]
        numbers = np.array(list_one, dtype=np.float32)
        target = 9
        # should return higher in a draw
        val = fastnearest(numbers, target)
        self.assertEqual(val, 7)

    def test_if_empty(self):
        list = np.array([], dtype = np.float32)
        target = 9
        val = fastnearest(list, target)
        self.assertIsNone(val)

class FastWithinAbsTol(unittest.TestCase):
    def test_fast_within_abs_total(self):
        nums = np.array([100.0, 100.5, 101.0, 101.5, 102.0, 102.5, 103.0, 103.5, 104.0], dtype=np.float32)
        target = 103.5
        result = fastwithin_abstol(nums, target, 0.5)
        # should hit index 7, append 8 and 6
        self.assertEqual(result, [7, 8, 6])

    def test_falling_off_edge(self):
        nums = np.array([100.0, 100.5, 101.0, 101.5, 102.0, 102.5, 103.0, 103.5, 104.0], dtype=np.float32)
        target = 104
        result = fastwithin_abstol(nums, target, 0.05)
        self.assertEqual(result, [8])

    def test_is_empty(self):
        nums = []
        target = 100
        result = fastwithin_abstol(nums, target, 50)
        self.assertTrue(np.array_equal(result, nums))


class FastPeakDetect(unittest.TestCase):
    #if not 2d, should exit out, and if list is empty
    def test_data_is_empty(self):

        nums = np.empty((0,2))
        result = fastpeakdetect(nums)
        self.assertTrue(np.array_equal(result, nums))

#find max intensity then traverse towards the edges til you find a lower value
# class FastCalcFWHM(unittest.TestCase):
#     def test_is_valid(self):
#         expected = (0.5, [2.0, 2.5])
#         peak = 8
#         matrix = [[2.0,3.5,2],
#                    [4,8.0,10.0],
#                    [2.5,3.5,4.5]]
#         result = fastcalc_FWHM(matrix, peak)
#         self.assertEqual(expected, result)
#
#



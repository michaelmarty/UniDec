import unittest

import numpy as np
import pandas as pd

from unidec.modules import matchtools, peakstructure


def make_peaks(masses, heights=None, integrals=None):
    if heights is None:
        heights = np.ones(len(masses))
    peaks = peakstructure.Peaks()
    peaks.add_peaks(np.column_stack((masses, heights)))
    if integrals is not None:
        for peak, integral in zip(peaks.peaks, integrals):
            peak.integral = integral
    return peaks


class TestModificationParsing(unittest.TestCase):
    def test_fixed_variable_and_global_modifications(self):
        fixed = pd.DataFrame([
            {"Mass": 10.0, "Name": "A", "Number": 2},
            {"Mass": 5.0, "Name": "B"},
        ])
        mass, name = matchtools.parse_fmoddf(fixed)
        self.assertEqual(mass, 25.0)
        self.assertEqual(name, "+2.0[A]+1[B]")

        variable = pd.DataFrame({"Mass": [0.0, 3.0], "Name": ["Unmodified", "Ox"]})
        masses, names = matchtools.parse_vmoddf(variable)
        np.testing.assert_array_equal(masses, [0.0, 3.0])
        np.testing.assert_array_equal(names, ["Unmodified", "Ox"])

        for value, expected in (("7.5", 7.5), (np.nan, 0), ("bad", 0)):
            with self.subTest(value=value):
                self.assertEqual(matchtools.parse_global_fixed_mod({"Global Fixed Mod": value}), expected)

    def test_fixed_modifications_apply_to_selected_sequence(self):
        fixed = pd.DataFrame([{"Mass": 10.0, "Name": "Tag", "Number": 2}])
        row = {
            "Sequence 1": 100.0,
            "Sequence 2": 200.0,
            "Apply Fixed Mods": "Seq1",
            "Correct": "Seq1+Seq2",
        }
        masses, labels = matchtools.calc_pairs(row, fmoddf=fixed)
        np.testing.assert_array_equal(masses, [320.0])
        np.testing.assert_array_equal(labels, ["Correct"])


class TestScientificMatching(unittest.TestCase):
    @staticmethod
    def _match(mass, tolerance=1.0, row=None, favor="Closest"):
        if row is None:
            row = {"Correct": 100.0}
        peaks = make_peaks([mass], [2.0])
        result = matchtools.UPP_check_peaks(row.copy(), peaks, tolerance, favor=favor)
        return result, peaks.peaks[0]

    def test_tolerance_is_strict_at_boundary(self):
        cases = ((100.999, "Correct", 100.0), (101.0, "Unknown", 0), (101.001, "Unknown", 0))
        for mass, label, matched_mass in cases:
            with self.subTest(mass=mass):
                _, peak = self._match(mass)
                self.assertEqual(peak.label, label)
                self.assertEqual(peak.match, matched_mass)

    def test_competing_candidates_honor_closest_and_favored_label(self):
        row = {"Correct": 100.0, "Incorrect": 100.5}
        _, closest = self._match(100.4, row=row)
        self.assertEqual(closest.label, "Incorrect")
        self.assertEqual(closest.match, 100.5)

        _, favored = self._match(100.4, row=row, favor="Correct")
        self.assertEqual(favored.label, "Correct")
        self.assertEqual(favored.match, 100.0)
        self.assertEqual(favored.numberalts, 2)

    def test_no_matches_and_empty_peak_collection(self):
        result, peak = self._match(200.0)
        self.assertEqual(peak.label, "Unknown")
        self.assertEqual(result["unmatched %"], 100.0)

        empty = make_peaks([])
        result = matchtools.UPP_check_peaks({"Correct": 100.0}, empty, 1.0)
        self.assertEqual(result["Total Height"], 0)
        self.assertEqual(result["correct %"], 0)
        self.assertEqual(result["incorrect %"], 0)
        self.assertEqual(result["unmatched %"], 0)

    def test_variable_and_global_modifications_shift_match_mass(self):
        variable = pd.DataFrame({"Mass": [0.0, 3.0], "Name": ["", "Ox"]})
        peaks = make_peaks([110.0], [1.0])
        row = {"Correct": 100.0, "Global Fixed Mod": 7.0}
        matchtools.UPP_check_peaks(row, peaks, 0.1, vmoddf=variable)
        self.assertEqual(peaks.peaks[0].match, 110.0)
        self.assertEqual(peaks.peaks[0].label, "Correct+Ox")


class TestBispecificAndDar(unittest.TestCase):
    def test_bispecific_classification_hand_calculation(self):
        row = {
            "BsAb (Correct)": "Seq1",
            "LC1 Mispair (Incorrect)": "Seq2",
            "LC2 Mispair (Incorrect)": "Seq3",
            "bsab (correct) Height": 2.0,
            "lc1 mispair (incorrect) Height": 1.0,
            "lc2 mispair (incorrect) Height": 1.0,
        }
        result = matchtools.calc_bispecific_correct(row)
        self.assertEqual(result["BsAb Pairing Calculated (%)"], 25.0)
        self.assertEqual(result["Light Chain Scrambled (%)"], 25.0)

    def test_dar_uses_height_or_integral_and_returns_minus_one_without_matches(self):
        peaks = make_peaks(
            [100.0, 110.0, 120.0],
            heights=[1.0, 2.0, 1.0],
            integrals=[1.0, 1.0, 2.0],
        )
        self.assertEqual(matchtools.dar_calc(peaks, 100.0, 10.0, 0, 3, 0.1), 1.0)
        self.assertEqual(matchtools.dar_calc(peaks, 100.0, 10.0, 0, 3, 0.1, integrate=True), 1.25)
        self.assertEqual(matchtools.dar_calc(peaks, 1000.0, 10.0, 0, 3, 0.1), -1)
        self.assertEqual(matchtools.dar_calc(make_peaks([]), 100.0, 10.0, 0, 3, 0.1), -1)


if __name__ == "__main__":
    unittest.main()

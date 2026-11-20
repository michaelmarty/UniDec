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

    def test_sequence_chemistry_uses_average_disulfide_and_pyroglutamate_masses(self):
        """Selected cysteines lose H; N-terminal E/Q lose H2O/NH3 case-insensitively."""
        for sequence, kwargs, expected in (
            ("AC", {"all_cyst_ox": True}, 191.22),
            ("ACC", {"all_cyst_ox": True}, 293.36),
            ("EAAA", {"pyroglu": True}, 342.35),
            ("eaaa", {"pyroglu": True}, 342.35),
            ("QAAA", {"pyroglu": True}, 342.35),
            ("qaaa", {"pyroglu": True}, 342.35),
        ):
            with self.subTest(sequence=sequence):
                self.assertAlmostEqual(matchtools.calc_pep_mass(sequence, **kwargs), expected, delta=0.01)

        masses, labels = matchtools.calc_pairs({
            "Sequence 1": "qac",
            "Disulfides Oxidized": "Seq1",
            "PyroGlu": True,
            "Correct": "Seq1",
        })
        np.testing.assert_allclose(masses, [302.33], atol=0.01, rtol=0)
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

    def test_reviewed_upp_sequence_modification_and_quantitation_contract(self):
        row = {
            "Sequence 1": "QAC",
            "Sequence 2": "MCK",
            "Sequence 3": "GAS",
            "Disulfides Oxidized": "Seq1+Seq2",
            "PyroGlu": True,
            "Apply Fixed Mods": "Seq1",
            "Correct": "Seq1+Seq2",
            "Incorrect": "2*Seq1",
            "Ignore": "Seq3",
        }
        fixed = pd.DataFrame([{"Mass": 57.0513, "Name": "Carbamidomethyl", "Number": 1}])
        variable = pd.DataFrame({"Mass": [0.0, 15.9994], "Name": ["", "Ox"]})

        candidate_masses, candidate_labels = matchtools.calc_pairs(row, fmoddf=fixed)
        np.testing.assert_array_equal(candidate_labels, ["Correct", "Incorrect", "Ignore"])
        np.testing.assert_allclose(candidate_masses, [738.8913, 718.7626, 233.22], atol=0.02, rtol=0)

        peak_masses = [738.8913, 734.7620, 233.22, 900.0]
        heights = [6.0, 2.0, 4.0, 2.0]
        integrals = [3.0, 4.0, 5.0, 8.0]
        height_peaks = make_peaks(peak_masses, heights, integrals)
        height_result = matchtools.UPP_check_peaks(
            row.copy(), height_peaks, 0.02, vmoddf=variable, fmoddf=fixed
        )
        self.assertEqual([p.label for p in height_peaks.peaks],
                         ["Correct", "Incorrect+Ox", "Ignore", "Unknown"])
        self.assertAlmostEqual(height_result["Total Height"], 10.0)
        self.assertAlmostEqual(height_result["correct %"], 60.0, delta=0.01)
        self.assertAlmostEqual(height_result["incorrect %"], 20.0, delta=0.01)
        self.assertAlmostEqual(height_result["unmatched %"], 20.0, delta=0.01)
        self.assertAlmostEqual(height_result["correct % Matched Only"], 75.0, delta=0.01)

        integral_peaks = make_peaks(peak_masses, heights, integrals)
        integral_result = matchtools.UPP_check_peaks(
            row.copy(), integral_peaks, 0.02, vmoddf=variable, fmoddf=fixed, integrate=True
        )
        self.assertAlmostEqual(integral_result["Total Height"], 15.0)
        self.assertAlmostEqual(integral_result["correct %"], 20.0, delta=0.01)
        self.assertAlmostEqual(integral_result["incorrect %"], 26.6666667, delta=0.01)
        self.assertAlmostEqual(integral_result["unmatched %"], 53.3333333, delta=0.01)
        self.assertAlmostEqual(integral_result["correct % Matched Only"], 42.8571429, delta=0.01)


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

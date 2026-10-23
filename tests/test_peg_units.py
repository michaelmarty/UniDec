"""Regression checks for the neutral-mass PEG600 screen."""

import contextlib
import importlib.util
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "PEG_Units.py"
spec = importlib.util.spec_from_file_location("peg_units", SCRIPT)
peg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(peg)


class TestPEGUnits(unittest.TestCase):
    def spectrum(self, masses, intensities=None):
        if intensities is None:
            intensities = [100] * len(masses)
        return pd.DataFrame({"mass": masses, "intensity": intensities})

    def test_free_peg600_neutral_masses_and_repeat_counts(self):
        masses = peg.WATER_MASS + np.arange(10, 18) * peg.REPEAT_MASS
        _, series, members = peg.find_peg_series(self.spectrum(masses))
        self.assertEqual(len(series), 1)
        self.assertTrue(series.iloc[0].PEG_candidate)
        self.assertEqual(members.free_PEG_n_if_mass_matches.tolist(), list(range(10, 18)))
        self.assertEqual(members.relative_repeat_step.tolist(), list(range(8)))
        self.assertNotIn("charge", members.columns)

    def test_unknown_end_group_still_detected_without_forced_n(self):
        masses = 436.228 + np.arange(7) * peg.REPEAT_MASS
        _, series, members = peg.find_peg_series(self.spectrum(masses))
        self.assertTrue(series.iloc[0].PEG_candidate)
        self.assertTrue(members.free_PEG_n_if_mass_matches.isna().all())

    def test_high_mass_peg_not_called_peg600(self):
        masses = 1492.75 + np.arange(10) * peg.REPEAT_MASS
        _, series, _ = peg.find_peg_series(self.spectrum(masses))
        self.assertEqual(len(series), 1)
        self.assertFalse(series.PEG_candidate.any())

    def test_same_high_mass_series_is_a_peg2000_candidate(self):
        masses = 1492.75 + np.arange(30) * peg.REPEAT_MASS
        _, series, members = peg.find_peg_series(self.spectrum(masses),
                                                mass_window=peg.PEG2000_MASS_WINDOW)
        self.assertTrue(series.iloc[0].PEG_candidate)
        self.assertGreater(series.iloc[0].window_peak_count, 20)
        self.assertTrue(members.free_PEG_n_if_mass_matches.isna().all())

    def test_isolated_and_two_peak_matches_are_not_a_distribution(self):
        _, series, members = peg.find_peg_series(self.spectrum([500, 544.02621474849, 600.5]))
        self.assertTrue(series.empty)
        self.assertTrue(members.empty)

    def test_neutral_axis_does_not_search_22_da_charge_spacing(self):
        masses = [500, 522.013107374245, 544.02621474849]
        _, series, _ = peg.find_peg_series(self.spectrum(masses))
        self.assertTrue(series.empty)

    def test_drift_not_accepted_and_subseries_not_duplicated(self):
        masses = 500 + np.arange(6) * (peg.REPEAT_MASS + 0.006)
        _, series, _ = peg.find_peg_series(self.spectrum(masses))
        self.assertTrue(series.empty)
        _, series, _ = peg.find_peg_series(self.spectrum(500 + np.arange(6) * peg.REPEAT_MASS))
        self.assertEqual(series.peak_count.tolist(), [6])

    def test_three_members_must_be_in_screening_window(self):
        _, series, _ = peg.find_peg_series(self.spectrum(740 + np.arange(5) * peg.REPEAT_MASS))
        self.assertEqual(series.iloc[0].window_peak_count, 2)
        self.assertFalse(series.iloc[0].PEG_candidate)

    def test_intensity_cutoff_is_adjustable(self):
        masses = [436.228, 480.25421474849, 524.28042949698, 1000]
        spectrum = self.spectrum(masses, [2, 2, 2, 10000])
        _, series, _ = peg.find_peg_series(spectrum)
        self.assertTrue(series.empty)
        _, series, _ = peg.find_peg_series(spectrum, min_relative_intensity=0.0001)
        self.assertTrue(series.PEG_candidate.any())

    def test_freestyle_and_plain_csv_formats(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.csv"
            for delimiter in [",", "\t", ";"]:
                path.write_text("SPECTRUM - MS\nexample.raw\nData points: 3\n"
                                + delimiter.join(["Mass", "Intensity"]) + "\n"
                                + delimiter.join(["500", "100"]) + "\n"
                                + delimiter.join(["550", "200"]) + "\n", encoding="utf-8-sig")
                self.assertEqual(peg.read_spectrum(path).mass.tolist(), [500, 550])
            path.write_text("500,100\n550,200\n", encoding="utf-8")
            self.assertEqual(len(peg.read_spectrum(path)), 2)
            path.write_text("Mass,Intensity\n500,nan\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid rows"):
                peg.read_spectrum(path)

    def test_editor_run_uses_saved_input_and_writes_screen_outputs(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.csv"
            self.spectrum(436.228 + np.arange(7) * peg.REPEAT_MASS).to_csv(path, index=False)
            with patch.object(peg, "INPUT_FILE", str(path)), patch("sys.argv", [str(SCRIPT), "--peg", "600"]):
                with contextlib.redirect_stdout(io.StringIO()):
                    peg.main()
            output = Path(directory) / "spectrum_PEG600_screen"
            candidates = pd.read_csv(output / "PEG600_candidates.csv")
            self.assertEqual(len(candidates), 1)
            for filename in ["PEG600_members.csv", "PEG_all_series.csv", "PEG_all_series_members.csv",
                             "PEG600_window_peaks.csv", "PEG600_report.txt", "PEG600_screen.png"]:
                self.assertGreater((output / filename).stat().st_size, 0)

    def test_no_candidates_still_saves_plot(self):
        spectrum = self.spectrum([500])
        _, series, members = peg.find_peg_series(spectrum)
        fig = peg.plot_screen(spectrum, series, members, peg.PEG600_MASS_WINDOW)
        peg.plt.close(fig)

    def test_editor_default_peg2000_writes_grade_specific_outputs(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.csv"
            masses = peg.WATER_MASS + np.arange(43, 48) * peg.REPEAT_MASS
            self.spectrum(masses, [10, 20, 100, 30, 10]).to_csv(path, index=False)
            with patch.object(peg, "INPUT_FILE", str(path)), patch("sys.argv", [str(SCRIPT)]):
                with contextlib.redirect_stdout(io.StringIO()):
                    peg.main()
            output = Path(directory) / "spectrum_PEG2000_screen"
            candidates = pd.read_csv(output / "PEG2000_candidates.csv")
            members = pd.read_csv(output / "PEG2000_members.csv")
            self.assertEqual(candidates.peak_count.tolist(), [5])
            self.assertAlmostEqual(candidates.iloc[0].window_mode_mass, masses[2])
            self.assertEqual(members.free_PEG_n_if_mass_matches.tolist(), list(range(43, 48)))
            self.assertTrue((output / "PEG2000_screen.png").exists())
            self.assertIn("PEG2000 screen", (output / "PEG2000_report.txt").read_text())


if __name__ == "__main__":
    unittest.main()

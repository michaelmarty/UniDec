import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis/entropy.py"
SPEC = importlib.util.spec_from_file_location("adduct_entropy", SCRIPT)
entropy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(entropy)
import SI_7_8 as strength


class TestEntropyByInstrument(unittest.TestCase):
    def tearDown(self):
        entropy.plt.close("all")

    def source(self):
        return pd.DataFrame([
            {"Dataset": dataset, "Ontology": ontology, "Polarity": polarity,
             "Adduct type": adduct, "Metabolite name": f"{ontology} {i}|chains",
             "MS/MS spectrum": "100:1" if i % 2 == 0 else "100:1 200:1 300:1",
             "Overall Match Status": accepted, "Validation Status": ""}
            for dataset in entropy.datasets for ontology in ["PC", "LPC", "EtherPC"]
            for polarity in ["positive", "negative"] for adduct in ["[M+H]+", "[M+Na]+"]
            for accepted in [True, False] for i in range(3)
        ])

    def test_calls_requested_calculator_on_normalized_merged_spectrum(self):
        with patch.object(entropy, "calc_spectral_entropy", wraps=entropy.calc_spectral_entropy) as calculate:
            value, peaks, status = entropy.spectral_entropy("100:1 100:2 200:3 300:0")
        self.assertEqual(status, "OK")
        self.assertEqual(peaks, 2)
        self.assertAlmostEqual(value, np.log(2), places=6)
        np.testing.assert_allclose(calculate.call_args.args[0], [[100, .5], [200, .5]])
        self.assertAlmostEqual(entropy.spectral_entropy("100:1e308 200:1e308")[0], np.log(2), places=6)
        self.assertEqual(entropy.spectral_entropy("100:5"), (0, 1, "OK"))

    def test_invalid_spectra_are_missing_not_zero(self):
        for spectrum in [None, np.nan, "", "100:0", "100:-1", "100:nan", "inf:1", "bad", "100:1:2"]:
            with self.subTest(spectrum=spectrum):
                value, _, status = entropy.spectral_entropy(spectrum)
                self.assertTrue(np.isnan(value))
                self.assertNotEqual(status, "OK")

    def test_grouping_and_population_match_si_7_8_without_class_collapse(self):
        audit = entropy.calculate_entropies(self.source())
        duplicate = audit.iloc[[0]].copy()
        duplicate[entropy.ENTROPY] = 2
        duplicate["Metabolite name"] = " pc 0 |other"
        provisional = audit.iloc[[0]].copy()
        provisional["Metabolite name"] = "PC provisional"
        provisional["Validation Status"] = "provisional"
        audit = pd.concat([audit, duplicate, provisional], ignore_index=True)
        for all_ids, excluded in [(False, False), (True, True), (True, False)]:
            actual = entropy.prepare_entropy(audit, all_ids, excluded)
            reference = strength.prepare_strength(audit.assign(**{"Class Strength": audit[entropy.ENTROPY]}),
                                                 all_ids=all_ids, excluded_only=excluded, allow_empty=True)
            columns = ["Plot Class", "Ontology", "Polarity", "Adduct type", "Lipid ID", "Dataset",
                       "Feature rows", "Accepted rows", "Excluded rows"]
            pd.testing.assert_frame_equal(actual[columns], reference[columns])
            np.testing.assert_allclose(actual[entropy.ENTROPY], reference["Class Strength"])
            self.assertEqual(set(actual.Ontology), {"PC", "LPC", "EtherPC"})
        accepted = entropy.prepare_entropy(audit)
        row = accepted.loc[accepted["Dataset"].eq("OTOT") & accepted["Ontology"].eq("PC")
                           & accepted["Polarity"].eq("positive") & accepted["Adduct type"].eq("[M+H]+")
                           & accepted["Lipid Name"].eq("PC 0")].iloc[0]
        self.assertEqual(row[entropy.ENTROPY], 1)

    def test_panel_count_rule_and_valid_zero_entropies(self):
        points = entropy.prepare_entropy(entropy.calculate_entropies(self.source()))
        retained, _ = entropy.comparable_points(points, min_lipids=3, remove_zero_groups=False)
        self.assertEqual(len(retained), len(points))
        none, _ = entropy.comparable_points(points, min_lipids=4, remove_zero_groups=False)
        self.assertTrue(none.empty)
        zero = points.copy()
        zero[entropy.ENTROPY] = 0
        retained, _ = entropy.comparable_points(zero, min_lipids=3, remove_zero_groups=False)
        self.assertEqual(len(retained), len(zero))
        missing = points.loc[points.Dataset.ne(entropy.datasets[-1])]
        retained, _ = entropy.comparable_points(missing, min_lipids=3, remove_zero_groups=False)
        self.assertTrue(retained.empty)

    def test_raw_statistics_and_legacy_strength_defaults_agree(self):
        points = entropy.prepare_entropy(entropy.calculate_entropies(self.source()))
        points[strength.NORMALIZED] = points[entropy.ENTROPY]
        original = strength.instrument_statistics(points, all_pairs=True)
        raw = strength.instrument_statistics(points, all_pairs=True, value_column=entropy.ENTROPY)
        pd.testing.assert_frame_equal(raw, original)
        self.assertEqual(len(raw), 3 * 2 * 2 * 10)
        self.assertTrue(raw["Mean difference (dataset - reference)"].abs().lt(1e-12).all())
        # The new optional raw mode must not change default normalized plot axes.
        legacy = strength.plot_strength(points)
        for fig in legacy.values():
            for ax in fig.axes:
                if ax.axison:
                    np.testing.assert_allclose(ax.get_yticks(), [0, .5, 1])
        points[entropy.ENTROPY] *= 5
        raw_figures = strength.plot_strength(points, value_column=entropy.ENTROPY,
                                             ylabel="Entropy (nats)", normalized_axis=False)
        self.assertTrue(all(fig.axes[0].get_ylim()[1] > 5 for fig in raw_figures.values()))

    def test_cli_exports_both_populations(self):
        with TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "input.csv"
            self.source().to_csv(source, index=False)
            with patch("matplotlib.figure.Figure.savefig") as save:
                entropy.main(["--input", str(source), "--output-dir", str(directory / "out"),
                              "--min-lipids", "3", "--no-show"])
            self.assertEqual(save.call_count, 8)
            for population in ["Accepted_IDs", "Excluded_IDs"]:
                path = directory / "out" / population / f"{entropy.OUTPUT_BASENAME}_Pairwise_Statistics.csv"
                self.assertEqual(len(pd.read_csv(path)), 120)
            self.assertEqual(entropy.plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Analysis"
SPEC = importlib.util.spec_from_file_location(
    "entropy_sharing", SCRIPT_DIR / "SI_Spectral_Entropy_Sharing_Correlation.py"
)
entropy = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(SCRIPT_DIR), *sys.path]):
    SPEC.loader.exec_module(entropy)


class TestSpectralEntropySharing(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_entropy_definition_and_intensity_scale(self):
        for spectrum, expected, peaks in [
            ("100:5", 0, 1),
            ("100:1 200:1", np.log(2), 2),
            ("100:1e308 200:1e308", np.log(2), 2),
            ("100:1 100:2 200:3 300:0", np.log(2), 2),
            ("100:1 200:3", -0.25 * np.log(0.25) - 0.75 * np.log(0.75), 2),
        ]:
            with self.subTest(spectrum=spectrum):
                result, count, status = entropy.spectral_entropy(spectrum)
                self.assertAlmostEqual(result, expected)
                self.assertEqual(count, peaks)
                self.assertEqual(status, "OK")

    def test_missing_and_invalid_spectra_are_not_zero(self):
        for spectrum in [None, np.nan, "", "100:0", "100:-1", "100:nan", "inf:1", "bad", "100:1:2"]:
            with self.subTest(spectrum=spectrum):
                value, _, status = entropy.spectral_entropy(spectrum)
                self.assertTrue(np.isnan(value))
                self.assertNotEqual(status, "OK")

    def fixture(self):
        first, second = entropy.datasets[:2]
        rows = [
            (first, "PC", "PC A|adduct1", "100:1", "positive", True),
            (first, "PC", "PC A|adduct2", "100:1", "positive", True),
            (second, "PC", "PC A|feature", "100:1 200:1", "positive", True),
            (first, "PC", "PC B", None, "positive", True),
            (first, "PE", "PE A", "100:1 200:1 300:1 400:1", "positive", True),
            (first, "PC", "PC B", "100:1", "negative", True),
            (second, "PC", "PC B", "100:1", "positive", False),
        ]
        return pd.DataFrame(rows, columns=[
            "Dataset", "Ontology", "Metabolite name", "MS/MS spectrum", "Polarity", "Overall Match Status"
        ])

    def test_missing_entropy_keeps_sharing_denominator_and_polarities_separate(self):
        source = entropy.polarity_subsets(self.fixture())["positive"]
        pooled, tables, audit, lipids = entropy.entropy_sharing_tables(source, "positive")
        original = entropy.class_sharing_table(source)
        pd.testing.assert_frame_equal(pooled[original.columns], original.rename_axis(None, axis=1))
        for dataset, original in entropy.dataset_sharing_tables(source).items():
            pd.testing.assert_frame_equal(tables[dataset][original.columns], original)
        self.assertEqual(pooled.loc["PC", "Pooled unique lipids"], 2)
        self.assertEqual(pooled.loc["PC", "Shared lipids (%)"], 50)
        self.assertEqual(pooled.loc["PC", "Accepted ID rows"], 4)
        self.assertEqual(pooled.loc["PC", "Rows with valid entropy"], 3)
        self.assertEqual(pooled.loc["PC", "Entropy row coverage (%)"], 75)
        self.assertEqual(pooled.loc["PC", entropy.ENTROPY], 0)
        self.assertEqual(len(audit), 5)
        self.assertEqual(len(lipids), 4)
        absent = tables[entropy.datasets[-1]]
        self.assertTrue(absent[entropy.ENTROPY].isna().all())
        self.assertTrue(absent["Shared lipids (%)"].isna().all())
        self.assertTrue(absent["Accepted ID rows"].eq(0).all())
        self.assertTrue(absent["Entropy row coverage (%)"].isna().all())

    def test_class_medians_use_all_rows_without_collapsing_lipids(self):
        source = entropy.polarity_subsets(self.fixture())["positive"]
        extra = source.iloc[[0]].copy()
        extra["Metabolite name"] = "PC C"
        extra["MS/MS spectrum"] = "100:1 200:1 300:1 400:1"
        pooled, tables, _, _ = entropy.entropy_sharing_tables(pd.concat([source, extra]), "positive")
        # Pooled valid rows: 0, 0, ln(2), ln(4); dataset 1: 0, 0, ln(4).
        self.assertAlmostEqual(pooled.loc["PC", entropy.ENTROPY], np.log(2) / 2)
        self.assertEqual(tables[entropy.datasets[0]].loc["PC", entropy.ENTROPY], 0)
        self.assertEqual(pooled.loc["PC", "Accepted ID rows"], 5)
        self.assertEqual(tables[entropy.datasets[0]].loc["PC", "Accepted ID rows"], 4)

    def test_all_missing_class_has_no_entropy_and_zero_valid_row_coverage(self):
        source = entropy.polarity_subsets(self.fixture())["positive"].copy()
        source["MS/MS spectrum"] = None
        pooled, _, _, _ = entropy.entropy_sharing_tables(source, "positive")
        self.assertTrue(pooled[entropy.ENTROPY].isna().all())
        self.assertTrue(pooled["Rows with valid entropy"].eq(0).all())
        self.assertTrue(pooled["Entropy row coverage (%)"].eq(0).all())

    def test_correlation_includes_zero_and_excludes_missing(self):
        table = pd.DataFrame({entropy.ENTROPY: [0, 1, 2, np.nan], "Shared lipids (%)": [10, 20, 30, 40]})
        result = entropy.entropy_correlation(table).iloc[0]
        self.assertEqual(result["Classes"], 3)
        self.assertAlmostEqual(result["Spearman rho"], 1)
        self.assertIn("fewer than three", entropy.entropy_correlation(table.iloc[:2]).iloc[0]["Status"])
        table["Shared lipids (%)"] = 100
        self.assertIn("constant", entropy.entropy_correlation(table).iloc[0]["Status"])

    def test_plot_handles_absent_datasets_and_both_polarities(self):
        analyses = {}
        for polarity, source in entropy.polarity_subsets(self.fixture()).items():
            pooled, tables, _, _ = entropy.entropy_sharing_tables(source, polarity)
            statistics = pd.concat([
                entropy.entropy_correlation(pooled).assign(Dataset="Pooled"),
                *[entropy.entropy_correlation(table).assign(Dataset=dataset) for dataset, table in tables.items()],
            ])
            analyses[polarity] = pooled, tables, statistics
        self.assertEqual(len(entropy.plot_entropy_correlations(analyses).axes), 1)
        self.assertEqual(len(entropy.plot_entropy_correlations(analyses, by_dataset=True).axes), 6)


if __name__ == "__main__":
    unittest.main()

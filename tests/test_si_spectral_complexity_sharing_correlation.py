import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Analysis"
SPEC = importlib.util.spec_from_file_location(
    "complexity_sharing", SCRIPT_DIR / "SI_Spectral_Complexity_Sharing_Correlation.py"
)
complexity = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(SCRIPT_DIR), *sys.path]):
    SPEC.loader.exec_module(complexity)


class TestSpectralComplexitySharing(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_inclusive_threshold_and_intensity_scaling(self):
        for spectrum, expected, positive in [
            ("100:100 200:1 300:0.99 400:0", 2, 3),
            ("100:10000 200:100 300:99", 2, 3),
            ("100:1e308 200:1e306 300:1e305", 2, 3),
            ("100:2", 1, 1),
        ]:
            with self.subTest(spectrum=spectrum):
                self.assertEqual(complexity.spectral_complexity(spectrum), (expected, positive, "OK"))

    def test_duplicates_merge_before_threshold_and_base_peak(self):
        # Merging raises the base peak from 100 to 200: intensity 1.5 falls below 1%.
        self.assertEqual(complexity.spectral_complexity("100:100 100:100 200:1.5"), (1, 2, "OK"))
        # Two sub-threshold entries at the same m/z become one retained peak.
        self.assertEqual(complexity.spectral_complexity("100:100 200:0.6 200:0.6"), (2, 2, "OK"))

    def test_unusable_spectra_remain_missing(self):
        for spectrum in [None, np.nan, "", "100:0", "100:-1", "100:nan", "inf:1", "bad", "100:1:2"]:
            with self.subTest(spectrum=spectrum):
                value, _, status = complexity.spectral_complexity(spectrum)
                self.assertTrue(np.isnan(value))
                self.assertNotEqual(status, "OK")

    def fixture(self):
        first, second = complexity.datasets[:2]
        return pd.DataFrame([
            (first, "PC", "PC A|adduct1", "100:1", "positive", True),
            (first, "PC", "PC A|adduct2", "100:1", "positive", True),
            (second, "PC", "PC A|feature", "100:1 200:1", "positive", True),
            (first, "PC", "PC B", None, "positive", True),
            (first, "PC", "PC C", "100:1 200:1 300:1 400:1", "positive", True),
            (first, "PC", "PC B", "100:1", "negative", True),
            (second, "PC", "PC B", "100:1", "positive", False),
        ], columns=["Dataset", "Ontology", "Metabolite name", "MS/MS spectrum", "Polarity", "Overall Match Status"])

    def test_medians_and_sharing_denominators(self):
        selected = complexity.polarity_subsets(self.fixture())["positive"]
        pooled, tables, audit, lipids = complexity.complexity_sharing_tables(selected, "positive")
        original = complexity.class_sharing_table(selected)
        pd.testing.assert_frame_equal(pooled[original.columns], original.rename_axis(None, axis=1))
        for dataset, original in complexity.dataset_sharing_tables(selected).items():
            pd.testing.assert_frame_equal(tables[dataset][original.columns], original)
        # Pooled valid rows: 1, 1, 2, 4; B stays missing. Repeats contribute.
        self.assertEqual(pooled.loc["PC", complexity.COMPLEXITY], 1.5)
        self.assertAlmostEqual(pooled.loc["PC", "Shared lipids (%)"], 100 / 3)
        self.assertEqual(pooled.loc["PC", "Accepted ID rows"], 5)
        self.assertEqual(pooled.loc["PC", "Rows with valid complexity"], 4)
        self.assertEqual(pooled.loc["PC", "Complexity row coverage (%)"], 80)
        self.assertEqual(tables[complexity.datasets[0]].loc["PC", complexity.COMPLEXITY], 1)
        self.assertTrue(tables[complexity.datasets[-1]][complexity.COMPLEXITY].isna().all())
        self.assertTrue(tables[complexity.datasets[-1]]["Accepted ID rows"].eq(0).all())
        self.assertTrue(tables[complexity.datasets[-1]]["Complexity row coverage (%)"].isna().all())
        self.assertEqual(len(audit), 5)
        self.assertEqual(len(lipids), 4)

    def test_all_missing_class_has_no_complexity_and_zero_valid_row_coverage(self):
        source = complexity.polarity_subsets(self.fixture())["positive"].copy()
        source["MS/MS spectrum"] = None
        pooled, _, _, _ = complexity.complexity_sharing_tables(source, "positive")
        self.assertTrue(pooled[complexity.COMPLEXITY].isna().all())
        self.assertTrue(pooled["Rows with valid complexity"].eq(0).all())
        self.assertTrue(pooled["Complexity row coverage (%)"].eq(0).all())

    def test_correlation_omits_missing_and_reports_undefined(self):
        table = pd.DataFrame({complexity.COMPLEXITY: [1, 2, 3, np.nan], "Shared lipids (%)": [30, 20, 10, 40]})
        result = complexity.complexity_correlation(table).iloc[0]
        self.assertEqual(result["Classes"], 3)
        self.assertAlmostEqual(result["Spearman rho"], -1)
        self.assertIn("fewer than three", complexity.complexity_correlation(table.iloc[:2]).iloc[0]["Status"])
        table["Shared lipids (%)"] = 100
        self.assertIn("constant", complexity.complexity_correlation(table).iloc[0]["Status"])

    def test_cli_exports_both_polarities_and_figures(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "input.csv"
            self.fixture().to_csv(source, index=False)
            # Draw real panels, but omit costly raster/PDF/SVG rendering in this test.
            with patch("matplotlib.figure.Figure.savefig") as save:
                complexity.main(["--input", str(source), "--output-dir", str(folder), "--no-show"])
            self.assertEqual(save.call_count, 4)
            self.assertEqual(len(list(folder.glob("*Statistics.csv"))), 3)
            for polarity in ["Positive", "Negative"]:
                self.assertTrue((folder / f"{complexity.OUTPUT_BASENAME}_{polarity}_Complexity_Audit.csv").exists())
            exported = pd.read_csv(folder / f"{complexity.OUTPUT_BASENAME}_Positive_Classes.csv")
            self.assertEqual(exported.loc[0, complexity.COMPLEXITY], 1.5)
            self.assertEqual(exported.loc[0, "Accepted ID rows"], 5)
            stats = pd.read_csv(folder / f"{complexity.OUTPUT_BASENAME}_Statistics.csv")
            self.assertTrue(stats["Spectral summary basis"].str.contains("all accepted ID rows").all())
            self.assertEqual(plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

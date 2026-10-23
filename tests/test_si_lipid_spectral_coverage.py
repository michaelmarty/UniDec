import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Analysis"
SPEC = importlib.util.spec_from_file_location("lipid_coverage", SCRIPT_DIR / "SI_Lipid_Spectral_Coverage.py")
coverage = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(SCRIPT_DIR), *sys.path]):
    SPEC.loader.exec_module(coverage)


class TestLipidSpectralCoverage(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def source(self):
        return pd.DataFrame([
            (" OTOT\u00a0", "Cer_AP", "Cer A|adduct", "100:1", "Positive", True),
            ("OTOT", "Cer_NS", "Cer A|feature", "100:1", "positive", True),
            ("OTIT", "Cer_NS", "Cer A", "100:1 200:1", "positive", True),
            ("ITIT", "Cer_NS", "Cer A", None, "positive", True),
            ("Stellar", "Cer_NS", "Cer A", "100:1", "negative", True),
            ("Stellarboth", "Cer_NS", "Cer A", "100:1", "positive", False),
            ("Other", "Cer_NS", "Cer A", "100:1", "positive", True),
            ("OTOT", "PC", "PC A", "100:1 200:1", "positive", True),
            ("OTOT", "PC", "PC missing", None, "positive", True),
            ("OTIT", "PC", "PC missing", None, "negative", True),
        ], columns=["Dataset", "Ontology", "Metabolite name", "MS/MS spectrum", "Polarity", "Overall Match Status"])

    def test_unique_identity_and_same_polarity_dataset_coverage(self):
        lipids, audit = coverage.lipid_coverage_tables(self.source())
        self.assertEqual(len(audit), 8)
        self.assertEqual(len(lipids), 5)
        self.assertFalse(lipids.duplicated(coverage.KEYS).any())
        cer = lipids.loc[lipids["Lipid Class"].eq("Cer")].set_index("Polarity")
        self.assertEqual(cer.loc["positive", coverage.COVERAGE], 3)
        self.assertEqual(cer.loc["negative", coverage.COVERAGE], 1)
        self.assertEqual(cer.loc["positive", "Accepted ID rows"], 4)
        self.assertEqual(cer.loc["positive", "Valid entropy spectra"], 3)
        self.assertEqual(cer.loc["positive", "Detected datasets"], "OTOT, OTIT, ITIT")
        # Row median of 0, 0, ln(2), rather than median of dataset medians.
        self.assertEqual(cer.loc["positive", coverage.ENTROPY], 0)
        self.assertEqual(cer.loc["positive", coverage.COMPLEXITY], 1)

    def test_missing_spectra_retain_coverage_but_are_not_tested(self):
        lipids, _ = coverage.lipid_coverage_tables(self.source())
        missing = lipids.loc[lipids["Lipid Name"].eq("PC missing")]
        self.assertTrue(missing[coverage.COVERAGE].eq(1).all())
        self.assertTrue(missing[coverage.ENTROPY].isna().all())
        self.assertFalse(missing["Included in correlations"].any())
        classes = coverage.class_coverage_summary(lipids)
        pc = classes.loc[("positive", "PC")]
        self.assertEqual(pc["Accepted lipid identities"], 2)
        self.assertEqual(pc["Lipids with valid spectra"], 1)
        self.assertEqual(pc["Lipids in 1 datasets"], 2)

    def test_statistics_use_true_integer_coverage_and_same_lipids(self):
        lipids = pd.DataFrame({"Polarity": ["positive"] * 4,
                               coverage.ENTROPY: [0, 1, 2, np.nan],
                               coverage.COMPLEXITY: [1, 2, 3, np.nan],
                               coverage.COVERAGE: [5, 3, 1, 2],
                               "Included in correlations": [True, True, True, False]})
        results = coverage.coverage_statistics(lipids)
        self.assertTrue(results["Lipids tested"].eq(3).all())
        self.assertTrue(results["Accepted lipid identities"].eq(4).all())
        np.testing.assert_allclose(results["Spearman rho"], [-1, -1])
        lipids[coverage.COVERAGE] = 1
        self.assertTrue(coverage.coverage_statistics(lipids)["Status"].str.contains("constant").all())

    def test_cli_exports_two_polarity_panels_per_metric(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "input.csv"
            self.source().to_csv(path, index=False)
            with patch("matplotlib.figure.Figure.savefig") as save:
                coverage.main(["--input", str(path), "--no-show"])
            self.assertEqual(save.call_count, 4)
            for suffix in ["Lipids", "Row_Audit", "Classes", "Statistics", "Pairwise_Statistics"]:
                self.assertTrue((root / f"{coverage.OUTPUT_BASENAME}_{suffix}.csv").exists())
            lipids, _ = coverage.lipid_coverage_tables(self.source())
            fig = coverage.plot_coverage(lipids, coverage.pairwise_statistics(lipids), coverage.ENTROPY)
            self.assertEqual(len(fig.axes), 2)
            self.assertEqual(fig.axes[0].get_yticks().tolist(), [1, 2, 3, 4, 5])

    def test_four_planned_comparisons_and_holm_per_metric(self):
        lipids = pd.DataFrame({"Polarity": ["positive"] * 6,
                               coverage.ENTROPY: [1, 2, 3, 4, 5, np.nan],
                               coverage.COMPLEXITY: [1, 2, 3, 4, 5, np.nan],
                               coverage.COVERAGE: [1, 2, 3, 4, 5, 5],
                               "Included in correlations": [True] * 5 + [False]})
        with patch.object(coverage, "mannwhitneyu", side_effect=[(1, p) for p in [0.01, 0.04, 0.03, 0.2] * 2]) as test:
            result = coverage.pairwise_statistics(lipids)
        self.assertEqual(test.call_count, 8)
        self.assertEqual(result["Comparison"].tolist(), ["4 vs 5", "3 vs 5", "2 vs 5", "1 vs 5"] * 2)
        np.testing.assert_allclose(result["Holm-adjusted p-value"], [0.04, 0.09, 0.09, 0.2] * 2)
        self.assertTrue(result["Reference lipids"].eq(1).all())
        self.assertEqual(test.call_args.kwargs["alternative"], "two-sided")

    def test_pairwise_missing_reference_and_tied_distributions(self):
        lipids, _ = coverage.lipid_coverage_tables(self.source())
        result = coverage.pairwise_statistics(lipids)
        self.assertTrue(result["Holm-adjusted p-value"].isna().all())
        self.assertTrue(result["Status"].str.contains("empty").all())
        tied = pd.DataFrame({"Polarity": ["positive"] * 10,
                             coverage.ENTROPY: [1] * 10, coverage.COMPLEXITY: [1] * 10,
                             coverage.COVERAGE: [1, 2, 3, 4, 5] * 2,
                             "Included in correlations": [True] * 10})
        self.assertTrue(coverage.pairwise_statistics(tied)["Holm-adjusted p-value"].eq(1).all())


if __name__ == "__main__":
    unittest.main()

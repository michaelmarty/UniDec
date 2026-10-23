import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Analysis"
SPEC = importlib.util.spec_from_file_location("class_size_coverage", SCRIPT_DIR / "SI_Class_Size_Coverage.py")
coverage = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(SCRIPT_DIR), *sys.path]):
    SPEC.loader.exec_module(coverage)


class TestClassSizeCoverage(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def source(self):
        return pd.DataFrame([
            ("OTOT", "Cer_AP", "Cer A|one", True, "positive"),
            ("OTOT", "Cer_NS", "Cer A|two", True, "positive"),
            ("OTIT", "Cer_NS", "Cer A", True, "positive"),
            ("OTOT", "PC", "PC A", True, "positive"),
            ("OTOT", "PC", "PC B", True, "positive"),
            ("OTOT", "TG", "TG A", True, "positive"),
            ("OTIT", "TG", "TG A", True, "positive"),
            ("OTOT", "TG", "TG B", True, "positive"),
            ("OTOT", "TG", "TG C", True, "positive"),
            ("ITIT", "PC", "PC A", True, "negative"),
            ("ITIT", "PC", "PC B", False, "positive"),
            ("Other", "PC", "PC B", True, "positive"),
        ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status", "Polarity"])

    def test_unique_counts_exact_coverage_and_duplicate_handling(self):
        tables, counts, statistics = coverage.class_size_tables(self.source())
        positive = tables["positive"]
        self.assertEqual(positive.loc["Cer", coverage.SIZE], 1)
        self.assertEqual(positive.loc["Cer", "Accepted ID rows"], 3)
        self.assertEqual(positive.loc["Cer", coverage.MEAN], 2)
        self.assertEqual(positive.loc["PC", coverage.MEAN], 1)
        self.assertAlmostEqual(positive.loc["TG", coverage.MEAN], 4 / 3)
        self.assertEqual(positive.loc["TG", "Median dataset coverage"], 1)
        self.assertEqual(counts["positive"].loc["PC", 1], 2)
        self.assertEqual(tables["negative"].loc["PC", coverage.SIZE], 1)
        result = statistics.set_index("Polarity").loc["positive"]
        self.assertEqual(result["Classes"], 3)
        self.assertAlmostEqual(result["Spearman rho"], -0.5)

    def test_cli_exports_and_handles_single_polarity(self):
        source = self.source().loc[lambda frame: frame["Polarity"].eq("positive")]
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source.to_csv(root / "input.csv", index=False)
            with patch("matplotlib.figure.Figure.savefig") as save:
                coverage.main(["--input", str(root / "input.csv"), "--no-show"])
            self.assertEqual(save.call_count, 4)
            self.assertTrue((root / f"{coverage.OUTPUT_BASENAME}_Positive_Classes.csv").exists())
            self.assertTrue((root / f"{coverage.OUTPUT_BASENAME}_Statistics.csv").exists())
            self.assertEqual(plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

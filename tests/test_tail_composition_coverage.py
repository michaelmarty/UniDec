import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestTailCompositionCoverage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location("composition_coverage", directory / "SI_8_9_10.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]):
            spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def fixture(self):
        rows = [
            ("OTOT", "PE", "shared", "16:0", 2),
            ("OTIT", "PE", "shared", "16:0", 2),
            ("OTOT", "FA", "unique", "16:0", 1),
            ("OTOT", "FA", "second", "18:0", 1),
            ("OTIT", "FA", "third", "20:0", 1),
            ("OTOT", "SM", "backbone", "LCB 18:1", 1),
        ]
        occurrences = pd.DataFrame(rows, columns=["Dataset", "Ontology", "Species", "Tail", "Tail copies"])
        lipids = occurrences[["Dataset", "Ontology", "Species"]].drop_duplicates()
        return lipids, occurrences

    def test_exact_denominators_multiplicity_and_shared_lipid_counting(self):
        lipids, occurrences = self.fixture()
        paired, results = self.script.comparison_tables(lipids, occurrences, ["OTOT", "OTIT", "Stellar"])
        values = paired.set_index("Tail")
        self.assertEqual(values.loc["16:0", "S7 OTOT (%)"], 60)
        self.assertAlmostEqual(values.loc["16:0", "Mean S7 composition (%)"], (60 + 200 / 3) / 2)
        self.assertAlmostEqual(values.loc["16:0", "Figure 5A dataset-specific (%)"], 100 / 3)
        self.assertEqual(values.loc["16:0", "Pooled tail copies"], 3)
        self.assertEqual(values.loc["18:0", "S7 OTIT (%)"], 0)
        self.assertEqual(values.loc["18:0", "Mean S7 composition (%)"], 10)
        self.assertNotIn("LCB 18:1", paired.Tail.tolist())
        self.assertTrue(values["S7 Stellar (%)"].isna().all())
        self.assertTrue(values["Nonempty datasets"].eq(2).all())
        self.assertTrue(np.isfinite(results.iloc[0]["Pearson r"]))
        self.assertEqual(results.iloc[0]["Tails"], 3)

    def test_empty_display_and_constant_coverage_are_unavailable(self):
        lipids, occurrences = self.fixture()
        paired, results = self.script.comparison_tables(lipids, occurrences, ["OTOT", "OTIT"], 100)
        self.assertTrue(paired.empty)
        self.assertTrue(results["Pearson r"].isna().all())
        self.script.plot_comparison(paired, results).canvas.draw()
        occurrences = occurrences.assign(Dataset="OTOT").drop_duplicates(["Dataset", "Ontology", "Species", "Tail"])
        lipids = occurrences[["Dataset", "Ontology", "Species"]].drop_duplicates()
        paired, results = self.script.comparison_tables(lipids, occurrences, ["OTOT", "OTIT"])
        self.assertIn("constant", results.iloc[0].Status)
        self.assertTrue(results["Pearson nominal p"].isna().all())


if __name__ == "__main__":
    unittest.main()

import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestTailComposition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location("tail_dotplots", directory / "SI_5_6.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]):
            spec.loader.exec_module(cls.script)


    def tearDown(self):
        self.script.plt.close("all")


    def fixture(self):
        occurrences = pd.DataFrame([
            {"Dataset": "OTOT", "Ontology": "PE", "Species": "PE 16:0_16:0", "Tail": "16:0", "Tail copies": 2},
            {"Dataset": "OTOT", "Ontology": "SM", "Species": "SM LCB 18:1_18:0", "Tail": "LCB 18:1", "Tail copies": 1},
            {"Dataset": "OTOT", "Ontology": "SM", "Species": "SM LCB 18:1_18:0", "Tail": "18:0", "Tail copies": 1},
            {"Dataset": "OTIT", "Ontology": "PE", "Species": "PE 18:0_18:0", "Tail": "18:0", "Tail copies": 2},
        ])
        lipids = occurrences[["Dataset", "Ontology", "Species"]].drop_duplicates()
        return lipids, occurrences


    def test_counts_denominators_hidden_tails_and_missing_datasets(self):
        lipids, occurrences = self.fixture()
        table, totals, percentages = self.script.composition_tables(
            lipids, occurrences, ["OTOT", "OTIT", "Stellar"])
        self.assertEqual(totals.to_dict(), {"OTOT": 4, "OTIT": 2, "Stellar": 0})
        counts = table.set_index(["Dataset", "Tail"])
        self.assertEqual(counts.loc[("OTOT", "16:0"), "Tail occurrences"], 2)
        self.assertEqual(counts.loc[("OTOT", "16:0"), "Percentage"], 50)
        self.assertEqual(counts.loc[("OTIT", "16:0"), "Percentage"], 0)
        self.assertFalse(counts.loc[("OTOT", "LCB 18:1"), "Displayed"])
        self.assertNotIn("LCB 18:1", percentages.index)
        self.assertEqual(percentages.loc["18:0", "OTOT"], 25)
        self.assertTrue(percentages["Stellar"].isna().all())
        np.testing.assert_allclose(table.groupby("Dataset")["Percentage"].sum().loc[["OTOT", "OTIT"]], 100)

    def test_filtering_deduplication_and_repeated_chains(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "ids.csv"
            pd.DataFrame({"Dataset": ["OTOT"] * 4, "Ontology": ["PE"] * 4,
                          "Metabolite name": ["PE 16:0_16:0", "PE 16:0_16:0", "PE 18:1_18:1", "PE 34:1"],
                          "Overall Match Status": [True, True, False, True]}).to_csv(path, index=False)
            lipids, occurrences, audit = self.script.prepare_lipids(path, ["OTOT"])
        self.assertEqual(len(lipids), 1)
        self.assertEqual(occurrences["Tail copies"].tolist(), [2])
        self.assertTrue(audit["Exclusion reason"].str.contains("Duplicate").any())


    def test_correlations_use_paired_tails_and_handle_undefined_values(self):
        percentages = pd.DataFrame({
            "OTOT": [1., 2., 3., 4.], "OTIT": [2., 4., 6., np.nan],
            "ITIT": [4., 3., 2., 1.], "Stellar": [5., 5., 5., 5.],
            "Stellarboth": [1., np.nan, np.nan, 4.],
        }, index=["14:0", "16:0", "18:0", "18:1"])
        original = percentages.copy()
        _, table = self.script.plot_tail_correlations(percentages)
        rows = table.set_index(["Dataset A", "Dataset B"])
        self.assertEqual(len(table), 10)
        self.assertAlmostEqual(rows.loc[("OTOT", "OTIT"), "Pearson r"], 1)
        self.assertEqual(rows.loc[("OTOT", "OTIT"), "Tails"], 3)
        self.assertAlmostEqual(rows.loc[("OTOT", "ITIT"), "Pearson r"], -1)
        self.assertTrue(np.isnan(rows.loc[("OTOT", "Stellar"), "Pearson r"]))
        self.assertIn("constant", rows.loc[("OTOT", "Stellar"), "Status"])
        self.assertIn("fewer than 3", rows.loc[("OTOT", "Stellarboth"), "Status"])
        pd.testing.assert_frame_equal(percentages, original)


    def test_cli_exports_only_heatmap_and_correlations(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "ids.csv"
            output = Path(directory) / "output"
            pd.DataFrame({"Dataset": ["OTOT", "OTIT"], "Ontology": ["PE", "PE"],
                          "Metabolite name": ["PE 16:0_16:0", "PE 18:0_18:0"],
                          "Overall Match Status": [True, True]}).to_csv(path, index=False)
            self.script.main(["--input", str(path), "--output-dir", str(output),
                              "--datasets", "OTOT", "OTIT", "--no-show"])
            counts = pd.read_csv(output / "SI_Tail_Dotplots_Tail_Composition.csv")
            self.assertEqual(counts["Total tail occurrences"].tolist(), [2, 2, 2, 2])
            expected = {f"SI_Tail_Dotplots_{name}.{ext}" for name in ["Overview", "Correlations"]
                        for ext in ["png", "pdf"]}
            expected.update(f"SI_Tail_Dotplots_{name}.csv" for name in
                            ["Tail_Composition", "Denominators", "Correlations"])
            self.assertEqual({file.name for file in output.iterdir()}, expected)

    def test_all_ten_correlations_fit_one_manuscript_page(self):
        percentages = pd.DataFrame({d: [1., 2., 3., 4.] for d in self.script.DATASETS},
                                   index=["14:0", "16:0", "18:0", "18:1"])
        fig, table = self.script.plot_tail_correlations(percentages)
        np.testing.assert_allclose(fig.get_size_inches(), [6.0, 6.25])
        self.assertEqual(len(table), 10)
        active = [ax for ax in fig.axes if ax.axison]
        self.assertEqual(len(active), 10)
        for ax in active:
            self.assertEqual(ax.get_xlim(), active[0].get_xlim())
            self.assertEqual(ax.get_ylim(), active[0].get_ylim())
        fig.canvas.draw()
        bounds = fig.get_tightbbox(fig.canvas.get_renderer())
        self.assertGreaterEqual(bounds.x0, -0.001)
        self.assertGreaterEqual(bounds.y0, -0.001)
        self.assertLessEqual(bounds.x1, 6.001)
        self.assertLessEqual(bounds.y1, 6.251)


if __name__ == "__main__":
    unittest.main()

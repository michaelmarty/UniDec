import contextlib
import importlib.util
import io
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import pandas as pd


class TestSIClassHeatmap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                / "DDA Analysis" / "SI_Class_Heatmap.py")
        spec = importlib.util.spec_from_file_location("si_class_heatmap", path)
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(path.parent), *sys.path]):
            spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def source(self):
        return pd.DataFrame([
            (" otot\u00a0", "PC", "PC 34:1|PC 16:0_18:1", True),
            ("OTOT", "PC", "PC 34:1|PC 16:1_18:0", "TRUE"),
            ("OTOT", "PC", "PC 36:2", True),
            ("OTOT", "EtherPC", "PC O-34:1", True),
            ("OTOT", "PE", "PE 36:2", True),
            ("OTIT", "LPC", "LPC 18:0", True),
            ("OTIT", "TG", "TG 54:3", False),
            ("Stellar", "PC", "PC 34:1", False),
            ("OTIT", "", "Missing class", True),
            ("OTIT", "PC", None, True),
            ("Unknown", "PC", "PC 34:1", True),
        ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status"])

    def test_class_counts_include_summed_ids_and_group_subclasses(self):
        lipids = self.script.prepare_lipids(self.source(), self.script.DATASETS)
        counts, totals, percentages = self.script.class_tables(lipids, self.script.DATASETS)
        self.assertEqual(totals.tolist(), [4, 1, 0, 0, 0])
        self.assertEqual(counts.at["PC", "OTOT"], 3)
        self.assertEqual(counts.at["PC", "OTIT"], 1)
        self.assertEqual(counts.at["PE", "OTOT"], 1)
        self.assertEqual(counts.at["PE", "OTIT"], 0)
        self.assertEqual(list(counts.index), ["PC", "PE"])
        self.assertTrue(counts["Stellar"].isna().all())
        self.assertTrue(percentages["Stellar"].isna().all())
        self.assertAlmostEqual(percentages.at["PC", "OTOT"], 75)
        self.assertAlmostEqual(percentages["OTOT"].sum(), 100)

    def test_selection_keeps_full_denominators_and_handles_empty_display(self):
        lipids = self.script.prepare_lipids(self.source(), self.script.DATASETS)
        counts, totals, percentages = self.script.class_tables(lipids, self.script.DATASETS, 50)
        self.assertNotIn("PE", counts.index)
        self.assertEqual(totals["OTOT"], 4)
        self.assertAlmostEqual(percentages["OTOT"].sum(), 75)
        counts, totals, percentages = self.script.class_tables(lipids, ["OTOT"], 100)
        self.assertTrue(counts.empty)
        fig = self.script.plot_class_heatmap(counts, percentages, totals, "Accepted")
        self.assertTrue(any("No classes" in text.get_text() for text in fig.axes[0].texts))

    def test_grouping_and_duplicate_identity_match_figure3(self):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        with patch.object(sys, "path", [str(directory), *sys.path]):
            from Figure3 import class_groups, prepare_coverage

        classes = list(class_groups) + ["LPC", "EtherPC", "OxPC", "EtherLPE", "TG"]
        source = pd.DataFrame([
            ("OTOT", lipid_class, f"{lipid_class} 34:1", True)
            for lipid_class in classes
        ] + [
            ("OTOT", "Cer_AP", "Cer 34:1|first assignment", True),
            ("OTOT", "Cer_NS", "Cer 34:1|second assignment", True),
            ("OTIT", "Cer_NS", "Cer 34:1", True),
        ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status"])
        lipids = self.script.prepare_lipids(source, self.script.DATASETS)
        counts, totals, _ = self.script.class_tables(lipids, self.script.DATASETS)
        coverage, _ = prepare_coverage(source)
        pd.testing.assert_series_equal(
            counts.sum(axis=1).astype("int64"),
            coverage.mul(coverage.columns, axis=1).sum(axis=1),
            check_names=False,
        )
        self.assertEqual(totals["OTOT"], len(classes) + 1)
        self.assertEqual(totals["OTIT"], 1)
        self.assertEqual(len(lipids.loc[lipids["Lipid Match Key"].eq("Cer 34:1")]), 2)

    def test_all_ids_and_invalid_inputs(self):
        frame = self.source().drop(columns="Overall Match Status")
        with self.assertRaisesRegex(ValueError, "Overall Match Status"):
            self.script.prepare_lipids(frame, self.script.DATASETS)
        lipids = self.script.prepare_lipids(frame, self.script.DATASETS, all_ids=True)
        self.assertIn("TG", set(lipids["Ontology"]))
        for threshold in [-1, 101, float("nan")]:
            with self.assertRaises(ValueError):
                self.script.class_tables(lipids, self.script.DATASETS, threshold)
        with self.assertRaisesRegex(ValueError, "No eligible"):
            self.script.prepare_lipids(frame.iloc[:0], self.script.DATASETS, all_ids=True)

    def test_panels_show_counts_first_and_percentages_without_counts(self):
        lipids = self.script.prepare_lipids(self.source(), self.script.DATASETS)
        counts, totals, percentages = self.script.class_tables(lipids, self.script.DATASETS)
        fig = self.script.plot_class_heatmap(counts, percentages, totals, "Accepted")
        first, second = fig.axes[:2]
        self.assertIn("A", [text.get_text() for text in first.texts])
        self.assertIn("B", [text.get_text() for text in second.texts])
        self.assertEqual(first.images[0].get_array()[0, 0], 3)
        self.assertAlmostEqual(second.images[0].get_array()[0, 0], 75)
        self.assertIn("75.0%", [text.get_text() for text in second.texts])
        self.assertFalse(any("n=" in text.get_text() for text in second.texts))

    def test_main_exports_only_combined_heatmap(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            self.source().to_csv(source, index=False)
            with contextlib.redirect_stdout(io.StringIO()):
                self.script.main(["--input", str(source), "--output-dir", str(root), "--no-show"])
            for extension in ["png", "pdf"]:
                self.assertGreater((root / f"SI_ClassHeatmap.{extension}").stat().st_size, 0)
            self.assertEqual({p.name for p in root.iterdir()},
                             {"input.csv", "SI_ClassHeatmap.png", "SI_ClassHeatmap.pdf"})
            self.assertEqual(self.script.plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

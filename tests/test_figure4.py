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
import numpy as np
import pandas as pd
from PIL import Image


class TestFigure5(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_path = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                       / "DDA Paper Figures" / "Figure5.py")
        spec = importlib.util.spec_from_file_location("figure4_script", script_path)
        cls.script = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def prepare(self, rows, all_ids=False):
        source = pd.DataFrame(rows, columns=[
            "Dataset", "Ontology", "Metabolite name", "Overall Match Status",
        ])
        with patch.object(self.script.pd, "read_excel", return_value=source):
            return self.script.prepare_lipids(Path("input.xlsx"), self.script.DATASETS, all_ids)

    def tables(self, rows, threshold=1):
        lipids, occurrences, audit = self.prepare(rows)
        overall = self.script.frequency_summary(lipids, occurrences, ["Dataset"])
        counts, denominators, percentages = self.script.tail_prevalence_tables(
            overall, self.script.DATASETS, threshold
        )
        return counts, denominators, percentages, occurrences, audit

    def test_distinct_tails_classes_and_compositions_count_once_per_method(self):
        counts, denominators, percentages, occurrences, audit = self.tables([
            ("OTOT", "PC", "PC 34:1|PC 16:0_18:1", True),
            (" otot\u00a0", "PC", "PC 18:1/16:0", "TRUE"),
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "TG", "TG 16:0_18:1_18:1", True),
            ("OTOT", "TG", "TG 16:0_16:0_18:1", True),
            ("ITIT", "PC", "PC 16:0_18:1", True),
            ("OTIT", "PC", "PC 34:1", True),
            ("Stellar", "PC", "PC 24:0_26:0", False),
            ("Stellarboth", "CE", "CE 30:0", True),
            ("Unselected", "PC", "PC 24:0_26:0", True),
        ])
        self.assertEqual(denominators.to_dict(), {
            "OTOT": 4, "OTIT": 0, "ITIT": 1, "Stellar": 0, "Stellarboth": 1,
        })
        self.assertEqual(counts.at["16:0", "OTOT"], 4)
        self.assertEqual(counts.at["18:1", "OTOT"], 3)
        self.assertEqual(percentages.at["18:1", "OTOT"], 75)
        self.assertEqual(percentages.at["16:0", "OTOT"], 100)
        self.assertEqual(percentages["OTOT"].sum(), 175)
        repeated = occurrences.loc[occurrences["Species"].eq("PE 16:0_16:0")]
        self.assertEqual(len(repeated), 1)
        self.assertEqual(repeated.iloc[0]["Tail copies"], 2)
        self.assertEqual(audit["Exclusion reason"].eq("Duplicate chain composition within dataset").sum(), 1)
        self.assertNotIn("24:0", counts.index)
        instrument_percentages = self.script.instrument_tail_tables(occurrences, self.script.DATASETS)[1]
        self.assertEqual(instrument_percentages.at["18:1", "OTOT"], 100)
        coverage = self.script.tail_coverage_tables(occurrences, self.script.DATASETS, ["18:1"])[0]
        self.assertEqual(coverage.at["18:1", 1], 2)
        self.assertEqual(coverage.at["18:1", 2], 1)

    def test_absent_tails_are_zero_and_methods_without_resolved_lipids_are_missing(self):
        counts, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("ITIT", "CE", "CE 30:0", True),
            ("OTIT", "PC", "PC 34:1", True),
        ])
        self.assertEqual(counts.at["16:0", "ITIT"], 0)
        self.assertEqual(percentages.at["16:0", "ITIT"], 0)
        for dataset in ["OTIT", "Stellar", "Stellarboth"]:
            self.assertEqual(denominators[dataset], 0)
            self.assertTrue(counts[dataset].isna().all())
            self.assertTrue(percentages[dataset].isna().all())
        self.assertEqual(list(percentages.columns), self.script.DATASETS)

    def test_selection_uses_lipid_prevalence_and_never_renormalizes(self):
        rows = [
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("OTOT", "PC", "PC 20:0_22:0", True),
        ]
        counts, denominators, percentages, _, _ = self.tables(rows, 30)
        self.assertEqual(set(counts.index), {"16:0", "18:1", "20:0", "22:0"})
        self.assertAlmostEqual(percentages.at["18:1", "OTOT"], 100 / 3)
        selected_counts, selected_denominators, selected_percentages, _, _ = self.tables(rows, 60)
        self.assertEqual(list(selected_counts.index), ["16:0"])
        pd.testing.assert_series_equal(denominators, selected_denominators)
        self.assertAlmostEqual(selected_percentages.at["16:0", "OTOT"], 200 / 3)
        empty_counts, empty_denominators, empty_percentages, _, _ = self.tables(rows, 100)
        self.assertTrue(empty_counts.empty)
        self.assertTrue(empty_percentages.empty)
        self.assertEqual(empty_denominators["OTOT"], 3)

    def test_structural_exclusions_preserve_denominators_and_modified_and_p_tails(self):
        counts, denominators, percentages, _, _ = self.tables([
            ("OTOT", "EtherPC", "PC O-16:0_18:1", True),
            ("OTOT", "EtherPE", "PE P-16:0_20:1", True),
            ("OTOT", "Cer", "Cer 18:1;O2/24:0", True),
            ("OTOT", "OxPC", "PC 16:0_18:1;O2", True),
            ("OTOT", "EtherPC", "PC O-10:0_O-12:0", True),
            ("ITIT", "EtherPC", "PC O-10:0_O-12:0", True),
        ])
        self.assertEqual(list(counts.index), ["16:0", "18:1", "18:1;O2", "20:1", "24:0", "P-16:0"])
        self.assertEqual(denominators["OTOT"], 5)
        self.assertEqual(denominators["ITIT"], 1)
        self.assertTrue(percentages["OTOT"].eq(20).all())
        self.assertTrue(percentages["ITIT"].eq(0).all())
        self.assertFalse(any(tail.startswith(("LCB ", "O-")) for tail in counts.index))

    def test_parser_excludes_sum_annotations_by_chain_count_and_retains_long_chains(self):
        for lipid_class, name in [("PC", "PC 34:1"), ("TG", "TG 54:3"),
                                  ("CL", "CL 72:4"), ("EtherPC", "PC O-34:1")]:
            _, chains, reason = self.script.parse_lipid_tails(name, lipid_class)
            self.assertEqual(chains, ())
            self.assertEqual(reason, "Unresolved total or incomplete chain annotation")
        counts, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PC", "PC 34:1", True),
            ("OTOT", "PC", "PC 34:1_4:0", True),
            ("OTOT", "FA", "FA 34:1", True),
            ("OTOT", "LPC", "LPC 34:1/0:0", True),
        ])
        self.assertEqual(list(counts.index), ["4:0", "34:1"])
        self.assertEqual(denominators["OTOT"], 3)
        self.assertEqual(counts.at["34:1", "OTOT"], 3)
        self.assertEqual(percentages.at["34:1", "OTOT"], 100)

    def test_all_ids_requires_explicit_flag_for_unvalidated_csv(self):
        source = pd.DataFrame([("OTOT", "PE", "PE 16:0_16:0")],
                              columns=["Dataset", "Ontology", "Metabolite name"])
        with patch.object(self.script.pd, "read_csv", return_value=source):
            with self.assertRaisesRegex(ValueError, "--all-ids"):
                self.script.prepare_lipids(Path("raw.csv"), self.script.DATASETS)
            lipids, _, _ = self.script.prepare_lipids(Path("raw.csv"), self.script.DATASETS, True)
        self.assertEqual(len(lipids), 1)

    def test_heatmap_has_one_scale_percent_annotations_and_method_denominators(self):
        _, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("ITIT", "PE", "PE 20:0_20:0", True),
        ])
        with patch.object(self.script, "save_figure") as save:
            fig = self.script.plot_tail_prevalence(percentages, denominators, Path("output"), "Accepted", 1)
        save.assert_called_once_with(fig, Path("output"), "Figure5")
        ax = fig.axes[0]
        self.assertEqual((ax.images[0].norm.vmin, ax.images[0].norm.vmax), (0, 100))
        mask = np.ma.getmaskarray(ax.images[0].get_array())
        self.assertFalse(mask[:, 0].any())
        self.assertTrue(mask[:, 1].all())
        labels = [text.get_text() for text in ax.texts]
        self.assertIn("0.0%", labels)
        self.assertIn("100.0%", labels)
        self.assertEqual(labels.count("N/A"), len(percentages) * 3)
        self.assertEqual([label.get_text() for label in ax.get_xticklabels()], [
            f"{self.script.DATASET_LABELS[d]}\nN = {denominators[d]}" for d in self.script.DATASETS
        ])
        self.assertEqual([label.get_text() for label in ax.get_yticklabels()], list(percentages.index))

    def test_empty_selection_still_shows_all_method_labels(self):
        _, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("OTOT", "PC", "PC 20:0_22:0", True),
        ], 100)
        with patch.object(self.script, "save_figure"):
            fig = self.script.plot_tail_prevalence(percentages, denominators, Path("output"), "Accepted", 100)
        self.assertEqual(len(fig.axes[0].get_xticklabels()), 5)
        self.assertIn("No eligible tails meet the threshold", [text.get_text() for text in fig.axes[0].texts])

    def test_main_exports_heatmap_tables_and_separate_overlap_figure(self):
        rows = [
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("ITIT", "PE", "PE 20:0_20:0", True),
        ]
        with TemporaryDirectory() as directory:
            output_dir = Path(directory)
            source = pd.DataFrame(rows, columns=[
                "Dataset", "Ontology", "Metabolite name", "Overall Match Status",
            ])
            input_path = output_dir / "input.csv"
            source.to_csv(input_path, index=False)
            with patch.object(sys, "argv", ["script", "--input", str(input_path),
                                             "--output-dir", directory, "--no-show"]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.script.main()
            for name in ["Figure5", "Whole_Lipid_Overlap_by_Tail"]:
                for extension in ["png", "pdf", "svg"]:
                    self.assertGreater((output_dir / f"{name}.{extension}").stat().st_size, 0)
            with Image.open(output_dir / "Figure5.png") as image:
                for dpi in image.info["dpi"]:
                    self.assertAlmostEqual(dpi, 600, places=1)
            counts = pd.read_csv(output_dir / "Figure5_Tail_Prevalence_Counts.csv", index_col="Tail")
            denominators = pd.read_csv(output_dir / "Figure5_Tail_Prevalence_Denominators.csv", index_col="Dataset")
            percentages = pd.read_csv(output_dir / "Figure5_Tail_Prevalence_Percentages.csv", index_col="Tail")
            self.assertEqual(counts.at["16:0", "OTOT"], 2)
            self.assertEqual(denominators.at["OTOT", "Resolved lipids"], 2)
            self.assertEqual(denominators.at["OTIT", "Resolved lipids"], 0)
            self.assertEqual(percentages.at["18:1", "OTOT"], 50)
            self.assertEqual(percentages.at["16:0", "ITIT"], 0)
            self.assertTrue(percentages["OTIT"].isna().all())
            self.assertEqual(denominators["Acquisition Method"].tolist(),
                             list(self.script.DATASET_LABELS.values()))
            self.assertFalse(self.script.plt.get_fignums())

    def test_main_retains_denominators_when_no_tails_are_displayed(self):
        for rows, threshold, expected_denominator in [
            ([("OTOT", "EtherPC", "PC O-10:0_O-12:0", True)], 1, 1),
            ([("OTOT", "PC", "PC 16:0_18:1", True),
              ("OTOT", "PC", "PC 20:0_22:0", True)], 100, 2),
        ]:
            with self.subTest(threshold=threshold), TemporaryDirectory() as directory:
                output_dir = Path(directory)
                prepared = self.prepare(rows)
                with patch.object(self.script, "prepare_lipids", return_value=prepared), \
                        patch.object(self.script, "save_figure"), \
                        patch.object(sys, "argv", ["script", "--output-dir", directory,
                                                   "--min-tail-percent", str(threshold), "--no-show"]), \
                        contextlib.redirect_stdout(io.StringIO()):
                    self.script.main()
                denominators = pd.read_csv(output_dir / "Figure5_Tail_Prevalence_Denominators.csv", index_col="Dataset")
                self.assertEqual(denominators.at["OTOT", "Resolved lipids"], expected_denominator)
                self.assertEqual(len(denominators), 5)
                for suffix in ["Counts", "Percentages"]:
                    table = pd.read_csv(output_dir / f"Figure5_Tail_Prevalence_{suffix}.csv", index_col="Tail")
                    self.assertTrue(table.empty)
                    self.assertEqual(list(table.columns), self.script.DATASETS)


if __name__ == "__main__":
    unittest.main()

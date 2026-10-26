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


class TestSITailHeatmap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_dir = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                      / "DDA Paper Figures")
        for attribute, filename in [("script", "SI_5.py"), ("figure4", "Figure5.py")]:
            spec = importlib.util.spec_from_file_location(attribute, script_dir / filename)
            module = importlib.util.module_from_spec(spec)
            with patch.object(sys, "path", [str(script_dir), *sys.path]):
                spec.loader.exec_module(module)
            setattr(cls, attribute, module)

    def tearDown(self):
        self.script.plt.close("all")

    def prepare(self, rows):
        source = pd.DataFrame(rows, columns=[
            "Dataset", "Ontology", "Metabolite name", "Overall Match Status",
        ])
        with patch.object(self.script.pd, "read_excel", return_value=source):
            return self.script.prepare_lipids(Path("input.xlsx"), self.script.DATASETS)

    def tables(self, rows, threshold=1):
        lipids, occurrences, audit = self.prepare(rows)
        summary = self.script.frequency_summary(lipids, occurrences, ["Dataset"])
        counts, denominators, percentages = self.script.tail_prevalence_tables(
            summary, self.script.DATASETS, threshold
        )
        return counts, denominators, percentages, occurrences, summary

    def test_repeated_tails_count_each_copy_after_lipid_deduplication(self):
        counts, denominators, percentages, _, summary = self.tables([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            (" otot\u00a0", "PE", "PE 16:0/16:0", "TRUE"),
            ("OTOT", "TG", "TG 16:0_18:1_18:1", True),
            ("OTOT", "CL", "CL 16:0_16:0_16:0_16:0", True),
            ("OTOT", "PC", "PC 34:1", True),
            ("Stellar", "PE", "PE 22:0_22:0", False),
        ])
        self.assertEqual(denominators["OTOT"], 9)
        self.assertEqual(counts.at["16:0", "OTOT"], 7)
        self.assertEqual(counts.at["18:1", "OTOT"], 2)
        self.assertAlmostEqual(percentages.at["16:0", "OTOT"], 700 / 9)
        self.assertAlmostEqual(percentages["OTOT"].sum(), 100)
        pd.testing.assert_frame_equal(
            percentages, self.figure4.tail_composition_table(summary, self.script.DATASETS)
        )

    def test_selection_uses_tail_copy_percentages_without_renormalizing(self):
        rows = [
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "TG", "TG 18:1_18:1_18:1", True),
        ]
        counts, denominators, percentages, _, _ = self.tables(rows, threshold=50)
        self.assertEqual(list(counts.index), ["18:1"])
        self.assertEqual(counts.at["18:1", "OTOT"], 3)
        self.assertEqual(denominators["OTOT"], 5)
        self.assertEqual(percentages.at["18:1", "OTOT"], 60)
        empty_counts, empty_denominators, empty_percentages, _, _ = self.tables(rows, threshold=100)
        self.assertTrue(empty_counts.empty)
        self.assertTrue(empty_percentages.empty)
        self.assertEqual(empty_denominators["OTOT"], 5)

    def test_hidden_tails_remain_in_total_tail_copy_denominator(self):
        counts, denominators, percentages, _, summary = self.tables([
            ("OTOT", "EtherPC", "PC O-16:0_18:1", True),
            ("OTOT", "Cer", "Cer 18:1;O2/24:0", True),
            ("OTOT", "EtherPE", "PE P-16:0_18:1", True),
            ("OTOT", "OxPC", "PC 22:3;O_18:1", True),
        ])
        self.assertEqual(denominators["OTOT"], 8)
        self.assertEqual(list(counts.index), ["18:1", "22:3;O", "24:0"])
        self.assertEqual(percentages.at["18:1", "OTOT"], 37.5)
        self.assertEqual(percentages.at["22:3;O", "OTOT"], 12.5)
        self.assertEqual(percentages["OTOT"].sum(), 62.5)
        reference = self.figure4.tail_composition_table(summary, self.script.DATASETS)
        pd.testing.assert_frame_equal(percentages, reference, check_exact=True)

    def test_tail_selection_and_percentages_match_figure4_at_every_cutoff(self):
        rows = [
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("OTOT", "EtherPE", "PE P-16:0_18:1", True),
            ("ITIT", "OxPC", "PC 22:3;O_18:1", True),
            ("ITIT", "EtherPC", "PC O-16:0_18:1", True),
            ("Stellar", "Cer", "Cer 18:1;O2/24:0", True),
        ]
        lipids, occurrences, audit = self.prepare(rows)
        source = pd.DataFrame(rows, columns=[
            "Dataset", "Ontology", "Metabolite name", "Overall Match Status",
        ])
        with patch.object(self.figure4.pd, "read_excel", return_value=source):
            reference_inputs = self.figure4.prepare_lipids(Path("input.xlsx"), self.script.DATASETS)
        for actual, expected in zip((lipids, occurrences, audit), reference_inputs):
            pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        summary = self.script.frequency_summary(lipids, occurrences, ["Dataset"])
        reference_summary = self.figure4.frequency_summary(*reference_inputs[:2], ["Dataset"])
        pd.testing.assert_frame_equal(summary, reference_summary, check_exact=True)
        for threshold in [0, 1, 5, 50, 100]:
            with self.subTest(threshold=threshold):
                _, _, percentages = self.script.tail_prevalence_tables(summary, self.script.DATASETS, threshold)
                reference = self.figure4.tail_composition_table(reference_summary, self.script.DATASETS, threshold)
                pd.testing.assert_frame_equal(percentages, reference, check_exact=True)

    def test_instrument_detection_exports_use_figure4_tail_catalog(self):
        _, occurrences, _ = self.prepare([
            ("OTOT", "EtherPE", "PE P-16:0_18:1", True),
            ("ITIT", "OxPC", "PC 22:3;O_18:1", True),
            ("Stellar", "Cer", "Cer 18:1;O2/24:0", True),
        ])
        actual = self.script.instrument_tail_tables(occurrences, self.script.DATASETS)
        expected = self.figure4.instrument_tail_tables(occurrences, self.script.DATASETS)
        for actual_table, expected_table in zip(actual, expected):
            pd.testing.assert_frame_equal(actual_table, expected_table, check_exact=True)

    def test_absent_tails_are_zero_and_missing_datasets_are_na(self):
        counts, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("ITIT", "TG", "TG 18:1_18:1_18:1", True),
        ])
        self.assertEqual(counts.at["16:0", "ITIT"], 0)
        self.assertEqual(percentages.at["16:0", "ITIT"], 0)
        self.assertEqual(percentages.at["18:1", "ITIT"], 100)
        self.assertEqual(denominators["ITIT"], 3)
        for dataset in ["OTIT", "Stellar", "Stellarboth"]:
            self.assertEqual(denominators[dataset], 0)
            self.assertTrue(counts[dataset].isna().all())
            self.assertTrue(percentages[dataset].isna().all())

    def test_pooled_coverage_and_overlap_count_shared_tail_copies_once(self):
        _, occurrences, _ = self.prepare([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("ITIT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("ITIT", "PC", "PC 18:1_16:0", True),
            ("ITIT", "TG", "TG 16:0_16:0_16:0", True),
        ])
        duplicated = pd.concat([occurrences, occurrences.iloc[[0]]], ignore_index=True)
        for tails in [["16:0", "18:1"], ["20:0"], []]:
            with self.subTest(tails=tails):
                actual = self.script.tail_coverage_tables(duplicated, self.script.DATASETS, tails)
                expected = self.figure4.tail_coverage_tables(duplicated, self.script.DATASETS, tails)
                for actual_table, expected_table in zip(actual, expected):
                    pd.testing.assert_frame_equal(actual_table, expected_table)
                actual = self.script.tail_overlap_tables(duplicated, self.script.DATASETS, tails)
                expected = self.figure4.tail_overlap_tables(duplicated, self.script.DATASETS, tails)
                for actual_table, expected_table in zip(actual, expected):
                    pd.testing.assert_frame_equal(actual_table, expected_table)

    def test_heatmap_labels_and_annotations_describe_tail_copies(self):
        _, denominators, percentages, _, _ = self.tables([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "PC", "PC 16:0_18:1", True),
        ])
        with patch.object(self.script, "save_figure") as save:
            fig = self.script.plot_tail_prevalence(percentages, denominators, Path("output"), "Accepted", 1)
        save.assert_called_once_with(fig, Path("output"), "SI_TailHeatmap")
        ax = fig.axes[0]
        self.assertEqual((ax.images[0].norm.vmin, ax.images[0].norm.vmax), (0, 100))
        self.assertIn("75.0%", [text.get_text() for text in ax.texts])
        self.assertIn("25.0%", [text.get_text() for text in ax.texts])
        self.assertEqual(fig.axes[1].get_ylabel(), "Accepted Tail Copies (%)")
        self.assertEqual(ax.get_title(), "Tail Composition of Accepted Lipids by Dataset")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), "OT/OT\nTail Copies = 4")

    def test_main_exports_tail_copy_counts_denominators_and_percentages(self):
        prepared = self.prepare([
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("OTOT", "PC", "PC 16:0_18:1", True),
        ])
        with TemporaryDirectory() as directory:
            output_dir = Path(directory)
            with patch.object(self.script, "prepare_lipids", return_value=prepared), \
                    patch.object(self.script, "save_figure"), \
                    patch.object(sys, "argv", ["script", "--output-dir", directory, "--no-show"]), \
                    contextlib.redirect_stdout(io.StringIO()) as output:
                self.script.main()
            denominators = pd.read_csv(output_dir / "Tail_Prevalence_Denominators.csv", index_col="Dataset")
            counts = pd.read_csv(output_dir / "Tail_Prevalence_Counts.csv", index_col="Tail")
            percentages = pd.read_csv(output_dir / "Tail_Prevalence_Percentages.csv", index_col="Tail")
            self.assertEqual(denominators.at["OTOT", "Total tail occurrences"], 4)
            self.assertEqual(counts.at["16:0", "OTOT"], 3)
            self.assertEqual(percentages.at["16:0", "OTOT"], 75)
            self.assertTrue(percentages["OTIT"].isna().all())
            self.assertIn("Within-method tail-copy composition", output.getvalue())
            self.assertIn("Tail copies (%)", output.getvalue())


if __name__ == "__main__":
    unittest.main()

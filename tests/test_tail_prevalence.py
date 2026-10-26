import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd


class TestTailPrevalence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_dir = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Paper Scripts"
        spec = importlib.util.spec_from_file_location(
            "tail_prevalence_script", script_dir / "Tail_Prevalence_by_Class.py"
        )
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(script_dir), *sys.path]):
            spec.loader.exec_module(cls.script)

    def load_lipids(self, rows):
        frame = pd.DataFrame(rows, columns=[
            "Dataset", "Ontology", "Metabolite name", "Overall Match Status"
        ])
        with patch.object(self.script.pd, "read_excel", return_value=frame):
            return self.script.accepted_lipids(Path("input.xlsx"))

    def test_repeated_tail_counts_once_in_each_lipid(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 36:2|PC 18:1_18:1", True),
            ("Stellar", "TG", "TG 54:3|TG 18:1_18:1_18:1", True),
        ])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(counts.at["18:1", "OTOT"], 1)
        self.assertEqual(counts.at["18:1", "Stellar"], 1)
        self.assertEqual(counts.to_numpy().sum(), 2)
        self.assertEqual(pooled.at["18:1"], 2)

    def test_features_and_shared_lipids_deduplicate_at_the_right_level(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 34:1|PC 16:0_18:1", True),
            (" otot\u00a0", "PC", "PC 18:1/16:0", "TRUE"),
            ("itit", "PC", "PC 34:1|PC 18:1_16:0", True),
        ])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(len(lipids), 1)
        for tail in ["16:0", "18:1"]:
            with self.subTest(tail=tail):
                self.assertEqual(counts.at[tail, "OTOT"], 1)
                self.assertEqual(counts.at[tail, "ITIT"], 1)
                self.assertEqual(counts.loc[tail].sum(), 2)
                self.assertEqual(pooled.at[tail], 1)

    def test_unresolved_and_rejected_lipids_do_not_add_tails(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 34:1", True),
            ("OTOT", "PC", "PC 26:0_28:0", False),
            ("Stellar", "PC", "PC 16:0_18:1", True),
        ])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(set(counts.index), {"16:0", "18:1"})
        self.assertEqual(counts["OTOT"].sum(), 0)
        self.assertEqual(counts["Stellar"].sum(), 2)
        self.assertEqual(pooled.sum(), 2)
        self.assertEqual(self.script.class_counts(lipids).at["PC"], 2)

    def test_tail_modifications_and_backbones_remain_distinct(self):
        lipids = self.load_lipids([
            ("Stellar", "PC", "PC 16:0_18:1", True),
            ("Stellar", "EtherPC", "PC O-16:0_18:1", True),
            ("Stellar", "Cer", "Cer 18:1;O2/18:1", True),
        ])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(set(counts.index), {"16:0", "18:1", "O-16:0", "LCB 18:1;O2"})
        self.assertEqual(pooled.at["18:1"], 3)
        for tail in ["16:0", "O-16:0", "LCB 18:1;O2"]:
            self.assertEqual(counts.at[tail, "Stellar"], 1)
        self.assertEqual(set(self.script.figure4_tails(lipids)), {"16:0", "18:1"})

    def test_figure4_selection_uses_copies_but_plot_counts_use_lipids(self):
        rows = [
            ("OTOT", "PC", f"PC {carbon}:{db}_18:1", True)
            for carbon in range(14, 35) for db in range(5)
        ]
        rows.append(("OTOT", "CL", "CL 14:0_14:0_14:0_14:0", True))
        lipids = self.load_lipids(rows)
        self.assertEqual(set(self.script.figure4_tails(lipids)), {"14:0", "18:1"})
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(counts.at["14:0", "OTOT"], 2)
        self.assertEqual(pooled.at["14:0"], 2)

    def test_different_classes_count_as_different_lipids(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("OTOT", "PE", "PE 16:0_18:1", True),
        ])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertEqual(counts.at["16:0", "OTOT"], 2)
        self.assertEqual(pooled.at["18:1"], 2)

    def test_all_unresolved_retains_class_counts_and_empty_tail_counts(self):
        lipids = self.load_lipids([("OTOT", "PC", "PC 34:1", True)])
        counts, pooled = self.script.tail_counts(lipids)
        self.assertTrue(counts.empty)
        self.assertTrue(pooled.empty)
        self.assertTrue(self.script.figure4_tails(lipids).empty)
        self.assertEqual(list(counts.columns), self.script.DATASETS)
        self.assertEqual(self.script.class_counts(lipids).at["PC"], 1)

    def test_method_prevalence_uses_unique_resolved_lipids_as_denominators(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 34:1|PC 16:0_18:1", True),
            ("otot", "PC", "PC 18:1/16:0", True),
            ("OTOT", "PC", "PC 16:0_16:0", True),
            ("OTOT", "PE", "PE 16:0_16:0", True),
            ("ITIT", "PC", "PC 16:0_18:1", True),
            ("OTIT", "PC", "PC 34:1", True),
            ("Stellar", "PC", "PC 26:0_28:0", False),
            ("Stellarboth", "CE", "CE 30:0", True),
        ])
        counts, denominators, percentages = self.script.tail_method_tables(lipids)
        self.assertEqual(denominators.to_dict(), {
            "OTOT": 3, "OTIT": 0, "ITIT": 1, "Stellar": 0, "Stellarboth": 1,
        })
        self.assertEqual(counts.at["16:0", "OTOT"], 3)
        self.assertEqual(counts.at["18:1", "OTOT"], 1)
        self.assertAlmostEqual(percentages.at["18:1", "OTOT"], 100 / 3)
        self.assertEqual(percentages.at["16:0", "OTOT"], 100)
        self.assertEqual(percentages.at["18:1", "ITIT"], 100)
        self.assertGreater(percentages["ITIT"].sum(), 100)
        self.assertEqual(counts.at["16:0", "Stellarboth"], 0)
        self.assertEqual(percentages.at["16:0", "Stellarboth"], 0)
        for dataset in ["OTIT", "Stellar"]:
            self.assertTrue(counts[dataset].isna().all())
            self.assertTrue(percentages[dataset].isna().all())
        self.assertEqual(list(counts.columns), self.script.DATASETS)

    def test_parser_excludes_totals_by_class_chain_count_not_chain_length(self):
        for lipid_class, annotation in [
            ("PC", "PC 34:1"), ("PE", "PE 32:0"), ("DG", "DG 34:1"),
            ("TG", "TG 54:3"), ("CL", "CL 72:4"), ("EtherPC", "PC O-34:1"),
            ("TG", "TG 16:0_18:1"), ("PC", "PC 0:0_16:0"),
        ]:
            with self.subTest(annotation=annotation):
                _, chains, reason = self.script.parse_lipid_tails(annotation, lipid_class)
                self.assertEqual(chains, ())
                self.assertEqual(reason, "Unresolved total or incomplete chain annotation")
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 34:1", True),
            ("OTOT", "PC", "PC 34:1_4:0", True),
            ("OTOT", "FA", "FA 34:1", True),
            ("OTOT", "LPC", "LPC 34:1/0:0", True),
        ])
        counts, denominators, percentages = self.script.tail_method_tables(lipids)
        self.assertEqual(list(counts.index), ["4:0", "34:1"])
        self.assertEqual(denominators["OTOT"], 3)
        self.assertEqual(counts.at["34:1", "OTOT"], 3)
        self.assertEqual(percentages.at["34:1", "OTOT"], 100)

    def test_method_prevalence_preserves_tail_multiplicity_in_lipid_identity(self):
        lipids = self.load_lipids([
            ("OTOT", "TG", "TG 16:0_18:0_18:0", True),
            ("OTOT", "TG", "TG 18:0/16:0/18:0", True),
            ("OTOT", "TG", "TG 16:0_16:0_18:0", True),
        ])
        counts, denominators, percentages = self.script.tail_method_tables(lipids)
        self.assertEqual(denominators["OTOT"], 2)
        self.assertEqual(counts.at["16:0", "OTOT"], 2)
        self.assertEqual(counts.at["18:0", "OTOT"], 2)
        self.assertEqual(percentages["OTOT"].sum(), 200)

    def test_non_acyl_chains_are_omitted_without_changing_denominators(self):
        lipids = self.load_lipids([
            ("OTOT", "EtherPC", "PC O-16:0_18:1", True),
            ("OTOT", "EtherPE", "PE P-16:0_20:1", True),
            ("OTOT", "Cer", "Cer 18:1;O2/24:0", True),
            ("OTOT", "EtherPC", "PC O-10:0_O-12:0", True),
            ("ITIT", "EtherPC", "PC O-10:0_O-12:0", True),
        ])
        counts, denominators, percentages = self.script.tail_method_tables(lipids)
        self.assertEqual(set(counts.index), {"18:1", "20:1", "24:0"})
        self.assertEqual(denominators["OTOT"], 4)
        self.assertEqual(denominators["ITIT"], 1)
        self.assertTrue(percentages["OTOT"].eq(25).all())
        self.assertTrue(counts["ITIT"].eq(0).all())
        self.assertTrue(percentages["ITIT"].eq(0).all())

    def test_method_prevalence_retains_rare_and_modified_tails(self):
        rows = [
            ("OTOT", "PC", f"PC {carbon}:{db}_18:1", True)
            for carbon in range(14, 35) for db in range(5)
        ]
        rows.append(("OTOT", "PC", "PC 16:0_18:1;O2", True))
        lipids = self.load_lipids(rows)
        counts, denominators, percentages = self.script.tail_method_tables(lipids)
        self.assertEqual(denominators["OTOT"], 106)
        self.assertEqual(counts.at["34:4", "OTOT"], 1)
        self.assertEqual(counts.at["18:1;O2", "OTOT"], 1)
        self.assertLess(percentages.at["34:4", "OTOT"], 1)
        self.assertNotIn("34:4", self.script.figure4_tails(lipids))

    def test_heatmap_labels_percentages_zeros_and_missing_data_on_one_scale(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("ITIT", "PE", "PE 20:0_20:0", True),
        ])
        _, denominators, percentages = self.script.tail_method_tables(lipids)
        fig = self.script.plot_tail_method_heatmap(percentages, denominators)
        try:
            ax = fig.axes[0]
            heatmap = ax.images[0]
            self.assertEqual((heatmap.norm.vmin, heatmap.norm.vmax), (0, 100))
            mask = np.ma.getmaskarray(heatmap.get_array())
            self.assertFalse(mask[:, 0].any())
            self.assertTrue(mask[:, 1].all())
            texts = [item.get_text() for item in ax.texts]
            self.assertIn("0.0%", texts)
            self.assertIn("100.0%", texts)
            self.assertEqual(texts.count("NA"), len(percentages) * 3)
            self.assertEqual(len(texts), len(percentages) * 5)
            self.assertEqual([label.get_text() for label in ax.get_yticklabels()],
                             list(percentages.index))
            labels = [label.get_text() for label in ax.get_xticklabels()]
            self.assertEqual(len(labels), 5)
            self.assertIn("N = 1", labels[0])
            self.assertIn("N = 0", labels[1])
        finally:
            self.script.plt.close(fig)

    def test_empty_heatmap_retains_all_method_denominators_in_csv(self):
        for annotation, lipid_class, expected_denominator in [
            ("PC 34:1", "PC", 0), ("PC O-10:0_O-12:0", "EtherPC", 1),
        ]:
            with self.subTest(annotation=annotation):
                lipids = self.load_lipids([("OTOT", lipid_class, annotation, True)])
                counts, denominators, percentages = self.script.tail_method_tables(lipids)
                self.assertTrue(counts.empty)
                self.assertTrue(percentages.empty)
                summary = self.script.tail_method_summary(counts, denominators, percentages)
                self.assertEqual(summary["Dataset"].tolist(), self.script.DATASETS)
                self.assertEqual(summary.iloc[0]["Resolved lipid IDs"], expected_denominator)
                self.assertTrue(summary["Lipids with tail (%)"].isna().all())
                self.assertTrue(summary["Lipids containing tail"].isna().all())
                fig = self.script.plot_tail_method_heatmap(percentages, denominators)
                try:
                    self.assertEqual(len(fig.axes[0].get_xticklabels()), 5)
                finally:
                    self.script.plt.close(fig)

    def test_main_exports_new_heatmap_and_preserves_existing_figure_outputs(self):
        lipids = self.load_lipids([
            ("OTOT", "PC", "PC 16:0_18:1", True),
            ("ITIT", "PE", "PE 20:0_20:0", True),
        ])
        with TemporaryDirectory() as directory:
            output_dir = Path(directory)
            with patch.object(self.script, "accepted_lipids", return_value=lipids), \
                    patch.object(sys, "argv", ["script", "--output-dir", directory, "--no-show"]):
                self.script.main()
            for filename in [
                "Pooled_Accepted_Lipid_IDs.csv", "Accepted_Lipids_by_Class.csv",
                "Tail_Prevalence_by_Dataset.csv", "Tail_Prevalence_by_Dataset.pdf",
                "Tail_Prevalence_by_Dataset_Page1.png", "Tail_Prevalence_by_Method_Heatmap.csv",
                "Tail_Prevalence_by_Method_Heatmap.pdf", "Tail_Prevalence_by_Method_Heatmap.png",
            ]:
                self.assertGreater((output_dir / filename).stat().st_size, 0)
            summary = pd.read_csv(output_dir / "Tail_Prevalence_by_Method_Heatmap.csv")
            self.assertEqual(len(summary), 15)
            row = summary.loc[summary["Tail"].eq("16:0") & summary["Dataset"].eq("OTOT")].iloc[0]
            self.assertEqual(row["Lipids containing tail"], 1)
            self.assertEqual(row["Resolved lipid IDs"], 1)
            self.assertEqual(row["Lipids with tail (%)"], 100)
            missing = summary.loc[summary["Dataset"].eq("OTIT")]
            self.assertTrue(missing["Lipids with tail (%)"].isna().all())
            self.assertTrue(missing["Lipids containing tail"].isna().all())
            self.assertTrue(missing["Data status"].eq("No resolved lipid IDs").all())
            self.assertFalse(self.script.plt.get_fignums())


if __name__ == "__main__":
    unittest.main()

import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import pandas as pd


class TestFigure5TailCopies(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_path = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                       / "DDA Paper Figures" / "Figure5.py")
        spec = importlib.util.spec_from_file_location("figure4_tail_copies", script_path)
        cls.script = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.script)

    def prepare(self, rows, datasets):
        source = pd.DataFrame(rows, columns=["Dataset", "Ontology", "Metabolite name"])
        with patch.object(self.script.pd, "read_csv", return_value=source):
            return self.script.prepare_lipids(Path("input.csv"), datasets, all_ids=True)

    def test_identical_tails_each_contribute_one_count(self):
        for lipid_class, name, copies in [
            ("PE", "PE 16:0_16:0", 2),
            ("TG", "TG 16:0_16:0_16:0", 3),
            ("CL", "CL 16:0_16:0_16:0_16:0", 4),
        ]:
            with self.subTest(lipid=name):
                _, occurrences, _ = self.prepare([("OTOT", lipid_class, name)], ["OTOT"])
                counts, percentages, _ = self.script.tail_coverage_tables(
                    occurrences, ["OTOT"], ["16:0"]
                )
                overlap, _ = self.script.tail_overlap_tables(occurrences, ["OTOT"], ["16:0"])
                self.assertEqual(counts.at["16:0", 1], copies)
                self.assertEqual(percentages.at["16:0", 1], 100)
                self.assertEqual(overlap.at["OT/OT Only", "16:0"], copies)

    def test_shared_lipids_and_duplicate_features_do_not_repeat_their_copies(self):
        datasets = ["OTOT", "ITIT", "Stellar"]
        lipids, occurrences, audit = self.prepare([
            ("OTOT", "PE", "PE 16:0_16:0"),
            ("OTOT", "PE", "PE 16:0/16:0"),
            ("ITIT", "PE", "PE 16:0_16:0"),
            ("OTOT", "PC", "PC 16:0_18:1"),
            ("ITIT", "PC", "PC 18:1/16:0"),
            ("ITIT", "TG", "TG 16:0_16:0_16:0"),
        ], datasets)
        self.assertEqual(len(lipids), 5)
        self.assertEqual(audit["Exclusion reason"].eq("Duplicate chain composition within dataset").sum(), 1)
        # Also exercise duplicate occurrence rows at the table boundary.
        duplicated = pd.concat([occurrences, occurrences.iloc[[0]]], ignore_index=True)
        counts, percentages, membership = self.script.tail_coverage_tables(
            duplicated, datasets, ["16:0", "18:1"]
        )
        self.assertEqual(counts.loc["16:0"].tolist(), [3, 3, 0])
        self.assertEqual(counts.loc["18:1"].tolist(), [0, 1, 0])
        self.assertEqual(percentages.loc["16:0"].tolist(), [50, 50, 0])
        self.assertEqual(percentages.loc["18:1"].tolist(), [0, 100, 0])
        self.assertEqual(membership["Tail copies"].sum(), 7)
        overlap, overlap_membership = self.script.tail_overlap_tables(
            duplicated, datasets, ["16:0", "18:1"]
        )
        self.assertEqual(overlap.loc["2 Datasets"].to_dict(), {"16:0": 3, "18:1": 1})
        self.assertEqual(overlap.at["IT/IT Only", "16:0"], 3)
        self.assertEqual(overlap.loc["OT/OT Only"].sum(), 0)
        self.assertEqual(overlap.to_numpy().sum(), 7)
        self.assertEqual(overlap_membership["Tail copies"].sum(), 7)

    def test_ether_tails_are_hidden_and_companion_acyl_tails_still_count(self):
        datasets = ["OTOT", "ITIT"]
        lipids, occurrences, _ = self.prepare([
            ("OTOT", "EtherPE", "PE P-16:0_18:1"),
            ("ITIT", "EtherPE", "PE P-16:0_18:1"),
            ("OTOT", "EtherPC", "PC O-16:0_18:1"),
            ("OTOT", "Cer", "Cer 18:1;O2/24:0"),
            ("OTOT", "PE", "PE 16:0_16:0"),
            ("ITIT", "EtherLPC", "LPC P-16:0"),
        ], datasets)
        self.assertIn("P-16:0", occurrences["Tail"].values)
        overall = self.script.frequency_summary(lipids, occurrences, ["Dataset"])
        selected = self.script.tail_composition_table(overall, datasets)
        self.assertEqual(list(selected.index), ["16:0", "18:1", "24:0"])
        # Hidden chains remain in the full tail-copy denominator for selection.
        self.assertEqual(selected.at["18:1", "OTOT"], 25)
        self.assertAlmostEqual(selected.at["18:1", "ITIT"], 100 / 3)
        coverage, _, _ = self.script.tail_coverage_tables(occurrences, datasets, selected.index)
        overlap, _ = self.script.tail_overlap_tables(occurrences, datasets, selected.index)
        self.assertEqual(coverage.loc["18:1"].tolist(), [1, 1])
        self.assertEqual(overlap.at["2 Datasets", "18:1"], 1)
        self.assertEqual(overlap.at["OT/OT Only", "18:1"], 1)
        self.assertEqual(overlap.at["OT/OT Only", "16:0"], 2)
        counts, _, metadata = self.script.instrument_tail_tables(occurrences, datasets)
        self.assertEqual(list(counts.index), ["16:0", "18:1", "24:0"])
        self.assertEqual(counts.loc["18:1"].tolist(), [2, 1])
        self.assertTrue(metadata["Total pooled lipids"].eq(5).all())

    def test_absent_tails_and_empty_selections(self):
        datasets = ["OTOT", "ITIT"]
        _, occurrences, _ = self.prepare([("OTOT", "PE", "PE 16:0_16:0")], datasets)
        for tails in [["20:0"], []]:
            with self.subTest(tails=tails):
                counts, percentages, membership = self.script.tail_coverage_tables(
                    occurrences, datasets, tails
                )
                overlap, overlap_membership = self.script.tail_overlap_tables(occurrences, datasets, tails)
                self.assertEqual(list(counts.index), tails)
                self.assertEqual(list(counts.columns), [1, 2])
                self.assertEqual(counts.to_numpy().sum(), 0)
                self.assertTrue(percentages.isna().all().all())
                self.assertEqual(list(overlap.columns), tails)
                self.assertEqual(overlap.to_numpy().sum(), 0)
                self.assertTrue(membership.empty)
                self.assertTrue(overlap_membership.empty)


if __name__ == "__main__":
    unittest.main()

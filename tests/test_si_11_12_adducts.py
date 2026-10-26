import importlib.util
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestAdductPanels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location("si_11_12", directory / "SI_6_7.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]), patch.dict(os.environ, {"MPLBACKEND": "Agg"}):
            spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def source(self):
        return pd.DataFrame([
            {"Dataset": "OTOT", "Ontology": "TG", "Polarity": "positive",
             "Metabolite name": lipid, "Adduct type": adduct,
             "Overall Match Status": accepted, "Class Strength": strength}
            for lipid, adduct, strength, accepted in [
                ("TG 50:1", "[M+Na]+", 8, True),
                ("TG 50:1", "[M+Na]+", 12, True),
                ("TG 52:1", "[M+Na]+", 6, True),
                ("TG 50:1", "[M+NH4]+", 2, True),
                ("TG 52:1", "[M+NH4]+", 4, True),
                ("TG 50:1", "[M+NH4]+", 20, False),
            ]])

    def test_duplicates_and_normalization_are_separate_within_each_adduct(self):
        points = self.script.prepare_strength(self.source(), all_ids=False)
        self.assertEqual(len(points), 4)
        sodium = points.loc[points["Adduct type"].eq("[M+Na]+")].iloc[0]
        self.assertEqual(sodium["Class Strength"], 10)
        self.assertEqual(sodium["Feature rows"], 2)
        self.assertEqual(sodium["Normalization maximum"], 10)
        ammonium = points.loc[points["Adduct type"].eq("[M+NH4]+")]
        self.assertTrue(ammonium["Normalization maximum"].eq(4).all())
        np.testing.assert_allclose(sorted(ammonium[self.script.NORMALIZED]), [0.5, 1.0])
        self.assertTrue(self.script.normalization_audit(points)["Pass"].all())
        self.assertEqual(set(points["Plot Class"]), {"TG\n[M+Na]+", "TG\n[M+NH4]+"})

    def test_excluded_reference_keeps_status_and_adduct_separate(self):
        source = self.source()
        source.loc[source["Overall Match Status"].eq(True), "Class Strength"] = 100
        points = self.script.prepare_strength(source, excluded_only=True)
        self.assertEqual(len(points), 1)
        self.assertEqual(points.iloc[0]["Class Strength"], 20)
        self.assertEqual(points.iloc[0]["Normalization maximum"], 20)
        self.assertEqual(points.iloc[0]["Adduct type"], "[M+NH4]+")
        self.assertEqual(points.iloc[0][self.script.NORMALIZED], 1.0)

    def test_missing_adduct_is_excluded(self):
        source = self.source()
        source.loc[0, "Adduct type"] = None
        points = self.script.prepare_strength(source, all_ids=False)
        self.assertEqual(points["Feature rows"].sum(), 4)
        with self.assertRaisesRegex(ValueError, "Adduct type"):
            self.script.prepare_strength(source.drop(columns="Adduct type"))

    def test_panel_pagination_preserves_every_adduct(self):
        points = self.script.prepare_strength(self.source(), all_ids=False).iloc[[0]]
        points = pd.concat([points.assign(**{"Plot Class": f"Class{i}\n[M+Na]+"})
                            for i in range(13)], ignore_index=True)
        figures = self.script.plot_strength(points)
        self.assertEqual(set(figures), {"positive_Page_1", "positive_Page_2", "negative"})
        titles = [ax.get_title() for key, fig in figures.items() if key.startswith("positive")
                  for ax in fig.axes if ax.get_title()]
        self.assertEqual(len(titles), 13)
        self.assertEqual(set(titles), set(points["Plot Class"]))

    def test_export_has_separate_adduct_panels_and_comparisons(self):
        source = pd.concat([self.source().assign(Dataset=d) for d in self.script.datasets],
                           ignore_index=True)
        reference = self.script.prepare_strength(source)
        with TemporaryDirectory() as directory:
            output = Path(directory)
            self.script.export_population(source, output, output / "input.csv", False, True,
                                          reference=reference, min_lipids=2)
            prefix = self.script.OUTPUT_BASENAME
            points = self.script.prepare_strength(source, all_ids=False)
            stats = self.script.instrument_statistics(points, all_pairs=True)
            self.assertEqual(len(stats), 20)
            self.assertEqual(stats["Plot Class"].nunique(), 2)
            self.assertTrue(self.script.normalization_audit(points)["Pass"].all())
            self.assertEqual({p.suffix for p in output.iterdir()}, {".png", ".pdf"})
            for mode in ["Positive", "Negative"]:
                for extension in ["png", "pdf"]:
                    self.assertGreater((output / f"{prefix}_{mode}.{extension}").stat().st_size, 0)

    def network_source(self):
        return pd.DataFrame([
            {"Dataset": "OTOT", "Ontology": "TG", "Polarity": "Positive",
             "Metabolite name": f"TG {i}", "Adduct type": adduct,
             "Average Mz": 800.0, "MS/MS spectrum": "100:100 200:50", "Tolerance": 0.05,
             "Overall Match Status": True, "Class Strength": 999.0}
            for adduct in ["[M+Na]+", "[M+NH4]+"] for i in range(3)
        ])

    def test_networks_split_before_strength_and_ignore_pooled_scores(self):
        source = self.network_source()
        original = source.copy(deep=True)
        scored, audit = self.script.recompute_adduct_networks(source)
        # Six identical spectra in a pooled network would have strength 5;
        # each three-node adduct network instead has strength 2.
        np.testing.assert_allclose(scored["Class Strength"], 2.0)
        self.assertEqual(len(audit), 2)
        self.assertTrue(audit["Edges"].eq(3).all())
        self.assertTrue(audit["Cross-adduct edges"].eq(0).all())
        pd.testing.assert_frame_equal(source, original)

    def test_missing_spectra_or_small_network_never_reuses_old_strength(self):
        source = self.network_source()
        source.loc[0, "MS/MS spectrum"] = None
        scored, audit = self.script.recompute_adduct_networks(source)
        self.assertTrue(scored.loc[scored["Adduct type"].eq("[M+Na]+"), "Class Strength"].isna().all())
        self.assertTrue(audit["Status"].str.startswith("Skipped").any())
        with self.assertRaisesRegex(ValueError, "MS/MS spectrum"):
            self.script.recompute_adduct_networks(source.drop(columns="MS/MS spectrum"))

    def test_mixed_network_tolerance_is_rejected(self):
        source = self.network_source()
        source.loc[0, "Tolerance"] = 0.7
        with self.assertRaisesRegex(ValueError, "Tolerance"):
            self.script.recompute_adduct_networks(source)

    def test_dataset_tolerance_reaches_fragment_and_neutral_loss_matching(self):
        source = self.network_source().loc[lambda f: f["Adduct type"].eq("[M+Na]+")].copy()
        source["MS/MS spectrum"] = ["100:100", "100.4:100", "200:100"]
        source = pd.concat([source.assign(Dataset=d, Tolerance=tol)
                            for d, tol in [("OTOT", 0.05), ("Stellar", 0.7)]], ignore_index=True)
        scored, audit = self.script.recompute_adduct_networks(source)
        np.testing.assert_allclose(scored.loc[scored["Dataset"].eq("OTOT"), "Class Strength"], 0)
        np.testing.assert_allclose(scored.loc[scored["Dataset"].eq("Stellar"), "Class Strength"], [1, 1, 0])
        self.assertEqual(audit["Tolerance"].tolist(), [0.05, 0.7])

    def test_coverage_drops_panels_with_missing_datasets(self):
        source = pd.concat([self.source().assign(Dataset=d) for d in self.script.datasets],
                           ignore_index=True)
        source = source.loc[~(source["Dataset"].eq("OTOT") & source["Adduct type"].eq("[M+NH4]+"))]
        points = self.script.prepare_strength(source, all_ids=False)
        selected, coverage = self.script.comparable_points(points, min_lipids=2)
        self.assertEqual(set(selected["Adduct type"]), {"[M+Na]+"})
        missing = coverage.loc[coverage["Plot Class"].str.contains("NH4", regex=False)]
        self.assertFalse(missing["Retained"].any())
        self.assertTrue(missing["Status"].str.contains("No eligible lipids in: OTOT", regex=False).all())
        self.assertTrue(coverage["Maximum lipids across datasets"].eq(2).all())
        single_dataset = points.loc[points["Dataset"].eq("Stellar")]
        kept, _ = self.script.comparable_points(single_dataset, min_lipids=2)
        self.assertTrue(kept.empty)
        empty, _ = self.script.comparable_points(points, min_lipids=5)
        self.assertTrue(empty.empty)
        with self.assertRaises(ValueError):
            self.script.comparable_points(points, min_lipids=1)

    def test_zero_strength_groups_stay_zero_and_are_omitted(self):
        points = self.script.prepare_strength(self.source().assign(**{"Class Strength": 0}),
                                              all_ids=False)
        self.assertTrue(points[self.script.NORMALIZED].eq(0).all())
        self.assertTrue(self.script.normalization_audit(points)["Pass"].all())
        selected, _ = self.script.comparable_points(points, min_lipids=2)
        self.assertTrue(selected.empty)

    def cutoff_source(self, last_count=10):
        rows = []
        for dataset in self.script.datasets:
            count = last_count if dataset == self.script.datasets[-1] else 5
            for index in range(count):
                rows.append({"Dataset": dataset, "Ontology": "PC", "Polarity": "positive",
                             "Metabolite name": f"PC {index}:0", "Adduct type": "[M+H]+",
                             "Overall Match Status": True, "Class Strength": index + 1})
        return pd.DataFrame(rows)

    def test_cutoff_10_uses_largest_dataset_and_includes_exact_boundary(self):
        self.assertEqual(self.script.MIN_LIPIDS, 10)
        for count in [11, 10, 9, 0]:
            with self.subTest(last_dataset_count=count):
                points = self.script.prepare_strength(self.cutoff_source(count), all_ids=False)
                kept, coverage = self.script.comparable_points(points)
                self.assertEqual(len(coverage), 5)
                self.assertTrue(coverage["Maximum lipids across datasets"].eq(max(5, count)).all())
                if count >= 10:
                    self.assertEqual(len(kept), 20 + count)
                    self.assertTrue(coverage["Retained"].all())
                else:
                    self.assertTrue(kept.empty)
                    self.assertFalse(coverage["Retained"].any())

    def test_cutoff_counts_unique_lipids_after_status_and_provisional_filters(self):
        for column, value in [("Overall Match Status", False), ("Validation Status", "Provisional")]:
            with self.subTest(column=column):
                source = self.cutoff_source()
                source.loc[29, column] = value
                # Repeated annotations must not inflate nine identities to ten.
                repeated = source.iloc[[20]].copy()
                repeated["Metabolite name"] += "|duplicate"
                points = self.script.prepare_strength(pd.concat([source, repeated], ignore_index=True),
                                                      all_ids=False)
                kept, coverage = self.script.comparable_points(points)
                self.assertTrue(kept.empty)
                self.assertEqual(coverage.loc[coverage["Dataset"].eq(self.script.datasets[-1]),
                                               "Lipids"].iloc[0], 9)

    def test_all_zero_largest_dataset_does_not_qualify_panel(self):
        source = self.cutoff_source()
        source.loc[source["Dataset"].eq(self.script.datasets[-1]), "Class Strength"] = 0
        points = self.script.prepare_strength(source, all_ids=False)
        kept, coverage = self.script.comparable_points(points)
        self.assertTrue(kept.empty)
        self.assertTrue(coverage["Maximum lipids across datasets"].eq(5).all())
        self.assertEqual(coverage["Zero-strength lipids removed"].sum(), 10)

    def test_zero_dataset_excludes_panel_even_when_another_reaches_ten(self):
        for mode in ["missing", "all-zero", "singleton"]:
            with self.subTest(mode=mode):
                source = self.cutoff_source()
                first_dataset = source["Dataset"].eq(self.script.datasets[0])
                if mode == "missing":
                    source = source.loc[~first_dataset]
                elif mode == "all-zero":
                    source.loc[first_dataset, "Class Strength"] = 0
                else:
                    source = source.loc[~first_dataset | source.index.isin([0])]
                points = self.script.prepare_strength(source, all_ids=False)
                kept, coverage = self.script.comparable_points(points)
                self.assertTrue(coverage["Maximum lipids across datasets"].eq(10).all())
                if mode == "singleton":
                    self.assertEqual(len(kept), 26)
                    self.assertTrue(coverage["Retained"].all())
                else:
                    self.assertTrue(kept.empty)
                    self.assertFalse(coverage["Retained"].any())
                    self.assertTrue(coverage["Status"].str.contains("No eligible lipids in: OTOT", regex=False).all())

    def test_zero_original_class_is_removed_as_its_own_panel(self):
        source = pd.concat([self.source().assign(Dataset=d) for d in self.script.datasets],
                           ignore_index=True)
        zero = source.assign(Ontology="EtherTG", **{"Class Strength": 0})
        points = self.script.prepare_strength(pd.concat([source, zero], ignore_index=True), all_ids=False)
        selected, coverage = self.script.comparable_points(points, min_lipids=2)
        self.assertEqual(set(selected["Ontology"]), {"TG"})
        self.assertTrue(self.script.normalization_audit(selected)["Normalized maximum"].eq(1).all())
        self.assertTrue(coverage["Zero-strength lipids removed"].gt(0).any())

    def test_original_classes_stay_separate_in_points_summaries_and_statistics(self):
        classes = [prefix + parent for parent in ["PC", "PE", "PG", "PI", "PS", "PA"]
                   for prefix in ["", "L", "Ether", "EtherL", "Ox"]]
        classes += ["NAPE", "LNAPE", "NAPS", "LNAPS"]
        source = pd.concat([self.source().assign(Dataset=d, Ontology=lipid_class)
                            for d in ["OTOT", "Stellar"]
                            for lipid_class in classes], ignore_index=True)
        points = self.script.prepare_strength(source, all_ids=False)
        self.assertEqual(set(points["Class Group"]), set(classes))
        for lipid_class in classes:
            selected = points.loc[points["Plot Class"].str.startswith(lipid_class + "\n")]
            self.assertEqual(set(selected["Ontology"]), {lipid_class})
        expected_panels = {lipid_class + "\n" + adduct
                           for lipid_class in classes
                           for adduct in ["[M+Na]+", "[M+NH4]+"]}
        self.assertEqual(set(self.script.strength_summary(points)["Plot Class"]), expected_panels)
        self.assertEqual(set(self.script.instrument_statistics(points, all_pairs=True)["Plot Class"]),
                         expected_panels)

    def test_main_recomputes_networks_and_exports_audits(self):
        source = pd.concat([self.network_source().assign(Dataset=d) for d in self.script.datasets],
                           ignore_index=True)
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source.to_csv(root / "input.csv", index=False)
            output = root / "out"
            with patch.object(self.script, "export_population", wraps=self.script.export_population) as export:
                self.script.main(["--input", str(root / "input.csv"), "--output-dir", str(output),
                                  "--accepted-only", "--min-lipids", "2", "--no-show"])
            prefix = self.script.OUTPUT_BASENAME
            scored = export.call_args.args[0]
            np.testing.assert_allclose(scored["Class Strength"], 2.0)
            self.assertFalse((output / f"{prefix}_Recomputed_Scores.csv").exists())
            points = self.script.prepare_strength(scored, all_ids=False)
            self.assertTrue(points[self.script.NORMALIZED].eq(1).all())
            files = [p for p in output.rglob("*") if p.is_file()]
            self.assertEqual(len(files), 4)
            self.assertEqual({p.suffix for p in files}, {".png", ".pdf"})

    def test_default_workbook_regenerates_scores_without_intermediate_file(self):
        self.assertEqual(self.script.INPUT_FILE.name, "Combined_Fragment_Assignment_Colored.xlsx")
        source = pd.concat([self.network_source().assign(Dataset=d) for d in self.script.datasets],
                           ignore_index=True).drop(columns="Class Strength")
        with TemporaryDirectory() as directory:
            workbook = Path(directory) / self.script.INPUT_FILE.name
            source.to_excel(workbook, sheet_name="Sheet1", index=False)
            with patch.object(self.script, "INPUT_FILE", workbook), \
                    patch.object(self.script, "export_population") as export:
                self.script.main(["--accepted-only", "--no-show"])
            np.testing.assert_allclose(export.call_args.args[0]["Class Strength"], 2.0)
            self.assertFalse(list(Path(directory).rglob("*Scores.csv")))
            self.assertFalse((Path(directory) / "Figure6_Class_Scores.csv").exists())


if __name__ == "__main__":
    unittest.main()

"""Standalone equal-IS processing of the standard Skyline transition export."""

import ast
import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location(
    "skyline_equal_is", ROOT / "Scripts/Annika/skyline_normalization_equal_is.py")
script = importlib.util.module_from_spec(spec)
spec.loader.exec_module(script)


class TestSkylineEqualIS(unittest.TestCase):
    groups = {"c1": "Control", "c2": "Control", "s1": "Sample", "s2": "Sample"}
    standards = {"PE": "IS PE standard"}

    def source(self):
        rows = []
        for replicate, factor, spike in (("c1", 1, 1), ("c2", 1, 2),
                                          ("s1", 2, 1), ("s2", 2, 2)):
            for adduct, denominator, ratio in (("[M-H]", 10, 2), ("[M+H]", 100, 4)):
                for molecule, area in (("PE analyte", denominator * ratio * factor * spike),
                                       ("IS PE standard", denominator * spike)):
                    for fragment, fraction in (("precursor", 10), ("T1", .25), ("H1", .75)):
                        rows.append({"Molecule": molecule, "Molecule List Name": "PE",
                                     "Precursor Adduct": adduct, "Replicate Name": replicate,
                                     "Fragment Ion": fragment, "Area": area * fraction,
                                     "Retention Time": 5., "Precursor Mz": 700.})
        return pd.DataFrame(rows).sample(frac=1, random_state=12)

    def test_products_sums_match_existing_plotter(self):
        # Load only the pure extraction functions, avoiding the plotter's GUI
        # backend, optional plotting dependencies, and module-level IS table.
        path = ROOT / "public/UniDec/PublicScripts/Lipids/SkylineResultsPlotter.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        names = {"clean_classes", "filter_transition_mode", "sum_transitions"}
        subset = ast.Module(body=[node for node in tree.body
                                 if isinstance(node, ast.FunctionDef) and node.name in names],
                            type_ignores=[])
        namespace = {"pd": pd, "np": np}
        exec(compile(subset, str(path), "exec"), namespace)
        source = self.source()
        expected, replicates = namespace["sum_transitions"](
            source, mode="Products", drop_IS=False, normalize_IS=False, normalize_TIC=False)
        actual = script.prepare_areas(source).pivot(
            index=["Molecule", "Adduct"], columns="Replicate Name", values="Area")
        expected = expected.set_index(["Molecule", "Adduct"])[replicates]
        pd.testing.assert_frame_equal(actual.sort_index().sort_index(axis=1),
                                      expected.sort_index().sort_index(axis=1),
                                      check_names=False)

    def test_adduct_ratios_are_averaged_after_normalization(self):
        normalized = script.normalize_areas(self.source(), self.groups, self.standards)
        values = script.molecule_values(normalized)
        # Correct result is mean(2, 4) = 3, rather than sum(20, 400)/sum(10, 100).
        np.testing.assert_allclose(values.loc[("PE analyte", "PE"), ["c1", "c2"]], 3)
        np.testing.assert_allclose(values.loc[("PE analyte", "PE"), ["s1", "s2"]], 6)
        _, folds = script.summarize(values, self.groups, "Control")
        self.assertEqual(folds.loc[0, "FoldChange"], 2)
        self.assertEqual(folds.loc[0, "Log2FC"], 1)

    def test_transition_modes(self):
        expected = {"Products": 20, "Precursors": 200, "Tails": 5, "Heads": 15, "All": 220}
        for mode, area in expected.items():
            with self.subTest(mode=mode):
                prepared = script.prepare_areas(self.source(), mode)
                row = prepared[(prepared.Molecule == "PE analyte") &
                               (prepared.Adduct == "[M-H]") & (prepared["Replicate Name"] == "c1")]
                self.assertEqual(row.Area.iloc[0], area)

    def test_missing_standard_and_measurements(self):
        source = self.source()
        mask = (source.Molecule == "IS PE standard") & (source["Fragment Ion"] != "precursor")
        source.loc[mask, "Area"] = 0
        with self.assertRaisesRegex(ValueError, "zero IS"):
            script.normalize_areas(source, self.groups, self.standards)
        source = self.source()
        source.loc[source.Molecule == "PE analyte", "Area"] = np.nan
        normalized = script.normalize_areas(source, self.groups, self.standards)
        self.assertTrue(normalized.Normalized_Area.isna().all())

    def test_labelled_adduct_alias_and_class_from_export(self):
        source = self.source()
        source.loc[(source.Molecule == "IS PE standard") &
                   (source["Precursor Adduct"] == "[M+H]"), "Precursor Adduct"] = "[M7H2+H]"
        # Molecule List Name is authoritative even for names without a class prefix.
        source["Molecule"] = source.Molecule.replace({"PE analyte": "analyte"})
        normalized = script.normalize_areas(source, self.groups, self.standards)
        self.assertEqual(set(normalized.Class), {"PE"})
        self.assertEqual(set(normalized.Adduct), {"[M-H]", "[M+H]"})

    def test_pc_preferred_adduct(self):
        source = self.source().replace({"Molecule List Name": {"PE": "PC"},
                                        "Precursor Adduct": {"[M-H]": "[M+CH3COO]"}})
        normalized = script.normalize_areas(source, self.groups, {"PC": "IS PE standard"})
        self.assertEqual(set(normalized.Adduct), {"[M+CH3COO]"})

    def test_direct_export_to_workbook_and_total_normalization(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "Molecule Transition Results.csv"
            self.source().to_csv(path, index=False, sep="\t")
            script.run(path, root / "out", self.groups, "Control", self.standards)
            folds = pd.read_csv(root / "out/lipid_fold_changes.csv")
            self.assertEqual(folds.loc[0, "FoldChange"], 2)
            classes = pd.read_csv(root / "out/class_fold_changes.csv")
            self.assertEqual(classes.loc[0, "FoldChange"], 2)
            with pd.ExcelFile(root / "out/normalized_abundances.xlsx") as book:
                self.assertEqual(len(book.sheet_names), 6)
                settings = pd.read_excel(book, sheet_name="Settings").set_index("Setting")
                self.assertEqual(settings.loc["Transition mode", "Value"], "Products")
            norm = script.normalize_areas(self.source(), self.groups, self.standards, True)
            np.testing.assert_allclose(script.molecule_values(norm).sum(), 1)

    def test_molecule_level_input_still_supported(self):
        source = self.source()
        source = source[source["Precursor Adduct"] == "[M-H]"]
        source = source[source["Fragment Ion"] == "T1"].drop(
            columns=["Fragment Ion", "Precursor Adduct", "Molecule List Name"])
        norm = script.normalize_areas(source, self.groups, self.standards)
        self.assertEqual(len(norm), 4)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            script.normalize_areas(pd.concat([source, source.iloc[[0]]]), self.groups, self.standards)

    def test_volcano_statistics_and_missing_values(self):
        from scipy import stats

        index = pd.MultiIndex.from_tuples(
            [("PE A", "PE"), ("PG B", "PG"), ("CL zero", "CL"), ("PE missing", "PE")],
            names=["Molecule", "Class"])
        values = pd.DataFrame({"s1": [2., 1., 2., 1.], "s2": [4., 2., 2., np.nan],
                               "s3": [6., 3., 2., np.nan], "c1": [1., 3., 0., 1.],
                               "c2": [2., 4., 0., 2.], "c3": [3., 5., 0., 3.]}, index=index)
        samples, controls = ["s1", "s2", "s3"], ["c1", "c2", "c3"]
        result = script.lipid_volcano_statistics(values, samples, controls).set_index("Molecule")
        self.assertEqual(result.loc["PE A", "Log2FC"], 1.)
        for name in ("PE A", "PG B"):
            row = values.xs(name, level="Molecule").iloc[0]
            p = stats.ttest_ind(row[samples], row[controls], equal_var=True).pvalue
            self.assertAlmostEqual(result.loc[name, "PValue"], p)
        self.assertTrue(np.isnan(result.loc["CL zero", "Log2FC"]))
        self.assertEqual(result.loc["CL zero", "PValue"], 0.)
        self.assertTrue(np.isnan(result.loc["PE missing", "PValue"]))
        ordered = result.dropna(subset=["PValue"]).sort_values("PValue")
        p = ordered.PValue.to_numpy()
        q = np.minimum.accumulate((p * len(p) / np.arange(1, len(p)+1))[::-1])[::-1].clip(0, 1)
        np.testing.assert_allclose(ordered.BH_QValue, q)

    def test_independent_reference_statistics(self):
        from scipy import stats

        samples, controls = ["s1", "s2", "s3"], ["c1", "c2", "c3"]
        values = pd.DataFrame({"c1": [1., 2.], "c2": [2., 3.], "c3": [3., 4.],
                               "s1": [2., 2.2], "s2": [4., 3.1], "s3": [6., 4.4]},
                              index=pd.Index(["CL", "PE"], name="Class"))
        result = script.independent_class_statistics(values, samples, controls).set_index("Class")
        for cls in values.index:
            sample = values.loc[cls, samples].to_numpy()
            control = values.loc[cls, controls].to_numpy()
            expected_n = stats.ttest_ind(sample, control, equal_var=True)
            self.assertAlmostEqual(result.loc[cls, "FigN_PValue"], expected_n.pvalue)
            self.assertAlmostEqual(result.loc[cls, "Sample_CI95"],
                                   stats.t.ppf(.975, 4) * stats.sem(sample))
            fc = sample.mean() / control.mean()
            sem = fc * np.sqrt((stats.sem(sample) / sample.mean())**2 +
                               (stats.sem(control) / control.mean())**2)
            self.assertAlmostEqual(result.loc[cls, "FigM_Sample_mean"], sample.mean())
            self.assertAlmostEqual(result.loc[cls, "FigM_Control_mean"], control.mean())
            self.assertAlmostEqual(result.loc[cls, "FoldChange"], fc)
            self.assertAlmostEqual(result.loc[cls, "FC_CI95"], stats.norm.ppf(.975) * sem)
            self.assertAlmostEqual(result.loc[cls, "FigM_PValue"],
                                   2 * stats.norm.sf(abs((fc - 1) / sem)))
        for figure in ("FigM", "FigN"):
            ordered = result.sort_values(f"{figure}_PValue")
            p = ordered[f"{figure}_PValue"].to_numpy()
            q = np.minimum.accumulate((p * 2 / np.arange(1, 3))[::-1])[::-1].clip(0, 1)
            np.testing.assert_allclose(ordered[f"{figure}_BH_QValue"], q)
        # Independent groups do not depend on correspondence/order of replicates.
        reversed_result = script.independent_class_statistics(values, samples[::-1], controls)
        np.testing.assert_allclose(reversed_result.FigN_PValue, result.FigN_PValue)
        # Unequal independent group sizes are supported.
        unequal = script.independent_class_statistics(values, samples[:2], controls)
        self.assertTrue((unequal.Sample_N == 2).all())
        self.assertTrue((unequal.Control_N == 3).all())
        with self.assertRaisesRegex(ValueError, "different replicates"):
            script.independent_class_statistics(values, samples, samples)

    def test_class_fold_changes_preserve_global_abundance_increase(self):
        # Doubling every class must remain visible in abundance mode. Only
        # explicitly total-normalized inputs should remove the global increase.
        samples, controls = ["s1", "s2", "s3"], ["c1", "c2", "c3"]
        values = pd.DataFrame({"c1": [1., 2.], "c2": [2., 4.], "c3": [3., 6.],
                               "s1": [2., 4.], "s2": [4., 8.], "s3": [6., 12.]},
                              index=pd.Index(["CL", "PE"], name="Class"))
        abundance = script.independent_class_statistics(values, samples, controls)
        np.testing.assert_allclose(abundance.FoldChange, 2.)
        composition = script.independent_class_statistics(values / values.sum(), samples, controls)
        np.testing.assert_allclose(composition.FoldChange, 1.)


if __name__ == "__main__":
    unittest.main()

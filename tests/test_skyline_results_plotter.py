"""Manual research-script checks; no GUI backend or local data required."""

import ast
from copy import deepcopy
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests


def load_quantification_functions():
    path = Path(__file__).resolve().parents[1] / "PublicScripts/Lipids/SkylineResultsPlotter.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = {"filter_transition_mode", "clean_classes", "sum_transitions", "normalize_tic",
             "normalize_is", "check_standard_variation", "dereplicate_adducts", "get_is",
             "calc_is", "get_area_matched_is", "get_is_conc_from_isdf",
             "compare_pipeline", "class_calcs", "stats_calcs", "safedivide", "split_names",
             "detection_table", "comparison_stats", "filter_normdf", "make_volcano_plot",
             "heatmap_plot", "make_class_plots", "select_internal_standards"}
    constants = {"preferred_quant_adducts", "area_match_is_classes", "default_is_mapper",
                 "IS_AREA_VARIATION_WARNING_FOLD", "PC_HEADGROUP_MZ",
                 "PC_HEADGROUP_TOLERANCE", "QUANTIFICATION_NOTE_COLUMNS"}
    nodes = [node for node in tree.body if
             isinstance(node, ast.FunctionDef) and node.name in names or
             isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
             and node.targets[0].id in constants]
    namespace = {"pd": pd, "np": np, "default_isdf": pd.DataFrame(), "os": os,
                 "deepcopy": deepcopy, "stats": stats, "multipletests": multipletests}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


class TestSkylineResultsPlotter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.quant = load_quantification_functions()

    def test_reviewed_lpi_standard_is_shared_between_groups(self):
        selected = "IS LPI(+[2]H5) 17:0"
        rejected = "IS LPI(+[2]H5) 15:0"
        reps = ["a1", "a2", "b1", "b2"]
        rows = []
        for i, rep in enumerate(reps):
            for name, area, rt in [("LPI analyte", 20. + 10. * i, 1.),
                                    (selected, 10., 5.)]:
                rows.append({"Molecule": name, "Molecule List Name": "LPI", "Replicate Name": rep,
                             "Precursor Adduct": "[M-H]", "Precursor Mz": 700., "Product Mz": 100.,
                             "Fragment Ion": "product", "Area": area, "Retention Time": rt})
            if rep.startswith("a"):
                rows.append(dict(rows[-1], Molecule=rejected, Area=100., **{"Retention Time": 1.}))
        source = pd.DataFrame(rows)
        audit = []
        filtered = self.quant["select_internal_standards"](source, {"LPI": selected}, reps, audit)
        self.assertNotIn(rejected, set(filtered["Molecule"]))
        self.assertEqual(audit[0]["Rows removed"], 2)
        standards = pd.DataFrame({"Class": ["LPI"], "Lipid": [selected.replace("IS ", "")], "Concentration (uM)": [1.]})
        result, _ = self.quant["sum_transitions"](
            filtered, normalize_TIC=False, customisdf=standards, is_replicate_groups=[reps[:2], reps[2:]])
        row = result.set_index("Molecule").loc["LPI analyte"]
        np.testing.assert_allclose(row[reps].to_numpy(dtype=float), [2., 3., 4., 5.])
        missing = source.loc[~(source.Molecule.eq(selected) & source["Replicate Name"].eq("b2"))]
        with self.assertRaisesRegex(ValueError, "b2"):
            self.quant["select_internal_standards"](missing, {"LPI": selected}, reps)

    def test_lpa_concentration_matches_catalogue_with_leading_whitespace(self):
        row = pd.Series({"Molecule": "IS LPA 18:1(+[2]H9)", "Molecule List Name": "LPA"})
        catalogue = pd.DataFrame({"Class": ["LPA"], "Lipid": [" IS LPA 18:1(+[2]H9) "],
                                  "Concentration (uM)": [0.5]})
        self.assertEqual(self.quant["get_is_conc_from_isdf"](row, catalogue, "Concentration (uM)"), 0.5)

    def test_missing_entries_remain_unavailable_and_distinct_from_failures(self):
        reps = ["a1", "a2", "b1", "b2"]
        quantified = pd.DataFrame([
            ["drop", "PE", 1., 2., np.nan, np.nan],
            ["increase", "PG", np.nan, np.nan, 3., 4.],
            ["failure", "PE", 2., 3., np.nan, np.nan],
            ["shared", "PE", 1., 2., 4., 6.]],
            columns=["Molecule", "Molecule List Name"] + reps)
        rows = []
        for name, cls, group in [("drop", "PE", reps[:2]), ("increase", "PG", reps[2:]),
                                  ("failure", "PE", reps), ("shared", "PE", reps)]:
            rows.extend({"Molecule": name, "Molecule List Name": cls, "Replicate Name": r} for r in group)
        annotated = self.quant["detection_table"](quantified, pd.DataFrame(rows), reps)
        result, _ = self.quant["comparison_stats"](annotated, reps[:2], reps[2:])
        result = result.set_index("Molecule")
        self.assertEqual(result.loc["drop", "Detection b1"], "Not measured")
        self.assertTrue(np.isnan(result.loc["drop", "b1"]))
        self.assertEqual(result.loc["drop", "Detection change"], "Set1 detected; Set2 unavailable")
        self.assertTrue(np.isnan(result.loc["drop", "Log2FC"]))
        self.assertTrue(np.isnan(result.loc["drop", "p-value"]))
        self.assertEqual(result.loc["increase", "Detection a1"], "Not measured")
        self.assertTrue(np.isnan(result.loc["increase", "a1"]))
        self.assertEqual(result.loc["failure", "Detection b1"], "Quantification failed")
        self.assertTrue(np.isnan(result.loc["failure", "b1"]))
        self.assertEqual(result.loc["shared", "Comparison status"], "Tested")
        selected = self.quant["filter_normdf"](result.reset_index(), fold_range=100.)
        self.assertTrue({"drop", "increase", "failure"}.issubset(set(selected["Molecule"])))

    def test_all_nd_comparisons_and_exclusive_classes_are_retained(self):
        source = pd.DataFrame({"Molecule": ["drop", "increase"], "Molecule List Name": ["PE", "PG"],
                               "a1": [1., 0.], "a2": [2., 0.], "b1": [0., 4.], "b2": [0., 6.]})
        result, _ = self.quant["comparison_stats"](source, ["a1", "a2"], ["b1", "b2"])
        self.assertEqual(len(result), 2)
        self.assertTrue(result["p-value"].isna().all())
        classes, _ = self.quant["class_calcs"](source, ["a1", "a2"], ["b1", "b2"], retain_unavailable=True)
        self.assertEqual(len(classes), 2)
        self.assertEqual(classes.loc[classes["Molecule List Name"].eq("PE"), "Set1_Mean"].iloc[0], 1.)

    def test_zero_tests_share_fdr_family_with_positive_measurements(self):
        a, b = ["a1", "a2", "a3"], ["b1", "b2", "b3"]
        source = pd.DataFrame([
            ["drop", 10., 11., 12., 0., 0., 0.],
            ["low", 20., 22., 24., .02, .022, .024],
            ["partial", 40., 44., 48., 0., .04, .05]], columns=["Molecule"] + a + b)
        result, _ = self.quant["comparison_stats"](source, a, b, include_nondetections=True)
        raw = [stats.ttest_ind(row[a].to_numpy(dtype=float), row[b].to_numpy(dtype=float)).pvalue
               for _, row in source.iterrows()]
        np.testing.assert_allclose(result["Raw p-value"], raw)
        np.testing.assert_allclose(result["p-value"], multipletests(raw, method="fdr_bh")[1])
        self.assertTrue(np.isneginf(result.loc[0, "Log2FC"]))
        self.assertTrue(result["Comparison status"].eq("Tested").all())
        self.assertEqual(result.loc[0, "Test basis"], "Includes zero non-detections")
        self.assertAlmostEqual(result.loc[1, "Log2FC"], np.log2(.001))
        self.assertTrue(np.isnan(result.loc[0, "CI95_FC"]))
        paired, _ = self.quant["comparison_stats"](source, a, b, paired=True, include_nondetections=True)
        self.assertTrue(np.isneginf(paired.loc[0, "Log2FC"]))
        self.assertTrue(np.isfinite(paired.loc[0, "p-value"]))

    def test_missing_set2_mode_retains_drops_low_values_and_quantification_failures(self):
        a, b = ["a1", "a2", "a3"], ["b1", "b2", "b3"]
        with TemporaryDirectory() as temporary:
            first, second = Path(temporary) / "control.csv", Path(temporary) / "analyte.csv"
            for path, reps, entries in [
                    (first, a, {"drop": [10., 11., 12.], "low": [20., 22., 24.],
                                "failure": [30., 32., 34.], "control_failure": [np.nan] * 3}),
                    (second, b, {"low": [.02, .022, .024], "failure": [np.nan] * 3,
                                 "set2_only": [5., 6., 7.]})]:
                rows = []
                for name, areas in entries.items():
                    for rep, area in zip(reps, areas):
                        rows.append({"Molecule": name, "Molecule List Name": "PE",
                                     "Replicate Name": rep, "Precursor Adduct": "[M-H]",
                                     "Precursor Mz": 700., "Product Mz": 100., "Fragment Ion": "product",
                                     "Area": area, "Retention Time": 5.})
                pd.DataFrame(rows).to_csv(path, index=False)
            result, _ = self.quant["compare_pipeline"](
                first, a, b, file2=second,
                normalize_IS=False, normalize_TIC=False, paired=False, fold_range=1.,
                plot_results=False, write_output=True)
            result = result.set_index("Molecule")
            self.assertEqual(set(result.index), {"drop", "low", "failure", "set2_only", "control_failure"})
            self.assertEqual(result.loc["drop", "Set2_AssignedZeros"], 3)
            self.assertTrue(result.loc["drop", "Major difference"])
            self.assertTrue(result.loc["drop", "Significant difference"])
            self.assertTrue(np.isneginf(result.loc["drop", "Log2FC"]))
            self.assertEqual(result.loc["low", "file2::b1"], .02)
            self.assertEqual(result.loc["low", "Set2_AssignedZeros"], 0)
            self.assertEqual(result.loc["failure", "Detection file2::b1"], "Quantification failed")
            self.assertTrue(np.isnan(result.loc["failure", "file2::b1"]))
            self.assertTrue(np.isnan(result.loc["set2_only", "file1::a1"]))
            self.assertEqual(result.loc["control_failure", "file2::b1"], 0.)
            self.assertTrue(np.isnan(result.loc["control_failure", "file1::a1"]))
            changes = pd.read_csv(Path(temporary) / "control_detection_changes.csv")
            self.assertIn("drop", set(changes.Molecule))
            with pd.ExcelFile(Path(temporary) / "control_stats.xlsx") as workbook:
                exported = pd.read_excel(workbook, "Molecule Stats").set_index("Molecule")
                self.assertEqual(exported.loc["drop", "Test basis"], "Includes zero non-detections")

    def test_volcano_places_nondetections_at_boundary_without_changing_fold_changes(self):
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
        source = pd.DataFrame({"Molecule": ["drop", "increase", "low"],
                               "Molecule List Name": ["PE"] * 3,
                               "Log2FC": [-np.inf, np.inf, -20.], "Log10p": [4., 3., 2.]})
        with patch.dict(self.quant, {"plt": plt, "get_color_for_class": lambda cls: "blue"}):
            fig, axes = plt.subplots(1, 2)
            try:
                self.quant["make_volcano_plot"](source, ax=axes[0], show_detection_summary=False,
                                                show_group_label=False, xlims=[-6., 6.])
                offsets = axes[0].collections[-1].get_offsets()
                np.testing.assert_allclose(offsets, [[-5.85, 4.], [5.85, 3.], [-5.85, 2.]])
                self.assertEqual(len(axes[0].get_legend().texts), 2)
                self.quant["heatmap_plot"](source, ["Log2FC"], ax=axes[1], normtype="zerocentered")
                self.assertTrue(np.isneginf(source.loc[0, "Log2FC"]))
                self.assertTrue(np.isposinf(source.loc[1, "Log2FC"]))
                self.assertEqual(source.loc[2, "Log2FC"], -20.)
                fig.canvas.draw()
            finally:
                plt.close(fig)

    def test_pipeline_exports_missing_group_lipids_without_fabricated_statistics(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / 'b2.csv'
            rows = []
            groups = [["a1", "a2"], ["b1", "b2"]]
            for group, name, cls in zip(groups, ["drop", "increase"], ["PE", "PG"]):
                for molecule, lipid_class in [(name, cls), ("shared", "PE")]:
                    for i, r in enumerate(group):
                        rows.append({"Molecule": molecule, "Molecule List Name": lipid_class,
                                     "Replicate Name": r, "Precursor Adduct": "[M-H]",
                                     "Precursor Mz": 700., "Product Mz": 100., "Fragment Ion": "product",
                                     "Area": float(i + 1), "Retention Time": 5.})
            pd.DataFrame(rows).to_csv(path, index=False)
            result, classes = self.quant["compare_pipeline"](
                path, groups[0], groups[1], normalize_IS=False, normalize_TIC=False, paired=False,
                missing_set2_as_zero=False, plot_results=False, write_output=True)
            self.assertEqual(set(result["Molecule"]), {"drop", "increase", "shared"})
            drop = result.set_index("Molecule").loc["drop"]
            self.assertTrue(np.isnan(drop["b1"]))
            self.assertTrue(np.isnan(drop["p-value"]))
            changes = pd.read_csv(Path(temporary) / 'b2_detection_changes.csv')
            self.assertEqual(set(changes["Molecule"]), {"drop", "increase"})
            with pd.ExcelFile(Path(temporary) / 'b2_stats.xlsx') as workbook:
                self.assertIn("Detection Changes", workbook.sheet_names)
                self.assertEqual(len(pd.read_excel(workbook, "Molecule Stats")), 3)

    def test_plots_keep_original_fold_changes_and_mask_unavailable_cells(self):
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
        frame = pd.DataFrame({"Molecule": ["tested", "drop"], "Molecule List Name": ["PE", "PE"],
                              "Log2FC": [20., np.nan], "Log10p": [3., np.nan],
                              "Detection change": ["Detected in both", "Set1 only; Set2 ND"],
                              "a": [1., 2.], "b": [np.nan, 0.]})
        with patch.dict(self.quant, {"plt": plt, "get_color_for_class": lambda cls: "blue"}):
            fig, axes = plt.subplots(1, 3)
            try:
                self.quant["make_volcano_plot"](frame, ax=axes[0])
                self.assertEqual(frame.loc[0, "Log2FC"], 20.)
                self.assertTrue(any("Set1 only" in label.get_text() for label in axes[0].texts))
                self.quant["heatmap_plot"](frame, ["a", "b"], ax=axes[1])
                self.assertTrue(axes[1].images[0].get_array().mask[0, 1])
                self.quant["make_class_plots"](pd.DataFrame({"FC": [np.nan], "CI95_FC": [np.nan]}), ax=axes[2])
                fig.canvas.draw()
            finally:
                plt.close(fig)

    def test_product_filter_excludes_renamed_precursors_and_preserves_fragments(self):
        source = pd.DataFrame({
            "Fragment Ion": ["precursor", "Ion [700/700]", "HG", "same mass different charge", "missing"],
            "Precursor Mz": [700.] * 5, "Product Mz": [700., 700.0001, 184., 700., np.nan],
            "Precursor Charge": [1] * 5, "Product Charge": [1, 1, 1, 2, 1]})
        products = self.quant["filter_transition_mode"](source, "Products")
        self.assertEqual(products.index.tolist(), [2, 3, 4])
        self.assertEqual(len(self.quant["filter_transition_mode"](source, "All")), 5)
        legacy = source[["Fragment Ion"]]
        self.assertEqual(len(self.quant["filter_transition_mode"](legacy, "Products")), 4)
        shifted = source.iloc[[1]].copy()
        shifted["Product Mz"] = 710.
        self.assertTrue(self.quant["filter_transition_mode"](shifted, "Products").empty)

    def test_reimported_project_product_sums_match_original(self):
        rows = []
        for replicate in ["original", "reimported"]:
            for name, product_mz, area in [("precursor", 700., 1000.), ("HG", 184., 20.)]:
                rows.append({"Molecule": "PI analyte", "Molecule List Name": "PI",
                             "Precursor Adduct": "[M-H]", "Replicate Name": replicate,
                             "Precursor Mz": 700., "Product Mz": product_mz,
                             "Fragment Ion": name, "Area": area, "Retention Time": 5.})
        rows.append(dict(rows[2], **{"Fragment Ion": "Ion [700/700]"}))
        result, _ = self.quant["sum_transitions"](
            pd.DataFrame(rows), normalize_IS=False, drop_IS=False, normalize_TIC=False)
        self.assertEqual(result["original"].iloc[0], 20.)
        self.assertEqual(result["reimported"].iloc[0], 20.)

    def test_pc_products_use_matching_headgroup_only_in_all_projects(self):
        rows = []
        for replicate in ["original", "reimported"]:
            for cls in ["PC", "EtherPC"]:
                for molecule, adduct, head, area in [
                        (f"{cls} analyte", "[M+H]", 184.0738, 20.),
                        (f"IS {cls} d5", "[M5H2+H]", 184.0733, 10.)]:
                    # Glycerol-d5 does not shift the phosphocholine headgroup.
                    # Reject the old, incorrectly shifted target even if intense.
                    for mass, value in [(head, area), (head + .02, 50.), (189.1047, 50000.), (400., 900.)]:
                        rows.append({"Molecule": molecule, "Molecule List Name": cls,
                                     "Precursor Adduct": adduct, "Replicate Name": replicate,
                                     "Precursor Mz": 800., "Product Mz": mass,
                                     "Fragment Ion": "HG" if replicate == "original" else "Ion [...]",
                                     "Area": value, "Retention Time": 5.})
        source = pd.DataFrame(rows)
        products = self.quant["filter_transition_mode"](source, "Products")
        self.assertEqual(len(products), 8)
        areas, _ = self.quant["sum_transitions"](
            source, normalize_IS=False, drop_IS=False, normalize_TIC=False)
        for _, row in areas.iterrows():
            expected = 10. if row["Molecule"].startswith("IS ") else 20.
            self.assertEqual(row["original"], expected)
            self.assertEqual(row["reimported"], expected)
        # Acetate and other classes continue to use all real product fragments.
        acetate = source.copy()
        acetate["Precursor Adduct"] = "[M+CH3COO]"
        self.assertEqual(len(self.quant["filter_transition_mode"](acetate, "Products")), len(acetate))
        pe = source.copy()
        pe["Molecule List Name"] = "PE"
        self.assertEqual(len(self.quant["filter_transition_mode"](pe, "Products")), len(pe))
        # Already-filtered inputs use All inside compare_pipeline; do not re-sum excluded peaks.
        all_areas, _ = self.quant["sum_transitions"](
            products, mode="All", normalize_IS=False, drop_IS=False, normalize_TIC=False)
        pd.testing.assert_frame_equal(areas, all_areas)
        # The full extraction/IS-normalization path must use headgroup/headgroup,
        # rather than dividing a full analyte fragment sum by one IS fragment.
        pc = source.loc[source["Molecule List Name"].eq("PC")]
        normalized, _ = self.quant["sum_transitions"](
            pc, normalize_IS=True, drop_IS=True, normalize_TIC=False, customisdf=None)
        np.testing.assert_allclose(normalized[["original", "reimported"]].to_numpy(), [[2., 2.]])

    def test_incorrectly_shifted_pc_standard_is_not_reinterpreted(self):
        source = pd.DataFrame({
            "Molecule": ["PC analyte", "IS PC glycerol-d5"],
            "Molecule List Name": ["PC", "PC"],
            "Precursor Adduct": ["[M+H]", "[M5H2+H]"],
            "Fragment Ion": ["HG(PC,184)", "HG(PC,184)"],
            "Product Mz": [184.0733, 189.1047]})
        products = self.quant["filter_transition_mode"](source, "Products")
        self.assertEqual(products["Molecule"].tolist(), ["PC analyte"])

    def test_pc_headgroup_requires_measured_product_mass(self):
        source = pd.DataFrame({"Molecule": ["PC analyte"], "Molecule List Name": ["PC"],
                               "Precursor Adduct": ["[M+H]"], "Fragment Ion": ["HG"]})
        with self.assertRaisesRegex(ValueError, "requires Product Mz"):
            self.quant["filter_transition_mode"](source, "Products")

    def test_pc_uses_matching_is_and_prefers_acetate_per_replicate(self):
        for lipid_class in ["PC", "EtherPC"]:
            with self.subTest(lipid_class=lipid_class):
                source = pd.DataFrame([
                    {"Molecule": "analyte", "Adduct": "[M+H]", "r1": 10., "r2": 14.},
                    {"Molecule": "analyte", "Adduct": "[M+CH3COO]", "r1": 40., "r2": 0.},
                    {"Molecule": "analyte", "Adduct": "[M+Na]", "r1": 999., "r2": 999.},
                    {"Molecule": "IS PC standard", "Adduct": "[M+H]", "r1": 2., "r2": 2.},
                    {"Molecule": "IS PC standard", "Adduct": "[M+CH3COO]", "r1": 10., "r2": 10.}])
                source["Molecule List Name"] = lipid_class
                source.loc[source["Molecule"].str.startswith("IS "), "Molecule List Name"] = "PC"
                source["RT"] = 5.
                result = self.quant["normalize_is"](source, ["r1", "r2"], "unused", isdf=None)
                row = result.loc[result["Molecule"].eq("analyte")].iloc[0]
                self.assertEqual(row["r1"], 4.)
                self.assertEqual(row["r2"], 7.)
                self.assertEqual(row["Quant Adduct r1"], "[M+CH3COO]")
                self.assertEqual(row["Quant Adduct r2"], "[M+H]")

    def test_missing_acetate_standard_uses_protonated_pc(self):
        source = pd.DataFrame([
            {"Molecule": "PC analyte", "Adduct": "[M+H]", "r": 20.},
            {"Molecule": "PC analyte", "Adduct": "[M+CH3COO]", "r": 40.},
            {"Molecule": "IS PC standard", "Adduct": "[M+H]", "r": 2.}])
        source["Molecule List Name"] = "PC"
        source["RT"] = 5.
        result = self.quant["normalize_is"](source, ["r"], "unused", isdf=None)
        row = result.loc[result["Molecule"].eq("PC analyte")].iloc[0]
        self.assertEqual(row["r"], 10.)
        self.assertEqual(row["Quant Adduct r"], "[M+H]")

    def test_standard_variation_warns_and_only_invalid_denominators_are_rejected(self):
        source = pd.DataFrame([
            ["IS LPC 15", "[M+H]", 100., 10., 10., 10.],
            ["IS LPC 17", "[M+H]", 10., 11., 100., 110.],
            ["IS LPC 19", "[M+H]", 10., 10., 10., 10.],
            ["IS LPC 15", "[M+CH3COO]", 10., 10., 10., 10.],
            ["LPC analyte", "[M+H]", 100., 1., 100., 1.]],
            columns=["Molecule", "Adduct", "a1", "a2", "b1", "b2"])
        result = self.quant["check_standard_variation"](source, [["a1", "a2"], ["b1", "b2"]])
        self.assertEqual(result.index.tolist(), [0, 1, 2, 3, 4])
        self.assertEqual(len(self.quant["check_standard_variation"](source, None)), 5)
        invalid = source.copy()
        invalid.loc[1, "a1"] = np.nan
        self.assertNotIn(1, self.quant["check_standard_variation"](invalid, [["a1", "a2"]]).index)

    def test_other_classes_still_average_valid_adducts(self):
        source = pd.DataFrame({"Molecule": ["PE analyte"] * 2, "Molecule List Name": ["PE"] * 2,
                               "Adduct": ["[M-H]", "[M+H]"], "r": [2., 4.]})
        self.assertEqual(self.quant["dereplicate_adducts"](source, ["r"])["r"].iloc[0], 3.)

    def test_failed_pc_quantification_does_not_become_zero(self):
        source = pd.DataFrame({"Molecule": ["PC analyte"] * 2, "Molecule List Name": ["PC"] * 2,
                               "Adduct": ["[M+CH3COO]", "[M+H]"], "r": [0., np.nan]})
        self.assertTrue(np.isnan(self.quant["dereplicate_adducts"](source, ["r"])["r"].iloc[0]))

    def test_variable_lpc_standard_is_retained_during_normalization(self):
        source = pd.DataFrame([
            ["LPC analyte", 2., 20., 20.],
            ["IS LPC 15", 2., 100., 10.],
            ["IS LPC 17", 3., 5., 5.],
            ["IS LPC 19", 5., 10., 10.]], columns=["Molecule", "RT", "r1", "r2"])
        source["Molecule List Name"] = "LPC"
        source["Adduct"] = "[M+H]"
        audit = []
        result = self.quant["normalize_is"](source, ["r1", "r2"], "unused", isdf=None,
                                             replicate_groups=[["r1", "r2"]], exclusion_audit=audit)
        self.assertIn("IS LPC 15", result["Molecule"].values)
        self.assertEqual(audit[0]["Molecule"], "IS LPC 15")
        self.assertIn("10.00-fold", audit[0]["Reason"])
        self.assertEqual(audit[0]["Rows removed"], 0)
        self.assertEqual(audit[0]["Stage"], "IS QC warning")
        row = result.loc[result["Molecule"].eq("LPC analyte")].iloc[0]
        np.testing.assert_allclose(row[["r1", "r2"]].to_numpy(dtype=float), [.2, 2.])

    def test_pc_standards_use_individual_concentrations_without_area_ratio_rejection(self):
        source = pd.DataFrame([
            ["PC low analyte", 20., 40.], ["PC high analyte", 1000., 2000.],
            ["PC very high analyte", 100000., 200000.],
            ["IS PC low standard", 10., 10.], ["IS PC high standard", 1000., 1000.]],
            columns=["Molecule", "a", "b"])
        source["Molecule List Name"] = "PC"
        source["Adduct"] = "[M+H]"
        source["RT"] = 5.
        catalogue = pd.DataFrame({"Class": ["PC", "PC"],
                                  "Lipid": ["PC low standard", "PC high standard"],
                                  "Concentration (uM)": [1., 10.]})
        result = self.quant["normalize_is"](source, ["a", "b"], "Concentration (uM)",
                                             isdf=catalogue, replicate_groups=[["a", "b"]])
        result = result.set_index("Molecule")
        for name, expected in [("PC low analyte", [2., 4.]), ("PC high analyte", [10., 20.]),
                                ("PC very high analyte", [1000., 2000.])]:
            np.testing.assert_allclose(result.loc[name, ["a", "b"]].to_numpy(dtype=float), expected)

    def test_pipeline_keeps_formerly_excluded_lipids_in_totals_and_exports(self):
        rows = []
        replicates = ["a1", "a2", "b1", "b2"]
        for name, cls, areas in [("PE ordinary", "PE", [1., 2., 3., 4.]),
                                  ("PG ordinary", "PG", [4., 3., 2., 1.]),
                                  ("PE 36:0|PE 18:0_18:0", "PE", [999.] * 4)]:
            for replicate, area in zip(replicates, areas):
                rows.append({"Molecule": name, "Molecule List Name": cls, "Replicate Name": replicate,
                             "Precursor Adduct": "[M-H]", "Precursor Mz": 700., "Product Mz": 100.,
                             "Fragment Ion": "product", "Area": area, "Retention Time": 5.})
        original_cwd = Path.cwd()
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "results.csv"
            pd.DataFrame(rows).to_csv(path, index=False)
            try:
                result, _ = self.quant["compare_pipeline"](
                    str(path), replicates[:2], replicates[2:], normalize_IS=False,
                    normalize_TIC=True, paired=False, plot_results=False, write_output=True)
                self.assertEqual(set(result["Molecule"]), {"PE ordinary", "PG ordinary", "PE 36:0|PE 18:0_18:0"})
                pe = result.loc[result["Molecule"].eq("PE ordinary")].iloc[0]
                np.testing.assert_allclose(pe[replicates].to_numpy(dtype=float), np.array([1., 2., 3., 4.]) / 1004.)
                notes = pd.read_csv(path.with_name("results_quantification_notes.csv"))
                self.assertTrue(notes["Rows removed"].eq(0).all())
                with pd.ExcelFile(path.with_name("results_stats.xlsx")) as workbook:
                    self.assertIn("Quantification Notes", workbook.sheet_names)
                    exported = pd.read_excel(workbook, "Molecule Stats")
                    self.assertIn("PE 36:0|PE 18:0_18:0", set(exported.Molecule))
            finally:
                os.chdir(original_cwd)

    def test_two_csvs_select_each_group_only_from_its_assigned_file(self):
        with TemporaryDirectory() as temporary:
            first = Path(temporary) / "b2.csv"
            second = Path(temporary) / "edta.csv"
            for path, selected, other, areas in [
                    (first, ["a1", "a2"], ["b1", "b2"], [10., 20.]),
                    (second, ["b1", "b2"], ["a1", "a2"], [40., 60.])]:
                rows = []
                for rep, area in zip(selected + other, areas + [999., 999.]):
                    for name, cls, value in [("PE ordinary", "PE", area),
                                              ("PG ordinary", "PG", area + 5.)]:
                        rows.append({"Molecule": name, "Molecule List Name": cls,
                                     "Replicate Name": rep, "Precursor Adduct": "[M-H]",
                                     "Precursor Mz": 700., "Product Mz": 100., "Fragment Ion": "product",
                                     "Area": value, "Retention Time": 5.})
                pd.DataFrame(rows).to_csv(path, index=False)
            result, _ = self.quant["compare_pipeline"](
                first, {"file1": ["a1", "a2"]}, {"file2": ["b1", "b2"]}, file2=second,
                normalize_IS=False, normalize_TIC=False, paired=False, plot_results=False, write_output=False)
            row = result.set_index("Molecule").loc["PE ordinary"]
            np.testing.assert_allclose(row[["file1::a1", "file1::a2", "file2::b1", "file2::b2"]].to_numpy(dtype=float), [10., 20., 40., 60.])
            with self.assertRaisesRegex(ValueError, "not found in CSV: missing_sample"):
                self.quant["compare_pipeline"](first, {"file1": ["a1", "a2"]}, {"file2": ["missing_sample"]},
                                               file2=second, plot_results=False, write_output=False)

    def test_groups_can_select_either_or_both_files_and_keep_custom_labels(self):
        with TemporaryDirectory() as temporary:
            file1, file2 = Path(temporary) / "first.csv", Path(temporary) / "second.csv"
            for path, areas in [(file1, [10., 30., 50., 70.]), (file2, [20., 40., 60., 80.])]:
                rows = [{"Molecule": "PE ordinary", "Molecule List Name": "PE",
                         "Replicate Name": rep, "Precursor Adduct": "[M-H]",
                         "Precursor Mz": 700., "Product Mz": 100., "Fragment Ion": "product",
                         "Area": area, "Retention Time": 5.}
                        for rep, area in zip(["r1", "r2", "r3", "r4"], areas)]
                pd.DataFrame(rows).to_csv(path, index=False)
            choices = [
                ({"file1": ["r1"], "file2": ["r1"]},
                 {"file1": ["r2"], "file2": ["r2"]}, 15., 35.),
                ({"file2": ["r1", "r2"]}, {"file1": ["r1", "r2"]}, 30., 20.),
                ({"file1": ["r1", "r2"]}, {"file1": ["r3", "r4"]}, 20., 60.),
                ({"file2": ["r1", "r2"]}, {"file2": ["r3", "r4"]}, 30., 70.)]
            for set1, set2, control_mean, analyte_mean in choices:
                with self.subTest(set1=set1, set2=set2):
                    result, _ = self.quant["compare_pipeline"](
                        file1, set1, set2, file2=file2, normalize_IS=False, normalize_TIC=False,
                        paired=False, plot_results=False, write_output=True,
                        set1_name="My control", set2_name="My analyte")
                    row = result.iloc[0]
                    self.assertEqual(row["Set1_Mean"], control_mean)
                    self.assertEqual(row["Set2_Mean"], analyte_mean)
                    self.assertEqual(row["Set1 name"], "My control")
                    self.assertEqual(row["Set2 name"], "My analyte")
                    with pd.ExcelFile(Path(temporary) / "first_stats.xlsx") as book:
                        selection = pd.read_excel(book, "Sample Selection")
                        self.assertEqual(len(selection), 4)
                        self.assertEqual(set(selection["Group label"]), {"My control", "My analyte"})
            with self.assertRaisesRegex(ValueError, "occurs in both CSVs"):
                self.quant["compare_pipeline"](file1, ["r1"], ["r2"], file2=file2,
                                               plot_results=False, write_output=False)
            with self.assertRaisesRegex(ValueError, "selected more than once"):
                self.quant["compare_pipeline"](file1, {"file1": ["r1"]}, {"file1": ["r1"]},
                                               file2=file2, plot_results=False, write_output=False)

    def test_single_csv_selects_comparison_samples_before_is_normalization(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "results.csv"
            reps = ["a1", "a2", "b1", "b2"]
            rows = []
            for replicate, area in zip(reps + ["unused"], [10., 20., 40., 60., 999.]):
                for name, cls, value in [("PE ordinary", "PE", area),
                                          ("IS PE standard", "PE", 10. if replicate != "unused" else 0.),
                                          ("PG ordinary", "PG", 20. + area / 10),
                                          ("IS PG standard", "PG", 10. if replicate != "unused" else 0.)]:
                    rows.append({"Molecule": name, "Molecule List Name": cls, "Replicate Name": replicate,
                                 "Precursor Adduct": "[M-H]", "Precursor Mz": 700., "Product Mz": 100.,
                                 "Fragment Ion": "product", "Area": value, "Retention Time": 5.})
            pd.DataFrame(rows).to_csv(path, index=False)
            standards = pd.DataFrame({"Lipid": ["PE standard", "PG standard"], "Class": ["PE", "PG"],
                                      "Concentration (uM)": [1., 1.]})
            original_cwd = Path.cwd()
            with patch.dict(self.quant, {"default_isdf": standards}):
                result, _ = self.quant["compare_pipeline"](
                    path, reps[:2], reps[2:], normalize_IS=True, normalize_TIC=False,
                    paired=False, plot_results=False, write_output=False)
            self.assertEqual(Path.cwd(), original_cwd)
            row = result.loc[result.Molecule.eq("PE ordinary")].iloc[0]
            np.testing.assert_allclose(row[reps].to_numpy(dtype=float), [1., 2., 4., 6.])
            self.assertNotIn("unused", result.columns)
            with self.assertRaisesRegex(ValueError, "not found in CSV: missing_sample"):
                self.quant["compare_pipeline"](path, reps[:2], ["missing_sample"],
                                               plot_results=False, write_output=False)


if __name__ == "__main__":
    unittest.main()

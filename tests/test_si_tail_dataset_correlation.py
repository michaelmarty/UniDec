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
from scipy.special import logit


class TestSITailDatasetCorrelation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_dir = Path(__file__).resolve().parents[3] / "Scripts" / "Annika" / "DDA Analysis"
        spec = importlib.util.spec_from_file_location("tail_dataset_correlation", script_dir / "SI_Tail_Dataset_Correlation.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(script_dir), *sys.path]):
            spec.loader.exec_module(cls.script)
        cls.datasets = ["A", "B", "C", "D", "E"]

    def tearDown(self):
        self.script.plt.close("all")

    def synthetic(self, interaction=False, confounding=False):
        # All nonempty identification patterns: selected union, known paired outcomes.
        lipids, occurrences = [], []
        identity = 0
        for cls in ["PC", "PE"]:
            for tail in [0, 1]:
                for pattern in range(1, 32):
                    weight = 1 + 2 * bool(pattern & 1) + bool(pattern & 2)
                    if interaction and tail and pattern & 1:
                        weight *= 3
                    if confounding:
                        weight = (1 + 4 * bool(pattern & (1 if cls == "PC" else 2)))
                        weight *= 4 if (tail == (cls == "PC")) else 1
                    for _ in range(weight):
                        species = f"lipid {identity}"
                        identity += 1
                        for m, method in enumerate(self.datasets):
                            if pattern & (1 << m):
                                lipids.append({"Dataset": method, "Ontology": cls, "Species": species})
                                occurrences.append({"Dataset": method, "Ontology": cls, "Species": species,
                                                    "Tail": "16:0" if tail else "18:0"})
        return pd.DataFrame(lipids), pd.DataFrame(occurrences)

    def test_grid_filtering_deduplication_tail_copies_and_all_ten_pairs(self):
        rows = [
            ("A", "PC", "PC 16:0_16:0", True), ("A", "PC", "PC 16:0/16:0", True),
            ("B", "PC", "PC 16:0_16:0", True), ("A", "FA", "FA 18:0", True),
            ("C", "FA", "FA 20:0", True), ("D", "FA", "FA 22:0", True),
            ("A", "FA", "FA 24:0", False), ("A", "PC", "PC 34:1", True)]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            pd.DataFrame(rows, columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status"]).to_csv(path, index=False)
            lipids, occurrences, _ = self.script.prepare_lipids(path, self.datasets)
        grid = self.script.pooled_observations(lipids, self.datasets)
        self.assertEqual(len(grid), 4 * 5)
        self.assertEqual(grid.Identified.sum(), 5)
        self.assertEqual(grid.loc[grid.Method.eq("E"), "Identified"].sum(), 0)
        frame = self.script.tail_observations(grid, occurrences, "16:0", self.datasets)
        self.assertEqual(frame.Tail_present.sum(), 5)
        self.assertEqual(frame.loc[frame.Tail_present.eq(0), "Lipid"].nunique(), 3)
        omnibus, pairs, diag, counts = self.script.interaction_tables(lipids, occurrences, self.datasets)
        self.assertTrue(pairs.groupby("Tail").size().eq(10).all())
        self.assertTrue(pairs["p-value"].isna().all())
        self.assertTrue(diag.Sparse.all())
        self.assertNotIn("24:0", set(omnibus.Tail))
        self.assertTrue(counts.groupby("Tail").Lipids.sum().eq(20).all())
        with self.assertRaisesRegex(ValueError, "at least two"):
            self.script.pooled_observations(lipids, ["A"])
        with self.assertRaisesRegex(ValueError, "No resolved"):
            self.script.pooled_observations(lipids, ["Y", "Z"])

    def test_unequal_overall_coverage_does_not_imply_interaction(self):
        lipids, occurrences = self.synthetic()
        omnibus, pairs, diag, _ = self.script.interaction_tables(lipids, occurrences, self.datasets)
        self.assertTrue(diag.Testable.all(), diag.Status.to_string())
        self.assertTrue(omnibus["p-value"].gt(.99).all())
        np.testing.assert_allclose(pairs["Log ROR (A / B)"], 0, atol=1e-8)
        grid = self.script.pooled_observations(lipids, self.datasets)
        self.assertGreater(grid.loc[grid.Method.eq("A"), "Identified"].mean(),
                           grid.loc[grid.Method.eq("E"), "Identified"].mean())

    def test_interaction_estimate_ci_direction_and_reference_invariance(self):
        lipids, occurrences = self.synthetic(interaction=True)
        omnibus, pairs, diag, _ = self.script.interaction_tables(lipids, occurrences, self.datasets)
        self.assertTrue(diag.Testable.all(), diag.Status.to_string())
        row = pairs.loc[pairs.Tail.eq("16:0") & pairs["Method pair"].eq("A vs B")].iloc[0]
        grid = self.script.tail_observations(self.script.pooled_observations(lipids, self.datasets), occurrences, "16:0", self.datasets)
        rates = grid.groupby(["Method", "Tail_present"], observed=True).Identified.mean()
        expected = (logit(rates["A", 1]) - logit(rates["A", 0])) - (logit(rates["B", 1]) - logit(rates["B", 0]))
        self.assertAlmostEqual(row["Log ROR (A / B)"], expected)
        self.assertGreater(row["ROR (A / B)"], 1)
        self.assertLess(row["BH q-value"], .05)
        self.assertTrue(omnibus["BH q-value"].lt(.05).all())
        _, reverse, _, _ = self.script.interaction_tables(lipids, occurrences, list(reversed(self.datasets)))
        reverse = reverse.loc[reverse.Tail.eq("16:0") & reverse["Method pair"].eq("B vs A")].iloc[0]
        self.assertAlmostEqual(reverse["Log ROR (A / B)"], -expected)
        self.assertAlmostEqual(reverse["95% CI lower"], 1 / row["95% CI upper"])
        self.assertAlmostEqual(reverse["p-value"], row["p-value"])
        # Independent calculation of cluster sandwich covariance for the saturated model.
        result, effects, diagnostic = self.script.fit_interaction(grid, self.datasets)
        x, y = result.model.exog, result.model.endog
        mu = result.fittedvalues.to_numpy()
        bread = np.linalg.inv(x.T @ ((mu * (1-mu))[:, None] * x))
        scores = pd.DataFrame(x * (y-mu)[:, None]).groupby(grid.Lipid.to_numpy()).sum().to_numpy()
        cov = bread @ scores.T @ scores @ bread
        contrast = effects[0] - effects[1]
        self.assertAlmostEqual(row["Robust SE"], np.sqrt(contrast @ cov @ contrast))

    def test_class_adjustment_removes_constructed_class_confounding(self):
        lipids, occurrences = self.synthetic(confounding=True)
        _, primary, _, _ = self.script.interaction_tables(lipids, occurrences, self.datasets)
        _, adjusted, diag, _ = self.script.interaction_tables(lipids, occurrences, self.datasets, adjusted=True)
        self.assertTrue(diag.Testable.all(), diag.Status.to_string())
        self.assertGreater(primary["Log ROR (A / B)"].abs().max(), .5)
        np.testing.assert_allclose(adjusted["Log ROR (A / B)"], 0, atol=1e-7)
        self.assertTrue(diag.Formula.str.contains(r"C\(Method\) \* C\(Ontology\)").all())

    def test_nonestimable_separation_singular_covariance_and_nonconvergence(self):
        lipids, occurrences = self.synthetic()
        grid = self.script.tail_observations(self.script.pooled_observations(lipids, self.datasets), occurrences, "16:0", self.datasets)
        no_absent = grid.assign(Tail_present=1)
        self.assertFalse(self.script.fit_interaction(no_absent, self.datasets)[2]["Testable"])
        confounded = grid.assign(Ontology=np.where(grid.Tail_present.eq(1), "PC", "PE"))
        self.assertIn("Rank-deficient", self.script.fit_interaction(confounded, self.datasets, True)[2]["Status"])
        self.assertTrue(self.script.has_separation(np.array([[1,0],[1,1]]), np.array([0,1])))
        separated = grid.copy()
        separated.loc[separated.Method.eq("A"), "Identified"] = 1
        self.assertTrue(self.script.fit_interaction(separated, self.datasets)[2]["Separation"])
        real, _, _ = self.script.fit_interaction(grid, self.datasets)
        with patch.object(self.script.GEE, "fit", return_value=real):
            with patch.object(real, "cov_params", return_value=np.zeros((10,10))):
                self.assertIn("Singular", self.script.fit_interaction(grid, self.datasets)[2]["Status"])
            real.converged = False
            self.assertIn("Nonconvergence", self.script.fit_interaction(grid, self.datasets)[2]["Status"])

    def test_nonestimable_interactions_and_model_warnings_are_not_significant(self):
        lipids, occurrences = self.synthetic()
        grid = self.script.tail_observations(self.script.pooled_observations(lipids, self.datasets), occurrences, "16:0", self.datasets)
        result, effects, diagnostic = self.script.fit_interaction(grid, self.datasets)
        effects[1] = effects[0]
        with patch.object(self.script, "fit_interaction", return_value=(result, effects, diagnostic)):
            omnibus, pairs, _, _ = self.script.interaction_tables(lipids, occurrences, self.datasets)
        self.assertTrue(omnibus["p-value"].isna().all())
        failed = pairs.loc[pairs["Method pair"].eq("A vs B")]
        self.assertFalse(failed.Testable.any())
        self.assertTrue(failed["BH q-value"].isna().all())
        def warn_fit(*args, **kwargs):
            import warnings
            warnings.warn("Simulated unstable fit", RuntimeWarning)
            return result
        with patch.object(self.script.GEE, "fit", side_effect=warn_fit):
            self.assertFalse(self.script.fit_interaction(grid, self.datasets)[2]["Testable"])

    def test_polarity_recovery_preserves_both_labels_before_deduplication(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            pd.DataFrame([
                ("A", "FA", "FA 16:0", True, "Positive"),
                ("A", "FA", "FA 16:0", True, "Negative"),
                ("B", "FA", "FA 18:0", True, "Positive"),
                ("B", "FA", "FA 20:0", False, "unknown"),
            ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status", "Polarity"]).to_csv(path,index=False)
            lipids, _, audit = self.script.prepare_lipids(path, self.datasets)
            catalogs, note = self.script.polarity_catalogs(path, audit)
            self.assertEqual(len(lipids), 2)
            self.assertEqual(len(catalogs["positive"][0]), 2)
            self.assertEqual(len(catalogs["negative"][0]), 1)
            source = pd.read_csv(path)
            source.loc[0, "Polarity"] = "unknown"
            source.to_csv(path,index=False)
            self.assertEqual(self.script.polarity_catalogs(path, audit)[0], {})
            source.drop(columns="Polarity").to_csv(path,index=False)
            self.assertEqual(self.script.polarity_catalogs(path, audit)[0], {})

    def test_bh_families_and_display_independent_exports(self):
        np.testing.assert_allclose(self.script.bh_adjust([.04, np.nan, .01, .03]), [.04, np.nan, .03, .04], equal_nan=True)
        lipids, occurrences = self.synthetic(interaction=True)
        audit = pd.DataFrame({"Chains": [("16:0",)], "Exclusion reason": [""]})
        tables = []
        with TemporaryDirectory() as directory:
            for threshold in [0, 100]:
                output = Path(directory) / str(threshold)
                with patch.object(self.script, "prepare_lipids", return_value=(lipids, occurrences, audit)), \
                     patch.object(self.script, "polarity_catalogs", return_value=({}, "No polarity")), \
                     patch.object(self.script, "save_figure"), \
                     patch.object(sys, "argv", ["script", "--datasets", *self.datasets, "--output-dir", str(output), "--min-tail-percent", str(threshold), "--no-show"]), \
                     patch.object(self.script, "frequency_summary", return_value=pd.DataFrame()), \
                     patch.object(self.script, "tail_prevalence_tables", return_value=(None,None,pd.DataFrame(index=["16:0"] if threshold==0 else []))), \
                     contextlib.redirect_stdout(io.StringIO()):
                    self.script.main()
                # Recompute the same in-memory tables to check display-independent inference.
                analyses = [self.script.interaction_tables(lipids, occurrences, self.datasets, name, adjusted)
                            for name, adjusted in [("Primary", False), ("Class adjusted", True)]]
                table = pd.concat([a[1] for a in analyses], ignore_index=True)
                joint = pd.concat([a[0] for a in analyses], ignore_index=True)
                tables.append(table)
                for result in [table, joint]:
                    for _, group in result.groupby("Analysis"):
                        np.testing.assert_allclose(group["BH q-value"], self.script.bh_adjust(group["p-value"]), equal_nan=True)
                self.assertEqual(list(output.iterdir()), [])  # Figure saving is mocked.
            pd.testing.assert_frame_equal(*tables)

    def test_heatmap_zero_centering_signed_log_scale_and_masking(self):
        table = pd.DataFrame({"Tail": ["16:0", "18:0"], "Method pair": ["A vs B"]*2,
                              "Log ROR (A / B)": [2., -1.], "BH q-value": [.01, np.nan], "Testable": [True, False]})
        with patch.object(self.script, "save_figure"):
            fig, differences, _ = self.script.plot_comparisons(table, ["16:0", "18:0"], ["A", "B"], Path("output"))
            fig.canvas.draw()
            self.assertEqual((fig.axes[0].images[0].norm.vmin, fig.axes[0].images[0].norm.vmax), (-2, 2))
            self.assertEqual([text.get_text() for text in fig.axes[0].texts], ["2.0*", "N/A"])
            self.assertTrue(np.isnan(differences.loc["18:0", "A vs B"]))
            fig, _, _ = self.script.plot_comparisons(table, [], ["A", "B"], Path("output"))
            fig.canvas.draw()


if __name__ == "__main__":
    unittest.main()

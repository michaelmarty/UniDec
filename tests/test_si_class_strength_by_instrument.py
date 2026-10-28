import contextlib
import importlib.util
import io
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


class TestSIClassStrengthByInstrument(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location(
            "si_class_strength_by_instrument", directory / "SI_Class_Strength_By_Instrument.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]), patch.dict(os.environ, {"MPLBACKEND": "Agg"}):
            spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def source(self):
        return pd.DataFrame([
            {"Dataset": dataset, "Ontology": "PC", "Polarity": polarity,
             "Metabolite name": f"PC {i}|chains", "Overall Match Status": accepted,
             "Validation Status": "", "Class Strength": score}
            for polarity in ["positive", "negative"]
            for dataset, scores in [("OTOT", [0.7, 0.9, 1.0]), ("Stellar", [0.2, 0.4, 0.6])]
            for i, (score, accepted) in enumerate(zip(scores, [True, True, False]))
        ])

    def test_normalizes_raw_strength_using_population_specific_maxima(self):
        source = self.source()
        source["Max-Normalized Class Strength"] = 0  # Pooled input scores are ignored.
        points = self.script.prepare_strength(source)
        self.assertEqual(len(points), 12)
        for _, row in points.iterrows():
            expected = source.loc[source["Dataset"].eq(row["Dataset"])
                                  & source["Polarity"].eq(row["Polarity"])
                                  & source["Metabolite name"].str.startswith(row["Lipid Name"] + "|"),
                                  "Class Strength"].iloc[0]
            maximum = source.loc[source["Dataset"].eq(row["Dataset"])
                                 & source["Polarity"].eq(row["Polarity"]), "Class Strength"].max()
            expected /= maximum
            self.assertEqual(row[self.script.NORMALIZED], expected)
        accepted = self.script.prepare_strength(source, all_ids=False)
        self.assertEqual(len(accepted), 8)
        self.assertTrue(accepted["Excluded rows"].eq(0).all())
        excluded = self.script.prepare_strength(source, excluded_only=True)
        for _, row in excluded.iterrows():
            original = points.loc[points["Dataset"].eq(row["Dataset"])
                                  & points["Polarity"].eq(row["Polarity"])
                                  & points["Lipid ID"].eq(row["Lipid ID"])].iloc[0]
            self.assertEqual(row[self.script.NORMALIZED], original[self.script.NORMALIZED])
        for (dataset, polarity), group in accepted.groupby(["Dataset", "Polarity"]):
            maximum = source.loc[source["Dataset"].eq(dataset) & source["Polarity"].eq(polarity)
                                 & source["Overall Match Status"], "Class Strength"].max()
            np.testing.assert_allclose(group[self.script.NORMALIZED], group["Class Strength"] / maximum)
            self.assertEqual(group[self.script.NORMALIZED].max(), 1)
        changed = source.copy()
        changed.loc[changed["Dataset"].eq("Stellar"), "Class Strength"] *= 100
        pd.testing.assert_series_equal(self.script.prepare_strength(changed)[self.script.NORMALIZED],
                                       points[self.script.NORMALIZED])

    def test_duplicate_names_collapse_to_median_and_preserve_status_counts(self):
        source = self.source().iloc[[0]].copy()
        duplicate = source.assign(**{"Dataset": " otot\u00a0", "Ontology": "pc",
                                    "Metabolite name": " pc  0 |other chains",
                                    "Class Strength": 0.3, "Overall Match Status": False})
        points = self.script.prepare_strength(pd.concat([source, duplicate], ignore_index=True))
        self.assertEqual(len(points), 1)
        self.assertAlmostEqual(points.iloc[0]["Class Strength"], 0.5)
        self.assertAlmostEqual(points.iloc[0]["Normalization maximum"], 0.5)
        self.assertEqual(points.iloc[0][self.script.NORMALIZED], 1)
        self.assertEqual(points.iloc[0]["Feature rows"], 2)
        self.assertEqual(points.iloc[0]["Accepted rows"], 1)
        self.assertEqual(points.iloc[0]["Excluded rows"], 1)

    def test_invalid_unknown_provisional_and_unrepresented_classes_are_excluded(self):
        base = self.source().iloc[[0]].copy()
        rows = [base]
        for score in [np.nan, np.inf, -1, "bad"]:
            rows.append(base.assign(**{"Class Strength": score}))
        rows += [base.assign(Dataset="Unknown"), base.assign(**{"Validation Status": " provisional "}),
                 base.assign(**{"Overall Match Status": ""}),
                 base.assign(Ontology="FA", **{"Overall Match Status": False})]
        points = self.script.prepare_strength(pd.concat(rows, ignore_index=True))
        self.assertEqual(len(points), 1)
        self.assertEqual(points.iloc[0]["Feature rows"], 1)
        with self.assertRaisesRegex(ValueError, "Class Strength"):
            self.script.prepare_strength(base.drop(columns="Class Strength"))
        with self.assertRaisesRegex(ValueError, "No eligible"):
            self.script.prepare_strength(base.assign(**{"Class Strength": np.nan}))

    def test_summary_and_descriptive_differences_and_missing_instruments(self):
        points = self.script.prepare_strength(self.source())
        summary = self.script.strength_summary(points)
        comparisons = self.script.compare_with_otot(summary)
        stellar = comparisons.loc[comparisons["Dataset"].eq("Stellar")]
        np.testing.assert_allclose(stellar["Median"], 0.4 / 0.6)
        np.testing.assert_allclose(stellar["OT/OT Median"], 0.9)
        np.testing.assert_allclose(stellar["Median difference (dataset - OT/OT)"], 0.4 / 0.6 - 0.9)
        missing = comparisons.loc[comparisons["Dataset"].eq("OTIT")]
        self.assertTrue(missing["Median"].isna().all())
        self.assertTrue(missing["Status"].str.startswith("Unavailable").all())
        self.assertEqual(len(summary), 10)
        self.assertFalse(any("p-value" in c for c in comparisons.columns))

    def test_plots_show_all_scores_and_keep_polarities_separate(self):
        points = self.script.prepare_strength(self.source())
        figures = self.script.plot_strength(points)
        self.assertEqual(set(figures), {"positive", "negative"})
        for figure in figures.values():
            axis = figure.axes[0]
            self.assertEqual(sum(len(c.get_offsets()) for c in axis.collections), 6)
            self.assertIsNone(axis.spines["left"].get_bounds())
            self.assertEqual(axis.get_ylim(), (-0.03, 1.18))
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            title = axis.title.get_window_extent(renderer)
            self.assertGreater(title.y0, axis.transData.transform((0, 1))[1])
            for text in axis.texts:
                self.assertFalse(title.overlaps(text.get_window_extent(renderer)))
            figure.canvas.draw()
        figures = self.script.plot_strength(points.loc[points["Polarity"].eq("positive")])
        self.assertTrue(any(t.get_text() == "No eligible scores"
                            for axis in figures["negative"].axes for t in axis.texts))

    def test_main_exports_only_class_strength_analysis(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            self.source().to_csv(root / "input.csv", index=False)
            with contextlib.redirect_stdout(io.StringIO()):
                self.script.main(["--input", str(root / "input.csv"), "--no-show"])
            base = root / "Class_Strength_By_Instrument"
            output = base / "Excluded_IDs"
            expected = ["Points.csv", "Summary.csv", "OTOT_Comparison.csv", "Statistics.csv", "Pairwise_Statistics.csv", "Normalization.csv", "Methods.txt",
                        "Positive.png", "Positive.pdf", "Negative.png", "Negative.pdf"]
            self.assertEqual(len(list(output.iterdir())), len(expected))
            for suffix in expected:
                self.assertGreater((output / f"SI_Class_Strength_By_Instrument_{suffix}").stat().st_size, 0)
            self.assertEqual(self.script.plt.get_fignums(), [])
            points = pd.read_csv(output / "SI_Class_Strength_By_Instrument_Points.csv")
            self.assertEqual(len(points), 4)
            self.assertTrue(points["Accepted rows"].eq(0).all())
            accepted = pd.read_csv(base / "Accepted_IDs/SI_Class_Strength_By_Instrument_Points.csv")
            self.assertEqual(len(accepted), 8)
            self.assertTrue(accepted["Excluded rows"].eq(0).all())
            self.assertFalse((base / "All_IDs").exists())
            pairwise = pd.read_csv(output / "SI_Class_Strength_By_Instrument_Pairwise_Statistics.csv")
            self.assertEqual(len(pairwise), 20)
            otot = pd.read_csv(output / "SI_Class_Strength_By_Instrument_Statistics.csv")
            np.testing.assert_allclose(
                otot["Holm-adjusted p-value"],
                pairwise.loc[pairwise["Reference Dataset"].eq("OTOT"), "Holm-adjusted p-value"],
                equal_nan=True)
            accepted_norm = pd.read_csv(base / "Accepted_IDs/SI_Class_Strength_By_Instrument_Normalization.csv")
            excluded_norm = pd.read_csv(output / "SI_Class_Strength_By_Instrument_Normalization.csv")
            self.assertTrue(accepted_norm["Normalization population"].eq("Accepted only").all())
            self.assertTrue(excluded_norm["Normalization population"].eq("Accepted + excluded").all())
            np.testing.assert_allclose(accepted_norm.loc[accepted_norm["Dataset"].eq("OTOT"), "Normalization maximum"], 0.9)
            np.testing.assert_allclose(excluded_norm.loc[excluded_norm["Dataset"].eq("OTOT"), "Normalization maximum"], 1.0)

    def test_single_accepted_coq_uses_accepted_maximum_even_with_stronger_excluded_feature(self):
        source = self.source().iloc[[0]].assign(Ontology="CoQ", **{
            "Metabolite name": "CoQ 10", "Class Strength": 5})
        stronger = source.assign(**{"Metabolite name": "CoQ 9", "Class Strength": 10,
                                    "Overall Match Status": False})
        source = pd.concat([source, stronger], ignore_index=True)
        accepted = self.script.prepare_strength(source, all_ids=False)
        excluded = self.script.prepare_strength(source, excluded_only=True)
        self.assertEqual(accepted[self.script.NORMALIZED].tolist(), [1.0])
        self.assertEqual(accepted["Normalization maximum"].tolist(), [5])
        self.assertEqual(excluded["Normalization maximum"].tolist(), [10])
        self.assertEqual(excluded[self.script.NORMALIZED].tolist(), [1.0])

    def test_excluded_uses_both_statuses_and_accepted_zero_group_stays_zero(self):
        source = self.source().iloc[[0]].assign(**{"Class Strength": 10})
        weaker = source.assign(**{"Metabolite name": "PC 1", "Class Strength": 5,
                                  "Overall Match Status": False})
        frame = pd.concat([source, weaker], ignore_index=True)
        excluded = self.script.prepare_strength(frame, excluded_only=True)
        self.assertEqual(excluded[self.script.NORMALIZED].tolist(), [0.5])
        self.assertEqual(excluded["Normalization maximum"].tolist(), [10])
        frame.loc[0, "Class Strength"] = 0
        accepted = self.script.prepare_strength(frame, all_ids=False)
        self.assertEqual(accepted[self.script.NORMALIZED].tolist(), [0.0])
        self.assertEqual(accepted["Normalization maximum"].tolist(), [0])

    def test_coq_cl_and_cer_duplicates_are_collapsed_before_maximum(self):
        base = self.source().iloc[[0]]
        examples = [("CoQ", "positive", [0.809732, 6.189240]),
                    ("CL", "negative", [36.717130, 47.270897]),
                    ("Cer", "negative", [24.683578, 40.820958])]
        for ontology, polarity, strengths in examples:
            with self.subTest(ontology=ontology):
                source = pd.concat([
                    base.assign(Ontology=ontology, Polarity=polarity,
                                **{"Metabolite name": ontology + " species", "Class Strength": strength})
                    for strength in strengths], ignore_index=True)
                source["Adduct type"] = ["[M-H]-", "[M+CH3COO]-"]
                accepted = self.script.prepare_strength(source, all_ids=False)
                self.assertEqual(len(accepted), 1)
                self.assertAlmostEqual(accepted.iloc[0]["Class Strength"], np.median(strengths))
                self.assertAlmostEqual(accepted.iloc[0]["Normalization maximum"], np.median(strengths))
                self.assertEqual(accepted.iloc[0][self.script.NORMALIZED], 1)

    def test_largest_lipid_median_sets_maximum_not_largest_feature(self):
        base = self.source().iloc[[0]]
        source = pd.concat([
            base.assign(**{"Metabolite name": name, "Class Strength": score})
            for name, score in [("PC 0", 2), ("PC 0", 10), ("PC 1", 8)]
        ], ignore_index=True)
        accepted = self.script.prepare_strength(source, all_ids=False)
        np.testing.assert_allclose(accepted["Class Strength"], [6, 8])
        np.testing.assert_allclose(accepted[self.script.NORMALIZED], [0.75, 1])
        np.testing.assert_allclose(accepted["Normalization maximum"], [8, 8])

    def test_excluded_reference_uses_status_specific_lipid_medians(self):
        base = self.source().iloc[[0]]
        source = pd.concat([
            base.assign(**{"Class Strength": score, "Overall Match Status": accepted})
            for score, accepted in [(2, True), (10, True), (1, False), (3, False)]
        ], ignore_index=True)
        accepted = self.script.prepare_strength(source, all_ids=False)
        excluded = self.script.prepare_strength(source, excluded_only=True)
        self.assertEqual(accepted["Class Strength"].tolist(), [6])
        self.assertEqual(accepted[self.script.NORMALIZED].tolist(), [1])
        self.assertEqual(excluded["Class Strength"].tolist(), [2])
        self.assertEqual(excluded["Normalization maximum"].tolist(), [6])
        self.assertEqual(excluded[self.script.NORMALIZED].tolist(), [1 / 3])
        source.loc[~source["Overall Match Status"], "Class Strength"] = [5, 9]
        excluded = self.script.prepare_strength(source, excluded_only=True)
        self.assertEqual(excluded["Normalization maximum"].tolist(), [7])
        self.assertEqual(excluded[self.script.NORMALIZED].tolist(), [1])

    def test_mixed_status_duplicates_have_separate_medians(self):
        source = self.source().iloc[[0]].copy()
        excluded = source.assign(**{"Class Strength": 0.1, "Overall Match Status": False})
        source = pd.concat([source, excluded], ignore_index=True)
        accepted = self.script.prepare_strength(source, all_ids=False)
        excluded = self.script.prepare_strength(source, excluded_only=True)
        self.assertEqual(accepted["Lipid ID"].tolist(), excluded["Lipid ID"].tolist())
        self.assertAlmostEqual(accepted.iloc[0][self.script.NORMALIZED], 1)
        self.assertAlmostEqual(excluded.iloc[0][self.script.NORMALIZED], 0.1 / 0.7)
        self.assertEqual(excluded.iloc[0]["Accepted rows"], 0)
        self.assertEqual(excluded.iloc[0]["Feature rows"], 1)

    def test_original_classes_and_polarities_use_separate_maxima_and_zero_is_valid(self):
        source = self.source().iloc[[0]].copy()
        rows = [source.assign(**{"Class Strength": 10}),
                source.assign(**{"Class Strength": 20, "Metabolite name": "PC 1"}),
                source.assign(Ontology="LPC", **{"Class Strength": 3}),
                source.assign(Polarity="negative", **{"Class Strength": 0})]
        points = self.script.prepare_strength(pd.concat(rows, ignore_index=True))
        pc = points.loc[points["Ontology"].eq("PC") & points["Polarity"].eq("positive")]
        np.testing.assert_allclose(pc[self.script.NORMALIZED], [0.5, 1])
        self.assertEqual(points.loc[points["Ontology"].eq("LPC"), self.script.NORMALIZED].iloc[0], 1)
        self.assertEqual(points.loc[points["Polarity"].eq("negative"), self.script.NORMALIZED].iloc[0], 0)

    def test_empty_excluded_population_omits_panels_and_keeps_missing_statistics(self):
        source = self.source().assign(**{"Overall Match Status": True})
        reference = self.script.prepare_strength(source)
        excluded = self.script.prepare_strength(source, excluded_only=True, allow_empty=True)
        self.assertTrue(excluded.empty)
        summary = self.script.strength_summary(excluded, reference)
        statistics = self.script.instrument_statistics(excluded, reference)
        self.assertTrue(summary["Lipids"].eq(0).all())
        self.assertTrue(summary["Mean"].isna().all())
        self.assertTrue(statistics["p-value"].isna().all())
        accepted_figures = self.script.plot_strength(reference, all_ids=False, reference=reference)
        excluded_figures = self.script.plot_strength(excluded, excluded_only=True, reference=reference)
        for mode in ["positive", "negative"]:
            self.assertIn("PC", [ax.get_title() for ax in accepted_figures[mode].axes])
            self.assertTrue(all(not ax.get_title() for ax in excluded_figures[mode].axes))
            self.assertTrue(any(t.get_text() == "No eligible scores"
                                for ax in excluded_figures[mode].axes for t in ax.texts))
            self.assertIn("Excluded only", excluded_figures[mode]._suptitle.get_text())
            excluded_figures[mode].canvas.draw()

    def test_zero_panels_omitted_but_zero_observations_in_nonzero_class_remain(self):
        points = self.script.prepare_strength(self.source())
        zero = points.assign(**{"Plot Class": "Zero", self.script.NORMALIZED: 0.0})
        points.loc[points.loc[points["Polarity"].eq("positive")].index[0], self.script.NORMALIZED] = 0.0
        combined = pd.concat([points, zero], ignore_index=True)
        figures = self.script.plot_strength(combined, reference=combined)
        for figure in figures.values():
            self.assertIn("PC", [ax.get_title() for ax in figure.axes])
            self.assertNotIn("Zero", [ax.get_title() for ax in figure.axes])
        axis = figures["positive"].axes[0]
        offsets = np.concatenate([collection.get_offsets() for collection in axis.collections])
        self.assertTrue((offsets[:, 1] == 0).any())

    def test_brackets_connect_only_significant_pairs(self):
        points = self.statistical_fixture()
        comparisons = pd.DataFrame([
            {"Plot Class": "PC", "Polarity": "positive", "Reference Dataset": baseline,
             "Dataset": dataset, "Holm-adjusted p-value": p, "Significant": p < 0.05}
            for baseline, dataset, p in [("OTIT", "Stellar", 0.02),
                                         ("OTOT", "Stellar", 0.0001),
                                         ("ITIT", "Stellar", 0.7)]])
        figures = self.script.plot_strength(points, statistics=comparisons)
        axis = figures["positive"].axes[0]
        brackets = [line for line in axis.lines if np.min(line.get_ydata()) > 1.1]
        self.assertEqual(len(brackets), 2)
        self.assertEqual({tuple(line.get_xdata()) for line in brackets},
                         {(1, 1, 3, 3), (0, 0, 3, 3)})
        stars = [text.get_text() for text in axis.texts if "*" in text.get_text()]
        self.assertCountEqual(stars, ["*", "***"])
        annotations = [text for text in axis.texts if "*" in text.get_text()]
        for text, bracket in zip(annotations, brackets):
            x = bracket.get_xdata()
            self.assertEqual(text.get_position(), ((x[0] + x[-1]) / 2, max(bracket.get_ydata())))
            self.assertIsNotNone(text.get_bbox_patch())
        figures["positive"].canvas.draw()
        renderer = figures["positive"].canvas.get_renderer()
        title = axis.title.get_window_extent(renderer)
        for text in axis.texts:
            self.assertFalse(title.overlaps(text.get_window_extent(renderer)))
        self.assertNotIn("OT/OT", figures["positive"]._supxlabel.get_text())

    def statistical_fixture(self):
        rng = np.random.default_rng(82)
        shared = rng.normal(0, 0.04, 30)
        return pd.DataFrame([
            {"Plot Class": "PC", "Polarity": "positive", "Dataset": dataset,
             "Lipid ID": f"PC {i}", self.script.NORMALIZED: mean + shared[i] + rng.normal(0, 0.025)}
            for dataset, mean in [("OTOT", 0.8), ("Stellar", 0.5)] for i in range(30)
        ])

    def test_statistics_match_clustered_formula_reference_and_full_holm_family(self):
        from statsmodels.formula.api import ols
        from statsmodels.stats.multitest import multipletests
        points = self.statistical_fixture()
        reference_data = points.rename(columns={self.script.NORMALIZED: "score", "Lipid ID": "lipid"})
        reference = ols("score ~ 0 + C(Dataset)", reference_data).fit(
            cov_type="cluster", cov_kwds={"groups": reference_data["lipid"],
                                         "use_correction": True, "df_correction": True}, use_t=True)
        contrast = reference.t_test([-1, 1])
        stats = self.script.instrument_statistics(points)
        row = stats.loc[stats["Dataset"].eq("Stellar")].iloc[0]
        self.assertAlmostEqual(row["Mean difference (dataset - OT/OT)"], float(contrast.effect.item()))
        self.assertAlmostEqual(row["Standard error"], float(contrast.sd.item()))
        self.assertAlmostEqual(row["p-value"], float(contrast.pvalue.item()))
        self.assertEqual(row["DF"], 29)
        p = stats["p-value"].fillna(1).to_numpy()
        adjusted = multipletests(p, method="holm")[1]
        np.testing.assert_allclose(row["Holm-adjusted p-value"], adjusted[stats["Dataset"].eq("Stellar")][0])
        ci = contrast.conf_int()[0]
        np.testing.assert_allclose([row["CI95 low"], row["CI95 high"]], ci)
        self.assertTrue(row["Significant"])
        figures = self.script.plot_strength(points, statistics=stats)
        self.assertTrue(any(t.get_text() == "***" for ax in figures["positive"].axes for t in ax.texts))

    def test_direct_comparison_without_otot_matches_clustered_reference(self):
        from statsmodels.formula.api import ols
        from statsmodels.stats.multitest import multipletests
        points = self.statistical_fixture().replace({"Dataset": {"OTOT": "OTIT"}})
        data = points.rename(columns={self.script.NORMALIZED: "score", "Lipid ID": "lipid"})
        model = ols("score ~ 0 + C(Dataset)", data).fit(
            cov_type="cluster", cov_kwds={"groups": data["lipid"],
                                         "use_correction": True, "df_correction": True}, use_t=True)
        contrast = model.t_test([-1, 1])
        table = self.script.instrument_statistics(points, all_pairs=True)
        self.assertEqual(len(table), 10)
        row = table.loc[table["Reference Dataset"].eq("OTIT") & table["Dataset"].eq("Stellar")].iloc[0]
        self.assertAlmostEqual(row["Mean difference (dataset - reference)"], float(contrast.effect.item()))
        self.assertAlmostEqual(row["Standard error"], float(contrast.sd.item()))
        self.assertAlmostEqual(row["p-value"], float(contrast.pvalue.item()))
        np.testing.assert_allclose([row["CI95 low"], row["CI95 high"]], contrast.conf_int()[0])
        adjusted = multipletests(table["p-value"].fillna(1), method="holm")[1]
        np.testing.assert_allclose(table.loc[table["p-value"].notna(), "Holm-adjusted p-value"],
                                   adjusted[table["p-value"].notna()])
        self.assertTrue(table.loc[table["Reference Dataset"].eq("OTOT"), "p-value"].isna().all())

    def test_sparse_constant_and_no_shared_lipids_statistics(self):
        points = self.statistical_fixture()
        constant = points.assign(**{self.script.NORMALIZED: 0.5})
        stats = self.script.instrument_statistics(constant)
        self.assertTrue(stats["p-value"].isna().all())
        sparse = points.loc[points["Dataset"].eq("OTOT") | points.index.isin([30])]
        stats = self.script.instrument_statistics(sparse)
        self.assertTrue(stats["p-value"].isna().all())
        independent = points.copy()
        independent["Lipid ID"] += independent["Dataset"]
        stats = self.script.instrument_statistics(independent)
        row = stats.loc[stats["Dataset"].eq("Stellar")].iloc[0]
        self.assertEqual(row["DF"], 59)
        self.assertTrue(np.isfinite(row["p-value"]))


if __name__ == "__main__":
    unittest.main()

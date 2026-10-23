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
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd


class TestSIClassSharingCorrelation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_dir = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                      / "DDA Analysis")
        spec = importlib.util.spec_from_file_location(
            "si_class_sharing_correlation", script_dir / "SI_Class_Sharing_Correlation.py"
        )
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(script_dir), *sys.path]):
            spec.loader.exec_module(cls.script)

    def tearDown(self):
        self.script.plt.close("all")

    def source(self):
        return pd.DataFrame([
            ("OTOT", "Cer_AP", "Cer 34:1|resolved 1", True),
            ("OTOT", "Cer_NS", "Cer 34:1|resolved 2", True),
            ("OTIT", "Cer_NS", "Cer 34:1", "TRUE"),
            ("OTOT", "LPC", "LPC 18:0", True),
            ("OTOT", "EtherPC", "PC O-34:1", True),
            ("OTIT", "EtherPC", "PC O-34:1", True),
            ("OTOT", "TG", "TG 50:1", True),
            ("OTOT", "TG", "TG 52:1", True),
            ("OTOT", "TG", "TG 54:1", True),
            ("ITIT", "TG", "TG 50:1", False),
            ("Other", "TG", "TG 50:1", True),
        ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status"])

    def test_grouping_deduplication_and_coverage_match_figure4(self):
        table = self.script.class_sharing_table(self.source())
        expected = {"Cer": (1, 1, 100), "PC": (2, 1, 50), "TG": (3, 0, 0)}
        for name, (total, shared, percent) in expected.items():
            self.assertEqual(table.at[name, "Pooled unique lipids"], total)
            self.assertEqual(table.at[name, "Shared unique lipids"], shared)
            self.assertEqual(table.at[name, "Shared lipids (%)"], percent)
        counts, _ = self.script.prepare_coverage(self.source())
        np.testing.assert_array_equal(table.iloc[:, :5].to_numpy(), counts.to_numpy())
        stats = self.script.correlation_statistics(table).iloc[0]
        self.assertAlmostEqual(stats["Spearman rho"], -1)
        self.assertEqual(stats["Classes"], 3)

    def test_missing_annotations_and_no_eligible_ids(self):
        source = self.source()
        source.loc[len(source)] = ["OTOT", None, "Unknown", True]
        source.loc[len(source)] = ["OTOT", "PC", None, True]
        pd.testing.assert_frame_equal(self.script.class_sharing_table(source),
                                      self.script.class_sharing_table(self.source()))
        with self.assertRaisesRegex(ValueError, "No accepted"):
            self.script.class_sharing_table(source.assign(**{"Overall Match Status": False}))
        with self.assertRaisesRegex(ValueError, "Ontology"):
            self.script.class_sharing_table(source.drop(columns="Ontology"))

    def test_all_id_counts_retain_duplicate_features_and_adducts(self):
        source = self.source()
        source.loc[len(source)] = [" OTOT\u00a0", "Cer_NS", "Cer 34:1|adduct", True]
        source.loc[len(source)] = ["OTOT", None, "Unknown", True]
        source.loc[len(source)] = ["OTOT", "PC", None, True]
        pooled, tables = self.script.all_id_counts(source)
        self.assertEqual(pooled.to_dict(), {"Cer": 4, "PC": 3, "TG": 3})
        self.assertEqual(tables.at[("OTOT", "Cer")], 3)
        self.assertEqual(tables.at[("OTIT", "Cer")], 1)
        self.assertNotIn("Other", tables.index.get_level_values("Dataset"))
        self.assertNotIn("ITIT", tables.index.get_level_values("Dataset"))

    def test_dataset_sharing_uses_local_counts_and_global_membership(self):
        tables = self.script.dataset_sharing_tables(self.source())
        otot, otit = tables["OTOT"], tables["OTIT"]
        self.assertEqual(otot.at["PC", "Dataset unique lipids"], 2)
        self.assertEqual(otot.at["PC", "Shared unique lipids"], 1)
        self.assertEqual(otot.at["PC", "Shared lipids (%)"], 50)
        self.assertEqual(otit.at["PC", "Dataset unique lipids"], 1)
        self.assertEqual(otit.at["PC", "Shared lipids (%)"], 100)
        self.assertEqual(otot.at["Cer", "Dataset unique lipids"], 1)
        self.assertEqual(otit.at["TG", "Dataset unique lipids"], 0)
        self.assertTrue(pd.isna(otit.at["TG", "Shared lipids (%)"]))
        self.assertEqual(otot.at["TG", "Shared lipids (%)"], 0)
        stats = self.script.correlation_statistics(otot, "Dataset unique lipids").iloc[0]
        self.assertAlmostEqual(stats["Spearman rho"], -1)
        self.assertEqual(stats["Classes"], 3)
        # A method with no accepted IDs has no classes for correlation, not zeros.
        empty = tables["ITIT"]
        self.assertTrue(empty["Shared lipids (%)"].isna().all())
        stats = self.script.correlation_statistics(empty, "Dataset unique lipids")
        self.assertEqual(stats.iloc[0]["Classes"], 0)
        fig = self.script.plot_dataset_correlations(
            self.script.class_sharing_table(self.source()), tables,
            pd.concat([
                self.script.correlation_statistics(self.script.class_sharing_table(self.source())).assign(Dataset="Pooled"),
                *[self.script.correlation_statistics(t, "Dataset unique lipids").assign(Dataset=d)
                  for d, t in tables.items()],
            ]),
        )
        self.assertEqual(len(fig.axes), 6)
        self.assertEqual(len(fig.axes[3].collections[0].get_offsets()), 0)

    def test_undefined_correlations_are_explicit(self):
        table = self.script.class_sharing_table(self.source())
        for frame, reason in [
            (table.iloc[:2], "fewer than three"),
            (table.assign(**{"Shared lipids (%)": 50}), "constant"),
            (table.assign(**{"Pooled unique lipids": 1}), "constant"),
        ]:
            stats = self.script.correlation_statistics(frame)
            self.assertTrue(pd.isna(stats.iloc[0]["Spearman rho"]))
            self.assertIn(reason, stats.iloc[0]["Status"])
            fig = self.script.plot_correlation(frame, stats)
            self.assertTrue(any("Undefined" in label.get_text() for label in fig.axes[0].texts))

    def test_polarities_are_split_before_sharing_and_deduplication(self):
        source = pd.DataFrame([
            ("OTOT", "PC", "PC 34:1", True, " Positive\u00a0"),
            ("OTIT", "PC", "PC 34:1", True, "NEGATIVE"),
            ("OTIT", "PC", "PC 34:1|duplicate", True, "negative"),
            ("OTOT", "PC", "PC 36:2", True, "positive"),
            ("OTIT", "PC", "PC 36:2", True, "positive"),
            ("ITIT", "PC", "PC 34:1", True, None),
            ("Stellar", "PC", "PC 34:1", True, "unknown"),
        ], columns=["Dataset", "Ontology", "Metabolite name", "Overall Match Status", "Polarity"])
        with contextlib.redirect_stdout(io.StringIO()) as output:
            subsets = self.script.polarity_subsets(source)
        self.assertIn("Excluded 2", output.getvalue())
        positive = self.script.class_sharing_table(subsets["positive"])
        negative = self.script.class_sharing_table(subsets["negative"])
        self.assertEqual(positive.at["PC", "Pooled unique lipids"], 2)
        self.assertEqual(positive.at["PC", "Shared lipids (%)"], 50)
        self.assertEqual(negative.at["PC", "Pooled unique lipids"], 1)
        self.assertEqual(negative.at["PC", "Shared lipids (%)"], 0)
        tables = self.script.dataset_sharing_tables(subsets["positive"])
        self.assertEqual(tables["OTOT"].at["PC", "Shared lipids (%)"], 50)
        self.assertEqual(tables["OTIT"].at["PC", "Shared lipids (%)"], 100)

    def test_missing_polarity_column_and_empty_polarity(self):
        with self.assertRaisesRegex(ValueError, "Polarity"):
            self.script.polarity_subsets(self.source())
        with contextlib.redirect_stdout(io.StringIO()) as output:
            subsets = self.script.polarity_subsets(self.source().assign(Polarity="positive"))
        self.assertEqual(list(subsets), ["positive"])
        self.assertIn("Skipping negative", output.getvalue())
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "recognized polarity"):
            self.script.polarity_subsets(self.source().assign(Polarity=None))

    def test_overlay_keeps_both_series_and_separate_statistics(self):
        analyses = {}
        for polarity, source in [("positive", self.source()),
                                 ("negative", self.source().iloc[3:4])]:
            table = self.script.class_sharing_table(source)
            tables = self.script.dataset_sharing_tables(source)
            stats = pd.concat([
                self.script.correlation_statistics(table).assign(Dataset="Pooled"),
                *[self.script.correlation_statistics(t, "Dataset unique lipids").assign(Dataset=d)
                  for d, t in tables.items()],
            ])
            analyses[polarity] = table, tables, stats
        fig = self.script.plot_overlaid_correlations(analyses)
        ax = fig.axes[0]
        self.assertEqual(len(ax.collections), 2)
        np.testing.assert_allclose(ax.collections[1].get_facecolors()[0], to_rgba("#C43C39"))
        np.testing.assert_allclose(ax.collections[0].get_offsets()[:, 0], [1, 2, 3])
        self.assertEqual(ax.collections[1].get_offsets()[0, 0], 1)
        self.assertGreater(ax.get_xlim()[1], 3)
        self.assertEqual([t.get_text() for t in ax.get_legend().get_texts()], ["Positive", "Negative"])
        self.assertTrue(any("Positive: \u03c1 = -1.00" in t.get_text() for t in ax.texts))
        self.assertTrue(any("Negative: \u03c1 undefined" in t.get_text() for t in ax.texts))

    def test_main_exports_figure_and_calculation_tables(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pd.concat([self.source().assign(Polarity="positive"),
                       self.source().iloc[3:4].assign(Polarity="negative")]).to_csv(
                root / "input.csv", index=False
            )
            with contextlib.redirect_stdout(io.StringIO()):
                self.script.main(["--input", str(root / "input.csv"), "--no-show"])
            for polarity in ["Positive", "Negative"]:
                basename = f"SI_Class_Sharing_Correlation_{polarity}"
                for suffix in ["_Classes.csv", "_Statistics.csv"]:
                    self.assertGreater((root / f"{basename}{suffix}").stat().st_size, 0)
                for suffix in ["_Classes.csv"]:
                    self.assertGreater((root / f"{basename}_ByDataset{suffix}").stat().st_size, 0)
            for suffix in ["", "_ByDataset"]:
                for extension in ["png", "pdf"]:
                    self.assertGreater((root / f"SI_Class_Sharing_Correlation{suffix}.{extension}").stat().st_size, 0)
            stats = pd.read_csv(root / "SI_Class_Sharing_Correlation_Statistics.csv")
            # All three pooled class row counts are 3: correlation is undefined.
            self.assertTrue(pd.isna(stats.at[0, "Spearman rho"]))
            self.assertIn("constant", stats.at[0, "Status"])
            pooled = pd.read_csv(root / "SI_Class_Sharing_Correlation_Positive_Classes.csv").set_index("Lipid Class")
            self.assertEqual(pooled.at["Cer", "Pooled accepted ID rows"], 3)
            self.assertEqual(pooled.at["Cer", "Pooled unique lipids"], 1)
            self.assertEqual(pooled.at["Cer", "Shared lipids (%)"], 100)
            self.assertEqual(stats["Dataset"].tolist(), ["Pooled", *self.script.datasets] * 2)
            self.assertEqual(stats["Polarity"].tolist(), ["positive"] * 6 + ["negative"] * 6)
            by_dataset = pd.read_csv(root / "SI_Class_Sharing_Correlation_Positive_ByDataset_Classes.csv")
            pc = by_dataset.loc[by_dataset["Dataset"].eq("OTIT") & by_dataset["Lipid Class"].eq("PC")]
            self.assertEqual(pc.iloc[0]["Shared lipids (%)"], 100)
            self.assertEqual(pc.iloc[0]["Dataset accepted ID rows"], 1)
            self.assertEqual(self.script.plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

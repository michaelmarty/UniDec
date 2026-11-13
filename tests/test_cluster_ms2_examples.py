"""Regression checks for the standalone Annika MS2 example selector."""
import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestClusterMS2Examples(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location("cluster_ms2", directory / "MS2_scorecomparisons.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]):
            spec.loader.exec_module(cls.script)

    def source(self):
        return pd.DataFrame([
            {"Dataset": "OTOT", "Ontology": "PC", "Polarity": "positive",
             "Adduct type": "[M+H]+", "Metabolite name": name,
             "Overall Match Status": True, "Average Mz": 800.,
             "MS/MS spectrum": spectrum, "Class Strength": strength}
            for name, spectrum, strength in [
                ("PC A", "100:10 200:20", 9.),
                ("PC B", "100:20 200:10", 8.),
                ("PC C", "350:20 450:10", 0.)]
        ])

    def test_real_network_selects_distinct_center_and_disconnected_outlier(self):
        scored, _ = self.script.recompute_adduct_networks(self.source())
        examples, audit = self.script.select_examples(scored)
        self.assertEqual(audit.iloc[0]["Status"], "Selected")
        self.assertEqual(examples["Metabolite name"].tolist(), ["PC A", "PC C"])
        np.testing.assert_allclose(examples[self.script.NORMALIZED], [1., 0.])
        figure = self.script.plot_pair(examples)
        figure.canvas.draw()
        self.assertEqual(len(figure.axes), 2)
        self.script.plt.close(figure)

    def test_lowest_selected_without_cutoff_including_ties_and_zero_network(self):
        for strengths in [[9., 8., 7.], [0., 0., 0.]]:
            source = self.source()
            source["Class Strength"] = strengths
            examples, audit = self.script.select_examples(source)
            self.assertEqual(len(examples), 2)
            self.assertEqual(audit.iloc[0]["Status"], "Selected")
            self.assertEqual(examples.iloc[0]["Metabolite name"], "PC A")
            self.assertNotEqual(examples.iloc[0]["Lipid ID"], examples.iloc[1]["Lipid ID"])
            self.assertEqual(examples.iloc[1]["Class Strength"], min(strengths))

    def test_population_median_and_representative_feature(self):
        source = self.source()
        source = pd.concat([source, source.iloc[[0]].assign(**{
            "Metabolite name": "PC A|duplicate", "Class Strength": 7.}),
            source.iloc[[1]].assign(**{"Metabolite name": "PC excluded",
                                      "Class Strength": 20., "Overall Match Status": False}),
            source.iloc[[2]].assign(**{"Metabolite name": "PC provisional",
                                      "Class Strength": 100., "Validation Status": "provisional"})],
                           ignore_index=True)
        examples, _ = self.script.select_examples(source)
        central = examples.iloc[0]
        self.assertEqual(central["Metabolite name"], "PC A")
        self.assertEqual(central["Lipid median strength"], 8.)
        self.assertEqual(central["Input Row"], 0)
        examples, _ = self.script.select_examples(source, all_ids=True)
        self.assertEqual(examples.iloc[0]["Metabolite name"], "PC excluded")

    def test_invalid_spectra(self):
        for value in ["", "bad", "100:nan", "100:0", "100:-1", "inf:10", "-1:10"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.script.parse_spectrum(value)

    def network_source(self, lipid_class="PC", adduct="[M+H]+", count=10):
        return pd.concat([
            self.source().iloc[[0]].assign(**{
                "Dataset": dataset, "Ontology": lipid_class, "Adduct type": adduct,
                "Metabolite name": f"{lipid_class} {index}", "Class Strength": float(index)})
            for dataset in self.script.datasets for index in range(count)
        ], ignore_index=True)

    def test_all_network_classes_keep_classes_and_adducts_separate(self):
        source = pd.concat([self.network_source(),
                            self.network_source("PC", "[M+Na]+"),
                            self.network_source("LPC")], ignore_index=True)
        examples, audit = self.script.select_examples(source, network_classes=True)
        self.assertEqual(len(audit), 3 * len(self.script.datasets))
        self.assertEqual(len(examples), 2 * len(audit))
        for _, pair in examples.groupby(self.script.NETWORK_GROUPS):
            self.assertEqual(pair["Role"].tolist(), ["Central", "Outlier candidate"])
            np.testing.assert_allclose(pair[self.script.NORMALIZED], [1., 0.])

    def test_network_eligibility_matches_count_coverage_and_zero_rules(self):
        for source in [self.network_source(count=9),
                       self.network_source().loc[lambda x: x["Dataset"].ne("OTOT")],
                       self.network_source().assign(**{"Class Strength": 0.})]:
            with self.subTest(rows=len(source)):
                examples, _ = self.script.select_examples(source, network_classes=True)
                self.assertTrue(examples.empty)
                # Individual class plots still allow small or edgeless populations.
                examples, _ = self.script.select_examples(source)
                self.assertFalse(examples.empty)

    def test_network_threshold_only_needed_in_one_dataset(self):
        source = self.network_source()
        source = source.loc[source["Dataset"].eq("OTOT")
                            | source["Metabolite name"].isin(["PC 0", "PC 1"])]
        examples, _ = self.script.select_examples(source, network_classes=True)
        self.assertEqual(len(examples), 2 * len(self.script.datasets))

    def test_default_command_exports_separate_class_adduct_figures(self):
        source = pd.concat([self.network_source(), self.network_source("LPC"),
                            self.network_source("PC", "[M+Na]+")], ignore_index=True)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectra.csv"
            source.to_csv(path, index=False)
            with patch.object(self.script, "recompute_adduct_networks",
                              return_value=(source, pd.DataFrame())), \
                    patch.object(self.script, "plot_datasets") as plot, \
                    patch.object(self.script.plt, "close"), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.script.main(["--input", str(path), "--no-show"])
            self.assertEqual(plot.call_count, 3)
            self.assertEqual(plot.return_value.savefig.call_count, 6)
            for call in plot.call_args_list:
                group = call.args[0]
                self.assertEqual(len(group[["Ontology", "Adduct type", "Polarity"]].drop_duplicates()), 1)
            paths = [call.args[0] for call in plot.return_value.savefig.call_args_list]
            self.assertEqual(len(set(paths)), 6)

    def test_five_dataset_page_has_ordered_panels_and_no_clipped_labels(self):
        examples, _ = self.script.select_examples(self.source())
        combined = pd.concat([examples.assign(Dataset=dataset)
                              for dataset in self.script.datasets], ignore_index=True)
        for order in [self.script.datasets, self.script.datasets[:1]]:
            figure = self.script.plot_datasets(combined, order, "PC | [M+H]+ | positive")
            figure.canvas.draw()
            self.assertEqual(len(figure.axes), 2 * len(order))
            np.testing.assert_allclose(figure.get_size_inches(), [6.5, 7.5])
            renderer = figure.canvas.get_renderer()
            for index in range(len(order)):
                label = figure.axes[2 * index].texts[0]
                self.assertTrue(label.get_text().startswith(chr(65 + index) + "  "))
                bounds = label.get_window_extent(renderer)
                self.assertGreaterEqual(bounds.x0, 0)
                self.assertLessEqual(bounds.y1, figure.bbox.y1)
            self.script.plt.close(figure)


if __name__ == "__main__":
    unittest.main()

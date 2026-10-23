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


class TestNetworkDistributionAlignment(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[3] / "Scripts/Annika/DDA Analysis"
        spec = importlib.util.spec_from_file_location("aligned_networks", directory / "SI_8_9_10_11_12.py")
        cls.script = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(directory), *sys.path]), patch.dict(os.environ, {"MPLBACKEND": "Agg"}):
            spec.loader.exec_module(cls.script)

    def setUp(self):
        minimum = patch.object(self.script, "MIN_LIPIDS", 3)
        minimum.start()
        self.addCleanup(minimum.stop)

    def tearDown(self):
        self.script.plt.close("all")

    def source(self):
        rows = []
        for dataset, tolerance in [("OTOT", 0.05), ("Stellar", 0.7)]:
            for name, peak, accepted, validation in [
                ("PC 30:0|first", "100:100", True, "Accepted"),
                ("PC 30:0|second", "100.4:100", True, "Accepted"),
                ("PC 32:0", "100:100", True, "Accepted"),
                ("PC 34:0", "200:100", True, "Accepted"),
                ("PC 36:0", "100:100", False, "Failed"),
                ("PC 38:0", "100:100", True, " Provisional "),
                ("PC 40:0", "100:100", None, "Unknown"),
            ]:
                rows.append({"Dataset": dataset, "Ontology": "PC", "Polarity": "Positive",
                             "Adduct type": "[M+H]+", "Metabolite name": name, "MS/MS spectrum": peak,
                             "Average Mz": 500.0, "Overall Match Status": accepted,
                             "Validation Status": validation, "Tolerance": tolerance})
        return pd.DataFrame(rows)

    def context(self, source, confirmed_only=True):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.script.prepare_network_context(source, confirmed_only)

    def test_colors_exactly_match_distribution_lipid_medians_and_maxima(self):
        source = self.source()
        original = source.copy(deep=True)
        context = self.context(source)
        scored, points, raw_graphs = context
        expected = self.script.prepare_strength(scored, all_ids=False)
        pd.testing.assert_frame_equal(points, expected)
        networks = self.script.build_class_networks(source, "PC", ["OTOT", "Stellar"],
                                                     adduct="[M+H]+", context=context)
        self.assertEqual(len(networks[0][1]), 2)
        for graph in networks[0][1]:
            dataset = graph.graph["Dataset"]
            reference = expected.loc[expected["Dataset"].eq(dataset)].set_index("Lipid ID")
            self.assertEqual(graph.graph["Eligible Lipids"], 3)
            self.assertEqual(graph.number_of_nodes(), 4)  # Two features share one lipid median/color.
            for node, attrs in graph.nodes(data=True):
                point = reference.loc[attrs["Lipid ID"]]
                self.assertEqual(attrs[self.script.NORMALIZED], point[self.script.NORMALIZED])
                self.assertEqual(attrs["Lipid Median Class Strength"], point["Class Strength"])
                self.assertEqual(attrs["Normalization maximum"], point["Normalization maximum"])
            self.assertEqual(max(a[self.script.NORMALIZED] for _, a in graph.nodes(data=True)), 1)
        # Narrow-tolerance duplicate features have different original strengths,
        # but the same median-normalized color as the distribution point.
        graph = networks[0][1][0]
        duplicate = [attrs for _, attrs in graph.nodes(data=True) if attrs["Lipid ID"] == "pc|pc 30:0"]
        self.assertNotEqual(duplicate[0]["strength"], duplicate[1]["strength"])
        self.assertEqual(duplicate[0][self.script.NORMALIZED], duplicate[1][self.script.NORMALIZED])
        # All seven feature rows contributed to the original network, including
        # excluded/provisional/unknown rows. Only eligible accepted rows are displayed.
        self.assertTrue(all(graph.number_of_nodes() == 7 for graph in raw_graphs.values()))
        self.assertTrue(any(a["strength"] > graph.degree(n, weight="weight")
                            for n, a in graph.nodes(data=True)))
        pd.testing.assert_frame_equal(source, original)

    def test_adduct_and_original_class_maxima_do_not_mix(self):
        source = self.source()
        extra = source.copy()
        extra["Adduct type"] = "[M+Na]+"
        extra["MS/MS spectrum"] = "200:100"
        ether = source.copy()
        ether["Ontology"] = "EtherPC"  # Must remain a separate class in both figures.
        ether["MS/MS spectrum"] = "300:100"
        lyso = source.copy()
        lyso["Ontology"] = "LPC"
        source = pd.concat([source, extra, ether, lyso], ignore_index=True)
        context = self.context(source)
        points = context[1]
        graphs = self.script.build_class_networks(source, "PC", ["OTOT", "Stellar"],
                                                   adduct="[M+H]+", context=context)[0][1]
        for graph in graphs:
            self.assertEqual(graph.graph["Eligible Lipids"], 3)
            for _, attrs in graph.nodes(data=True):
                self.assertEqual(attrs["Adduct type"], "[M+H]+")
                self.assertEqual(attrs["Class"], "PC")
                expected = points.loc[points["Dataset"].eq(graph.graph["Dataset"])
                                      & points["Adduct type"].eq("[M+H]+")
                                      & points["Lipid ID"].eq(attrs["Lipid ID"])].iloc[0]
                self.assertEqual(attrs[self.script.NORMALIZED], expected[self.script.NORMALIZED])
            for a, b in graph.edges():
                self.assertEqual(graph.nodes[a]["Class"], graph.nodes[b]["Class"])
        for lipid_class in ["LPC", "EtherPC"]:
            own_graphs = self.script.build_class_networks(source, lipid_class, ["OTOT", "Stellar"],
                                                           adduct="[M+H]+", context=context)[0][1]
            for graph in own_graphs:
                self.assertEqual({attrs["Class"] for _, attrs in graph.nodes(data=True)}, {lipid_class})

    def test_threshold_counts_unique_eligible_lipids_after_selection(self):
        source = self.source()
        source.loc[source["Metabolite name"].eq("PC 34:0"),
                   "Overall Match Status"] = False
        # Other classes cannot supply the missing PC identity to satisfy the cutoff.
        source = pd.concat([source, self.source().assign(Ontology="LPC"),
                            self.source().assign(Ontology="EtherPC")], ignore_index=True)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.script.build_class_networks(source, "PC", ["OTOT", "Stellar"],
                                                               adduct="[M+H]+"), [])

    def test_all_phospholipid_families_select_only_the_requested_ontology(self):
        classes = [prefix + parent for parent in ["PC", "PE", "PG", "PI", "PS", "PA"]
                   for prefix in ["", "L", "Ether", "EtherL", "Ox"]]
        source = pd.concat([self.source().assign(Ontology=lipid_class) for lipid_class in classes],
                           ignore_index=True)
        context = self.context(source)
        for lipid_class in classes:
            with self.subTest(lipid_class=lipid_class):
                graphs = self.script.build_class_networks(source, lipid_class, ["OTOT", "Stellar"],
                                                           adduct="[M+H]+", context=context)[0][1]
                for graph in graphs:
                    self.assertEqual(graph.graph["Eligible Lipids"], 3)
                    self.assertEqual({attrs["Class"] for _, attrs in graph.nodes(data=True)}, {lipid_class})

    def test_zero_strength_groups_are_not_plotted(self):
        source = self.source().iloc[:4].copy()
        source["Metabolite name"] = ["PC 30:0", "PC 32:0", "PC 34:0", "PC 36:0"]
        source["MS/MS spectrum"] = ["100:100", "200:100", "300:100", "400:100"]
        context = self.context(source)
        self.assertTrue(context[1][self.script.NORMALIZED].eq(0).all())
        self.assertEqual(self.script.build_class_networks(source, "PC", ["OTOT"], context=context), [])

    def test_distribution_cutoff_allows_smaller_datasets_but_rejects_zero(self):
        self.assertEqual(self.script.DISTRIBUTION_MIN_LIPIDS, 10)
        for last_count in [10, 9, 3, 0]:
            rows = []
            for dataset in self.script.DATASETS:
                for index in range(last_count if dataset == self.script.DATASETS[-1] else 10):
                    rows.append({"Dataset": dataset, "Ontology": "PG", "Polarity": "negative",
                                 "Adduct type": "[M-H]-", "Metabolite name": f"PG {index}:0",
                                 "MS/MS spectrum": "100:100", "Average Mz": 500.0,
                                 "Overall Match Status": True})
            with self.subTest(last_count=last_count), \
                    patch.object(self.script, "MIN_LIPIDS", self.script.DISTRIBUTION_MIN_LIPIDS), \
                    contextlib.redirect_stdout(io.StringIO()):
                groups = self.script.build_class_networks(pd.DataFrame(rows), "PG", self.script.DATASETS,
                                                           adduct="[M-H]-")
                if last_count > 0:
                    self.assertEqual(len(groups[0][1]), 5)
                    self.assertEqual([graph.graph["Eligible Lipids"] for graph in groups[0][1]],
                                     [10, 10, 10, 10, last_count])
                else:
                    self.assertEqual(groups, [])

    def test_all_classes_exports_exactly_the_shared_filter_selections(self):
        source = pd.concat([self.source(), self.source().assign(Ontology="PS"),
                            self.source().loc[lambda f: f["Dataset"].eq("OTOT")].assign(Ontology="PA")],
                           ignore_index=True)
        with TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()), \
                patch.object(matplotlib.figure.Figure, "savefig") as save:
            self.script.plot_separate_networks(source, [], ["OTOT", "Stellar"],
                                                output_dir=directory, all_classes=True, close_figures=True)
            names = [call.args[0].name for call in save.call_args_list]
            self.assertTrue(all("_PC_" in name or "_PS_" in name for name in names))
            self.assertEqual(list(Path(directory).iterdir()), [])
            self.assertEqual(save.call_count, 4)
            self.assertEqual(self.script.plt.get_fignums(), [])

    def test_figures_use_fixed_0_1_scale_and_csv_input(self):
        source = self.source()
        with TemporaryDirectory() as directory, patch.object(self.script, "NETWORKS", [("PC", "[M+H]+")]), \
                patch.object(self.script.pd, "read_csv", return_value=source) as read, \
                patch.object(matplotlib.figure.Figure, "savefig") as save, \
                patch.object(self.script.plt, "show"), contextlib.redirect_stdout(io.StringIO()):
            self.script.main(["--input", str(Path(directory) / "scores.csv"), "--datasets", "OTOT", "Stellar"])
            read.assert_called_once()
            self.assertEqual(save.call_count, 2)
        figure = self.script.plt.figure(self.script.plt.get_fignums()[0])
        norms = [figure.axes[-1]._colorbar.norm]
        for axis in figure.axes[:2]:
            nodes = next(c for c in axis.collections if isinstance(c, matplotlib.collections.PathCollection))
            norms.append(nodes.norm)
            self.assertTrue(np.all((nodes.get_array() >= 0) & (nodes.get_array() <= 1)))
        self.assertTrue(all(norm is norms[0] for norm in norms))
        self.assertEqual((norms[0].vmin, norms[0].vmax), (0, 1))


if __name__ == "__main__":
    unittest.main()

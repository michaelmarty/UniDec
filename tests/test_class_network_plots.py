import contextlib
import importlib.util
import io
from pathlib import Path
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestClassNetworkPlots(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script_path = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                       / "DDA Analysis" / "SI_8_9.py")
        spec = importlib.util.spec_from_file_location("class_network_plots", script_path)
        cls.script = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(cls.script)
        except ModuleNotFoundError as exc:
            raise unittest.SkipTest(f"Lipid network dependency unavailable: {exc.name}") from exc
        from PublicScripts.Lipids import LipidFunctions
        cls.lipids = LipidFunctions

    def tearDown(self):
        self.script.plt.close("all")

    def setUp(self):
        # Small fixtures exercise plotting; threshold-specific tests use 15 below.
        minimum = patch.object(self.script, "MIN_LIPIDS", 3)
        minimum.start()
        self.addCleanup(minimum.stop)

    def frame(self, dataset="Stellar", polarity="Positive"):
        return pd.DataFrame({
            "Dataset": [dataset] * 3,
            "Polarity": [polarity] * 3,
            "Ontology": ["PC"] * 3,
            "Metabolite name": ["PC 30:0", "PC 32:0", "PC 34:0"],
            "Adduct type": ["[M+H]+"] * 3,
            "Average Mz": [500.0] * 3,
            "MS/MS spectrum": ["100:100", "100.4:100", "200:100"],
        })

    def plot(self, frame, datasets, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.script.plot_class_networks(frame, "PC", datasets, "strength", **kwargs)

    def test_tolerance_controls_fragment_and_neutral_loss_matches(self):
        for second_mz, second_peak, expected in [
            (500.0, "100.4:100", 1.0),  # Fragment and neutral loss match.
            (510.0, "100.4:100", 0.5),  # Only the fragment matches.
            (510.0, "110.4:100", 0.5),  # Only the neutral loss matches.
        ]:
            with self.subTest(precursor=second_mz, peak=second_peak):
                frame = self.frame().iloc[:2].copy()
                frame.loc[1, "Average Mz"] = second_mz
                frame.loc[1, "MS/MS spectrum"] = second_peak
                narrow = self.lipids.create_msms_score_matrix(frame, mz_tol=0.05)
                broad = self.lipids.create_msms_score_matrix(frame, mz_tol=0.7)
                self.assertAlmostEqual(narrow[0, 1], 0.0)
                self.assertAlmostEqual(broad[0, 1], expected)
                np.testing.assert_allclose(broad, broad.T)
                np.testing.assert_allclose(np.diag(broad), 1.0)

    def test_default_tolerance_and_exact_matches(self):
        frame = self.frame().iloc[:2].copy()
        frame.loc[1, "MS/MS spectrum"] = "100.04:100"
        self.assertAlmostEqual(self.lipids.create_msms_score_matrix(frame)[0, 1], 0.0)
        self.assertAlmostEqual(
            self.lipids.create_msms_score_matrix(frame, mz_tol=0.05)[0, 1], 1.0
        )
        frame.loc[1, "MS/MS spectrum"] = "100:100"
        self.assertAlmostEqual(self.lipids.create_msms_score_matrix(frame)[0, 1], 1.0)

    def test_workbook_tolerance_reaches_each_dataset_and_polarity(self):
        frames = []
        for polarity in ["Positive", "Negative"]:
            for dataset, tolerance in [("OTOT", "0.05"), ("Stellar", "0.7")]:
                frame = self.frame(dataset, polarity)
                frame["Tolerance"] = tolerance
                frames.append(frame)
        source = pd.concat(frames, ignore_index=True)
        original = source.copy(deep=True)
        figures = self.plot(source, ["OTOT", "Stellar"])
        self.assertEqual(len(figures), 2)
        for fig in figures:
            self.assertIn("3 lipids, 0 links", fig.axes[0].get_title())
            self.assertIn("3 lipids, 1 links", fig.axes[1].get_title())
        pd.testing.assert_frame_equal(source, original)

    def test_workbooks_without_tolerance_keep_class_default(self):
        frame = self.frame()
        frame.loc[1, "MS/MS spectrum"] = "100.08:100"
        figures = self.plot(frame, ["Stellar"])
        self.assertIn("3 lipids, 0 links", figures[0].axes[0].get_title())

    def test_adduct_filter_runs_before_network_and_tolerance_selection(self):
        protonated = self.frame()
        protonated["Tolerance"] = 0.05
        sodium = self.frame()
        sodium["Adduct type"] = "[M+Na]+"
        sodium["Tolerance"] = 0.7
        source = pd.concat([protonated, sodium], ignore_index=True)
        original = source.copy(deep=True)
        figures = self.plot(source, ["Stellar"], adduct=" [M+H]+ ")
        self.assertIn("3 lipids, 0 links", figures[0].axes[0].get_title())
        self.assertIn("[M+H]+", figures[0]._suptitle.get_text())
        pd.testing.assert_frame_equal(source, original)

    def test_missing_adduct_lists_available_values(self):
        with self.assertRaisesRegex(ValueError, r"Available adducts: \[M\+H\]\+"):
            self.plot(self.frame(), ["Stellar"], adduct="[M+Na]+")

    def test_abcd_selections_order_and_shared_color_scale(self):
        panels = [("PC", "[M+H]+", "Positive"), ("PC", "[M+Na]+", "Positive"),
                  ("PE", "[M-H]-", "Negative"), ("PI", "[M-H]-", "Negative")]
        frames = []
        for index, (lipid_class, adduct, polarity) in enumerate(panels):
            frame = self.frame(polarity=polarity)
            frame["Ontology"] = lipid_class
            frame["Adduct type"] = adduct
            frame["Tolerance"] = 0.05 if index % 2 == 0 else 0.7
            frames.append(frame)
        # Other datasets must not enter the four networks.
        frames.append(self.frame(dataset="OTOT"))
        source = pd.concat(frames, ignore_index=True)
        original = source.copy(deep=True)
        with contextlib.redirect_stdout(io.StringIO()):
            fig = self.script.plot_abcd_networks(source, [(c, a) for c, a, _ in panels], "Stellar")
        self.assertEqual(len(fig.axes), 5)  # Four panels and a shared colorbar.
        nodes = []
        for index, (ax, (lipid_class, adduct, polarity)) in enumerate(zip(fig.axes[:4], panels)):
            self.assertIn(f"{lipid_class} | {adduct} | {polarity}", ax.get_title())
            links = index % 2
            self.assertIn(f"3 lipids, {links} links", ax.get_title())
            self.assertIn("ABCD"[index], [text.get_text() for text in ax.texts])
            nodes.append(next(c for c in ax.collections
                              if isinstance(c, matplotlib.collections.PathCollection)))
        self.assertTrue(all(node.norm is nodes[0].norm for node in nodes))
        pd.testing.assert_frame_equal(source, original)

    def test_abcd_requires_four_panels_and_identifies_missing_panel(self):
        with self.assertRaisesRegex(ValueError, "exactly four"):
            self.script.plot_abcd_networks(self.frame(), [], "Stellar")
        panels = [("PC", "[M+H]+")] * 3 + [("SM", "[M+H]+")]
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, "Panel D"):
                self.script.plot_abcd_networks(self.frame(), panels, "Stellar")

    def test_main_abcd_plots_all_datasets_in_one_grid_with_separate_tolerances(self):
        frames = []
        for dataset, tolerance in [("OTOT", 0.05), ("Stellar", 0.7)]:
            frame = self.frame(dataset)
            frame["Tolerance"] = tolerance
            frames.append(frame)
        source = pd.concat(frames, ignore_index=True)
        with patch.object(self.script, "PANELS", [("PC", "[M+H]+")] * 4), \
                patch.object(self.script, "DATASETS", ["OTOT", "Stellar"]), \
                patch.object(self.script.pd, "read_excel", return_value=source), \
                patch.object(matplotlib.figure.Figure, "savefig") as save, \
                patch.object(self.script.plt, "show"), \
                contextlib.redirect_stdout(io.StringIO()):
            self.script.main(["--layout", "abcd"])
        figures = [self.script.plt.figure(n) for n in self.script.plt.get_fignums()]
        self.assertEqual(len(figures), 1)
        fig = figures[0]
        self.assertEqual(len(fig.axes), 13)  # Four headers, eight networks, one colorbar.
        headers = [ax for ax in fig.axes if any(t.get_text() in "ABCD" for t in ax.texts)]
        self.assertEqual(len(headers), 4)
        for index, ax in enumerate(headers):
            self.assertIn("ABCD"[index], [text.get_text() for text in ax.texts])
        fig.canvas.draw()
        positions = [ax.get_position() for ax in headers]
        self.assertLess(positions[0].x0, positions[1].x0)
        self.assertGreater(positions[0].y0, positions[2].y0)
        self.assertAlmostEqual(positions[0].y0, positions[1].y0)
        self.assertAlmostEqual(positions[2].y0, positions[3].y0)
        network_axes = [ax for ax in fig.axes if ax.get_title()]
        self.assertEqual(len(network_axes), 8)
        norms = []
        for index, ax in enumerate(network_axes):
            links = index % 2
            self.assertIn(f"3 lipids, {links} links", ax.get_title())
            self.assertTrue(ax.get_title().startswith(["OT/OT", "Stellar"][index % 2]))
            norms.append(next(c.norm for c in ax.collections
                              if isinstance(c, matplotlib.collections.PathCollection)))
        self.assertTrue(all(norm is norms[0] for norm in norms))
        self.assertEqual([call.args[0].name for call in save.call_args_list], [
            "Class_Network_ABCD_All_Datasets.png", "Class_Network_ABCD_All_Datasets.pdf",
        ])

    def test_default_layout_shows_all_five_instruments_for_each_pair(self):
        datasets = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]
        selections = [("PC", "[M+H]+"), ("PC", "[M+Na]+")]
        frames = []
        for _, adduct in selections:
            for index, dataset in enumerate(datasets):
                frame = self.frame(dataset)
                frame["Adduct type"] = adduct
                frame["Tolerance"] = 0.05 if index == 0 else 0.7
                frames.append(frame)
        source = pd.concat(frames, ignore_index=True)
        with patch.object(self.script, "NETWORKS", selections), \
                patch.object(self.script, "DATASETS", datasets), \
                patch.object(self.script.pd, "read_excel", return_value=source), \
                patch.object(matplotlib.figure.Figure, "savefig") as save, \
                patch.object(self.script.plt, "show"), \
                contextlib.redirect_stdout(io.StringIO()):
            self.script.main([])
        figures = [self.script.plt.figure(n) for n in self.script.plt.get_fignums()]
        self.assertEqual(len(figures), 2)
        for fig, (_, adduct) in zip(figures, selections):
            self.assertIn(adduct, fig._suptitle.get_text())
            norms = []
            for index, (ax, dataset) in enumerate(zip(fig.axes[:5], datasets)):
                self.assertIn(self.script.DATASET_DISPLAY_NAMES.get(dataset, dataset), ax.get_title())
                self.assertIn(f"3 lipids, {0 if index == 0 else 1} links", ax.get_title())
                norms.append(next(c.norm for c in ax.collections
                                  if isinstance(c, matplotlib.collections.PathCollection)))
            self.assertTrue(all(norm is norms[0] for norm in norms))
        self.assertEqual([call.args[0].name for call in save.call_args_list], [
            "Class_Network_PC_Positive_[M+H]+.png", "Class_Network_PC_Positive_[M+H]+.pdf",
            "Class_Network_PC_Positive_[M+Na]+.png", "Class_Network_PC_Positive_[M+Na]+.pdf",
        ])

    def test_color_scale_is_shared_across_different_class_figures_and_colorbars(self):
        frames = []
        for dataset in ["OTOT", "Stellar"]:
            connected = self.frame(dataset)
            connected["MS/MS spectrum"] = "100:100"
            isolated = self.frame(dataset)
            isolated["Ontology"] = "PE"
            isolated["MS/MS spectrum"] = ["100:100", "200:100", "300:100"]
            frames.extend([connected, isolated])
        with contextlib.redirect_stdout(io.StringIO()):
            figures = self.script.plot_separate_networks(
                pd.concat(frames, ignore_index=True), [("PC", "[M+H]+"), ("PE", "[M+H]+")],
                ["OTOT", "Stellar"])
        norms = []
        for fig in figures:
            for ax in fig.axes[:2]:
                nodes = next(c for c in ax.collections
                             if isinstance(c, matplotlib.collections.PathCollection))
                norms.append(nodes.norm)
            norms.append(fig.axes[-1]._colorbar.norm)
        self.assertTrue(all(norm is norms[0] for norm in norms))
        self.assertEqual((norms[0].vmin, norms[0].vmax), (0, 2))
        self.assertEqual(norms[0](1), 0.5)

    def sized_frame(self, dataset, count, lipid_class="PC"):
        frame = pd.concat([self.frame(dataset)] * max(1, (count + 2) // 3), ignore_index=True).iloc[:count].copy()
        frame["Ontology"] = lipid_class
        frame["Metabolite name"] = [f"{lipid_class} {index}:0" for index in range(count)]
        frame["MS/MS spectrum"] = "100:100"
        return frame

    def test_minimum_15_in_every_dataset_includes_boundary_and_excludes_small_or_missing(self):
        with patch.object(self.script, "MIN_LIPIDS", 15):
            for second_count in [15, 14, 0]:
                with self.subTest(second_count=second_count), contextlib.redirect_stdout(io.StringIO()):
                    source = pd.concat([self.sized_frame("OTOT", 15),
                                        self.sized_frame("Stellar", second_count)], ignore_index=True)
                    with patch.object(self.script, "class_network_analysis",
                                      wraps=self.script.class_network_analysis) as build:
                        networks = self.script.build_class_networks(source, "PC", ["OTOT", "Stellar"],
                                                                    adduct="[M+H]+")
                    if second_count == 15:
                        self.assertEqual([g.number_of_nodes() for g in networks[0][1]], [15, 15])
                        self.assertEqual(build.call_count, 2)
                    else:
                        self.assertEqual(networks, [])
                        build.assert_not_called()

    def test_minimum_applies_after_identification_and_usable_spectrum_filters(self):
        source = pd.concat([self.sized_frame("OTOT", 15), self.sized_frame("Stellar", 15)],
                           ignore_index=True)
        source["Overall Match Status"] = True
        for column, value in [("Overall Match Status", False), ("MS/MS spectrum", ""), ("Average Mz", np.nan)]:
            with self.subTest(column=column), patch.object(self.script, "MIN_LIPIDS", 15), \
                    contextlib.redirect_stdout(io.StringIO()):
                invalid = source.copy()
                invalid.loc[29, column] = value
                self.assertEqual(self.script.build_class_networks(invalid, "PC", ["OTOT", "Stellar"],
                                                                  confirmed_only=True), [])

    def test_excluded_choice_does_not_affect_shared_scale_and_all_excluded_is_empty(self):
        source = pd.concat([self.sized_frame("OTOT", 15), self.sized_frame("Stellar", 15),
                            self.sized_frame("OTOT", 30, "PE"), self.sized_frame("Stellar", 14, "PE")],
                           ignore_index=True)
        with patch.object(self.script, "MIN_LIPIDS", 15), contextlib.redirect_stdout(io.StringIO()):
            figures = self.script.plot_separate_networks(source, [("PC", "[M+H]+"), ("PE", "[M+H]+")],
                                                         ["OTOT", "Stellar"])
            excluded = self.script.plot_separate_networks(source, [("PE", "[M+H]+")], ["OTOT", "Stellar"])
        self.assertEqual(len(figures), 1)
        self.assertEqual(figures[0].axes[-1]._colorbar.norm.vmax, 14.5)  # Constant-strength range padding.
        self.assertEqual(excluded, [])

    def test_invalid_or_mixed_tolerances_are_rejected(self):
        for tolerances in [
            [np.nan] * 3, [np.inf] * 3, [-0.7] * 3, [0.0] * 3,
            ["invalid"] * 3, [0.7, np.nan, 0.7], [0.05, 0.7, 0.7],
        ]:
            with self.subTest(tolerances=tolerances):
                frame = self.frame()
                frame["Tolerance"] = tolerances
                with self.assertRaisesRegex(ValueError, "Stellar, Positive.*Tolerance"):
                    self.plot(frame, ["Stellar"])

    def test_confirmed_filter_runs_before_tolerance_selection(self):
        frame = self.frame()
        frame["Tolerance"] = 0.7
        frame["Overall Match Status"] = True
        rejected = frame.iloc[[0]].copy()
        rejected["Tolerance"] = np.nan
        rejected["Overall Match Status"] = False
        figures = self.plot(pd.concat([frame, rejected], ignore_index=True),
                            ["Stellar"], confirmed_only=True)
        self.assertIn("3 lipids, 1 links", figures[0].axes[0].get_title())


if __name__ == "__main__":
    unittest.main()

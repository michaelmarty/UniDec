import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestSIMSDialScatterPlots(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                  / "DDA Paper Scripts" / "SI_MSDial_Scatter_Plots.py")
        spec = importlib.util.spec_from_file_location("si_msdial_scatter", script)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def tearDown(self):
        self.module.plt.close("all")

    def fixture(self):
        return pd.DataFrame([
            ("OTOT", "Positive", "PC", True, "Accepted", 5, 700),
            (" otot\u00a0", "Negative", "HexCer_NS", False, "Failed", 6, 800),
            ("OTIT", "Positive", "NewClass", False, "Failed", 7, 900),
            ("Apex", "Positive", "PC", True, "Accepted", 5, 700),
            ("ITIT", "Positive", "PC", "TRUE", "Provisional", 8, 1000),
            ("ITIT", "Negative", "LPC", "True", "Accepted", 9, 1734),
            ("Stellar", "Negative", "PC", True, "Accepted", "bad", 800),
            ("Stellarboth", "Positive", "PC", True, "Accepted", 4, np.inf),
        ], columns=["Dataset", "Polarity", "Ontology", "Overall Match Status",
                    "Validation Status", "Average Rt(min)", "Average Mz"])

    def test_population_and_validation_filter(self):
        original = self.fixture()
        frame = self.module.prepare_lipids(original)
        self.assertEqual(len(frame), 5)
        self.assertEqual(int(frame["Accepted"].sum()), 2)
        self.assertFalse(frame.loc[4, "Accepted"])
        self.assertEqual(frame.loc[1, "Dataset"], "OTOT")
        self.assertEqual(frame.loc[1, "Plot Class"], "HexCer")
        self.assertEqual(frame.loc[2, "Plot Class"], "Other")
        self.assertEqual(len(original), 8)
        self.assertNotIn("Plot Class", original.columns)

    def test_pre_post_points_counts_and_shared_axes(self):
        frame = self.module.prepare_lipids(self.fixture())
        for polarity, before, after in [("positive", 3, 1), ("negative", 2, 1)]:
            for filtered, expected in [(False, before), (True, after)]:
                fig = self.module.plot_scatter_panels(frame, polarity, filtered, 20, 1800)
                self.assertEqual(len(fig.axes), 6)
                self.assertEqual(sum(len(ax.collections[0].get_offsets())
                                     for ax in fig.axes[:5]), expected)
                first_count = 0 if polarity == "negative" and filtered else 1
                self.assertEqual(fig.axes[0].texts[1].get_text(), f"n = {first_count}")
                for ax in fig.axes[:5]:
                    self.assertEqual(ax.get_xlim(), (0, 20))
                    self.assertEqual(ax.get_ylim(), (0, 1800))
                if polarity == "negative" and filtered:
                    np.testing.assert_allclose(fig.axes[2].collections[0].get_offsets(),
                                               [[9, 1734]])
                self.module.plt.close(fig)

    def test_four_figure_pairs_are_saved(self):
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "figures"
            with patch.object(self.module.pd, "read_excel", return_value=self.fixture()):
                self.module.main(["--input", "input.xlsx", "--output-dir", str(output), "--no-show"])
            self.assertEqual(len(list(output.glob("*.png"))), 4)
            self.assertEqual(len(list(output.glob("*.pdf"))), 4)
            self.assertEqual(self.module.plt.get_fignums(), [])


if __name__ == "__main__":
    unittest.main()

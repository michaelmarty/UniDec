import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestSIClassIntensity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                  / "DDA Paper Figures" / "MS1_by_Class.py")
        spec = importlib.util.spec_from_file_location("si_class_intensity", script)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def tearDown(self):
        self.module.plt.close("all")

    def fixture(self):
        rng = np.random.default_rng(17)
        rows = []
        for polarity in ["Positive", "Negative"]:
            for dataset_index, dataset in enumerate(self.module.DATASETS):
                for lipid in range(6):
                    intensity = 10 ** (2 + lipid * 0.1 + dataset_index * 0.2
                                       + rng.normal(0, 0.08))
                    rows.append((dataset, "PC", polarity, True, "Accepted",
                                 f"PC {lipid}|detail", f"100:{intensity}", f"90:{intensity}"))
        return pd.DataFrame(rows, columns=["Dataset", "Ontology", "Polarity",
                                          "Overall Match Status", "Validation Status",
                                          "Metabolite name", "MS1 isotopic spectrum", "MS/MS spectrum"])

    def test_accepted_only_parent_classes_and_original_labels(self):
        raw = self.fixture().iloc[:6].copy()
        raw["Ontology"] = ["OxPI", "EtherLPC", "Cer_AP", "PC", "PC", "PC"]
        raw["Overall Match Status"] = [True, "TRUE", 1, False, True, None]
        raw.loc[4, "Validation Status"] = " Provisional "
        frame = self.module.prepare_intensity(raw)
        self.assertEqual(frame["Plot Class"].tolist(), ["PI", "PC", "Cer"])
        self.assertEqual(frame["Ontology"].tolist(), ["OxPI", "EtherLPC", "Cer_AP"])
        self.assertTrue(frame["ID Status"].eq("Accepted").all())
        self.assertNotIn("Plot Class", raw)

    def test_exports_and_plot_use_accepted_ids(self):
        raw = self.fixture()
        excluded = raw.iloc[[0]].copy()
        excluded["Overall Match Status"] = False
        raw = pd.concat([raw, excluded], ignore_index=True)
        with TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source = folder / "input.csv"
            raw.to_csv(source, index=False)
            self.module.main(["--input", str(source), "--output-dir", str(folder), "--no-show"])
            basename = "MS1_by_Class"
            points = pd.read_csv(folder / f"{basename}_Points.csv")
            self.assertEqual(len(points), 60)
            self.assertTrue(points["ID Status"].eq("Accepted").all())
            self.assertEqual(list(folder.glob("*_Statistics.xlsx")), [])
            for extension in ["png", "pdf"]:
                self.assertGreater((folder / f"{basename}.{extension}").stat().st_size, 0)
        frame = self.module.prepare_intensity(raw)
        fig = self.module.plot_intensity(frame)
        self.assertEqual(sum(len(c.get_offsets()) for ax in fig.axes for c in ax.collections), 60)
        self.assertEqual(len(fig.legends[0].texts), 5)
        self.assertFalse(any("p=" in text.get_text() or "NT" in text.get_text()
                             for ax in fig.axes for text in ax.texts))

    def test_ms2_script_defaults_to_summed_fragment_signal(self):
        raw = self.fixture()
        raw["MS1 isotopic spectrum"] = "100:999"
        raw["MS/MS spectrum"] = "50:2 75:3"
        raw.loc[0, "MS/MS spectrum"] = "bad"
        raw.loc[1, "MS/MS spectrum"] = "50:0"
        raw.loc[2, "Overall Match Status"] = False
        script = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                  / "DDA Paper Figures" / "Figure3.py")
        with TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source = folder / "input.csv"
            raw.to_csv(source, index=False)
            subprocess.run([sys.executable, str(script), "--input", str(source),
                            "--output-dir", str(folder), "--no-show"],
                           check=True, capture_output=True, text=True,
                           env={**os.environ, "MPLBACKEND": "Agg"})
            points = pd.read_csv(folder / "Figure3_Points.csv")
            self.assertEqual(len(points), 57)
            self.assertTrue(points["Signal Intensity"].eq(5).all())
            for extension in ["png", "pdf"]:
                self.assertGreater((folder / f"Figure3.{extension}").stat().st_size, 0)
            self.assertFalse((folder / "MS1_by_Class.png").exists())


if __name__ == "__main__":
    unittest.main()

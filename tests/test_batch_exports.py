import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from unidec import batch
from unidec.modules import peakstructure


class ExportRecordingEngine:
    def __init__(self):
        self.pks = peakstructure.Peaks()
        self.html_str = ""
        self.reset_config()

    def reset_config(self):
        self.config = SimpleNamespace(smashflag=0, peakwindow=999.0)

    def open_file(self, path, **kwargs):
        del kwargs
        self.config.filename = path
        self.config.reportfile = str(Path(path).with_suffix(".html"))

    def autorun(self, **kwargs):
        del kwargs
        name = Path(self.config.filename).name
        values = [[101.0, 11.0]] if name == "first.txt" else [[202.0, 22.0], [203.0, 23.0]]
        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(values)

    def gen_html_report(self, **kwargs):
        del kwargs
        name = Path(self.config.filename).name
        masses = ", ".join(str(peak.mass) for peak in self.pks.peaks)
        self.html_str = (
            f"<h1>{name}</h1>\n<p class=\"sample\">{name}</p>\n"
            f"<table><tr><th>Masses</th><td>{masses}</td></tr></table>\n"
        )
        Path(self.config.reportfile).write_text(
            f"<html>{self.html_str}</html>", encoding="utf-8"
        )
        return self.config.reportfile


class TestBatchExports(unittest.TestCase):
    def test_xlsx_peak_and_html_exports_preserve_row_and_value_associations(self):
        """Two distinct fake results expose row swaps while keeping XLSX and HTML formatting real."""
        with tempfile.TemporaryDirectory(prefix="unidec-batch-export-") as directory:
            root = Path(directory)
            for name in ("first.txt", "second.txt"):
                (root / name).write_text("100 1\n", encoding="utf-8")
            source = root / "batch.xlsx"
            pd.DataFrame(
                [
                    {
                        "Sample name": "first.txt",
                        "Data Directory": directory,
                        "Expected Value": 11.0,
                    },
                    {
                        "Sample name": "second.txt",
                        "Data Directory": directory,
                        "Expected Value": 22.0,
                    },
                ]
            ).to_excel(source, index=False)
            frame = pd.read_excel(source)

            processor = batch.UniDecBatchProcessor()
            processor.eng = ExportRecordingEngine()
            processor.filename = str(source)
            processor.top_dir = directory
            previous_directory = os.getcwd()
            try:
                result = processor.run_df(
                    frame,
                    decon=True,
                    write_html=True,
                    write_xlsx=True,
                    write_peaks=True,
                    individual_html=True,
                )
            finally:
                os.chdir(previous_directory)

            results_path = root / "batch_results.xlsx"
            peaks_path = root / "batch_peaks.xlsx"
            report_path = root / "batch_report.html"
            self.assertTrue(results_path.is_file())
            self.assertTrue(peaks_path.is_file())
            self.assertTrue(report_path.is_file())

            exported = pd.read_excel(results_path)
            self.assertEqual(exported["Sample name"].tolist(), ["first.txt", "second.txt"])
            self.assertEqual(exported["Expected Value"].tolist(), [11, 22])
            self.assertEqual(exported["NumPeaks"].tolist(), [1, 2])
            self.assertEqual([Path(path).name for path in exported["Reports"]], ["first.html", "second.html"])
            self.assertEqual(result["NumPeaks"].tolist(), [1, 2])

            peaks = pd.read_excel(peaks_path)
            self.assertEqual(peaks["Mass"].tolist(), [101, 202, 203])
            self.assertEqual(peaks["FileNumber"].tolist(), [0, 1, 1])
            self.assertEqual([Path(path).name for path in peaks["FileName"]], ["first.txt", "second.txt", "second.txt"])

            html = report_path.read_text(encoding="utf-8")
            self.assertIn("first.txt", html)
            self.assertIn("second.txt", html)
            self.assertIn("101.0", html)
            self.assertIn("202.0, 203.0", html)
            self.assertEqual(html.count("<h1>"), 1)


if __name__ == "__main__":
    unittest.main()

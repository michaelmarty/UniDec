import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd

from unidec import batch
from unidec.modules import peakstructure


class TestBatchParameters(unittest.TestCase):
    @staticmethod
    def _engine(**values):
        defaults = dict(
            smashflag=7, peakwindow=-1, peakthresh=-1, masslb=-1, massub=-1,
            minmz=-1, maxmz=-1, startz=-1, endz=-1, massbins=-1, mzsig=-1,
            psfun=-1, integratelb=-1, integrateub=-1, datanorm=-1, subbuff=-1,
            subtype=-1, mzbins=-1, linflag=-1, smooth=-1, peaknorm=-1,
            kernel="unchanged", doubledec=False, smashfile="unchanged",
            import_smashfile=Mock(),
        )
        defaults.update(values)
        return SimpleNamespace(config=SimpleNamespace(**defaults))

    def test_every_numeric_column_accepts_numeric_strings(self):
        engine = self._engine()
        row = {
            "Config Peak Window": "10.5", "Config Peak Thres": "0.2",
            "Config Low Mass": "100", "Config High Mass": "200",
            "Config Low mz": "50", "Config High m/z": "250",
            "Config Low z": "2.0", "Config High z": "20.0",
            "Config Sample Mass Every": "1.5", "Config mz Peak FWHM": "0.8",
            "Config m/z Peak Shape": "2.0", "Config Integral Lower": "90",
            "Config Integral Upper": "210", "Config Data Norm": "1.0",
            "Config Data Subbuff": "5.0", "Config Data Binevery": "3.0",
            "Config Data Lintype": "4.0", "Config Data Smoothing": "2.0",
            "Config Peak Norm": "1.0",
        }
        batch.set_param_from_row(engine, row)
        expected = {
            "peakwindow": 10.5, "peakthresh": 0.2, "masslb": 100.0, "massub": 200.0,
            "minmz": 50.0, "maxmz": 250.0, "startz": 2, "endz": 20,
            "massbins": 1.5, "mzsig": 0.8, "psfun": 2, "integratelb": 90.0,
            "integrateub": 210.0, "datanorm": 1, "subbuff": 5, "subtype": 2,
            "mzbins": 3, "linflag": 4, "smooth": 2, "peaknorm": 1,
        }
        for field, value in expected.items():
            with self.subTest(field=field):
                self.assertEqual(getattr(engine.config, field), value)
        self.assertEqual(engine.config.kernel, "unchanged")

    def test_missing_nan_and_invalid_values_leave_fields_unchanged(self):
        engine = self._engine()
        batch.set_param_from_row(engine, {})
        self.assertEqual(engine.config.peakwindow, -1)
        self.assertEqual(engine.config.smashflag, 0)

        batch.set_param_from_row(engine, {"Config Peak Window": np.nan, "Config Peak Thres": "nan"})
        self.assertEqual(engine.config.peakwindow, -1)
        self.assertEqual(engine.config.peakthresh, -1)

        batch.set_param_from_row(engine, {"Config Peak Window": "bad", "Config Data Norm": "bad"})
        self.assertEqual(engine.config.peakwindow, -1)
        self.assertEqual(engine.config.datanorm, -1)

    def test_invalid_integral_limits_reset_to_unspecified(self):
        engine = self._engine(integratelb=10, integrateub=20)
        batch.set_param_from_row(
            engine,
            {"Config Integral Lower": "bad", "Config Integral Upper": "bad"},
        )
        self.assertEqual(engine.config.integratelb, "")
        self.assertEqual(engine.config.integrateub, "")

    def test_kernel_and_smash_files_resolve_relative_to_batch_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            kernel = root / "kernel_mass.txt"
            np.savetxt(kernel, [[100.0, 1.0], [101.0, 2.0], [102.0, 1.0]])
            smash = root / "smash.txt"
            np.savetxt(smash, [[100.0, 101.0]])
            engine = self._engine()

            batch.set_param_from_row(
                engine,
                {"DoubleDec Kernel File": kernel.name, "Config Smash File": smash.name},
                directory,
            )

            self.assertEqual(Path(engine.config.kernel), kernel)
            self.assertTrue(engine.config.doubledec)
            self.assertEqual(Path(engine.config.smashfile), smash)
            self.assertEqual(engine.config.smashflag, 1)
            engine.config.import_smashfile.assert_called_once_with()

    def test_kernel_source_file_resolves_generated_mass_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "kernel.txt"
            source.write_text("100 1\n101 2\n", encoding="utf-8")
            generated = root / "kernel_unidecfiles" / "kernel_mass.txt"
            generated.parent.mkdir()
            generated.write_text("100 1\n101 2\n", encoding="utf-8")
            self.assertEqual(Path(batch.find_kernel_file(str(source))), generated)
            self.assertIsNone(batch.find_kernel_file(str(root / "missing_mass.txt")))


class TestBatchPathsAndTime(unittest.TestCase):
    def test_time_range_defaults_valid_values_and_invalid_values(self):
        self.assertIsNone(batch.get_time_range({}))
        self.assertEqual(batch.get_time_range({"Start Time": "1.5"}), [1.5, 1e12])
        self.assertEqual(batch.get_time_range({"End Time": "4.5"}), [0.0, 4.5])
        self.assertIsNone(batch.get_time_range({"Start Time": "bad"}))

    def test_file_paths_and_duplicate_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.txt"
            path.write_text("100 1\n", encoding="utf-8")
            processor = batch.UniDecBatchProcessor()
            processor.rundf = pd.DataFrame([
                {"Sample name": "sample.txt", "Data Directory": directory},
                {"Sample name": "sample.txt", "Data Directory": directory},
            ])
            expected = os.path.abspath(path)
            self.assertEqual(processor.get_file_path(processor.rundf.iloc[0]), expected)
            self.assertEqual(processor.check_duplicate_filenames(), [expected])


class RecordingEngine:
    def __init__(self):
        self.autorun_peakwindows = []
        self.opened = []
        self.pks = peakstructure.Peaks()
        self.reset_config()

    def reset_config(self):
        self.config = SimpleNamespace(smashflag=0, peakwindow=999.0)

    def open_file(self, path, **kwargs):
        self.opened.append((path, kwargs))

    def autorun(self, **kwargs):
        self.autorun_peakwindows.append(self.config.peakwindow)
        self.pks = peakstructure.Peaks()
        self.pks.add_peaks([[100.0, 1.0]])


class TestBatchRowOrchestration(unittest.TestCase):
    def test_rows_reset_config_keep_result_association_and_handle_missing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("first.txt", "second.txt"):
                (root / name).write_text("100 1\n", encoding="utf-8")
            frame = pd.DataFrame([
                {"Sample name": "first.txt", "Data Directory": directory, "Config Peak Window": 10.0},
                {"Sample name": "second.txt", "Data Directory": directory},
                {"Sample name": "missing.txt", "Data Directory": directory},
            ])
            processor = batch.UniDecBatchProcessor()
            engine = RecordingEngine()
            processor.eng = engine

            result = processor.run_df(
                frame, decon=True, write_html=False, write_xlsx=False,
                write_peaks=False, individual_html=False,
            )

            self.assertEqual(engine.autorun_peakwindows, [10.0, 999.0])
            self.assertEqual([Path(item[0]).name for item in engine.opened], ["first.txt", "second.txt"])
            self.assertEqual(result["NumPeaks"].tolist(), [1, 1, 0])
            self.assertEqual([p.filenumber for p in processor.pks.peaks], [0, 1])
            self.assertEqual([Path(p.filename).name for p in processor.pks.peaks], ["first.txt", "second.txt"])


if __name__ == "__main__":
    unittest.main()

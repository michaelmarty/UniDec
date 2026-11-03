import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np

from unidec.modules.ChromEng import ChromEngine
from unidec.unichrom_session import UniChromSession


UNIDEC_ROOT = Path(__file__).resolve().parents[1]


class TestUniChromSession(unittest.TestCase):
    @staticmethod
    def _importer():
        importer = Mock()
        importer.get_polarity.return_value = "Positive"
        importer.get_tic.return_value = [[0.0, 10.0], [1.0, 20.0], [2.0, 5.0]]
        importer.get_bpc.return_value = [[0.0, 4.0], [1.0, 8.0], [2.0, 3.0]]
        importer.scans = [1, 2, 3]
        importer.get_scans_from_times.side_effect = lambda values: [
            int(values[0]) + 1, int(values[1]) + 1]
        importer.get_avg_scan.return_value = np.array([[500.0, 10.0], [501.0, 5.0]])
        importer.get_times_from_scans.side_effect = lambda values: (
            float(values[0] - 1),
            float((values[0] + values[1]) / 2 - 1),
            float(values[1] - 1),
        )
        return importer

    @patch("unidec.modules.ChromEng.ImporterFactory.create_importer")
    def test_open_file_routes_derived_project_to_managed_directory(self, create_importer):
        create_importer.return_value = self._importer()
        with tempfile.TemporaryDirectory(prefix="unichrom-input-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="unichrom-work-") as working_directory:
            source = Path(input_directory, "sample.mzML")
            source.touch()
            session = UniChromSession(
                engine=ChromEngine(),
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )

            response = session.open_chromatogram(source)

            self.assertTrue(response["ok"], response)
            state = response["state"]
            self.assertEqual(state["mode"], "unichrom")
            self.assertEqual(state["chromatogram"]["scan_count"], 3)
            self.assertEqual(state["chromatogram"]["tic"]["points"], 3)
            self.assertEqual(state["chromatogram"]["bpc"]["maximum_intensity"], 8.0)
            managed = Path(session.engine.config.hdf_file)
            self.assertTrue(managed.is_relative_to(Path(working_directory)))
            self.assertTrue(managed.is_file())
            self.assertFalse(Path(input_directory, "sample.hdf5").exists())
            self.assertFalse(Path(input_directory, "sample.hdf5_unidecfiles").exists())

            parameters = session.get_parameters(("chromatogram_type", "time_window", "dtsig"))
            self.assertTrue(parameters["ok"], parameters)
            self.assertIn("unichrom", parameters["schema"]["time_window"]["modes"])
            updated = session.set_parameters({"chromatogram_type": 1, "dtsig": 0})
            self.assertTrue(updated["ok"], updated)
            rejected_parameters = session.set_parameters({"time_start": 2, "time_end": 1})
            self.assertFalse(rejected_parameters["ok"])
            self.assertEqual(rejected_parameters["error"]["code"], "invalid_parameters")
            self.assertEqual(session.engine.config.chromatogram_type, 1)

            failed = session.run()
            self.assertEqual(failed["error"]["code"], "no_windows")
            operation_id = failed["error"]["operation_id"]
            manifest = session.get_run_manifest(operation_id)["manifest"]
            self.assertEqual(manifest["workload"], {"scan_count": 3, "window_count": 0})
            exported = session.export_regression_fixture(operation_id, "fixtures/no-windows.json")
            self.assertTrue(exported["ok"], exported)
            fixture = json.loads(Path(exported["path"]).read_text(encoding="utf-8"))
            self.assertEqual(fixture["source"]["mode"], "unichrom")
            self.assertEqual(fixture["source"]["workload"]["window_count"], 0)

            windows = session.create_time_windows(width=1.0, start=0.0, end=2.0)
            self.assertTrue(windows["ok"], windows)
            self.assertEqual(windows["state"]["chromatogram"]["window_count"], 2)
            self.assertEqual(session.engine.data.v1name, "timemid")
            self.assertEqual(session.engine.data.v2name, "scanstart")
            self.assertTrue(managed.is_file())
            self.assertFalse(Path(input_directory, "sample.hdf5").exists())

            rejected = session.create_time_windows(width=0)
            self.assertFalse(rejected["ok"])
            self.assertEqual(rejected["error"]["code"], "invalid_windows")
            self.assertEqual(len(session.engine.data.spectra), 2)

    @patch("unidec.modules.ChromEng.ImporterFactory.create_importer")
    def test_open_vendor_directory_and_reject_outside_root(self, create_importer):
        create_importer.return_value = self._importer()
        with tempfile.TemporaryDirectory(prefix="unichrom-input-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="unichrom-other-") as other_directory, \
                tempfile.TemporaryDirectory(prefix="unichrom-work-") as working_directory:
            vendor_directory = Path(input_directory, "sample.raw")
            vendor_directory.mkdir()
            outside = Path(other_directory, "outside.raw")
            outside.mkdir()
            session = UniChromSession(
                engine=ChromEngine(),
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )

            opened = session.open_chromatogram(vendor_directory)
            self.assertTrue(opened["ok"], opened)
            self.assertEqual(opened["state"]["source"]["kind"], "directory")

            rejected = session.open_chromatogram(outside)
            self.assertFalse(rejected["ok"])
            self.assertEqual(rejected["error"]["code"], "path_not_allowed")

    @patch("unidec.modules.ChromEng.ImporterFactory.create_importer")
    def test_managed_windows_complete_native_workflow(self, create_importer):
        importer = self._importer()
        importer.get_avg_scan.return_value = np.loadtxt(
            UNIDEC_ROOT / "unidec" / "bin" / "Example Data" / "ADH.txt")
        create_importer.return_value = importer
        with tempfile.TemporaryDirectory(prefix="unichrom-input-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="unichrom-work-") as working_directory:
            source = Path(input_directory, "workflow.mzML")
            source.write_bytes(b"test source sentinel")
            original = source.read_bytes()
            session = UniChromSession(
                engine=ChromEngine(),
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )
            self.assertTrue(session.open_chromatogram(source)["ok"])
            self.assertTrue(session.create_time_windows(1.0, 0.0, 2.0)["ok"])
            update = session.set_parameters({
                "startz": 5,
                "endz": 50,
                "masslb": 10000,
                "massub": 250000,
                "mzbins": 1,
                "dtsig": 0,
            })
            self.assertTrue(update["ok"], update)

            response = session.run(stage="all")

            self.assertTrue(response["ok"], response)
            self.assertTrue(all(response["state"]["available"].values()))
            results = session.get_results(max_peaks=5)
            self.assertTrue(results["ok"], results)
            self.assertEqual(results["results"]["window_count"], 2)
            self.assertGreater(results["results"]["peak_count"], 0)
            self.assertLessEqual(len(results["results"]["peaks"]), 5)
            self.assertEqual(source.read_bytes(), original)
            manifest = session.get_run_manifest(response["operation_id"])["manifest"]
            self.assertEqual(manifest["workload"]["window_count"], 2)
            self.assertEqual(manifest["results"]["window_count"], 2)
            self.assertLessEqual(len(manifest["results"]["peaks"]), 25)


if __name__ == "__main__":
    unittest.main()

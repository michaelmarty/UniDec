import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

from unidec.meta_session import MetaUniDecCollectionSession
from unidec.metaunidec.mudeng import MetaUniDec

from tests._test_support import copy_workflow_spectra


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class TestMetaUniDecCollectionSession(unittest.TestCase):
    @staticmethod
    def _write_collection(path):
        with h5py.File(path, "w") as hdf:
            hdf.require_group("config").attrs["metamode"] = -1
            dataset = hdf.require_group("ms_dataset")
            dataset.attrs["v1name"] = "Collision Voltage"
            dataset.attrs["v2name"] = "Replicate"
            for index, name in enumerate(("low", "high")):
                spectrum = dataset.require_group(str(index))
                values = np.array([[100.0, index + 1.0], [101.0, index + 2.0]])
                spectrum.create_dataset("raw_data", data=values)
                spectrum.create_dataset("processed_data", data=values)
                spectrum.attrs["name"] = name
                spectrum.attrs["Collision Voltage"] = 10 * (index + 1)
                spectrum.attrs["Replicate"] = index + 1

    def test_open_uses_managed_copy_and_returns_bounded_summary(self):
        with tempfile.TemporaryDirectory(prefix="meta-session-input-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="meta-session-work-") as working_directory:
            source = Path(input_directory, "collection.hdf5")
            self._write_collection(source)
            original_digest = file_digest(source)

            session = MetaUniDecCollectionSession(
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )
            response = session.open_collection(source)

            self.assertTrue(response["ok"], response)
            state = response["state"]
            self.assertEqual(state["mode"], "metaunidec")
            self.assertEqual(state["status"], "ready")
            self.assertEqual(state["collection"]["spectrum_count"], 2)
            self.assertEqual(state["collection"]["spectra"][0]["name"], "low")
            self.assertEqual(state["collection"]["spectra"][1]["variable_1"], 20)
            self.assertFalse(state["collection"]["spectra_truncated"])

            managed_project = Path(session.engine.config.hdf_file)
            self.assertTrue(managed_project.resolve().is_relative_to(Path(working_directory).resolve()))
            self.assertNotEqual(managed_project, source)
            self.assertEqual(file_digest(source), original_digest)
            self.assertFalse(Path(str(source) + "_unidecfiles").exists())
            self.assertTrue(Path(str(managed_project) + "_unidecfiles").is_dir())

            parameters = session.get_parameters(("startz", "endz"))
            self.assertTrue(parameters["ok"], parameters)
            self.assertIn("metaunidec", parameters["schema"]["startz"]["modes"])

            original_startz = parameters["values"]["startz"]
            update = session.set_parameters({"startz": 4, "endz": 40})
            self.assertTrue(update["ok"], update)
            self.assertEqual(session.get_parameters(("startz",))["values"]["startz"], 4)
            rejected = session.set_parameters({"startz": 41, "endz": 40})
            self.assertFalse(rejected["ok"])
            self.assertEqual(rejected["error"]["code"], "invalid_parameters")
            self.assertEqual(session.get_parameters(("startz",))["values"]["startz"], 4)
            self.assertNotEqual(original_startz, 4)
            self.assertEqual(file_digest(source), original_digest)

    def test_open_rejects_a_collection_outside_allowed_roots(self):
        with tempfile.TemporaryDirectory(prefix="meta-session-allowed-") as allowed_directory, \
                tempfile.TemporaryDirectory(prefix="meta-session-other-") as other_directory, \
                tempfile.TemporaryDirectory(prefix="meta-session-work-") as working_directory:
            source = Path(other_directory, "collection.hdf5")
            self._write_collection(source)
            session = MetaUniDecCollectionSession(
                working_directory=working_directory,
                allowed_input_roots=(allowed_directory,),
            )

            response = session.open_collection(source)

            self.assertFalse(response["ok"])
            self.assertEqual(response["error"]["code"], "path_not_allowed")

    def test_managed_collection_runs_all_and_returns_bounded_results(self):
        with tempfile.TemporaryDirectory(prefix="meta-session-run-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="meta-session-work-") as working_directory:
            source = Path(input_directory, "workflow.hdf5")
            creator = MetaUniDec()
            creator.data.new_file(str(source))
            for spectrum in copy_workflow_spectra(input_directory):
                creator.data.add_file(path=str(spectrum))
            creator.config.write_hdf5(str(source))
            original_digest = file_digest(source)

            session = MetaUniDecCollectionSession(
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )
            self.assertTrue(session.open_collection(source)["ok"])
            update = session.set_parameters({
                "startz": 5,
                "endz": 50,
                "masslb": 10000,
                "massub": 250000,
                "mzbins": 1,
            })
            self.assertTrue(update["ok"], update)

            response = session.run(stage="all")

            self.assertTrue(response["ok"], response)
            self.assertTrue(all(response["state"]["available"].values()))
            results = session.get_results(max_peaks=5)
            self.assertTrue(results["ok"], results)
            self.assertEqual(results["results"]["spectrum_count"], 2)
            self.assertGreater(results["results"]["mass_points"], 0)
            self.assertGreater(results["results"]["peak_count"], 0)
            self.assertLessEqual(len(results["results"]["peaks"]), 5)
            self.assertEqual(file_digest(source), original_digest)
            manifest = session.get_run_manifest(response["operation_id"])
            self.assertTrue(manifest["ok"], manifest)
            self.assertEqual(manifest["manifest"]["mode"], "metaunidec")
            self.assertEqual(manifest["manifest"]["workload"]["spectrum_count"], 2)
            self.assertLessEqual(len(manifest["manifest"]["results"]["peaks"]), 25)

            checkpoint = session.save_checkpoint("checkpoints/workflow.hdf5")
            self.assertTrue(checkpoint["ok"], checkpoint)
            checkpoint_path = Path(checkpoint["path"])
            self.assertTrue(checkpoint_path.is_file())
            invalid_checkpoint = session.save_checkpoint("checkpoints/workflow.zip")
            self.assertFalse(invalid_checkpoint["ok"])
            self.assertEqual(invalid_checkpoint["error"]["code"], "invalid_path")
            duplicate = session.save_checkpoint("checkpoints/workflow.hdf5")
            self.assertFalse(duplicate["ok"])
            self.assertEqual(duplicate["error"]["code"], "already_exists")

            restored = MetaUniDecCollectionSession(
                working_directory=Path(working_directory, "restored"),
                allowed_input_roots=(input_directory, working_directory),
            )
            loaded = restored.load_checkpoint(checkpoint_path)
            self.assertTrue(loaded["ok"], loaded)
            self.assertEqual(loaded["state"]["collection"]["spectrum_count"], 2)
            self.assertEqual(restored.get_parameters(("startz",))["values"]["startz"], 5)
            stale = restored.get_results()
            self.assertFalse(stale["ok"])
            self.assertEqual(stale["error"]["code"], "stale_results")
            rerun = restored.run(stage="all")
            self.assertTrue(rerun["ok"], rerun)
            self.assertGreater(restored.get_results()["results"]["peak_count"], 0)
            self.assertEqual(file_digest(source), original_digest)

            invalid = session.run(stage="process")
            self.assertFalse(invalid["ok"])
            self.assertEqual(invalid["error"]["code"], "invalid_stage")

    def test_failed_run_exports_bounded_collection_fixture(self):
        with tempfile.TemporaryDirectory(prefix="meta-session-input-") as input_directory, \
                tempfile.TemporaryDirectory(prefix="meta-session-work-") as working_directory:
            source = Path(input_directory, "collection.hdf5")
            self._write_collection(source)
            session = MetaUniDecCollectionSession(
                working_directory=working_directory,
                allowed_input_roots=(input_directory,),
            )
            self.assertTrue(session.open_collection(source)["ok"])
            with patch.object(session.engine, "process_data", side_effect=RuntimeError("x" * 2100)):
                failed = session.run()

            self.assertEqual(failed["error"]["code"], "run_failed")
            operation_id = failed["error"]["operation_id"]
            manifest = session.get_run_manifest(operation_id)["manifest"]
            self.assertEqual(manifest["status"], "error")
            self.assertEqual(manifest["workload"], {"spectrum_count": 2})
            exported = session.export_regression_fixture(operation_id, "fixtures/meta.json")
            self.assertTrue(exported["ok"], exported)
            fixture = json.loads(Path(exported["path"]).read_text(encoding="utf-8"))
            self.assertEqual(fixture["source"]["mode"], "metaunidec")
            self.assertEqual(fixture["source"]["workload"]["spectrum_count"], 2)
            self.assertEqual(len(fixture["expected_error"]["message"]), 2000)
            self.assertTrue(fixture["expected_error"]["message_truncated"])
            self.assertEqual(
                session.export_regression_fixture(operation_id, "fixtures/meta.json")
                ["error"]["code"],
                "already_exists",
            )


if __name__ == "__main__":
    unittest.main()

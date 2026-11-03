import hashlib
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from unidec.meta_session import MetaUniDecCollectionSession


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
            self.assertTrue(managed_project.is_relative_to(Path(working_directory)))
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


if __name__ == "__main__":
    unittest.main()

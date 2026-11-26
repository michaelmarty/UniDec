import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

from unidec.metaunidec.mudstruct import MetaDataSet


class TestMetaDataSetIO(unittest.TestCase):
    @staticmethod
    def _engine():
        config = SimpleNamespace(dtype=np.float32, datanorm=0, masslb=0.0, massub=1000.0)
        return SimpleNamespace(config=config)

    def test_create_reopen_remove_reindex_and_variable_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "dataset.hdf5")
            dataset = MetaDataSet(self._engine())
            dataset.new_file(str(path))
            dataset.v1name = "Collision Voltage"
            dataset.v2name = "Replicate"
            dataset.add_data(np.array([[100.0, 1.0], [101.0, 2.0]]), name="first", export=False)
            dataset.add_data(np.array([[200.0, 3.0], [201.0, 4.0]]), name="second", export=False)
            dataset.spectra[0].var1, dataset.spectra[0].var2 = 10.0, 1.0
            dataset.spectra[1].var1, dataset.spectra[1].var2 = 20.0, 2.0
            dataset.export_hdf5()

            with h5py.File(path, "r") as hdf:
                self.assertEqual(hdf["ms_dataset"].attrs["num"], 2)
                self.assertEqual(sorted(hdf["ms_dataset"].keys()), ["0", "1"])
                self.assertEqual(hdf["ms_dataset/0/raw_data"].shape, (2, 2))
                self.assertEqual(hdf["ms_dataset/1"].attrs["name"], "second")

            reopened = MetaDataSet(self._engine())
            reopened.import_hdf5(str(path))
            self.assertEqual([s.name for s in reopened.spectra], ["first", "second"])
            np.testing.assert_array_equal(reopened.var1, [10.0, 20.0])
            np.testing.assert_array_equal(reopened.var2, [1.0, 2.0])
            np.testing.assert_array_equal(reopened.spectra[1].rawdata, [[200.0, 3.0], [201.0, 4.0]])

            reopened.remove_data([0])
            reopened.import_hdf5(str(path))
            self.assertEqual(reopened.len, 1)
            self.assertEqual(len(reopened.spectra), 1)
            self.assertEqual(reopened.spectra[0].index, "0")
            self.assertEqual(reopened.spectra[0].name, "second")
            with h5py.File(path, "r") as hdf:
                self.assertEqual(sorted(hdf["ms_dataset"].keys()), ["0"])
                self.assertEqual(hdf["ms_dataset"].attrs["num"], 1)

    def test_fresh_spectrum_keeps_zero_mass_placeholder_without_warning(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "dataset.hdf5")
            engine = self._engine()
            engine.config.datanorm = 1
            dataset = MetaDataSet(engine)
            dataset.new_file(str(path))
            raw = np.array([[100.0, 1.0], [101.0, 2.0]])
            dataset.add_data(raw, name="fresh", export=False)
            dataset.add_data(raw, name="deconvolved", export=False)
            dataset.spectra[1].massdat = np.array([[300.0, 2.0], [310.0, 4.0]])
            dataset.export_hdf5()

            reopened = MetaDataSet(engine)
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always", RuntimeWarning)
                reopened.import_hdf5(str(path))

            self.assertEqual([str(w.message) for w in recorded if w.category is RuntimeWarning], [])
            np.testing.assert_array_equal(reopened.spectra[0].massdat, [[0.0, 0.0], [1000.0, 0.0]])
            np.testing.assert_array_equal(reopened.spectra[1].massdat, [[300.0, 0.5], [310.0, 1.0]])

    def test_missing_optional_spectrum_datasets_load_as_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "minimal.hdf5")
            with h5py.File(path, "w") as hdf:
                group = hdf.create_group("ms_dataset")
                group.attrs["num"] = 1
                spectrum = group.create_group("0")
                spectrum.create_dataset("raw_data", data=np.array([[100.0, 1.0]], dtype=np.float32))
                spectrum.attrs["name"] = "minimal"
                hdf.create_group("config")

            dataset = MetaDataSet(self._engine())
            dataset.import_hdf5(str(path))
            self.assertEqual(dataset.spectra[0].name, "minimal")
            np.testing.assert_array_equal(dataset.spectra[0].data2, [[100.0, 1.0]])
            self.assertEqual(dataset.spectra[0].massdat.size, 0)
            self.assertEqual(dataset.spectra[0].mzgrid.size, 0)
            self.assertEqual(dataset.spectra[0].massgrid.size, 0)

            with h5py.File(path, "r") as hdf:
                self.assertNotIn("v1name", hdf["ms_dataset"].attrs)
                self.assertNotIn("v2name", hdf["ms_dataset"].attrs)

    def test_missing_required_structure_rejects_collection_without_modifying_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "invalid.hdf5")
            with h5py.File(path, "w") as hdf:
                hdf.create_group("unrelated")

            dataset = MetaDataSet(self._engine())
            with self.assertRaisesRegex(ValueError, "missing the ms_dataset group"):
                dataset.import_hdf5(str(path))
            self.assertIsNone(dataset.filename)
            with h5py.File(path, "r") as hdf:
                self.assertEqual(list(hdf.keys()), ["unrelated"])

    def test_invalid_raw_data_rejects_entire_collection_in_full_and_speedy_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "invalid-raw-data.hdf5")
            dataset = MetaDataSet(self._engine())
            dataset.new_file(str(path))
            dataset.add_data(np.array([[100.0, 1.0]], dtype=np.float32), name="first", export=False)
            dataset.add_data(np.array([[200.0, 2.0]], dtype=np.float32), name="second", export=False)
            dataset.export_hdf5()
            dataset.import_hdf5(str(path))
            original_first = dataset.spectra[0]
            original_rawdata = original_first.rawdata.copy()

            with h5py.File(path, "a") as hdf:
                hdf["ms_dataset/0/raw_data"][:] = [[999.0, 9.0]]
                del hdf["ms_dataset/1/raw_data"]

            for speedy in (False, True):
                with self.subTest(speedy=speedy):
                    with self.assertRaisesRegex(ValueError, "spectrum 1.*raw_data"):
                        dataset.import_hdf5(str(path), speedy=speedy)
                    self.assertIs(dataset.spectra[0], original_first)
                    np.testing.assert_array_equal(dataset.spectra[0].rawdata, original_rawdata)

            with h5py.File(path, "a") as hdf:
                hdf["ms_dataset/1"].create_dataset("raw_data", data=np.ones((2, 3), dtype=np.float32))
            with self.assertRaisesRegex(ValueError, r"raw_data must have shape \(N, 2\)"):
                dataset.import_hdf5(str(path))

    def test_reused_object_tracks_external_membership_growth_and_shrinkage(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "external-edits.hdf5")
            dataset = MetaDataSet(self._engine())
            dataset.new_file(str(path))
            dataset.v1name = "Variable 1"
            dataset.add_data(np.array([[100.0, 1.0]]), name="first", export=False)
            dataset.spectra[0].var1 = 10.0
            dataset.export_hdf5()
            dataset.import_hdf5(str(path))

            with h5py.File(path, "a") as hdf:
                group = hdf["ms_dataset"]
                second = group.create_group("1")
                second.create_dataset("raw_data", data=np.array([[200.0, 2.0]], dtype=np.float32))
                second.create_dataset("processed_data", data=np.array([[200.0, 2.0]], dtype=np.float32))
                second.attrs["name"] = "second"
                second.attrs["Variable 1"] = 20.0
                group.attrs["num"] = 2

            dataset.import_hdf5(str(path))
            self.assertEqual(dataset.len, 2)
            self.assertEqual([int(s.index) for s in dataset.spectra], [0, 1])
            self.assertEqual([s.name for s in dataset.spectra], ["first", "second"])
            np.testing.assert_array_equal(dataset.var1, [10.0, 20.0])
            np.testing.assert_array_equal(dataset.spectra[1].rawdata, [[200.0, 2.0]])

            with h5py.File(path, "a") as hdf:
                del hdf["ms_dataset/0"]
                hdf.move("ms_dataset/1", "ms_dataset/0")
                hdf["ms_dataset"].attrs["num"] = 1

            dataset.import_hdf5(str(path))
            self.assertEqual(dataset.len, 1)
            self.assertEqual(len(dataset.spectra), 1)
            self.assertEqual(int(dataset.spectra[0].index), 0)
            self.assertEqual(dataset.spectra[0].name, "second")
            np.testing.assert_array_equal(dataset.var1, [20.0])
            np.testing.assert_array_equal(dataset.spectra[0].rawdata, [[200.0, 2.0]])


if __name__ == "__main__":
    unittest.main()

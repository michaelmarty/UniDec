import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.modules.IMEng import UniDecIM

from tests._test_support import copy_unidec_example


class TestBundledIMProcessing(unittest.TestCase):
    @staticmethod
    def _stage_fixture(directory):
        root = Path(directory)
        spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
        config_dir = root / "aqpz_unidecfiles"
        config_dir.mkdir()
        config = copy_unidec_example(
            config_dir, "IMMS", "aqpz_unidecfiles", "aqpz_conf.dat"
        )
        return spectrum, config

    def test_import_projection_matches_independent_raw_row_sums(self):
        """The importer may fill zero cells, but its m/z projection must preserve source intensities."""
        with tempfile.TemporaryDirectory(prefix="unidec-bundled-im-") as directory:
            spectrum, _ = self._stage_fixture(directory)
            source = np.loadtxt(spectrum)
            source_mz, inverse = np.unique(source[:, 0], return_inverse=True)
            expected_projection = np.bincount(inverse, weights=source[:, 2])

            engine = UniDecIM(ignore_args=True, silent=True)
            engine.open_file(spectrum.name, directory, silent=True)

            np.testing.assert_allclose(engine.data.rawdata[:, 0], source_mz)
            np.testing.assert_allclose(
                engine.data.rawdata[:, 1], expected_projection, rtol=1e-12, atol=1e-12
            )
            self.assertEqual(engine.data.rawdata3.shape[1], 3)
            self.assertEqual(
                len(engine.data.rawdata3),
                len(np.unique(engine.data.rawdata3[:, 0]))
                * len(np.unique(engine.data.rawdata3[:, 1])),
            )
            self.assertTrue(np.isfinite(engine.data.rawdata3).all())
            self.assertEqual(engine.config.imflag, 1)
            self.assertEqual(engine.config.mzbins, 4.0)
            self.assertEqual(Path(engine.config.outfname).parent, Path(directory) / "aqpz_unidecfiles")

    def test_real_im_preprocessing_preserves_raw_grid_and_builds_dense_processed_grid(self):
        """The bundled config crops both axes and creates one value for every processed grid cell."""
        with tempfile.TemporaryDirectory(prefix="unidec-bundled-im-process-") as directory:
            spectrum, _ = self._stage_fixture(directory)
            engine = UniDecIM(ignore_args=True, silent=True)
            engine.open_file(spectrum.name, directory, silent=True)
            raw3 = engine.data.rawdata3.copy()
            raw2 = engine.data.rawdata.copy()

            result = engine.process_data(silent=True)

            self.assertIsNone(result)
            self.assertEqual(engine.config.procflag, 1)
            np.testing.assert_array_equal(engine.data.rawdata3, raw3)
            np.testing.assert_array_equal(engine.data.rawdata, raw2)
            self.assertTrue(np.isfinite(engine.data.data3).all())
            processed_mz = np.unique(engine.data.data3[:, 0])
            processed_dt = np.unique(engine.data.data3[:, 1])
            self.assertEqual(len(engine.data.data3), len(processed_mz) * len(processed_dt))
            self.assertGreaterEqual(processed_mz[0], engine.config.minmz)
            self.assertLessEqual(processed_mz[-1], engine.config.maxmz)
            self.assertGreaterEqual(processed_dt[0], engine.config.mindt)
            self.assertLessEqual(processed_dt[-1], engine.config.maxdt)
            intensity = engine.data.data3[:, 2].reshape(len(processed_mz), len(processed_dt))
            expected_projection = np.column_stack((processed_mz, np.sum(intensity, axis=1)))
            self.assertEqual(expected_projection.shape, (len(processed_mz), 2))
            self.assertTrue(Path(engine.config.infname).is_file())


if __name__ == "__main__":
    unittest.main()

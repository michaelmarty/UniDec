import tempfile
import unittest
from pathlib import Path

import numpy as np

from UniDecImporter.ImporterFactory import ImporterFactory
from unidec.modules.CDEng import UniDecCD

from tests._test_support import copy_unidec_example


class TestBundledCDProcessing(unittest.TestCase):
    def test_dmt_and_saved_npz_contain_identical_events_in_different_orders(self):
        """DMT sorts by m/z while the cached NPZ preserves acquisition order."""
        with tempfile.TemporaryDirectory(prefix="unidec-bundled-cd-readers-") as directory:
            root = Path(directory)
            dmt_dir = root / "dmt"
            npz_dir = root / "npz"
            dmt_dir.mkdir()
            npz_dir.mkdir()
            dmt_path = copy_unidec_example(dmt_dir, "CDMS", "BgalGroEL.dmt")
            npz_path = copy_unidec_example(
                npz_dir, "CDMS", "BgalGroEL_unidecfiles", "BgalGroEL_rawdata.npz"
            )

            dmt_reader = ImporterFactory.create_importer(dmt_path)
            npz_reader = ImporterFactory.create_importer(npz_path)
            try:
                dmt = dmt_reader.get_cdms_data()
                saved = npz_reader.get_cdms_data()
            finally:
                dmt_reader.close()
                npz_reader.close()

            self.assertEqual(dmt.shape, (81227, 5))
            self.assertEqual(saved.shape, (81227, 5))
            self.assertTrue(np.isfinite(dmt).all())
            self.assertTrue(np.isfinite(saved).all())
            self.assertFalse(np.array_equal(dmt, saved))
            dmt_order = np.lexsort(tuple(dmt[:, i] for i in range(4, -1, -1)))
            saved_order = np.lexsort(tuple(saved[:, i] for i in range(4, -1, -1)))
            np.testing.assert_array_equal(dmt[dmt_order], saved[saved_order])
            np.testing.assert_array_equal(dmt[:, 3], np.full(len(dmt), 5.0))
            np.testing.assert_array_equal(dmt[:, 4], np.full(len(dmt), -1.0))

    def test_dmt_configured_processing_matches_independent_event_mask_and_charge_formula(self):
        """Configured strict m/z/charge filters leave exactly the independently selected events."""
        with tempfile.TemporaryDirectory(prefix="unidec-bundled-cd-process-") as directory:
            root = Path(directory)
            dmt_path = copy_unidec_example(root, "CDMS", "BgalGroEL.dmt")
            config_dir = root / "BgalGroEL_unidecfiles"
            config_dir.mkdir()
            copy_unidec_example(
                config_dir, "CDMS", "BgalGroEL_unidecfiles", "BgalGroEL_conf.dat"
            )

            engine = UniDecCD()
            try:
                engine.open_file(str(dmt_path), refresh=True)
                original = engine.darray.copy()
                slope = engine.config.CDslope
                charge = original[:, 1] / slope
                expected_mask = (
                    (original[:, 0] > engine.config.minmz)
                    & (original[:, 0] < engine.config.maxmz)
                    & (charge > engine.config.startz)
                    & (charge < engine.config.endz + 1)
                )

                self.assertEqual(engine.process_data(transform=True), 1)
            finally:
                if getattr(engine, "importer", None) is not None:
                    engine.importer.close()

            self.assertEqual(len(engine.farray), int(np.count_nonzero(expected_mask)))
            np.testing.assert_allclose(engine.farray[:, 0], original[expected_mask, 0])
            np.testing.assert_allclose(engine.zarray, engine.farray[:, 1] / slope)
            self.assertEqual(engine.harray.shape, (len(engine.ztab), len(engine.mz)))
            self.assertEqual(engine.data.data2.shape, (len(engine.mz), 2))
            self.assertEqual(engine.data.massdat.ndim, 2)
            self.assertEqual(engine.data.massdat.shape[1], 2)
            self.assertTrue(np.isfinite(engine.harray).all())
            self.assertTrue(np.isfinite(engine.data.massdat).all())


if __name__ == "__main__":
    unittest.main()

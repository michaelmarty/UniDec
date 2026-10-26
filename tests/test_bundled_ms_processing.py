import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.engine import UniDec

from tests._test_support import copy_unidec_example


class TestBundledMSProcessing(unittest.TestCase):
    cases = (
        ("ADHclean.txt", "ADHclean_unidecfiles", "ADHclean_conf.dat"),
        ("BSA.txt", "BSA_unidecfiles", "BSA_conf.dat"),
        ("POPC_Nanodiscs.txt", "POPC_Nanodiscs_unidecfiles", "POPC_Nanodiscs_conf.dat"),
        ("GroEL UniDec.txt", "GroEL UniDec_unidecfiles", "GroEL UniDec_conf.dat"),
    )

    def test_real_import_and_preprocessing_for_additional_packaged_spectra(self):
        """Each example uses its own config and produces finite, ordered, normalized processed data."""
        for filename, config_directory, config_name in self.cases:
            with self.subTest(spectrum=filename):
                with tempfile.TemporaryDirectory(prefix="unidec-bundled-ms-") as directory:
                    root = Path(directory)
                    spectrum = copy_unidec_example(root, filename)
                    staged_config_directory = root / config_directory
                    staged_config_directory.mkdir()
                    copy_unidec_example(
                        staged_config_directory, config_directory, config_name
                    )

                    engine = UniDec(ignore_args=True, silent=True)
                    engine.open_file(spectrum.name, directory, silent=True)
                    raw = engine.data.rawdata.copy()

                    result = engine.process_data(silent=True)

                    self.assertIsNone(result)
                    self.assertEqual(engine.config.procflag, 1)
                    np.testing.assert_array_equal(engine.data.rawdata, raw)
                    self.assertEqual(raw.ndim, 2)
                    self.assertEqual(raw.shape[1], 2)
                    self.assertGreater(len(raw), 0)
                    self.assertTrue(np.isfinite(raw).all())
                    self.assertEqual(engine.data.data2.ndim, 2)
                    self.assertEqual(engine.data.data2.shape[1], 2)
                    self.assertGreater(len(engine.data.data2), 0)
                    self.assertTrue(np.isfinite(engine.data.data2).all())
                    self.assertTrue(np.all(np.diff(engine.data.data2[:, 0]) > 0))
                    self.assertGreaterEqual(engine.data.data2[0, 0], engine.config.minmz)
                    self.assertLessEqual(engine.data.data2[-1, 0], engine.config.maxmz)
                    if engine.config.datanorm == 1:
                        self.assertAlmostEqual(np.max(engine.data.data2[:, 1]), 1.0, places=12)
                    self.assertTrue(Path(engine.config.infname).is_file())
                    self.assertEqual(Path(engine.config.infname).parent, staged_config_directory)


if __name__ == "__main__":
    unittest.main()

import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.engine import UniDec

from tests._test_support import copy_unidec_example


class TestEngineState(unittest.TestCase):
    def test_headless_state_restores_processed_data_and_nondefault_config(self):
        """A fresh engine restores values from the ZIP after the staged source tree is removed."""
        with tempfile.TemporaryDirectory(prefix="unidec-engine-state-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "ADHclean.txt")
            config_dir = root / "ADHclean_unidecfiles"
            config_dir.mkdir()
            copy_unidec_example(config_dir, "ADHclean_unidecfiles", "ADHclean_conf.dat")

            source = UniDec(ignore_args=True, silent=True)
            source.open_file(spectrum.name, directory, silent=True)
            source.config.startz = 4
            source.config.endz = 44
            source.config.suppression_topn = 3
            source.process_data(silent=True)
            expected = source.data.data2.copy()
            state = root / "adhclean_state.zip"
            source.save_state(str(state))

            spectrum.unlink()
            shutil.rmtree(config_dir)

            restored = UniDec(ignore_args=True, silent=True)
            self.assertTrue(restored.load_state(str(state)))

            self.assertEqual(restored.config.startz, 4)
            self.assertEqual(restored.config.endz, 44)
            self.assertEqual(restored.config.suppression_topn, 3)
            self.assertEqual(restored.config.procflag, 1)
            np.testing.assert_allclose(restored.data.data2, expected, rtol=1e-6, atol=1e-6)


if __name__ == "__main__":
    unittest.main()

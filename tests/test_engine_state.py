import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path, PureWindowsPath

import numpy as np

from unidec.engine import UniDec

from tests._test_support import copy_unidec_example


class TestEngineState(unittest.TestCase):
    def test_historical_8_2_1_state_restores_from_original_archive(self):
        fixture = (Path(__file__).parent / "fixtures" / "historical_8_2_1" /
                   "unidec-8.2.1-state.zip")
        low = np.array([[500, 2], [501, 4], [502, 8], [503, 6], [504, 2]], dtype=float)
        with tempfile.TemporaryDirectory(prefix="unidec-historical-state-") as directory:
            state = Path(directory) / fixture.name
            shutil.copy2(fixture, state)
            with zipfile.ZipFile(state) as archive:
                self.assertEqual(set(archive.namelist()), {
                    "legacy_ms_conf.dat", "legacy_ms_input.dat", "legacy_ms_rawdata.txt",
                })
                config = archive.read("legacy_ms_conf.dat").decode("utf-8")
                values = dict(line.split(maxsplit=1) for line in config.splitlines() if line.strip())
                self.assertEqual(values["version"], "8.2.0")
                self.assertFalse(PureWindowsPath(values["input"]).anchor)
                self.assertFalse(PureWindowsPath(values["output"]).anchor)

            restored = UniDec(ignore_args=True, silent=True)
            self.assertTrue(restored.load_state(str(state)))
            self.assertEqual(restored.config.startz, 4)
            self.assertEqual(restored.config.endz, 12)
            self.assertEqual(restored.config.procflag, 1)
            np.testing.assert_allclose(restored.data.rawdata, low, rtol=0, atol=1e-6)
            np.testing.assert_allclose(
                restored.data.data2,
                [[501, 0.5], [502, 1.0], [503, 0.75]],
                rtol=0, atol=1e-6,
            )

    def test_short_name_state_restores_from_npz_cache(self):
        with tempfile.TemporaryDirectory(prefix="unidec-short-state-") as directory:
            root = Path(directory)
            spectrum = root / "sample.txt"
            np.savetxt(spectrum, [[500.0, 1.0], [501.0, 2.0], [502.0, 1.0]])

            source = UniDec(ignore_args=True, silent=True)
            source.open_file(spectrum.name, directory, silent=True)
            source.config.startz = 4
            source.process_data(silent=True)
            expected = source.data.data2.copy()
            state = root / "sample_state.zip"
            source.save_state(str(state))

            spectrum.unlink()
            shutil.rmtree(root / "sample_unidecfiles")

            restored = UniDec(ignore_args=True, silent=True)
            self.assertTrue(restored.load_state(str(state)))
            self.assertEqual(restored.config.startz, 4)
            self.assertEqual(restored.config.procflag, 1)
            self.assertTrue((root / "sample_unidecfiles" / "rawdata.npz").is_file())
            np.testing.assert_allclose(restored.data.data2, expected)

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

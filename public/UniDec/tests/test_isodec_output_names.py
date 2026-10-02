"""IsoDecGUI's short output names and legacy config loading."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.engine import UniDec


class TestIsoDecOutputNames(unittest.TestCase):
    def test_short_names_and_legacy_configuration(self):
        with tempfile.TemporaryDirectory() as directory:
            spectrum = Path(directory) / "ca_etd.dat"
            np.savetxt(spectrum, [[500., 1.], [501., 2.], [502., 1.]])
            output = Path(directory) / "ca_etd_unidecfiles"

            legacy = UniDec(ignore_args=True)
            legacy.open_file(spectrum.name, directory)
            legacy.config.startz = 7
            legacy.config.masslist = np.array([1000., 2000.])
            legacy.export_config()
            self.assertTrue((output / "ca_etd_conf.dat").is_file())

            engine = UniDec(ignore_args=True)
            engine.open_file(spectrum.name, directory, simple_output=True)
            self.assertEqual(engine.config.startz, 7)
            np.testing.assert_array_equal(engine.config.masslist, [1000., 2000.])
            self.assertEqual(Path(engine.config.confname), output / "conf.dat")
            self.assertEqual(Path(engine.config.infname), output / "input.dat")
            self.assertTrue((output / "rawdata.txt").is_file())
            self.assertEqual(Path(engine.config.outfname), output)

            engine.export_config()
            self.assertTrue((output / "mfile.dat").is_file())
            engine.config.startz = 9
            engine.open_file(spectrum.name, directory, simple_output=True)
            self.assertEqual(engine.config.startz, 7)


if __name__ == "__main__":
    unittest.main()

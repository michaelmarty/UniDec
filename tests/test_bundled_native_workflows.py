import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.engine import UniDec
from unidec.modules.IMEng import UniDecIM

from tests._test_support import copy_unidec_example


class TestBundledNativeWorkflows(unittest.TestCase):
    def test_bsa_native_outputs_match_imported_axis_dimensions(self):
        """A second conventional spectrum exercises native execution without duplicating ADH assertions."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-ms-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            config_dir = root / "BSA_unidecfiles"
            config_dir.mkdir()
            copy_unidec_example(config_dir, "BSA_unidecfiles", "BSA_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            self.assertEqual(engine.run_unidec(silent=True), 0)

            zlen = engine.config.endz - engine.config.startz + 1
            mzlen = len(engine.data.data2)
            self.assertEqual(engine.data.massdat.ndim, 2)
            self.assertEqual(engine.data.massdat.shape[1], 2)
            self.assertEqual(engine.data.fitdat.shape, (mzlen,))
            self.assertEqual(engine.data.mzgrid.shape, (mzlen * zlen, 3))
            self.assertTrue(np.isfinite(engine.data.massdat).all())
            self.assertTrue(np.isfinite(engine.data.fitdat).all())
            self.assertTrue(np.isfinite(engine.data.mzgrid).all())
            engine.pick_peaks()
            self.assertGreater(len(engine.pks.peaks), 0)
            self.assertTrue(
                all(engine.config.masslb <= peak.mass <= engine.config.massub for peak in engine.pks.peaks)
            )

    def test_aqpz_im_native_outputs_match_imported_axis_dimensions(self):
        """The packaged executable must produce a complete, internally consistent IM result set."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-im-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
            config_dir = root / "aqpz_unidecfiles"
            config_dir.mkdir()
            copy_unidec_example(
                config_dir, "IMMS", "aqpz_unidecfiles", "aqpz_conf.dat"
            )

            engine = UniDecIM(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            self.assertEqual(engine.run_unidec(silent=True), 0)

            mzlen = len(np.unique(engine.data.data3[:, 0]))
            dtlen = len(np.unique(engine.data.data3[:, 1]))
            zlen = engine.config.endz - engine.config.startz + 1
            self.assertEqual(engine.data.massdat.ndim, 2)
            self.assertEqual(engine.data.massdat.shape[1], 2)
            self.assertEqual(engine.data.ccsdata.ndim, 2)
            self.assertEqual(engine.data.ccsdata.shape[1], 2)
            self.assertEqual(
                engine.data.massccs.shape,
                (len(engine.data.massdat), len(engine.data.ccsdata)),
            )
            self.assertEqual(engine.data.ccsz.shape, (zlen, len(engine.data.ccsdata)))
            self.assertEqual(engine.data.mztgrid.shape, (mzlen, dtlen, zlen))
            self.assertEqual(engine.data.mzgrid.shape, (mzlen * zlen, 3))
            for values in (
                engine.data.massdat,
                engine.data.ccsdata,
                engine.data.massccs,
                engine.data.ccsz,
                engine.data.mztgrid,
                engine.data.mzgrid,
            ):
                self.assertTrue(np.isfinite(values).all())

            engine.pick_peaks()
            self.assertGreater(len(engine.pks.peaks), 0)
            self.assertTrue(
                all(engine.config.masslb <= peak.mass <= engine.config.massub for peak in engine.pks.peaks)
            )


if __name__ == "__main__":
    unittest.main()

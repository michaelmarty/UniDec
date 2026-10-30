import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from unidec.engine import UniDec
from unidec.modules.IMEng import UniDecIM

from tests._test_support import copy_unidec_example


class TestBundledNativeWorkflows(unittest.TestCase):
    def test_variable_peak_width_text_and_hdf5_inputs(self):
        with tempfile.TemporaryDirectory(prefix="unidec-variable-pw-") as directory:
            spectrum = copy_unidec_example(Path(directory), "BSA.txt")
            engine = UniDec(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)
            engine.config.variablepw = 1
            engine.data.fwhmlist = np.full(len(engine.data.data2), 1.0)

            self.assertEqual(engine.run_unidec(silent=True), 0)
            np.testing.assert_array_equal(np.loadtxt(engine.config.fwhmfile), engine.data.fwhmlist)
            self.assertTrue(np.isfinite(engine.data.massdat).all())

            np.savetxt(engine.config.fwhmfile, [1.0])
            invalid_text = subprocess.run(
                [engine.config.UniDecPath, engine.config.confname], cwd=Path(engine.config.confname).parent,
                capture_output=True, text=True, timeout=60
            )
            self.assertEqual(invalid_text.returncode, 12)
            self.assertIn("does not match processed data length", invalid_text.stderr)

            hdf_path = Path(directory, "variable-pw.hdf5")
            engine.config.write_hdf5(str(hdf_path))
            with h5py.File(hdf_path, "a") as hdf:
                group = hdf.require_group("ms_data")
                group.create_dataset("processed_data", data=engine.data.data2)
                group.create_dataset("fwhm", data=engine.data.fwhmlist)

            result = subprocess.run(
                [engine.config.UniDecPath, str(hdf_path)], capture_output=True, text=True, timeout=60
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with h5py.File(hdf_path, "r") as hdf:
                self.assertTrue(np.isfinite(hdf["ms_data/mass_data"][:]).all())

            with h5py.File(hdf_path, "a") as hdf:
                del hdf["ms_data/fwhm"]
                hdf["ms_data"].create_dataset("fwhm", data=[1.0])
            invalid_hdf5 = subprocess.run(
                [engine.config.UniDecPath, str(hdf_path)], capture_output=True, text=True, timeout=60
            )
            self.assertEqual(invalid_hdf5.returncode, 12)
            self.assertIn("does not match processed data length", invalid_hdf5.stderr)

    def test_native_text_config_uses_exact_parameter_keys(self):
        with tempfile.TemporaryDirectory(prefix="unidec-native-config-keys-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            engine = UniDec(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            command = [engine.config.UniDecPath, engine.config.confname, "-nthreads", "1"]
            control = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(control.returncode, 0, control.stdout + control.stderr)
            expected = np.loadtxt(engine.config.massdatfile)

            with open(engine.config.confname, "a", encoding="utf-8") as config_file:
                config_file.write("prefix_startz 20\n")
                config_file.write("endz_suffix 30\n")
                config_file.write("numit_extra 1\n")

            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            actual = np.loadtxt(engine.config.massdatfile)

            np.testing.assert_array_equal(actual, expected)

    def test_native_resolves_relative_config_paths_from_config_directory(self):
        with tempfile.TemporaryDirectory(prefix="unidec-native-relative-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            engine = UniDec(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            result = subprocess.run(
                [engine.config.UniDecPath, engine.config.confname],
                cwd=root,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(Path(engine.config.udir, "mass.txt").is_file())

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

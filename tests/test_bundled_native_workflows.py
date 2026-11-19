import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py
import numpy as np

from unidec.engine import UniDec
from unidec.modules.IMEng import UniDecIM

from tests._test_support import copy_unidec_example


REVIEWED_CONFIGS = Path(__file__).parent / "fixtures" / "native_regression"


def _copy_reviewed_config(destination, name):
    target = Path(destination, name)
    shutil.copy2(REVIEWED_CONFIGS / name, target)
    return target


def _charge_summary(engine, peak):
    """Return the weighted and apex charge at the mass-grid row nearest a picked peak."""
    grid = np.asarray(engine.data.massgrid).reshape(len(engine.data.massdat), len(engine.data.ztab))
    row = grid[np.argmin(np.abs(engine.data.massdat[:, 0] - peak.mass))]
    return np.average(engine.data.ztab, weights=row), engine.data.ztab[np.argmax(row)]


def _ccs_summary(engine, peak):
    """Return the weighted and apex CCS at the mass-grid row nearest a picked peak."""
    row = engine.data.massccs[np.argmin(np.abs(engine.data.massdat[:, 0] - peak.mass))]
    axis = engine.data.ccsdata[:, 0]
    return np.average(axis, weights=row), axis[np.argmax(row)]


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

    def test_adh_reviewed_native_results(self):
        """Reviewed 2026-11-18: bundled ADH settings produce three known mass/charge species."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-adh-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "ADH.txt")
            config_dir = root / "ADH_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "ADH_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            if not os.path.isfile(engine.config.UniDecPath):
                self.skipTest(f"Packaged native executable is unavailable: {engine.config.UniDecPath}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)
            self.assertEqual(engine.run_unidec(silent=True), 0)
            engine.pick_peaks()

            self.assertAlmostEqual(engine.config.error, 0.95644, delta=0.005)
            peaks = sorted(engine.pks.peaks, key=lambda peak: peak.mass)
            expected = (
                (78560, 5.06, 560, 18.51, 19),
                (110930, 2.19, 340, 10.24, 9),
                (147980, 100, 800, 26.31, 26),
            )
            self.assertEqual(len(peaks), len(expected))
            for peak, (mass, height, fwhm, mean_charge, apex_charge) in zip(peaks, expected):
                weighted_charge, charge_apex = _charge_summary(engine, peak)
                self.assertAlmostEqual(peak.mass, mass, delta=20)
                self.assertAlmostEqual(peak.height, height, delta=0.5)
                self.assertAlmostEqual(peak.errorFWHM, fwhm, delta=20)
                self.assertAlmostEqual(weighted_charge, mean_charge, delta=0.1)
                self.assertEqual(charge_apex, apex_charge)

    def test_bsa_native_outputs_match_imported_axis_dimensions(self):
        """Reviewed BSA results and native output dimensions remain stable."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-ms-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            config_dir = root / "BSA_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "BSA_conf.dat")

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
            self.assertAlmostEqual(engine.config.error, 0.99728, delta=0.005)
            self.assertEqual(len(engine.pks.peaks), 1)
            peak = engine.pks.peaks[0]
            weighted_charge, charge_apex = _charge_summary(engine, peak)
            self.assertAlmostEqual(peak.mass, 66430, delta=20)
            self.assertAlmostEqual(peak.height, 100, delta=0.5)
            self.assertAlmostEqual(peak.errorFWHM, 30, delta=20)
            self.assertAlmostEqual(weighted_charge, 14.92, delta=0.1)
            self.assertEqual(charge_apex, 15)
            self.assertTrue(
                all(engine.config.masslb <= peak.mass <= engine.config.massub for peak in engine.pks.peaks)
            )

    def test_aqpz_im_native_outputs_match_imported_axis_dimensions(self):
        """The packaged executable must produce a complete, internally consistent IM result set."""
        # MKL FFT kernels differ by CPU instruction set and shift this iterative IM result,
        # so pin the reference-value run to the AVX2 kernels used for the expected numbers.
        env = mock.patch.dict(os.environ, {"MKL_ENABLE_INSTRUCTIONS": "AVX2"})
        env.start()
        self.addCleanup(env.stop)
        with tempfile.TemporaryDirectory(prefix="unidec-native-im-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
            config_dir = root / "aqpz_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "aqpz_conf.dat")

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
            self.assertAlmostEqual(engine.config.error, 0.87040, delta=0.005)
            peaks = sorted(engine.pks.peaks, key=lambda peak: peak.mass)
            expected = (
                (98849.92, 100, 13.40, 13, 5477.51, 5021.27),
                (99619.92, 86.47, 13.35, 13, 5472.59, 4921.27),
                (100389.92, 47.63, 13.22, 13, 5397.56, 4921.27),
                (105909.92, 10.02, 15.06, 15, 6213.10, 6081.27),
                (108139.92, 12.38, 14.09, 14, 5642.64, 5421.27),
            )
            self.assertEqual(len(peaks), len(expected))
            for peak, values in zip(peaks, expected):
                mass, height, mean_charge, apex_charge, mean_ccs, apex_ccs = values
                weighted_charge, charge_apex = _charge_summary(engine, peak)
                weighted_ccs, ccs_apex = _ccs_summary(engine, peak)
                self.assertAlmostEqual(peak.mass, mass, delta=20)
                self.assertAlmostEqual(peak.height, height, delta=1)
                self.assertAlmostEqual(weighted_charge, mean_charge, delta=0.1)
                self.assertEqual(charge_apex, apex_charge)
                self.assertAlmostEqual(weighted_ccs, mean_ccs, delta=25)
                self.assertAlmostEqual(ccs_apex, apex_ccs, delta=20)
            self.assertTrue(
                all(engine.config.masslb <= peak.mass <= engine.config.massub for peak in engine.pks.peaks)
            )


if __name__ == "__main__":
    unittest.main()

import ctypes
import os
import platform
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py
import numpy as np

from unidec.engine import UniDec
from unidec.modules import IM_functions, MassSpecBuilder
from unidec.modules.IMEng import UniDecIM
from unidec.modules.unidecstructure import UniDecConfig

from tests._test_support import copy_unidec_example


REVIEWED_CONFIGS = Path(__file__).parent / "fixtures" / "native_regression"
REQUIRE_NATIVE_WORKFLOWS = os.environ.get("UNIDEC_REQUIRE_NATIVE_WORKFLOWS") == "1"


def _require_native(condition, message):
    if condition:
        return
    if REQUIRE_NATIVE_WORKFLOWS:
        raise AssertionError(message)
    raise unittest.SkipTest(message)


def _require_native_executable(engine):
    _require_native(
        os.path.isfile(engine.config.UniDecPath),
        f"Packaged native executable is unavailable: {engine.config.UniDecPath}",
    )


def _load_native_library():
    name = {
        "Windows": "unideclib.dll",
        "Darwin": "libunideclib.dylib",
        "Linux": "libunideclib.so",
    }[platform.system()]
    return ctypes.CDLL(str(Path(__file__).resolve().parents[1] / "unidec" / "bin" / name))


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


def _linear_ccs_reference(mass, charge, drift_time, config):
    """Independent Mason-Schamp conversion for the reviewed drift-tube calibration."""
    elementary_charge = 1.60217657e-19
    boltzmann = 1.3806488e-23
    loschmidt = 2.6867774e25
    atomic_mass = 1.6605389e-27
    standard_temperature = 273.15
    temperature = config.temp + standard_temperature
    ccs_constant = (
        np.sqrt(18 * np.pi) / 16
        * elementary_charge / np.sqrt(boltzmann * temperature) / loschmidt
        * config.volt / config.driftlength ** 2
        * 760 / config.pressure
        * temperature / standard_temperature
        * 1e20
    )
    reduced_mass = atomic_mass * mass * config.gasmass / (mass + config.gasmass)
    return charge * (drift_time - config.to) * 0.001 * np.sqrt(1 / reduced_mass) * ccs_constant


class TestBundledNativeWorkflows(unittest.TestCase):
    def test_native_regularization_and_suppression_operations(self):
        """Native regularizers preserve their conservative operation-level contracts."""
        library = _load_native_library()
        float_pointer = ctypes.POINTER(ctypes.c_float)
        int_pointer = ctypes.POINTER(ctypes.c_int)
        char_pointer = ctypes.POINTER(ctypes.c_char)
        library.softargmax.argtypes = [float_pointer, ctypes.c_int, ctypes.c_int, ctypes.c_float]
        library.softargmax.restype = None
        library.point_smoothing.argtypes = [
            float_pointer, float_pointer, float_pointer, char_pointer,
            ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ]
        library.point_smoothing.restype = None
        library.apply_suppressions.argtypes = [
            float_pointer, float_pointer, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_int, int_pointer, ctypes.c_int,
            ctypes.c_float, ctypes.c_float,
        ]
        library.apply_suppressions.restype = None

        def floats(values):
            return np.ascontiguousarray(values, dtype=np.float32)

        charge_rows = floats([[1, 2, 4], [4, 3, 1]])
        original_sums = charge_rows.sum(axis=1)
        original_fractions = charge_rows.max(axis=1) / original_sums
        library.softargmax(charge_rows.ctypes.data_as(float_pointer), 2, 3, 1.0)
        np.testing.assert_allclose(charge_rows.sum(axis=1), original_sums, rtol=1e-6, atol=1e-7)
        self.assertTrue(np.all(charge_rows.max(axis=1) / charge_rows.sum(axis=1) > original_fractions))
        self.assertTrue(np.isfinite(charge_rows).all())
        self.assertTrue((charge_rows >= 0).all())

        unsmoothed = floats([[4, 0], [0, 4], [4, 0], [0, 4], [4, 0]])
        smoothed = unsmoothed.copy()
        scratch = np.empty_like(smoothed)
        sums = np.empty(smoothed.shape[1], dtype=np.float32)
        mask = np.ones(smoothed.shape, dtype=np.int8)
        library.point_smoothing(
            smoothed.ctypes.data_as(float_pointer), scratch.ctypes.data_as(float_pointer),
            sums.ctypes.data_as(float_pointer), mask.ctypes.data_as(char_pointer), 5, 2, 1,
        )
        self.assertLess(np.abs(np.diff(smoothed, axis=0)).sum(), np.abs(np.diff(unsmoothed, axis=0)).sum())
        self.assertTrue(np.isfinite(smoothed).all())
        self.assertTrue((smoothed >= 0).all())

        ztab = np.arange(2, 9, dtype=np.int32)

        def suppress(values, *, satellite=0, harmonic=0, top_n=0, top_x=0):
            grid = floats([values])
            temporary = np.empty_like(grid)
            library.apply_suppressions(
                grid.ctypes.data_as(float_pointer), temporary.ctypes.data_as(float_pointer),
                1, grid.shape[1], satellite, harmonic, ztab.ctypes.data_as(int_pointer),
                top_n, top_x, 0.0,
            )
            return grid[0]

        np.testing.assert_allclose(
            suppress([5, 5, 2, 1], top_n=1), [5, 5, 0, 0], rtol=1e-6, atol=1e-7
        )
        np.testing.assert_allclose(
            suppress([5, 2, 1.99, 0], top_x=0.4), [5, 2, 0, 0], rtol=1e-6, atol=1e-7
        )
        np.testing.assert_allclose(
            suppress([0, 1, 5, 2, 0, 0, 0], satellite=1), [0, 0, 5, 0, 0, 0, 0],
            rtol=1e-6, atol=1e-7,
        )
        np.testing.assert_allclose(
            suppress([3, 0, 7, 0, 0, 0, 0], harmonic=1), [0, 0, 7, 0, 0, 0, 0],
            rtol=1e-6, atol=1e-7,
        )

    def test_native_regularization_and_suppression_workflow(self):
        """Enabled regularizers reach the executable without moving the BSA mass assignment."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-regularization-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            config_dir = root / "BSA_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "BSA_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            def run(settings, iterations=100):
                engine.config.numit = iterations
                engine.config.beta = 0
                engine.config.psig = 0
                engine.config.suppression_topn = 0
                engine.config.suppression_topx = 0
                engine.config.suppression_percent = 0
                engine.config.suppression_startit = 3
                engine.config.suppression_harmonic = 0
                engine.config.suppression_satellite = 0
                for name, value in settings.items():
                    setattr(engine.config, name, value)
                self.assertEqual(engine.run_unidec(silent=True), 0)
                mass_data = engine.data.massdat.copy()
                mz_grid = engine.data.mzgrid[:, 2].reshape(len(engine.data.data2), -1).copy()
                self.assertTrue(np.isfinite(mass_data).all())
                self.assertTrue(np.isfinite(mz_grid).all())
                self.assertTrue((mass_data[:, 1] >= 0).all())
                self.assertTrue((mz_grid >= 0).all())
                return mass_data, mz_grid, mass_data[np.argmax(mass_data[:, 1]), 0]

            control = run({})
            beta = run({"beta": 50})
            point_smoothed = run({"psig": 10})
            for result in (beta, point_smoothed):
                self.assertAlmostEqual(result[2], control[2], delta=engine.config.massbins)
                self.assertGreater(np.max(np.abs(result[1] - control[1])), 1e-6)
            self.assertLess(
                np.abs(np.diff(point_smoothed[1], axis=0)).sum(),
                np.abs(np.diff(control[1], axis=0)).sum(),
            )

            suppressed = run({
                "suppression_topn": 1,
                "suppression_topx": 0.4,
                "suppression_satellite": 1,
                "suppression_harmonic": 1,
                "suppression_startit": 0,
            })
            self.assertAlmostEqual(suppressed[2], control[2], delta=engine.config.massbins)
            self.assertGreater(np.max(np.abs(suppressed[1] - control[1])), 1e-6)
            self.assertAlmostEqual(suppressed[1].sum(), control[1].sum(), delta=0.25 * control[1].sum())

            short_control = run({}, iterations=4)
            boundary_inactive = run({"suppression_topn": 1, "suppression_startit": 3}, iterations=4)
            boundary_active = run({"suppression_topn": 1, "suppression_startit": 2}, iterations=4)
            np.testing.assert_array_equal(boundary_inactive[0], short_control[0])
            np.testing.assert_array_equal(boundary_inactive[1], short_control[1])
            self.assertGreater(np.max(np.abs(boundary_active[1] - short_control[1])), 1e-6)
            self.assertAlmostEqual(boundary_active[2], short_control[2], delta=engine.config.massbins)

    def test_aqpz_linear_ccs_formula_and_round_trip(self):
        """The saved aqpz drift-tube calibration defines CCS in square angstroms."""
        config = UniDecConfig()
        config.config_import(str(REVIEWED_CONFIGS / "aqpz_conf.dat"))
        mass = np.array([90000.0, 100000.0, 105000.0])
        charge = np.array([10.0, 13.0, 18.0])
        drift_time = np.array([10.08, 12.5, 14.94])
        expected = _linear_ccs_reference(mass, charge, drift_time, config)

        actual = IM_functions.calc_linear_ccs(mass, charge, drift_time, config)
        round_trip = IM_functions.calc_linear_dt(mass, charge, actual, config)

        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=0.01)
        np.testing.assert_allclose(round_trip, drift_time, rtol=1e-12, atol=1e-12)

    def test_variable_peak_width_text_and_hdf5_inputs(self):
        with tempfile.TemporaryDirectory(prefix="unidec-variable-pw-") as directory:
            spectrum = copy_unidec_example(Path(directory), "BSA.txt")
            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
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

    def test_native_rejects_invalid_text_and_hdf5_spectra_before_writing_outputs(self):
        engine = UniDec(ignore_args=True, silent=True)
        _require_native_executable(engine)

        text_cases = {
            "empty": ("", "at least two"),
            "singleton": ("1000 1\n", "at least two"),
            "nonfinite-mz": ("nan 1\n1001 2\n", "non-finite"),
            "nonfinite-intensity": ("1000 1\n1001 inf\n", "non-finite"),
            "truncated": ("1000 1\n1001\n", "two numeric fields"),
        }
        with tempfile.TemporaryDirectory(prefix="unidec-native-invalid-") as directory:
            root = Path(directory)
            for name, (contents, message) in text_cases.items():
                with self.subTest(format="text", case=name):
                    case_dir = root / f"text-{name}"
                    case_dir.mkdir()
                    shutil.copy2(REVIEWED_CONFIGS / "BSA_conf.dat", case_dir / "BSA_conf.dat")
                    (case_dir / "BSA_input.dat").write_text(contents, encoding="ascii")

                    result = subprocess.run(
                        [engine.config.UniDecPath, str(case_dir / "BSA_conf.dat")],
                        cwd=case_dir,
                        capture_output=True,
                        text=True,
                        timeout=60,
                    )

                    self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
                    self.assertIn(message, result.stderr)
                    self.assertEqual(
                        {path.name for path in case_dir.iterdir()},
                        {"BSA_conf.dat", "BSA_input.dat"},
                    )

            hdf5_cases = {
                "empty": (np.empty((0, 2)), "at least two"),
                "singleton": (np.array([[1000.0, 1.0]]), "at least two"),
                "nonfinite": (np.array([[1000.0, 1.0], [1001.0, np.nan]]), "non-finite"),
            }
            for name, (data, message) in hdf5_cases.items():
                with self.subTest(format="hdf5", case=name):
                    hdf_path = root / f"{name}.hdf5"
                    config = UniDecConfig()
                    config.config_import(str(REVIEWED_CONFIGS / "BSA_conf.dat"))
                    config.write_hdf5(str(hdf_path))
                    with h5py.File(hdf_path, "a") as hdf:
                        hdf.require_group("ms_data").create_dataset("processed_data", data=data)
                        before = []
                        hdf.visit(before.append)

                    result = subprocess.run(
                        [engine.config.UniDecPath, str(hdf_path)],
                        capture_output=True,
                        text=True,
                        timeout=60,
                    )

                    self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
                    self.assertIn(message, result.stderr)
                    with h5py.File(hdf_path, "r") as hdf:
                        after = []
                        hdf.visit(after.append)
                    self.assertEqual(after, before)

    def test_native_text_config_uses_exact_parameter_keys(self):
        with tempfile.TemporaryDirectory(prefix="unidec-native-config-keys-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
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
            _require_native_executable(engine)
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

    def test_synthetic_native_mass_and_charge_reconstruction(self):
        """A known 50 kDa charge series reconstructs its independently defined mass and charges."""
        target_mass = 50_000.0
        adduct_mass = 1.007276467
        source_charges = np.arange(8, 13)
        spectrum, _ = MassSpecBuilder.make_mass_spectrum(
            [[target_mass, 8.5, 10.2, 1.0, 1.0]],
            zrange=(source_charges[0], source_charges[-1]),
            mzrange=(4000, 6501),
            mz_bin_size=0.5,
            adductmass=adduct_mass,
            psfun=0,
            noise=0,
            baseline=0,
        )
        source_mz = (target_mass + adduct_mass * source_charges) / source_charges
        source_heights = spectrum[
            [np.argmin(np.abs(spectrum[:, 0] - mz)) for mz in source_mz], 1
        ]
        self.assertFalse(np.allclose(source_heights, source_heights[::-1]))

        with tempfile.TemporaryDirectory(prefix="unidec-native-synthetic-") as directory:
            spectrum_path = Path(directory, "synthetic.txt")
            np.savetxt(spectrum_path, spectrum)

            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
            engine.open_file(spectrum_path.name, directory, silent=True)

            for name, value in {
                "smooth": 0,
                "subbuff": 0,
                "intthresh": 0,
                "reductionpercent": 0,
                "mzbins": 0,
                "linflag": 2,
                "startz": 1,
                "endz": 50,
                "masslb": -45_000,
                "massub": -55_000,
                "massbins": 10,
                "adductmass": adduct_mass,
                "numit": 100,
                "mzsig": 0.85,
                "psfun": 0,
                "zzsig": 1,
                "psig": 0,
                "beta": 0,
                "msig": 0,
                "molig": 0,
                "baselineflag": 0,
                "poolflag": 0,
                "suppression_topn": 0,
                "suppression_topx": 0,
                "suppression_percent": 0,
                "suppression_harmonic": 0,
                "suppression_satellite": 0,
                "peakwindow": 100,
                "peakthresh": 0.1,
            }.items():
                setattr(engine.config, name, value)

            engine.process_data(silent=True)
            self.assertEqual(engine.run_unidec(silent=True), 0)

            expected_mass_axis = np.arange(45_000, 55_000, 10, dtype=float)
            np.testing.assert_array_equal(engine.data.massdat[:, 0], expected_mass_axis)
            self.assertEqual(engine.data.massgrid.size, len(expected_mass_axis) * 50)
            self.assertEqual(engine.data.mzgrid.shape, (len(engine.data.data2) * 50, 3))
            for values in (engine.data.massdat, engine.data.massgrid, engine.data.mzgrid, engine.data.fitdat):
                self.assertTrue(np.isfinite(values).all())
            for intensities in (
                engine.data.massdat[:, 1], engine.data.massgrid,
                engine.data.mzgrid[:, 2], engine.data.fitdat,
            ):
                self.assertTrue((intensities >= 0).all())

            dominant_mass = engine.data.massdat[np.argmax(engine.data.massdat[:, 1]), 0]
            self.assertAlmostEqual(dominant_mass, target_mass, delta=engine.config.massbins)

            charge_grid = engine.data.mzgrid[:, 2].reshape(len(engine.data.data2), 50)
            charge_axis = np.arange(1, 51)
            for mz, charge in zip(source_mz, source_charges):
                row = np.argmin(np.abs(engine.data.data2[:, 0] - mz))
                self.assertEqual(charge_axis[np.argmax(charge_grid[row])], charge)

            observed = engine.data.data2[:, 1] / np.max(engine.data.data2[:, 1])
            fitted = engine.data.fitdat / np.max(engine.data.fitdat)
            self.assertLessEqual(np.sqrt(np.mean((observed - fitted) ** 2)), 0.05)

            engine.pick_peaks(calc_dscore=False)
            self.assertEqual(len(engine.pks.peaks), 1)
            self.assertAlmostEqual(engine.pks.peaks[0].mass, target_mass, delta=engine.config.massbins)

    def test_adh_reviewed_native_results(self):
        """Reviewed 2026-11-18: bundled ADH settings produce three known mass/charge species."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-adh-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "ADH.txt")
            config_dir = root / "ADH_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "ADH_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
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
            _require_native_executable(engine)
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

    def test_bsa_native_results_agree_across_thread_counts(self):
        """One- and two-thread runs must produce the same scientific result."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-threads-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            config_dir = root / "BSA_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "BSA_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            _require_native_executable(engine)
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)

            outputs = []
            for threads in (1, 2):
                result = subprocess.run(
                    [engine.config.UniDecPath, engine.config.confname, "-nthreads", str(threads)],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                output = result.stdout + result.stderr
                if "OpenMP is unavailable" in output:
                    _require_native(False, "Packaged native executable does not support OpenMP")
                self.assertEqual(result.returncode, 0, output)
                self.assertIn(f"Number of Threads: {threads}", output)
                engine.unidec_imports()
                engine.pick_peaks()
                outputs.append((
                    engine.data.massdat.copy(),
                    engine.data.fitdat.copy(),
                    engine.data.mzgrid.copy(),
                    np.sort([peak.mass for peak in engine.pks.peaks]),
                ))

            single, double = outputs
            for left, right in zip(single[:3], double[:3]):
                self.assertEqual(left.shape, right.shape)
            np.testing.assert_array_equal(single[0][:, 0], double[0][:, 0])
            np.testing.assert_array_equal(single[2][:, :2], double[2][:, :2])
            np.testing.assert_allclose(single[0][:, 1], double[0][:, 1], rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(single[1], double[1], rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(single[2][:, 2], double[2][:, 2], rtol=1e-5, atol=1e-6)
            self.assertEqual(len(single[3]), len(double[3]))
            np.testing.assert_allclose(single[3], double[3], rtol=0, atol=engine.config.massbins / 2)

    def test_aqpz_im_native_outputs_match_imported_axis_dimensions(self):
        """The packaged executable must produce a complete, internally consistent IM result set."""
        # MKL FFT kernels differ by CPU and shift this iterative IM result, so run with
        # conditional numerical reproducibility; the expected numbers are for this mode.
        env = mock.patch.dict(os.environ, {"MKL_CBWR": "COMPATIBLE"})
        env.start()
        self.addCleanup(env.stop)
        with tempfile.TemporaryDirectory(prefix="unidec-native-im-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
            config_dir = root / "aqpz_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "aqpz_conf.dat")

            engine = UniDecIM(ignore_args=True, silent=True)
            _require_native_executable(engine)
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
            self.assertAlmostEqual(engine.config.error, 0.881318, delta=0.005)
            peaks = sorted(engine.pks.peaks, key=lambda peak: peak.mass)
            expected = (
                (98849.92, 100, 13.47, 13, 5508.26, 5001.27),
                (99619.92, 91.56, 13.37, 13, 5465.03, 4921.27),
                (100389.92, 47.03, 13.31, 13, 5440.15, 4921.27),
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

    def test_aqpz_native_integration_preserves_independent_mass_ccs_centroid(self):
        """Forward integration places the native m/z-drift grid on the calibrated CCS grid."""
        with tempfile.TemporaryDirectory(prefix="unidec-native-ccs-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
            config_dir = root / "aqpz_unidecfiles"
            config_dir.mkdir()
            _copy_reviewed_config(config_dir, "aqpz_conf.dat")

            engine = UniDecIM(ignore_args=True, silent=True)
            _require_native_executable(engine)
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)
            engine.config.poolflag = 0

            self.assertEqual(engine.run_unidec(silent=True), 0)

            mz_axis = np.unique(engine.data.data3[:, 0])
            drift_axis = np.unique(engine.data.data3[:, 1])
            charge_axis = np.asarray(engine.data.ztab)
            mass = (
                (mz_axis[:, np.newaxis] - engine.config.adductmass)
                * charge_axis[np.newaxis, :]
            )
            ccs = _linear_ccs_reference(
                mass[:, np.newaxis, :],
                charge_axis[np.newaxis, np.newaxis, :],
                drift_axis[np.newaxis, :, np.newaxis],
                engine.config,
            )
            mass_grid = np.broadcast_to(mass[:, np.newaxis, :], engine.data.mztgrid.shape)
            output_mass_axis = engine.data.massdat[:, 0]
            output_ccs_axis = engine.data.ccsdata[:, 0]
            interior = (
                (mass_grid > output_mass_axis[0])
                & (mass_grid < output_mass_axis[-1])
                & (ccs > output_ccs_axis[0])
                & (ccs < output_ccs_axis[-1])
            )
            source_weights = np.where(interior, engine.data.mztgrid, 0)
            expected_total = np.sum(source_weights, dtype=np.float64)
            expected_mass = np.sum(source_weights * mass_grid) / expected_total
            expected_ccs = np.sum(source_weights * ccs) / expected_total

            transformed = np.asarray(engine.data.massccs)
            actual_total = np.sum(transformed, dtype=np.float64)
            actual_mass = np.average(output_mass_axis, weights=np.sum(transformed, axis=1))
            actual_ccs = np.average(output_ccs_axis, weights=np.sum(transformed, axis=0))

            self.assertAlmostEqual(actual_total, expected_total, delta=expected_total * 1e-5)
            self.assertAlmostEqual(
                actual_mass, expected_mass, delta=engine.config.massbins / 2
            )
            self.assertAlmostEqual(
                actual_ccs, expected_ccs, delta=engine.config.ccsbins / 2
            )


if __name__ == "__main__":
    unittest.main()

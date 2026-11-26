import ctypes
import math
import os
import platform
import subprocess
import tempfile
import unittest
from pathlib import Path

from unidec.modules.unidecstructure import UniDecConfig


class TestNativeBuild(unittest.TestCase):
    def test_platform_executable_loads(self):
        config = UniDecConfig()
        config.initialize_system_paths()
        executable = Path(config.UniDecPath)
        self.assertTrue(executable.is_file(), f"Missing native executable: {executable}")

        result = subprocess.run(
            [str(executable)],
            capture_output=True,
            text=True,
            timeout=30,
        )

        self.assertEqual(result.returncode, 88)
        self.assertIn("Universal Deconvolution", result.stdout)

    def test_platform_shared_library_loads(self):
        system = platform.system()
        library_name = {
            "Windows": "unideclib.dll",
            "Darwin": "libunideclib.dylib",
            "Linux": "libunideclib.so",
        }.get(system)
        if library_name is None:
            self.skipTest(f"No native library name configured for {system}")

        library = Path(__file__).resolve().parents[1] / "unidec" / "bin" / library_name
        self.assertTrue(library.is_file(), f"Missing native shared library: {library}")
        ctypes.CDLL(str(library))

    def test_native_text_loader_reads_variable_peak_width_settings(self):
        from unidec.modules.unidecwrapper import Config, unideclib

        with tempfile.TemporaryDirectory(prefix="unidec-native-config-") as directory:
            path = Path(directory, "config.dat")
            path.write_text(
                "input input.dat\n"
                "output output\n"
                "variablepw 2\n"
                "minratio 0.25\n"
                "fwhmfile fwhm.dat\n",
                encoding="utf-8",
            )

            config = Config()
            unideclib.SetDefaultConfig(ctypes.byref(config))
            unideclib.LoadConfig.argtypes = [Config, ctypes.c_char_p]
            unideclib.LoadConfig.restype = Config

            restored = unideclib.LoadConfig(config, os.fsencode(path))

            self.assertEqual(restored.variablepw, 2)
            self.assertAlmostEqual(restored.minratio, 0.25)
            self.assertEqual(Path(os.fsdecode(restored.fwhmfile)), path.parent / "fwhm.dat")

    def test_native_hdf5_loader_declares_persistence_settings(self):
        source = Path(__file__).resolve().parents[1] / "unidec" / "src" / "h5io.c"
        text = source.read_text(encoding="utf-8")

        self.assertIn('"variablepw", config.variablepw', text)
        self.assertIn('"minratio", config.minratio', text)
        legacy_psfunz = 'config.zpsfun = int_attr(file_id, "/config", "zpsfn", config.zpsfun);'
        canonical_psfunz = 'config.zpsfun = int_attr(file_id, "/config", "psfunz", config.zpsfun);'
        self.assertIn(legacy_psfunz, text)
        self.assertIn(canonical_psfunz, text)
        self.assertLess(text.index(legacy_psfunz), text.index(canonical_psfunz))
        input_source = source.with_name("udio.c").read_text(encoding="utf-8")
        self.assertIn('strjoin(config->dataset, "/fwhm", fwhmdat)', input_source)

    def test_native_background_windows_are_symmetric(self):
        source = Path(__file__).resolve().parents[1] / "unidec" / "src" / "UD_dataproc.c"
        text = source.read_text(encoding="utf-8")

        self.assertIn("int start = i - radius;", text)
        self.assertIn("int end = i + radius + 1;", text)
        self.assertIn("const int threshold = gaussian_radius * 2 + 1;", text)
        self.assertIn("int start = i - gaussian_radius;", text)
        self.assertIn("int end = i + gaussian_radius + 1;", text)
        self.assertIn("gaussian[j - i + gaussian_radius]", text)

    def test_non_speedy_peak_shape_rows_include_both_endpoints(self):
        from unidec.modules.unidecwrapper import Config, Decon, Input, unideclib

        mz = (ctypes.c_float * 6)(98, 100, 101, 102, 103, 105)
        intensity = (ctypes.c_float * 6)(0.1, 0.4, 1, 0.8, 0.3, 0.1)
        charges = (ctypes.c_int * 1)(2)
        config = Config()
        unideclib.SetDefaultConfig(ctypes.byref(config))
        config.lengthmz = len(mz)
        config.startz = config.endz = 2
        config.numz = 1
        config.numit = 1
        config.masslb, config.massub, config.massbins = 190, 220, 1
        config.linflag = 2
        config.mzsig = 2.35482  # Gaussian FWHM for sigma = 1 m/z.
        config.psthresh = 1
        config.psfun = 0
        config.silent = 1
        unideclib.PostImport(ctypes.byref(config))
        self.assertEqual(config.speedyflag, 0)

        decon = Decon()
        free_decon = unideclib.FreeDecon
        free_decon.argtypes = [Decon]
        free_decon.restype = None
        self.addCleanup(free_decon, decon)
        inputs = Input(dataMZ=mz, dataInt=intensity, nztab=charges)
        self.assertEqual(unideclib.run_unidec_core(config, inputs, ctypes.byref(decon), 0), 0)

        # Adjacent interior windows span m/z 100–102 and 101–103.
        # Each Gaussian row needs three slots, including its right endpoint.
        self.assertEqual(decon.maxlength, 3)
        flank = math.exp(-0.5)
        for row, start in ((2, 1), (3, 2)):
            with self.subTest(row=row):
                self.assertEqual((decon.starttab[row], decon.endtab[row]), (start, start + 2))
                offset = row * decon.maxlength
                for column, expected in enumerate((flank, 1, flank)):
                    self.assertAlmostEqual(decon.mzdist[offset + column], expected, delta=1e-6)

    def test_constant_resolution_modes_preserve_fwhm_in_mz(self):
        from types import SimpleNamespace

        import numpy as np
        from unidec import tools

        config = UniDecConfig()
        config.initialize_system_paths()
        executable = Path(os.environ.get("UNIDEC_TEST_EXECUTABLE", config.UniDecPath))
        fwhm = 6.0
        bin_size = 0.15
        raw = np.array([[150.0, 1.0], [1200.0, 1.0]])
        centers = np.array([200.0, 1000.0])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "conf.dat")
            for linflag in (1, 4):
                mz = tools.linearize(raw, bin_size, linflag)[:, 0]
                apexes = np.searchsorted(mz, centers)
                sticks = np.zeros_like(mz, dtype=np.float32)
                sticks[apexes] = 1.0
                np.savetxt(Path(directory, "input.dat"), np.column_stack((mz, sticks)))
                sticks.tofile(Path(directory, "result_grid.bin"))
                path.write_text(
                    "input input.dat\noutput result\nstartz 2\nendz 2\nadductmass 0\n"
                    f"psfun 0\nmzsig {fwhm}\nlinflag {linflag}\n",
                    encoding="utf-8",
                )
                peaks = SimpleNamespace(plen=len(apexes), peaks=[
                    SimpleNamespace(mztab=np.array([[mz[index], 1.0, index]]))
                    for index in apexes
                ])
                fallback = tools.makeconvspecies(
                    np.column_stack((mz, np.ones_like(mz))), peaks,
                    SimpleNamespace(linflag=linflag, mzbins=bin_size, mzsig=fwhm, psfun=0),
                )
                # Write both peaks at once with only one charge state. The
                # native output has one row per peak regardless of charge count.
                np.savetxt(Path(directory, "result_peaks.dat"),
                           np.column_stack((2 * mz[apexes], np.ones(len(apexes)))))
                result = subprocess.run(
                    [str(executable), str(path), "-conv"], cwd=directory,
                    capture_output=True, text=True, timeout=30,
                    env=dict(os.environ, OMP_NUM_THREADS="2"),
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                native = np.fromfile(Path(directory, "result_conv.bin"), dtype=np.float32)
                self.assertEqual(native.size, len(apexes) * len(mz))
                native = native.reshape(len(apexes), len(mz))
                for center, apex, native_profile, python_profile in zip(centers, apexes, native, fallback):
                    for backend, profile in (("native", native_profile), ("Python", python_profile)):
                        with self.subTest(linflag=linflag, backend=backend, center=center):
                            self.assertEqual(np.argmax(profile), apex)
                            half_height = profile[apex] / 2
                            above = np.flatnonzero(profile >= half_height)
                            left, right = above[0], above[-1]
                            left_mz = np.interp(half_height, profile[left - 1:left + 1], mz[left - 1:left + 1])
                            right_mz = np.interp(half_height, profile[[right + 1, right]], mz[[right + 1, right]])
                            self.assertAlmostEqual(right_mz - left_mz, fwhm, delta=0.1)

    def test_variable_gaussian_widths_are_fwhm_in_mz(self):
        import numpy as np

        config = UniDecConfig()
        config.initialize_system_paths()
        executable = Path(os.environ.get("UNIDEC_TEST_EXECUTABLE", config.UniDecPath))
        mz = np.arange(100.0, 201.0, 0.25)
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory, "spectrum")
            project.mkdir()
            path = project / "conf.dat"
            widths = np.where(mz < 150, 2.0, 8.0)
            # The broad peak at 150 must reach targets whose own widths are narrow.
            widths[(abs(mz - 150) <= 5) & (mz != 150)] = 0.25
            np.savetxt(project / "fwhm.dat", widths)
            for linflag in (2, 0):
                for variablepw, center, expected in ((1, 125, 2), (1, 150, 8), (1, 175, 8),
                                                     (2, 125, 2), (2, 150, 6), (2, 175, 6)):
                    with self.subTest(linflag=linflag, variablepw=variablepw, center=center):
                        sticks = (mz == center).astype(np.float32)
                        sticks.tofile(project / "result_grid.bin")
                        np.savetxt(project / "input.dat", np.column_stack((mz, sticks)))
                        np.savetxt(project / "result_peaks.dat", [[2 * center, 1.0]])
                        path.write_text(
                            "input input.dat\noutput result\nstartz 2\nendz 2\nadductmass 0\n"
                            "psfun 0\nmzsig 6\nfwhmfile fwhm.dat\npeakshapeinflate 3\n"
                            f"linflag {linflag}\nvariablepw {variablepw}\n",
                            encoding="utf-8",
                        )
                        result = subprocess.run(
                            [str(executable), str(path), "-conv"],
                            cwd=project if linflag == 2 else directory,
                            capture_output=True, text=True, timeout=30,
                            env=dict(os.environ, OMP_NUM_THREADS="2"),
                        )
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                        profile = np.fromfile(project / "result_conv.bin", dtype=np.float32)
                        self.assertEqual(profile.size, mz.size)
                        self.assertAlmostEqual(profile[mz == center].item(), 1)
                        for sign in (-1, 1):
                            half_mz = center + sign * expected / 2
                            self.assertAlmostEqual(profile[mz == half_mz].item(), 0.5, delta=1e-5)

    def test_zero_iterations_write_blank_outputs_when_no_grid_points_survive(self):
        import numpy as np

        config = UniDecConfig()
        config.initialize_system_paths()
        executable = Path(os.environ.get("UNIDEC_TEST_EXECUTABLE", config.UniDecPath))
        data = np.array([[198, 0], [199, 1], [200, 0], [201, 0.5], [202, 0]])
        # A single charge with the default charge smoothing has no neighbors
        # at either impulse. The fit and grids should still be writable.
        with tempfile.TemporaryDirectory(prefix="unidec-zero-iterations-") as directory:
            root = Path(directory)
            np.savetxt(root / "input.dat", data)
            path = root / "conf.dat"
            path.write_text(
                "input input.dat\noutput result\nnumit 0\nstartz 2\nendz 2\n"
                "adductmass 0\nlinflag 4\nmzsig 6\npoolflag 0\nrawflag 0\n"
                "masslb 390\nmassub 410\nmassbins 1\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [str(executable), str(path)], cwd=root, capture_output=True, text=True,
                timeout=30, env=dict(os.environ, OMP_NUM_THREADS="1"),
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for name, size in (("fitdat", len(data)), ("grid", len(data)), ("massgrid", 20)):
                with self.subTest(output=name):
                    values = np.fromfile(root / f"result_{name}.bin", dtype=np.float32)
                    np.testing.assert_array_equal(values, np.zeros(size))
            mass = np.loadtxt(root / "result_mass.txt")
            np.testing.assert_array_equal(mass[:, 0], np.arange(390, 410))
            np.testing.assert_array_equal(mass[:, 1], np.zeros(20))
            stats = dict(line.split(" = ") for line in (root / "result_error.txt").read_text().splitlines())
            self.assertEqual(int(stats["iterations"]), 0)
            self.assertAlmostEqual(float(stats["error"]), np.sum(data[:, 1] ** 2))

    def test_in_process_wrapper_rejects_invalid_fwhm_arrays(self):
        import numpy as np
        from types import SimpleNamespace
        from unidec.modules.unidecwrapper import Config, py_to_c_input

        config = Config()
        config.variablepw = 1
        python_config = SimpleNamespace(masslist=[])
        data = SimpleNamespace(data2=np.array([[100.0, 1.0], [101.0, 0.5]]), fwhmlist=[0.2])

        with self.assertRaisesRegex(ValueError, "one finite, positive FWHM value"):
            py_to_c_input(data, config, python_config)

        data.fwhmlist = [0.2, 0.3]
        restored = py_to_c_input(data, config, python_config)
        self.assertAlmostEqual(restored.fwhmlist[1], 0.3)


if __name__ == "__main__":
    unittest.main()

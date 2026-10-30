import ctypes
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
            self.assertEqual(bytes(restored.fwhmfile).rstrip(b"\0"), b"fwhm.dat")

    def test_native_hdf5_loader_declares_variable_peak_width_settings(self):
        source = Path(__file__).resolve().parents[1] / "unidec" / "src" / "h5io.c"
        text = source.read_text(encoding="utf-8")

        self.assertIn('"variablepw", config.variablepw', text)
        self.assertIn('"minratio", config.minratio', text)
        input_source = source.with_name("udio.c").read_text(encoding="utf-8")
        self.assertIn('strjoin(config->dataset, "/fwhm", fwhmdat)', input_source)

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

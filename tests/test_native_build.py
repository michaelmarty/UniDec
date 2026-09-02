import ctypes
import platform
import subprocess
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


if __name__ == "__main__":
    unittest.main()

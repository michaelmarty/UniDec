import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.engine import UniDec
from unidec.modules.DoubleDec import DoubleDec, cconv2, make_kernel
from tests._test_support import copy_unidec_example


class TestDoubleDecReconstruction(unittest.TestCase):
    @staticmethod
    def _problem():
        """Two separated impulses blurred by a normalized asymmetric kernel."""
        offsets = np.arange(64, dtype=float)
        axis = 1000.0 + offsets
        kernel_table = np.column_stack((
            axis,
            np.exp(-0.5 * ((offsets - 9) / 2) ** 2)
            + 0.25 * np.exp(-0.5 * ((offsets - 13) / 1.5) ** 2),
        ))
        kernel = make_kernel(kernel_table.copy())
        source = np.zeros(64)
        source[18] = 1.0
        source[39] = 0.6
        observed = cconv2(source, kernel)
        observed /= np.max(observed)
        return axis, kernel_table, kernel, source, observed

    def _assert_reconstruction_contract(self, reconstruction, kernel, source, observed):
        reconstruction = np.asarray(reconstruction, dtype=float)
        reconstruction /= np.max(reconstruction)
        for expected_index in np.flatnonzero(source):
            window = slice(expected_index - 2, expected_index + 3)
            recovered_index = expected_index - 2 + np.argmax(reconstruction[window])
            self.assertLessEqual(abs(recovered_index - expected_index), 1)

        expected_indexes = np.flatnonzero(source)
        expected_ratio = source[expected_indexes[1]] / source[expected_indexes[0]]
        recovered_ratio = reconstruction[expected_indexes[1]] / reconstruction[expected_indexes[0]]
        self.assertAlmostEqual(recovered_ratio, expected_ratio, delta=0.01)

        forward = cconv2(reconstruction, kernel)
        forward /= np.max(forward)
        forward_rmse = np.sqrt(np.mean((forward - observed) ** 2))
        self.assertLess(forward_rmse, 0.03)
        return reconstruction

    def test_python_reconstructs_asymmetric_two_peak_problem(self):
        _, _, kernel, source, observed = self._problem()

        reconstruction = DoubleDec().dd_core(observed.copy(), kernel.copy())

        self._assert_reconstruction_contract(reconstruction, kernel, source, observed)

    def test_native_matches_python_for_equal_sampling_kernel(self):
        with tempfile.TemporaryDirectory(prefix="unidec-doubledec-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "BSA.txt")
            config_dir = root / "BSA_unidecfiles"
            config_dir.mkdir()
            reviewed_config = (
                Path(__file__).parent / "fixtures" / "native_regression" / "BSA_conf.dat"
            )
            shutil.copy2(reviewed_config, config_dir / "BSA_conf.dat")

            engine = UniDec(ignore_args=True, silent=True)
            executable = Path(engine.config.UniDecPath)
            if not executable.is_file():
                self.skipTest(f"Packaged UniDec executable is unavailable: {executable}")
            engine.open_file(spectrum.name, directory, silent=True)
            engine.process_data(silent=True)
            engine.config.doubledec = False
            self.assertEqual(engine.run_unidec(silent=True), 0)
            baseline = np.asarray(engine.data.massdat).copy()

            offsets = np.arange(len(baseline), dtype=float)
            kernel_table = np.column_stack((
                baseline[:, 0],
                np.exp(-0.5 * ((offsets - 9) / 0.75) ** 2)
                + 0.1 * np.exp(-0.5 * ((offsets - 11) / 0.5) ** 2),
            ))
            kernel = make_kernel(kernel_table.copy())
            python_result = DoubleDec().dd_core(baseline[:, 1] / np.max(baseline[:, 1]), kernel.copy())
            kernel_path = root / "kernel.txt"
            np.savetxt(kernel_path, kernel_table)
            engine.config.doubledec = True
            engine.config.kernel = str(kernel_path)

            self.assertEqual(engine.run_unidec(silent=True), 0)
            native_mass = np.asarray(engine.data.massdat)

        np.testing.assert_array_equal(native_mass[:, 0], baseline[:, 0])
        native_result = native_mass[:, 1] / np.max(native_mass[:, 1])
        np.testing.assert_allclose(
            native_result, python_result, rtol=5e-3, atol=5e-4
        )


if __name__ == "__main__":
    unittest.main()

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from unidec.engine import UniDec


class TestRunUniDecFailures(unittest.TestCase):
    def setUp(self):
        self.engine = UniDec(silent=True)
        self.engine.config.procflag = 1
        self.engine.config.doubledec = False

    def test_preprocesses_before_successful_native_call(self):
        self.engine.config.procflag = 0

        with patch.object(self.engine, "process_data") as process_data, \
                patch.object(self.engine, "check_badness", return_value=0), \
                patch.object(self.engine, "export_config"), \
                patch("unidec.engine.ud.unidec_call", return_value=0) as native_call, \
                patch.object(self.engine, "unidec_imports") as imports:
            result = self.engine.run_unidec(silent=True, efficiency=True)

        self.assertEqual(result, 0)
        process_data.assert_called_once_with()
        native_call.assert_called_once_with(self.engine.config, silent=True)
        imports.assert_called_once_with(True)

    def test_invalid_config_prevents_native_call(self):
        with patch.object(self.engine, "check_badness", return_value=1), \
                patch("unidec.engine.ud.unidec_call") as native_call, \
                patch.object(self.engine, "unidec_imports") as imports:
            result = self.engine.run_unidec(silent=True)

        self.assertEqual(result, 1)
        native_call.assert_not_called()
        imports.assert_not_called()

    def test_native_error_is_propagated_and_stale_results_are_cleared(self):
        self.engine.data.massdat = np.ones((2, 2))
        self.engine.data.massgrid = np.ones(2)
        self.engine.data.fitdat = np.ones(2)
        self.engine.data.mzgrid = np.ones((2, 3))
        self.engine.pks.add_peaks([[100.0, 1.0]])

        with patch.object(self.engine, "check_badness", return_value=0), \
                patch.object(self.engine, "export_config"), \
                patch("unidec.engine.ud.unidec_call", return_value=7), \
                patch.object(self.engine, "unidec_imports") as imports:
            result = self.engine.run_unidec(silent=True)

        self.assertEqual(result, 7)
        imports.assert_not_called()
        self.assertEqual(len(self.engine.pks.peaks), 0)
        for value in (self.engine.data.massdat, self.engine.data.massgrid,
                      self.engine.data.fitdat, self.engine.data.mzgrid):
            self.assertEqual(np.asarray(value).size, 0)

    def test_missing_doubledec_kernel_aborts_before_export(self):
        with tempfile.TemporaryDirectory(prefix="unidec-kernel-") as directory:
            self.engine.config.doubledec = True
            self.engine.config.kernel = str(Path(directory, "missing_kernel.txt"))
            with patch.object(self.engine, "check_badness", return_value=0), \
                    patch.object(self.engine, "export_config") as export_config, \
                    patch("unidec.engine.ud.unidec_call") as native_call:
                result = self.engine.run_unidec(silent=True)

        self.assertEqual(result, 1)
        export_config.assert_not_called()
        native_call.assert_not_called()


class TestUniDecOutputImports(unittest.TestCase):
    def _engine_with_paths(self, directory):
        engine = UniDec(silent=True)
        root = Path(directory)
        engine.config.massdatfile = str(root / "mass.txt")
        engine.config.massgridfile = str(root / "massgrid.bin")
        engine.config.fitdatfile = str(root / "fit.bin")
        engine.config.mzgridfile = str(root / "mzgrid.bin")
        engine.config.errorfile = str(root / "missing_error.txt")
        engine.config.outfname = str(root / "output")
        engine.config.startz = 1
        engine.config.endz = 2
        engine.data.data2 = np.array([[100.0, 1.0], [101.0, 3.0]])
        return engine

    def test_empty_and_truncated_mass_outputs_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="unidec-output-") as directory:
            engine = self._engine_with_paths(directory)
            mass_path = Path(engine.config.massdatfile)
            for contents in ("", "100\n", "100 1\n101\n"):
                with self.subTest(contents=contents):
                    mass_path.write_text(contents, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, "at least two columns"):
                        engine.unidec_imports(efficiency=True)
                    self.assertEqual(engine.data.massdat.size, 0)

    def test_single_mass_row_is_preserved_as_a_two_dimensional_table(self):
        with tempfile.TemporaryDirectory(prefix="unidec-output-") as directory:
            engine = self._engine_with_paths(directory)
            np.savetxt(engine.config.massdatfile, [[100.0, 2.5]])

            engine.unidec_imports(efficiency=True)

            np.testing.assert_allclose(engine.data.massdat, [[100.0, 2.5]])
            self.assertEqual(engine.data.massdat.shape, (1, 2))
            self.assertEqual(engine.config.massdatnormtop, 2.5)

    def test_efficiency_mode_does_not_read_binary_grids(self):
        with tempfile.TemporaryDirectory(prefix="unidec-output-") as directory:
            engine = self._engine_with_paths(directory)
            np.savetxt(engine.config.massdatfile, [[100.0, 2.5]])

            with patch("unidec.engine.np.fromfile") as fromfile:
                engine.unidec_imports(efficiency=True)

            fromfile.assert_not_called()

    def test_normal_mode_reads_fit_mass_and_mz_grids(self):
        with tempfile.TemporaryDirectory(prefix="unidec-output-") as directory:
            engine = self._engine_with_paths(directory)
            np.savetxt(engine.config.massdatfile, [[100.0, 2.5]])
            np.asarray([1.0, 2.0], dtype=engine.config.dtype).tofile(engine.config.massgridfile)
            np.asarray([0.5, 0.75], dtype=engine.config.dtype).tofile(engine.config.fitdatfile)
            np.asarray([1.0, 2.0, 3.0, 4.0], dtype=engine.config.dtype).tofile(engine.config.mzgridfile)

            engine.unidec_imports(efficiency=False)

            np.testing.assert_allclose(engine.data.massgrid, [1.0, 2.0])
            np.testing.assert_allclose(engine.data.fitdat, [0.5, 0.75])
            self.assertEqual(np.asarray(engine.data.mzgrid).shape, (4, 3))


if __name__ == "__main__":
    unittest.main()

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import h5py

from unidec.modules.ChromEng import ChromEngine
from unidec.modules.unidecstructure import UniDecConfig


class TestUniChromWorkflow(unittest.TestCase):
    def test_unichrom_settings_round_trip_through_hdf5(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "config.hdf5")
            config = UniDecConfig()
            self.assertEqual(config.UClineardecon, 1)
            self.assertEqual(config.UCtype, 0)
            config.UClineardecon = 0
            config.UCtype = 1
            config.write_hdf5(str(path))
            with h5py.File(path, "r") as hdf:
                self.assertEqual(hdf["config"].attrs["UClineardecon"], 0)
                self.assertEqual(hdf["config"].attrs["UCtype"], 1)
            restored = UniDecConfig()
            restored.read_hdf5(str(path))
            self.assertEqual(restored.UClineardecon, 0)
            self.assertEqual(restored.UCtype, 1)

    def test_chromatogram_folder_gets_retention_time_metadata(self):
        engine = object.__new__(ChromEngine)
        engine.chromdat = Mock()
        engine.chromdat.get_avg_scan.return_value = "spectrum"
        engine.chromdat.get_times_from_scans.return_value = (1.0, 1.4, 1.8)

        self.assertEqual(engine.get_data_from_scans((10, 12)), "spectrum")

        self.assertEqual(engine.attrs["retention_time"], 1.4)

    @patch("unidec.modules.ChromEng.metaunidec_call", return_value=0)
    def test_positive_width_is_written_for_c_dispatch(self, call):
        engine = object.__new__(ChromEngine)
        engine.config = Mock()
        engine.config.dtsig = 2.5
        engine.pks = Mock()
        engine.data = Mock()
        engine.check_badness = Mock(return_value=0)
        engine.update_history = Mock()

        result = engine.run_unidec()

        self.assertEqual(result, 0)
        engine.config.write_hdf5.assert_called_once_with()
        call.assert_called_once_with(engine.config)
        engine.data.import_hdf5.assert_called_once_with()
        engine.update_history.assert_called_once_with()

    @patch("unidec.modules.ChromEng.metaunidec_call", return_value=0)
    def test_zero_width_uses_regular_metaunidec(self, call):
        engine = object.__new__(ChromEngine)
        engine.config = Mock()
        engine.config.dtsig = 0
        engine.pks = Mock()
        engine.data = Mock()
        engine.check_badness = Mock(return_value=0)
        engine.update_history = Mock()

        self.assertEqual(engine.run_unidec(), 0)
        call.assert_called_once_with(engine.config)
        engine.data.import_hdf5.assert_called_once_with()
        engine.update_history.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()

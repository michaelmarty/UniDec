import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import h5py

from unidec.modules.ChromEng import ChromEngine
from unidec.modules.unidecstructure import UniDecConfig


class TestUniChromWorkflow(unittest.TestCase):
    def test_linear_decon_setting_round_trips_through_hdf5(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "config.hdf5")
            config = UniDecConfig()
            self.assertEqual(config.UClineardecon, 1)
            config.UClineardecon = 0
            config.write_hdf5(str(path))
            with h5py.File(path, "r") as hdf:
                self.assertEqual(hdf["config"].attrs["UClineardecon"], 0)
            restored = UniDecConfig()
            restored.read_hdf5(str(path))
            self.assertEqual(restored.UClineardecon, 0)

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

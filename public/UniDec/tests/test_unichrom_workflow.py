import unittest
from unittest.mock import Mock, patch

from unidec.modules.ChromEng import ChromEngine


class TestUniChromWorkflow(unittest.TestCase):
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

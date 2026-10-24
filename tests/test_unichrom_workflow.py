import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, call, patch

import h5py
import numpy as np

from unidec.modules.ChromEng import ChromEngine
from unidec.modules.unidecstructure import UniDecConfig


class TestUniChromWorkflow(unittest.TestCase):
    def test_unichrom_settings_round_trip_through_hdf5(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "config.hdf5")
            config = UniDecConfig()
            self.assertEqual(config.UClineardecon, 1)
            self.assertEqual(config.UCtype, 0)
            self.assertEqual(config.chromatogram_type, 0)
            config.UClineardecon = 0
            config.UCtype = 1
            config.chromatogram_type = 1
            config.write_hdf5(str(path))
            with h5py.File(path, "r") as hdf:
                self.assertEqual(hdf["config"].attrs["UClineardecon"], 0)
                self.assertEqual(hdf["config"].attrs["UCtype"], 1)
                self.assertEqual(hdf["config"].attrs["chromatogram_type"], 1)
            restored = UniDecConfig()
            restored.read_hdf5(str(path))
            self.assertEqual(restored.UClineardecon, 0)
            self.assertEqual(restored.UCtype, 1)
            self.assertEqual(restored.chromatogram_type, 1)

            with h5py.File(path, "a") as hdf:
                del hdf["config"].attrs["chromatogram_type"]
            legacy = UniDecConfig()
            legacy.read_hdf5(str(path))
            self.assertEqual(legacy.chromatogram_type, 0)

    def test_selected_chromatogram_uses_configured_cache(self):
        engine = object.__new__(ChromEngine)
        engine.config = Mock(chromatogram_type=0)
        engine.ticdat = "tic"
        engine.bpcdat = "bpc"

        self.assertEqual(engine.get_selected_chromatogram(), ("tic", "Total Ion Chromatogram"))

        engine.config.chromatogram_type = 1
        self.assertEqual(engine.get_selected_chromatogram(), ("bpc", "Base Peak Chromatogram"))

    def test_peak_detection_uses_selected_bpc(self):
        engine = object.__new__(ChromEngine)
        engine.config = Mock(chromatogram_type=1, chrom_peak_width=2)
        engine.ticdat = np.array([[0.0, 100.0], [1.0, 1.0], [2.0, 100.0]])
        engine.bpcdat = np.array([[0.0, 2.0], [1.0, 8.0], [2.0, 2.0]])

        with patch("unidec.modules.ChromEng.ud.gsmooth", side_effect=lambda data, width: data), \
                patch("unidec.modules.ChromEng.ud.noise_level2", return_value=0.25), \
                patch("unidec.modules.ChromEng.ud.peakdetect_nonlinear", return_value=[]) as detect:
            self.assertEqual(engine.get_chrom_peaks(), ([], []))

        detected = detect.call_args.args[0]
        np.testing.assert_array_equal(detected[:, 0], engine.bpcdat[:, 0])
        np.testing.assert_array_equal(detected[:, 1], [0.0, 1.0, 0.0])
        self.assertEqual(detect.call_args.kwargs, {"window": 2, "threshold": 0.25})

    def test_time_parsing_uses_selected_bpc_axis(self):
        engine = object.__new__(ChromEngine)
        engine.config = Mock(chromatogram_type=1, time_window=2, time_start=None, time_end=None)
        engine.ticdat = np.array([[0.0, 1.0], [1.0, 1.0], [2.0, 1.0]])
        engine.bpcdat = np.array([[10.0, 1.0], [12.0, 1.0], [14.0, 1.0]])
        engine.data = Mock()
        engine.attrs = {"source": "bpc"}
        engine.get_data_from_times = Mock(return_value="spectrum")

        engine.add_regular_times()

        self.assertEqual(
            engine.get_data_from_times.call_args_list,
            [call(10.0, 12.0), call(12.0, 14.0)],
        )

    @patch("unidec.modules.ChromEng.ImporterFactory.create_importer")
    def test_opening_chromatogram_caches_tic_and_bpc(self, create_importer):
        tic = [[0.0, 10.0], [1.0, 20.0]]
        bpc = [[0.0, 4.0], [1.0, 7.0]]
        importer = create_importer.return_value
        importer.get_polarity.return_value = "Positive"
        importer.get_tic.return_value = tic
        importer.get_bpc.return_value = bpc
        importer.scans = [1, 2]

        engine = object.__new__(ChromEngine)
        engine.data = Mock()
        engine.setup_filenames = Mock()
        engine.clear = Mock()
        engine.open = Mock()
        engine.update_history = Mock()

        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory, "sample.raw"))
            self.assertFalse(engine.load_mzml(path, load_hdf5=False))

        create_importer.assert_called_once_with(path)
        importer.get_tic.assert_called_once_with()
        importer.get_bpc.assert_called_once_with()
        np.testing.assert_array_equal(engine.ticdat, tic)
        np.testing.assert_array_equal(engine.bpcdat, bpc)

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

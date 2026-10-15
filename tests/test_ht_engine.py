import unittest
from types import SimpleNamespace

import numpy as np
import scipy.fft as fft

from unidec.modules.HTEng import HTEng


class TestHTDemultiplexing(unittest.TestCase):
    @staticmethod
    def _engine(**config):
        engine = HTEng.__new__(HTEng)
        defaults = dict(
            demultiplexmode="HT", HToutputub=0, HToutputlb=0,
            HTtimepad=0, HTksmooth=0, FTapodize=0, FTflatten=0,
        )
        defaults.update(config)
        engine.config = SimpleNamespace(**defaults)
        engine.dtype = float
        return engine

    def test_ht_delta_kernel_preserves_signal_position_and_scale(self):
        engine = self._engine()
        engine.padindex = 0
        engine.shiftindex = 0
        engine.rollindex = 0
        engine.htkernel = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        engine.fftk = fft.rfft(engine.htkernel).conj()
        data = np.arange(6.0)

        output, returned_data = engine.htdecon(data)

        np.testing.assert_allclose(output, data, atol=1e-12)
        self.assertIs(returned_data, data)

    def test_mrs_pools_matching_cycle_positions(self):
        engine = self._engine()
        engine.fulltime = np.arange(8.0)
        engine.fullscans = np.arange(1, 9)
        engine.setup_mrs(cycleindex=2)
        data = np.arange(1.0, 9.0)

        output, returned_data = engine.mrs_decon(data)

        np.testing.assert_allclose(output, [1 + 3 + 5 + 7, 2 + 4 + 6 + 8])
        np.testing.assert_array_equal(engine.decontime, [0.0, 1.0])
        self.assertIs(returned_data, data)

    def test_ft_returns_known_impulse_spectrum_and_zero_trace(self):
        engine = self._engine()
        engine.decontime = np.arange(5.0)
        impulse = np.zeros(8)
        impulse[0] = 2.0
        output, returned_data = engine.ftdecon(impulse)
        np.testing.assert_allclose(output, np.full(5, 2.0))
        np.testing.assert_array_equal(returned_data, impulse)

        zeros, _ = engine.ftdecon(np.zeros(8))
        np.testing.assert_array_equal(zeros, np.zeros(5))

    def test_setup_ft_constructs_frequency_derived_time_axis(self):
        engine = self._engine()
        engine.fulltime = np.array([0.0, 1 / 60, 2 / 60, 3 / 60])
        engine.fullscans = np.arange(4)
        engine.setup_ft(FTstart=0.0, FTend=1.0, nzp=0)
        np.testing.assert_allclose(engine.decontime, [0.0, 750.0, 1500.0])
        np.testing.assert_array_equal(engine.deconscans, [0, 1, 2])

    def test_run_demultiplex_chops_time_axis_and_unknown_mode_is_identity(self):
        engine = self._engine(HToutputlb=1.5, HToutputub=4.5)
        engine.decontime = np.arange(6.0)
        engine.deconscans = np.arange(6)
        data = np.arange(6.0)
        output, returned_data = engine.run_demultiplex(data, mode="unknown", chop=True)
        np.testing.assert_array_equal(output, [2.0, 3.0, 4.0])
        np.testing.assert_array_equal(engine.decontime, [2.0, 3.0, 4.0])
        np.testing.assert_array_equal(engine.deconscans, [2, 3, 4])
        self.assertIs(returned_data, data)


if __name__ == "__main__":
    unittest.main()

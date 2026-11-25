import unittest
import zlib
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

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

    def test_bgal_ht_peak_matches_single_injection_controls(self):
        """Compare the Figure 2 multiplexed TIC with separately acquired runs."""
        fixture = Path(__file__).parent / "fixtures" / "ht_sec" / "bgal_5bit_tic.u8"
        counts = np.fromfile(fixture, dtype=np.uint8).astype(float)
        self.assertEqual(counts.size, 474)
        scan_bin_minutes = 41 / 474
        engine = self._engine(htbit=5, HTtimepad=10, HTtimeshift=10,
                              HTcycleindex=-1, HTksmooth=1)
        engine.fullscans = np.arange(1, counts.size + 1)
        engine.fulltime = engine.fullscans * scan_bin_minutes
        engine.padindex = engine.shiftindex = 0

        engine.setup_demultiplex()
        output, _ = engine.run_demultiplex(counts)

        self.assertEqual(output.shape, counts.shape)
        apex = int(np.argmax(output))
        # The median single-injection apex projects to output index 110.
        # One compressed scan bin covers the measured control apex spread.
        self.assertLessEqual(abs(apex - 110), 1)
        self.assertAlmostEqual(engine.decontime[apex], 9.597744360902256,
                               delta=scan_bin_minutes)

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

    def test_protein_mix_mrs_peak_matches_ensemble_control(self):
        fixture = Path(__file__).parent / "fixtures" / "mrs" / "protein_mix_two_injections.u8.zlib"
        counts = np.frombuffer(zlib.decompress(fixture.read_bytes()), dtype=np.uint8)
        scans_per_run = 588
        self.assertEqual(counts.size, 2 * scans_per_run)

        engine = self._engine(demultiplexmode="MRS")
        engine.fullscans = np.arange(1, counts.size + 1)
        first_scan_minutes = 0.0022250317
        last_scan_minutes = 5.5025247
        scan_period_minutes = (last_scan_minutes - first_scan_minutes) / (scans_per_run - 1)
        engine.fulltime = first_scan_minutes + (engine.fullscans - 1) * scan_period_minutes

        engine.setup_demultiplex(cycleindex=scans_per_run)
        output, _ = engine.run_demultiplex(counts)

        self.assertEqual(output.shape, (scans_per_run,))
        apex = int(np.argmax(output))
        self.assertEqual(apex, 385)
        self.assertAlmostEqual(output[apex], 94.0, places=9)
        # The RAW endpoints approximate the archived time axis within 0.00017 min.
        self.assertAlmostEqual(engine.decontime[apex], 3.6098639, delta=0.00017)
        # Allow one scan interval from each of the CD-MS and ensemble acquisitions.
        self.assertAlmostEqual(engine.decontime[apex], 3.6001669,
                               delta=scan_period_minutes + 0.00809)

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

    def test_gdh_ft_peak_matches_reviewed_4kv_drift_time(self):
        fixture = Path(__file__).parent / "fixtures" / "ft_im" / "gdh_32plus_counts.u8.zlib"
        counts = np.frombuffer(zlib.decompress(fixture.read_bytes()), dtype=np.uint8)
        self.assertEqual(counts.size, 166)
        step_hz = 3
        start_hz = 5
        engine = self._engine(demultiplexmode="FT", FTstart=start_hz,
                              FTend=start_hz + counts.size * step_hz,
                              HTtimepad=1, FTapodize=1, FTflatten=1)
        engine.fullscans = np.arange(1, counts.size + 1)
        engine.fulltime = engine.fullscans * (15 * 0.016044)

        engine.setup_demultiplex()
        output, _ = engine.run_demultiplex(counts)

        self.assertEqual(output.shape, (129,))
        apex = int(np.argmax(output))
        self.assertEqual(apex, 76)  # Direct Fourier sum of the measured trace.
        # Figure 2's 128.5 nm^2 apex corresponds to 98.23 ms at 4 kV.
        native_bin_ms = 1000 / (counts.size * step_hz)
        self.assertAlmostEqual(engine.decontime[apex], 98.23, delta=native_bin_ms)

    def test_absorption_ft_removes_dominant_tone_phase(self):
        """A cosine at bin 6 has DFT coefficient N*A/2 times exp(i*phase)."""
        n = 32
        sample = np.arange(n)
        for phase in (0.0, np.pi / 3, -np.pi / 2):
            for keepcomplex in (False, True):
                with self.subTest(phase=phase, keepcomplex=keepcomplex):
                    engine = self._engine()
                    data = (3 * np.cos(2 * np.pi * 6 * sample / n + phase)
                            + np.cos(2 * np.pi * 9 * sample / n + np.pi / 6))
                    original = data.copy()
                    expected = np.zeros(n // 2 + 1, dtype=complex)
                    expected[6] = 48.0
                    expected[9] = 16 * np.exp(1j * (np.pi / 6 - phase))

                    output, returned_data = engine.ftdecon(data, aFT=True, keepcomplex=keepcomplex)

                    if not keepcomplex:
                        expected = expected.real
                    np.testing.assert_allclose(output, expected, atol=1e-12)
                    np.testing.assert_array_equal(returned_data, original)
                    np.testing.assert_array_equal(data, original)

    def test_absorption_ft_normalization_and_magnitude_have_analytic_amplitudes(self):
        """Phase rotation affects absorption, while magnitude retains both amplitudes."""
        sample = np.arange(32)
        data = (3 * np.cos(2 * np.pi * 6 * sample / 32 + np.pi / 2)
                + np.cos(2 * np.pi * 9 * sample / 32))
        engine = self._engine()
        absorption, _ = engine.ftdecon(data, aFT=True, normalize=True)
        magnitude, _ = engine.ftdecon(data, normalize=True)
        expected_absorption = np.zeros(17)
        expected_absorption[6] = 1.0
        expected_magnitude = expected_absorption.copy()
        expected_magnitude[9] = 1 / 3
        np.testing.assert_allclose(absorption, expected_absorption, atol=1e-12)
        np.testing.assert_allclose(magnitude, expected_magnitude, atol=1e-12)

    def test_masked_ht_matches_direct_circular_correlations(self):
        """Prescribed masks remove known bits; direct dot products supply the oracle."""
        # The 3-bit sequence 1110100 is sampled once every two scans. Each
        # masked bit contributes zero rather than its original +1/-1 weight.
        signed = np.array([1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0])
        data = np.array([2.1, 0.3, 1.7, 4.6, 0.2, 3.9, 1.3, 0.1,
                         5.7, 2.2, 0.8, 1.4, 3.6, 0.4])

        def correlate(mask):
            kernel = np.zeros(14)
            weights = signed.copy()
            weights[list(mask)] = 0
            # These masks preserve the first positive bit, so the kernel
            # origin is unchanged and padding/shift conventions are irrelevant.
            kernel[::2] = weights
            return np.array([sum(data[(i + j) % 14] * kernel[j] for j in range(14))
                             for i in range(14)])

        cases = (
            ("rand", [np.array([1, 3]), np.array([6, 2])], ((1, 3), (6, 2))),
            ("con", [5, 2], ((5, 6), (2, 3))),
            ("alt", [1, 2], ((1, 3), (2, 4))),
        )
        for mode, draws, masks in cases:
            with self.subTest(mode=mode):
                engine = self._engine(htbit=3, HTtimeshift=-1, HTcycleindex=0)
                engine.fulltime = np.arange(14.0)
                engine.fullscans = np.arange(14)
                engine.padindex = engine.shiftindex = 0
                baseline = correlate(())
                # Keep sign decisions away from zero, where FFT roundoff can
                # change a strict-negative classification.
                for mask in masks:
                    self.assertGreater(np.min(np.abs(correlate(mask))), 0.01)
                negative_counts = sum((correlate(mask) < 0).astype(int) for mask in masks)
                expected = baseline / (1 + negative_counts)

                with patch("numpy.random.randint", side_effect=draws) as random_draw:
                    output, returned_data = engine.masked_demultiplex_ht(data, win=2, n=2, mode=mode)

                self.assertEqual(random_draw.call_count, 2)
                np.testing.assert_allclose(output, expected, atol=1e-12)
                np.testing.assert_array_equal(returned_data, data)

    def test_masked_ht_uses_literal_negative_classification(self):
        """Exact zero is neutral, but every representable negative increments the divisor."""
        engine = self._engine(htbit=3, HTtimeshift=-1, HTcycleindex=0)
        engine.fulltime = np.arange(4.0)
        engine.fullscans = np.arange(4)
        tiny_negative = np.nextafter(0.0, -1.0)
        transforms = (
            (np.full(4, 8.0), np.arange(4.0)),
            (np.array([0.0, tiny_negative, 0.0, -1.0]), np.arange(4.0)),
            (np.array([0.0, 0.0, tiny_negative, 1.0]), np.arange(4.0)),
        )

        with patch("numpy.random.randint", return_value=0), \
                patch.object(engine, "setup_ht"), patch.object(
                    engine, "htdecon", side_effect=transforms):
            output, returned_data = engine.masked_demultiplex_ht(
                np.arange(4.0), win=0, n=2, mode="con")

        np.testing.assert_array_equal(output, [8.0, 4.0, 4.0, 4.0])
        np.testing.assert_array_equal(returned_data, np.arange(4.0))

    def test_masked_ht_wraparound_selects_circular_sequence_positions(self):
        """Consecutive and alternating windows wrap on the seven-bit sequence."""
        original = np.array([1, 1, 1, 0, 1, 0, 0])
        cases = (
            ("con", 5, 3, (5, 6, 0)),
            ("alt", 6, 2, (6, 1)),
        )
        for mode, start, win, masked_indexes in cases:
            with self.subTest(mode=mode):
                engine = self._engine(htbit=3, HTtimeshift=-1, HTcycleindex=0)
                engine.fulltime = np.arange(7.0)
                engine.fullscans = np.arange(7)
                sequences = []

                def record_sequence(seq):
                    sequences.append(np.asarray(seq).copy())

                with patch("numpy.random.randint", return_value=start), \
                        patch.object(engine, "setup_ht", side_effect=record_sequence), \
                        patch.object(engine, "htdecon", return_value=(np.ones(7), np.ones(7))):
                    engine.masked_demultiplex_ht(np.ones(7), win=win, n=1, mode=mode)

                expected = original.copy()
                expected[list(masked_indexes)] = -1
                np.testing.assert_array_equal(sequences, [original, expected])

    def test_ht_padding_shift_and_masked_kernel_origin_match_direct_correlation(self):
        """Padding and a masked first bit preserve the explicitly shifted time origin."""
        engine = self._engine(
            htbit=3, HTtimeshift=1.0, HTcycleindex=0, HTtimepad=2.0,
        )
        engine.fulltime = np.arange(14.0)
        engine.fullscans = np.arange(14)
        engine.padindex = engine.shiftindex = 0
        engine.setup_ht(seq=np.array([-1, 1, 1, 0, 1, 0, 0]))
        data = np.array([0.2, 1.1, 0.4, 2.3, 1.7, 0.8, 3.2,
                         0.6, 1.4, 2.8, 0.1, 1.9, 0.7, 2.1])

        output, returned_data = engine.htdecon(data)

        segment = data[1:13]
        direct = np.array([
            sum(segment[(i + j) % len(segment)] * engine.htkernel[j]
                for j in range(len(segment)))
            for i in range(len(segment))
        ])
        expected = np.roll(np.concatenate((np.zeros(2), direct)), -3)
        self.assertEqual(engine.padindex, 2)
        self.assertEqual(engine.shiftindex, 1)
        self.assertEqual(engine.kernelroll, -1)
        self.assertEqual(engine.rollindex, -3)
        np.testing.assert_allclose(output, expected, atol=1e-12)
        self.assertIs(returned_data, data)

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

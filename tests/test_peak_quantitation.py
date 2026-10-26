import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from unittest.mock import Mock

import numpy as np

from unidec import tools
from unidec.engine import UniDec
from unidec.modules import MassFitter, peakstructure


class TestPeakIntegrationAndNormalization(unittest.TestCase):
    def test_rectangular_areas_and_overlapping_windows(self):
        data = np.column_stack((np.arange(7.0), np.ones(7)))
        peaks = peakstructure.Peaks()
        peaks.add_peaks([[2.0, 1.0], [4.0, 1.0]])
        peaks.integrate(data, lb=-1.0, ub=2.0)
        np.testing.assert_array_equal(peaks.areas, [2.0, 2.0])
        self.assertEqual(peaks.peaks[0].integralrange, [1.0, 4.0])
        self.assertEqual(peaks.peaks[1].integralrange, [3.0, 6.0])

    def test_all_peak_normalization_modes(self):
        expected = {
            0: ([1.0, 2.0], [2.0, 4.0]),
            1: ([50.0, 100.0], [50.0, 100.0]),
            2: ([100.0 / 3.0, 200.0 / 3.0], [100.0 / 3.0, 200.0 / 3.0]),
        }
        for mode, (heights, integrals) in expected.items():
            with self.subTest(mode=mode):
                engine = UniDec(silent=True)
                engine.pks = peakstructure.Peaks()
                engine.pks.add_peaks([[100.0, 1.0], [200.0, 2.0]])
                engine.pks.peaks[0].integral = 2.0
                engine.pks.peaks[1].integral = 4.0
                engine.config.peaknorm = mode
                engine.normalize_peaks()
                np.testing.assert_allclose([p.height for p in engine.pks.peaks], heights)
                np.testing.assert_allclose([p.integral for p in engine.pks.peaks], integrals)


class TestFwhmAndFitting(unittest.TestCase):
    def test_symmetric_triangular_peak_has_known_fwhm_and_centroid(self):
        data = np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 4.0], [3.0, 2.0], [4.0, 0.0]])
        peaks = peakstructure.Peaks()
        peaks.add_peaks([[2.0, 4.0]])
        tools.peaks_error_FWHM(peaks, data)
        peak = peaks.peaks[0]
        self.assertEqual(peak.errorFWHM, 2.0)
        self.assertEqual(peak.intervalFWHM, [1.0, 3.0])
        self.assertEqual(peak.centroid, 2.0)
        self.assertFalse(peak.badFWHM)

    def test_two_gaussian_fit_recovers_centers_widths_and_areas(self):
        x = np.linspace(-8.0, 8.0, 321)
        y = (
            tools.make_peak_shape(x, 0, 1.0, -2.0, norm_area=True) * 3.0
            + tools.make_peak_shape(x, 0, 1.5, 2.0, norm_area=True) * 2.0
        )
        fitter = MassFitter.MassFitter(
            np.column_stack((x, y)), np.array([[-1.8, 3.0], [2.2, 2.0]]), 0, "smallguess"
        )
        fit_data, parameters = fitter.perform_fit("sort")
        np.testing.assert_allclose(parameters, [[-2.0, 1.0, 3.0], [2.0, 1.5, 2.0]], atol=1e-5)
        np.testing.assert_allclose(fit_data, y, atol=1e-6)

    def test_charge_state_mass_uncertainty_has_known_weighted_standard_deviation(self):
        massdat = np.column_stack((np.arange(8.0, 13.0), np.ones(5)))
        charge_grid = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 0.0],
                [0.0, 3.0],
                [0.0, 0.0],
            ]
        )
        peaks = peakstructure.Peaks()
        peaks.add_peaks([[10.0, 1.0]])
        peaks.peaks[0].errorFWHM = 1.0

        tools.peaks_error_mean(
            peaks,
            charge_grid.ravel(),
            np.array([1.0, 2.0]),
            massdat,
            SimpleNamespace(peakwindow=10.0),
        )

        peak = peaks.peaks[0]
        np.testing.assert_array_equal(peak.errormasses, [9.0, 11.0])
        np.testing.assert_array_equal(peak.errorints, [1.0, 3.0])
        self.assertAlmostEqual(peak.errormean, np.sqrt(0.75), places=12)


class TestPeakFilteringAndExport(unittest.TestCase):
    def test_filter_cutoff_is_strict(self):
        engine = UniDec(silent=True)
        engine.pks = peakstructure.Peaks()
        engine.pks.add_peaks([[100.0, 1.0], [200.0, 2.0], [300.0, 3.0]])
        for peak, score in zip(engine.pks.peaks, [0.39, 0.40, 0.41]):
            peak.dscore = score
        engine.pick_peaks = Mock()
        engine.dscore = Mock()
        engine.setup_peaks = Mock()
        engine.filter_peaks(minscore=0.40)
        np.testing.assert_array_equal(engine.setup_peaks.call_args.args[0], [[300.0, 3.0]])

    def test_peak_parameter_export_columns_and_values(self):
        with tempfile.TemporaryDirectory() as directory:
            engine = UniDec(silent=True)
            engine.config.outfname = str(Path(directory, "result"))
            engine.config.batchflag = 1
            engine.config.peaknorm = 0
            engine.config.integratelb = -1.0
            engine.config.integrateub = 2.0
            engine.config.mzsig = 0.5
            engine.data.massdat = np.column_stack((np.arange(7.0), np.ones(7)))
            engine.data.ztab = np.array([1.0, 2.0])
            engine.pks = peakstructure.Peaks()
            engine.pks.add_peaks([[2.0, 3.0]])
            peak = engine.pks.peaks[0]
            peak.mztab = np.array([[10.0, 1.0], [11.0, 3.0]])
            peak.centroid = 2.0
            peak.errorFWHM = 2.0
            peak.errormean = 0.25
            peak.dscore = 0.8

            engine.export_params(silent=True)

            exported = np.atleast_2d(np.loadtxt(engine.config.outfname + "_peakparam.dat"))
            self.assertEqual(exported.shape, (1, 10))
            expected = [2.0, 0.875, 1.75, np.sqrt(0.1875), 3.0, 2.0, 2.0, 2.0, 0.25, 0.8]
            np.testing.assert_allclose(exported[0], expected, rtol=1e-5, atol=1e-5)

    def test_postfit_export_uses_fit_width_area_fraction_and_replicate_error(self):
        with tempfile.TemporaryDirectory() as directory:
            engine = UniDec(silent=True)
            engine.config.outfname = str(Path(directory, "postfit"))
            engine.config.batchflag = 1
            engine.config.peaknorm = 0
            engine.config.integratelb = -1.0
            engine.config.integrateub = 1.0
            engine.config.mzsig = 0.5
            engine.data.massdat = np.column_stack((np.arange(5.0), np.ones(5)))
            engine.data.ztab = np.array([1.0, 2.0])
            engine.pks = peakstructure.Peaks()
            engine.pks.add_peaks([[2.0, 3.0]])
            peak = engine.pks.peaks[0]
            peak.mztab = np.array([[10.0, 1.0], [11.0, 3.0]])
            peak.centroid = 2.0
            peak.errorFWHM = 1.5
            peak.errorreplicate = 0.4
            peak.dscore = 0.8
            engine.massfit = np.array([[2.0, 1.2, 5.0]])

            engine.export_params(e="PostFit", silent=True)

            exported = np.atleast_2d(np.loadtxt(engine.config.outfname + "_peakparam.dat"))
            self.assertEqual(exported.shape, (1, 11))
            expected = [
                2.0,
                0.875,
                1.75,
                np.sqrt(0.1875),
                1.0,
                1.2,
                1.0,
                2.0,
                1.5,
                0.4,
                0.8,
            ]
            np.testing.assert_allclose(exported[0], expected, rtol=1e-5, atol=1e-5)


class TestPeakScoreComponents(unittest.TestCase):
    def test_mzstack_uses_one_dimensional_fit_data_for_uniqueness_score(self):
        engine = UniDec(silent=True)
        engine.config = SimpleNamespace(massbins=1.0, adductmass=0.0, orbimode=0)
        engine.data = SimpleNamespace(
            data2=np.array([[4.0, 9.0], [5.0, 9.0], [6.0, 9.0]]),
            fitdat=np.array([1.0, 2.0, 3.0]),
            ztab=np.array([1.0]),
            mzgrid=np.array([[4.0, 1.0, 1.0], [5.0, 1.0, 2.0], [6.0, 1.0, 3.0]]),
        )
        engine.pks = peakstructure.Peaks()
        engine.pks.add_peaks([[5.0, 2.0]])
        peak = engine.pks.peaks[0]
        peak.intervalFWHM = [4.0, 6.0]

        with patch.object(engine, "get_mzstack", wraps=engine.get_mzstack) as get_mzstack:
            engine.pks_uscore(xfwhm=1)

        get_mzstack.assert_called_once_with(xfwhm=1, fitdat=True)
        np.testing.assert_allclose(np.asarray(peak.mzstack[0][:, 1], dtype=float), engine.data.fitdat)
        self.assertEqual(peak.uscore, 1.0)

        engine.pks_uscore(xfwhm=1, uscore_fitdat=False)
        np.testing.assert_allclose(np.asarray(peak.mzstack[0][:, 1], dtype=float), engine.data.data2[:, 1])

    def test_hand_computable_score_components_and_combined_score(self):
        engine = UniDec(silent=True)
        engine.config = SimpleNamespace(
            massbins=1.0,
            adductmass=1.0,
            orbimode=0,
            error=0.8,
        )
        engine.data = SimpleNamespace(
            ztab=np.array([1.0, 2.0]),
            massdat=np.array([[9.0, 0.5], [10.0, 2.0], [11.0, 0.5]]),
        )
        engine.pks = peakstructure.Peaks()
        engine.pks.add_peaks([[10.0, 2.0]])
        peak = engine.pks.peaks[0]
        peak.intervalFWHM = [9.0, 11.0]
        peak.badFWHM = False
        peak.zstack = np.array(
            [
                [9.0, 10.0],
                [1.0, 1.0],
                [2.0, 2.0],
            ]
        )
        peak.mzstack = np.array(
            [
                np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]]),
                np.array([[1.0, 2.0, 1.0], [2.0, 2.0, 1.0]]),
            ],
            dtype=object,
        )

        with (
            patch.object(engine, "get_zstack"),
            patch.object(engine, "get_mzstack"),
        ):
            engine.pks_mscore()
            engine.pks_uscore(power=1)
            engine.pks_csscore()
            engine.pks_fscore()
            self.assertEqual(peak.mscore, 1.0)
            self.assertEqual(peak.uscore, 0.75)
            self.assertEqual(peak.cs_score, 1.0)
            self.assertEqual(peak.fscore, 1.0)

            engine.dscore(power=1)

        self.assertEqual(peak.dscore, 0.75)
        self.assertEqual(peak.lscore, 0.75)
        self.assertAlmostEqual(engine.pks.uniscore, 0.6, places=12)

    def test_asymmetric_bad_fwhm_applies_minimum_penalty(self):
        engine = UniDec(silent=True)
        engine.config = SimpleNamespace(massbins=1.0)
        engine.data = SimpleNamespace(
            massdat=np.array(
                [
                    [6.0, 1.5],
                    [7.0, 1.5],
                    [8.0, 1.8],
                    [9.0, 1.9],
                    [10.0, 2.0],
                    [11.0, 1.0],
                ]
            )
        )
        engine.pks = peakstructure.Peaks()
        engine.pks.add_peaks([[10.0, 2.0]])
        peak = engine.pks.peaks[0]
        peak.intervalFWHM = [7.0, 11.0]
        peak.badFWHM = True

        engine.pks_fscore()

        self.assertEqual(peak.fscore, 0.5)


if __name__ == "__main__":
    unittest.main()

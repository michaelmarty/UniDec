import unittest
from types import SimpleNamespace

import numpy as np

from unidec.modules import CDEng


class TestCDCalibrationAndFiltering(unittest.TestCase):
    @staticmethod
    def _engine(**config):
        engine = CDEng.UniDecCD.__new__(CDEng.UniDecCD)
        defaults = dict(CDslope=10.0, subtype=1)
        defaults.update(config)
        engine.config = SimpleNamespace(**defaults)
        engine.noise = 2.0
        return engine

    def test_linear_and_quadratic_calibration_models(self):
        x = np.array([0.0, 1.0, 2.0])
        np.testing.assert_array_equal(CDEng.slopefunc(x, 3.0), [0.0, 3.0, 6.0])
        np.testing.assert_array_equal(CDEng.linearmodel(x, 3.0, 2.0), [2.0, 5.0, 8.0])
        np.testing.assert_array_equal(CDEng.quadraticmodel(x, 2.0, 3.0, 4.0), [4.0, 9.0, 18.0])

    def test_charge_conversion_positive_and_noise_scaled_modes(self):
        engine = self._engine()
        engine.farray = np.array([[100.0, 10.0], [101.0, 20.0], [102.0, 40.0]])
        engine.convert()
        np.testing.assert_allclose(engine.zarray, [1.0, 2.0, 4.0])

        engine.config.CDslope = -5.0
        engine.convert()
        np.testing.assert_allclose(engine.zarray, [1.0, 2.0, 4.0])
        np.testing.assert_allclose(engine.simp_convert(np.array([10.0, 20.0])), [1.0, 2.0])

    def test_mz_and_charge_filter_boundaries_are_exclusive(self):
        engine = self._engine()
        engine.farray = np.array([
            [100.0, 10.0], [101.0, 20.0], [102.0, 30.0], [103.0, 40.0]
        ])
        engine.filter_mz([100.0, 103.0])
        np.testing.assert_array_equal(engine.farray[:, 0], [101.0, 102.0])

        engine.farray = np.array([[100.0, 10.0], [101.0, 20.0], [102.0, 30.0]])
        engine.filter_z([1.0, 3.0])
        np.testing.assert_array_equal(engine.farray, [[101.0, 20.0]])


class TestCDHistogramAndTransform(unittest.TestCase):
    @staticmethod
    def _engine(**config):
        engine = CDEng.UniDecCD.__new__(CDEng.UniDecCD)
        defaults = dict(
            minmz=-1, maxmz=-1, startz=-1, endz=-1, datanorm=0,
            CDiitflag=False, adductmass=1.0, massbins=10.0,
            masslb=0.0, massub=1000.0, poolflag=0,
        )
        defaults.update(config)
        engine.config = SimpleNamespace(**defaults)
        engine.data = SimpleNamespace()
        engine.invinjtime = None
        return engine

    def test_tiny_histogram_counts_bin_edges_and_unequal_axes(self):
        engine = self._engine()
        x = np.array([100.0, 100.5, 101.0, 102.0, 102.5])
        y = np.array([1.0, 1.0, 2.0, 2.0, 2.5])
        engine.farray = np.column_stack((x, np.ones(len(x)), np.zeros((len(x), 3))))
        engine.zarray = y

        engine.histogram(
            mzbins=1.0, zbins=1.0, x=x, y=y,
            mzrange=[100.0, 103.0], zrange=[1.0, 3.0],
        )

        self.assertEqual(engine.harray.shape, (2, 3))
        np.testing.assert_array_equal(engine.mz, [100.0, 101.0, 102.0])
        np.testing.assert_array_equal(engine.ztab, [1.0, 2.0])
        np.testing.assert_array_equal(engine.harray, [[1.0, 1.0, 0.0], [0.0, 1.0, 2.0]])

    def test_histogram_empty_selection_returns_empty_marker(self):
        engine = self._engine()
        engine.farray = np.empty((0, 5))
        engine.zarray = np.array([])
        self.assertEqual(engine.histogram(x=np.array([]), y=np.array([])), 0)
        self.assertEqual(engine.harray, [])

    def test_mass_transform_uses_mz_minus_adduct_times_charge(self):
        engine = self._engine()
        engine.mz = np.array([11.0, 21.0])
        engine.ztab = np.array([1.0, 2.0])
        engine.X, engine.Y = np.meshgrid(engine.mz, engine.ztab, indexing="xy")
        engine.mass = (engine.X - engine.config.adductmass) * engine.Y
        np.testing.assert_array_equal(engine.mass, [[10.0, 20.0], [20.0, 40.0]])
        harray = np.array([[1.0, 2.0], [3.0, 4.0]])

        result = engine.transform(harray=harray)

        expected_massdat = np.array([
            [0.0, 0.0], [10.0, 1.0], [20.0, 5.0], [30.0, 0.0],
            [40.0, 4.0], [50.0, 0.0], [60.0, 0.0],
        ])
        np.testing.assert_array_equal(result.massdat, expected_massdat)
        self.assertEqual(engine.config.massdatnormtop, 5.0)


class TestCDRichardsonLucy(unittest.TestCase):
    @staticmethod
    def _engine(observed, **config):
        engine = CDEng.UniDecCD.__new__(CDEng.UniDecCD)
        defaults = dict(
            psig=0, beta=0, zzsig=0, msig=0, mzsig=1, csig=0,
            numit=200, datanorm=0, rawflag=1,
        )
        defaults.update(config)
        engine.config = SimpleNamespace(**defaults)
        engine.harray = observed.copy()
        engine.data = SimpleNamespace()
        return engine

    def test_one_iteration_with_asymmetric_blur_matches_hand_calculation(self):
        """A two-tap blur has an explicit forward sum and reversed adjoint."""
        observed = np.array([[2.0, 4.0, 2.0, 0.0, 0.0]])
        engine = self._engine(observed, numit=1)
        engine.kernel = np.array([[0.75, 0.25, 0.0, 0.0, 0.0]])
        engine.ckernel = np.array([[0.75, 0.0, 0.0, 0.0, 0.25]])
        # Forward blur is [3/2, 7/2, 5/2, 1/2, 0]; dividing observations
        # gives [4/3, 8/7, 4/5, 0, 0]. Apply the adjoint and multiply by D.
        expected = np.array([[18 / 7, 148 / 35, 6 / 5, 0.0, 0.0]])

        result = engine.decon_core()

        np.testing.assert_allclose(result, expected, atol=1e-12)
        expected_fit = 0.75 * expected[0] + 0.25 * np.roll(expected[0], 1)
        np.testing.assert_allclose(engine.data.fitdat, expected_fit, atol=1e-12)
        np.testing.assert_array_equal(observed, [[2.0, 4.0, 2.0, 0.0, 0.0]])

    def test_known_blurred_impulse_reconstructs_source_and_forward_fit(self):
        """A unit-sum three-tap blur spreads an eight-count impulse to 2,4,2."""
        source = np.zeros((3, 5))
        source[1, 2] = 8.0
        observed = (0.5 * source + 0.25 * np.roll(source, 1, axis=1)
                    + 0.25 * np.roll(source, -1, axis=1))
        engine = self._engine(observed)
        engine.kernel = np.zeros_like(source)
        engine.kernel[0, [0, 1, -1]] = [0.5, 0.25, 0.25]
        engine.ckernel = engine.kernel.copy()

        result = engine.decon_core()

        # The iterative solver stops at a relative-change threshold of 1e-4;
        # require recovery to within 1% of source amplitude, not exact inversion.
        np.testing.assert_allclose(result, source, atol=0.08, rtol=0)
        self.assertTrue(np.all(np.isfinite(result)))
        self.assertTrue(np.all(result >= 0))
        self.assertAlmostEqual(np.sum(result), 8.0, places=10)
        forward = (0.5 * result + 0.25 * np.roll(result, 1, axis=1)
                   + 0.25 * np.roll(result, -1, axis=1))
        np.testing.assert_allclose(forward, observed, atol=0.08, rtol=0)
        np.testing.assert_allclose(engine.data.fitdat, np.sum(forward, axis=0), atol=1e-12)

    def test_zero_width_kernel_preserves_singleton_axes_and_optional_normalization(self):
        """With blur and regularizers disabled, the operator is the identity."""
        for shape in ((1, 1), (1, 5), (3, 1), (3, 5)):
            for normalize in (0, 1):
                with self.subTest(shape=shape, normalize=normalize):
                    observed = np.arange(1.0, np.prod(shape) + 1).reshape(shape)
                    engine = self._engine(observed, mzsig=0, datanorm=normalize)
                    engine.mz = np.arange(shape[1], dtype=float)
                    engine.ztab = np.arange(shape[0], dtype=float)
                    expected = observed / np.max(observed) if normalize else observed
                    engine.data.data2 = np.column_stack((engine.mz, np.sum(expected, axis=0)))
                    engine.make_kernel()

                    result = engine.decon_core()

                    self.assertEqual(result.shape, shape)
                    np.testing.assert_allclose(result, expected, atol=1e-12)
                    np.testing.assert_allclose(engine.data.fitdat, np.sum(expected, axis=0), atol=1e-12)


if __name__ == "__main__":
    unittest.main()

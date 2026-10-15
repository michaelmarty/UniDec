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


if __name__ == "__main__":
    unittest.main()

import unittest
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np

from unidec.modules.HTEng import UniChromCDEng


class TestUniChromCDStackShapes(unittest.TestCase):
    @staticmethod
    def _engine(grid_shape=(3, 4), stack_shape=None):
        engine = UniChromCDEng.__new__(UniChromCDEng)
        engine.config = SimpleNamespace(adductmass=1.007276467)
        engine.ztab = np.arange(grid_shape[0], dtype=float) + 1
        engine.mz = np.arange(grid_shape[1], dtype=float) + 100
        engine.fulltime = np.array([0.0, 1.0])
        engine.fullscans = np.array([1, 2])
        if stack_shape is None:
            stack_shape = (2,) + grid_shape
        engine.fullhstack = np.zeros(stack_shape)
        engine.topharray = np.zeros(grid_shape)
        engine.X, engine.Y = np.meshgrid(engine.mz, engine.ztab, indexing="xy")
        engine.mass = (engine.X - engine.config.adductmass) * engine.Y
        return engine

    def test_rebuilds_transposed_derived_coordinate_grids(self):
        engine = self._engine(grid_shape=(3, 4))
        engine.X = engine.X.T
        engine.Y = engine.Y.T
        engine.mass = engine.mass.T
        engine.process_data_scans = Mock()

        stack = engine._ensure_full_histogram_stack()

        self.assertIs(stack, engine.fullhstack)
        self.assertEqual(engine.mass.shape, (3, 4))
        np.testing.assert_allclose(engine.mass, (engine.X - engine.config.adductmass) * engine.Y)
        engine.process_data_scans.assert_not_called()

    def test_rebuilds_cached_stack_when_spatial_axes_are_transposed(self):
        engine = self._engine(grid_shape=(3, 4), stack_shape=(2, 4, 3))

        def rebuild():
            engine.fullhstack = np.zeros((2, 3, 4))

        engine.process_data_scans = Mock(side_effect=rebuild)

        stack = engine._ensure_full_histogram_stack()

        self.assertEqual(stack.shape, (2, 3, 4))
        engine.process_data_scans.assert_called_once_with()

    def test_rejects_external_stack_with_wrong_spatial_shape(self):
        engine = self._engine(grid_shape=(3, 4))

        with self.assertRaisesRegex(ValueError, "charge/m/z shape must be"):
            engine._validate_histogram_stack(np.zeros((2, 4, 3)), name="Demultiplexed histogram stack")


if __name__ == "__main__":
    unittest.main()

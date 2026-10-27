import unittest
from types import SimpleNamespace
from unittest.mock import Mock, call

import numpy as np

from unidec.modules.HTEng import UniChromCDEng
from unidec.UniChromCD import UniChromCDApp


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

    def test_transform_stacks_rebuilds_stale_cached_stack(self):
        engine = self._engine(grid_shape=(3, 4), stack_shape=(2, 3, 9))
        engine.massaxis = np.arange(5, dtype=float)
        engine.fullhstack_ht = None
        engine.config.poolflag = 0
        engine.config.datanorm = 0

        def rebuild():
            engine.fullhstack = np.zeros((2, 3, 4))

        engine.process_data_scans = Mock(side_effect=rebuild)
        engine.transform_array = Mock(return_value=(np.zeros((2, 5)), np.zeros((2, 5, 3))))

        engine.transform_stacks()

        engine.process_data_scans.assert_called_once_with()
        self.assertIs(engine.transform_array.call_args.args[0], engine.fullhstack)

    def test_tic_rebuilds_processed_stack_when_histogram_filtering_is_enabled(self):
        for intthresh, reductionpercent in ((1, 0), (0, 50)):
            with self.subTest(intthresh=intthresh, reductionpercent=reductionpercent):
                engine = self._engine(grid_shape=(2, 2))
                engine.config.intthresh = intthresh
                engine.config.reductionpercent = reductionpercent
                engine.config.datanorm = 0
                engine.fullhstack = np.array([
                    [[1.0, 0.0], [2.0, 0.0]],
                    [[0.0, 3.0], [0.0, 4.0]],
                ])
                engine.process_data_scans = Mock()

                tic = engine.get_tic(normalize=engine.config.datanorm)

                engine.process_data_scans.assert_called_once_with()
                np.testing.assert_array_equal(tic, [[0.0, 3.0], [1.0, 7.0]])
                np.testing.assert_array_equal(engine.fulltic, [3.0, 7.0])

    def test_eic_demultiplex_transforms_full_stack_once_for_mass_eics(self):
        tic = SimpleNamespace(label="TIC", massmode=False)
        mass_eic = SimpleNamespace(
            label="Mass: 145000-152000", massmode=True,
            massrange=[145000, 152000], zrange=[20, 55], chromdat=None, decondat=None,
        )
        mz_eic = SimpleNamespace(
            label="m/z: 3000-4000", massmode=False, mzrange=[3000, 4000], zrange=[20, 55],
            sarray=[-1, -1, -1, -1], chromdat=None, decondat=None,
        )
        raw_mass = np.array([[0.0, 1.0]])
        demultiplexed_mass = np.array([[0.0, 2.0]])
        raw_mz = np.array([[0.0, 3.0]])
        demultiplexed_mz = np.array([[0.0, 4.0]])

        app = UniChromCDApp.__new__(UniChromCDApp)
        app.eng = SimpleNamespace(
            config=SimpleNamespace(confname="conf.dat", datanorm=0),
            cc=SimpleNamespace(chromatograms=[tic, mass_eic, mz_eic]),
            run_all_ht=Mock(), transform_stacks=Mock(),
            get_mass_eic=Mock(side_effect=[raw_mass, demultiplexed_mass]),
            eic_ht=Mock(return_value=(demultiplexed_mz, raw_mz)),
        )
        app.export_config = Mock()
        app.plot_chromatograms = Mock()

        app.on_run_eic_ht()

        app.eng.run_all_ht.assert_called_once_with()
        app.eng.transform_stacks.assert_called_once_with()
        self.assertEqual(app.eng.get_mass_eic.call_args_list, [
            call(mass_eic.massrange, mass_eic.zrange),
            call(mass_eic.massrange, mass_eic.zrange, ht=True),
        ])
        app.eng.eic_ht.assert_called_once_with(
            mz_eic.mzrange, mz_eic.zrange, normalize=0, sarray=mz_eic.sarray
        )
        self.assertIs(mass_eic.chromdat, raw_mass)
        self.assertIs(mass_eic.decondat, demultiplexed_mass)
        self.assertIs(mz_eic.chromdat, raw_mz)
        self.assertIs(mz_eic.decondat, demultiplexed_mz)


if __name__ == "__main__":
    unittest.main()

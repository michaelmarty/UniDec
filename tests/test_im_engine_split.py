import unittest
from unittest.mock import patch

import numpy as np

from unidec.engine import UniDec
from unidec.modules.IMEng import UniDec as ModuleUniDecIM
from unidec.modules.IMEng import UniDecIM


class TestEngineModeSplit(unittest.TestCase):
    def test_ms_engine_remains_in_ms_mode_after_reset(self):
        engine = UniDec(ignore_args=True)
        self.assertEqual(engine.config.imflag, 0)
        engine.config.imflag = 1
        engine.reset_config()
        self.assertEqual(engine.config.imflag, 0)

    def test_im_engine_remains_in_im_mode_after_reset(self):
        engine = UniDecIM(ignore_args=True)
        self.assertEqual(engine.config.imflag, 1)
        engine.config.imflag = 0
        engine.reset_config()
        self.assertEqual(engine.config.imflag, 1)

    def test_im_engine_reuses_common_unidec_workflow(self):
        self.assertTrue(issubclass(UniDecIM, UniDec))
        self.assertIs(ModuleUniDecIM, UniDecIM)
        self.assertIs(UniDecIM.run_unidec, UniDec.run_unidec)
        self.assertIs(UniDecIM.pick_peaks, UniDec.pick_peaks)

    def test_im_processing_keeps_the_two_dimensional_shape(self):
        engine = UniDecIM(ignore_args=True)
        engine.data.rawdata = np.array([[100.0, 3.0], [101.0, 7.0]])
        engine.data.rawdata3 = np.array(
            [[100.0, 1.0, 1.0], [100.0, 2.0, 2.0], [101.0, 1.0, 3.0], [101.0, 2.0, 4.0]]
        )
        engine.config.minmz = 100.0
        engine.config.maxmz = 101.0
        engine.config.mindt = 1.0
        engine.config.maxdt = 2.0
        mz = np.array([[100.0, 100.0], [101.0, 101.0]])
        dt = np.array([[1.0, 2.0], [1.0, 2.0]])
        intensity = np.array([[1.0, 2.0], [3.0, 4.0]])

        with (
            patch.object(engine, "export_config"),
            patch.object(engine, "check_badness", return_value=0),
            patch("unidec.modules.IMEng.IM_func.process_data_2d", return_value=(mz, dt, intensity)),
            patch("unidec.modules.IMEng.ud.dataexportbin"),
        ):
            engine.process_data(silent=True)

        np.testing.assert_allclose(
            engine.data.data3,
            np.array(
                [[100.0, 1.0, 1.0], [100.0, 2.0, 2.0], [101.0, 1.0, 3.0], [101.0, 2.0, 4.0]]
            ),
        )
        self.assertEqual(engine.config.procflag, 1)

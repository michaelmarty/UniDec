import os
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import wx

from unidec.GUniDec import UniDecApp
from unidec.UniDecIM import UniDecIMApp
from unidec.modules.unidec_presbase import UniDecPres

from _test_support import copy_unidec_example, find_importer_test_data, has_gui_display


@unittest.skipUnless(has_gui_display(), "wxPython requires a graphical display")
class TestUniDecWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory(prefix="unidec-workflow-")
        with patch.object(UniDecPres, "on_load_default", lambda self, *args, **kwargs: None):
            cls.app = UniDecApp(ignore_args=True)

    @classmethod
    def tearDownClass(cls):
        cls.app.view.Destroy()
        cls.app.wx_app.Yield()
        cls.app.wx_app.Destroy()
        cls.tempdir.cleanup()

    def setUp(self):
        self.app.on_reset(0)
        self.assertEqual(self.app.eng.config.imflag, 0)
        self.assertFalse(hasattr(self.app.view, "plot1im"))

    def test_engine_paths_exist(self):
        config = self.app.eng.config
        self.assertTrue(os.path.isfile(config.UniDecPath))
        self.assertTrue(os.path.isfile(config.cdcreaderpath))

    def test_suppression_cut_percent_is_not_exposed_in_gui(self):
        self.assertFalse(hasattr(self.app.view.controls, "ctlsuppressionpercent"))
        self.app.eng.config.suppression_percent = 0.25
        self.app.import_config()
        self.app.export_config()
        self.assertEqual(self.app.eng.config.suppression_percent, 0.25)

    def test_ms_process_deconvolve_pick_and_restore_state(self):
        spectrum = copy_unidec_example(self.tempdir.name, "ADH.txt")
        self.app.on_open_file(spectrum.name, str(spectrum.parent), clean=True)

        config = self.app.eng.config
        config.startz = 5
        config.endz = 20
        config.masslb = 10000
        config.massub = 200000
        config.mzbins = 1
        self.app.import_config()

        self.app.on_dataprep_button(0)
        self.app.on_unidec_button(0)
        self.app.on_pick_peaks(0)
        self.app.on_replot(0)

        self.assertGreater(len(self.app.eng.data.data2), 0)
        self.assertGreater(len(self.app.eng.data.massdat), 0)
        self.assertGreater(len(self.app.eng.pks.peaks), 0)
        self.assertTrue(np.isfinite(self.app.eng.data.massdat).all())

        state_path = os.path.join(self.tempdir.name, "adh_state.zip")
        self.app.on_save_state(0, state_path)
        self.assertTrue(os.path.isfile(state_path))
        self.app.on_load_state(0, state_path)
        self.assertGreater(len(self.app.eng.data.massdat), 0)


@unittest.skipUnless(has_gui_display(), "wxPython requires a graphical display")
class TestUniDecIMWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory(prefix="unidec-im-workflow-")
        with patch.object(UniDecPres, "on_load_default", lambda self, *args, **kwargs: None):
            cls.app = UniDecIMApp(ignore_args=True)

    @classmethod
    def tearDownClass(cls):
        cls.app.view.Destroy()
        cls.app.wx_app.Yield()
        cls.app.wx_app.Destroy()
        cls.tempdir.cleanup()

    def setUp(self):
        self.app.on_reset(0)
        self.assertEqual(self.app.eng.config.imflag, 1)
        self.assertTrue(hasattr(self.app.view, "plot1im"))

    def test_imms_process_deconvolve_and_pick(self):
        importer_data = find_importer_test_data()
        if importer_data is None:
            if os.environ.get("UNIDEC_REQUIRE_IMPORTER_TEST_DATA") == "1":
                self.fail("UNIDEC_IMPORTER_TEST_DATA does not contain the required importer fixtures")
            self.skipTest("Set UNIDEC_IMPORTER_TEST_DATA to a UniDecImporter TestData checkout")

        source = importer_data / "IMMS" / "test_watersimms_txt.txt"
        if not source.is_file():
            if os.environ.get("UNIDEC_REQUIRE_IMPORTER_TEST_DATA") == "1":
                self.fail(f"Required UniDecImporter IM-MS fixture is missing: {source}")
            self.skipTest(f"UniDecImporter IM-MS fixture is missing: {source}")
        self.assertGreater(source.stat().st_size, 1_000_000, "IM fixture appears to be a Git LFS pointer")

        spectrum = os.path.join(self.tempdir.name, source.name)
        shutil.copy2(source, spectrum)

        self.app.on_open_file(os.path.basename(spectrum), os.path.dirname(spectrum), clean=True)
        config = self.app.eng.config
        config.startz = 10
        config.endz = 18
        config.mzbins = 4
        self.app.import_config()

        self.app.on_dataprep_button(0)
        self.app.on_unidec_button(0)
        self.app.on_pick_peaks(0)

        self.assertGreater(len(self.app.eng.data.data3), 0)
        self.assertGreater(len(self.app.eng.data.massdat), 0)
        self.assertGreater(len(self.app.eng.data.ccsdata), 0)
        self.assertGreater(self.app.eng.data.massccs.size, 0)
        self.assertGreater(len(self.app.eng.pks.peaks), 0)
        self.assertTrue(np.isfinite(self.app.eng.config.error))

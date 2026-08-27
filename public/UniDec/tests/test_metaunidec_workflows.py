import atexit
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import wx

from unidec.MetaUniDec import UniDecApp
from unidec.modules.unidec_presbase import UniDecPres

from _test_support import copy_workflow_spectra, has_gui_display


@unittest.skipUnless(has_gui_display(), "wxPython requires a graphical display")
class TestMetaUniDecWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory(prefix="metaunidec-workflow-")
        with patch.object(UniDecPres, "on_load_default", lambda self, *args, **kwargs: None):
            cls.app = UniDecApp(ignore_args=True)

        cls.spectra = [str(path) for path in copy_workflow_spectra(cls.tempdir.name)]

    @classmethod
    def tearDownClass(cls):
        atexit.unregister(cls.app.repack_hdf5)
        cls.app.view.Destroy()
        cls.app.wx_app.Yield()
        cls.app.wx_app.Destroy()
        cls.tempdir.cleanup()

    def test_create_process_deconvolve_and_pick(self):
        hdf5_path = Path(self.tempdir.name, "workflow.hdf5")
        self.app.new_file(str(hdf5_path))
        self.app.add_files(self.spectra)

        config = self.app.eng.config
        config.startz = 5
        config.endz = 50
        config.masslb = 10000
        config.massub = 250000
        config.mzbins = 1
        self.app.import_config()

        self.app.on_dataprep_button(0)
        self.app.on_unidec_button(0)
        self.app.on_pick_peaks(0)
        self.app.on_replot(0)

        self.assertTrue(hdf5_path.is_file())
        self.assertEqual(len(self.app.eng.data.spectra), 2)
        self.assertGreater(len(self.app.eng.data.massdat), 0)
        self.assertGreater(len(self.app.eng.pks.peaks), 0)
        self.assertTrue(np.isfinite(self.app.eng.data.massdat).all())
        self.app.on_ignore([0])
        self.app.on_isolate([1])
        self.app.on_repopulate()
        self.assertEqual(len(self.app.eng.data.spectra), 2)

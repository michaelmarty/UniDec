import atexit
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import wx

from unidec.MetaUniDec import UniDecApp
from unidec.modules.unidec_presbase import UniDecPres

from _test_support import copy_workflow_spectra, has_gui_display


UNIDEC_ROOT = Path(__file__).resolve().parents[1]


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

    def test_suppression_controls_round_trip_to_hdf5(self):
        controls = self.app.view.controls
        config = self.app.eng.config
        config.suppression_topn = 3
        config.suppression_topx = 0.2
        config.suppression_satellite = 1
        config.suppression_harmonic = 0
        config.suppression_startit = 6
        controls.import_config_to_gui()
        self.assertEqual(controls.ctlsuppressiontopn.GetValue(), "3")
        self.assertEqual(controls.ctlsuppressiontopx.GetValue(), "0.2")
        self.assertEqual(controls.ctlsuppressionsatellite.GetValue(), "1")
        self.assertFalse(controls.ctlsuppressionharmonic.GetValue())
        self.assertEqual(controls.ctlsuppressionstartit.GetValue(), "6")

        controls.ctlsuppressiontopn.SetValue("4")
        controls.ctlsuppressiontopx.SetValue("0.15")
        controls.ctlsuppressionsatellite.SetValue("2")
        controls.ctlsuppressionharmonic.SetValue(True)
        controls.ctlsuppressionstartit.SetValue("7")
        controls.export_gui_to_config()

        self.assertEqual(config.suppression_topn, 4)
        self.assertEqual(config.suppression_topx, 0.15)
        self.assertEqual(config.suppression_satellite, 2)
        self.assertEqual(config.suppression_harmonic, 1)
        self.assertEqual(config.suppression_startit, 7)

        hdf5_path = Path(self.tempdir.name, "suppression-config.hdf5")
        config.write_hdf5(str(hdf5_path))
        with h5py.File(hdf5_path, "r") as hdf:
            attrs = hdf["config"].attrs
            self.assertEqual(attrs["suppression_topn"], 4)
            self.assertEqual(attrs["suppression_topx"], 0.15)
            self.assertEqual(attrs["suppression_satellite"], 2)
            self.assertEqual(attrs["suppression_harmonic"], 1)
            self.assertEqual(attrs["suppression_startit"], 7)

        config.suppression_topn = 0
        config.suppression_topx = 0
        config.suppression_satellite = 0
        config.suppression_harmonic = 0
        config.suppression_startit = 3
        controls.import_config_to_gui()

    def test_c_hdf5_loader_reads_suppression_settings(self):
        source = Path(UNIDEC_ROOT, "unidec", "src", "h5io.c").read_text(encoding="utf-8")
        for setting in (
            "suppression_topn",
            "suppression_topx",
            "suppression_percent",
            "suppression_startit",
            "suppression_harmonic",
            "suppression_satellite",
        ):
            self.assertIn(f'config.{setting} =', source)
            self.assertIn(f'"{setting}"', source)

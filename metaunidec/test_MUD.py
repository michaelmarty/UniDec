import unittest
from mudpres import UniDecApp
import os
import wx
import numpy as np

global topdir, testdir
topdir = os.getcwd()
topdir = os.path.dirname(topdir)
testdir = os.path.join(topdir, "TestSpectra")


class TestGUI(unittest.TestCase):
    @classmethod
    def setUp(cls):
        os.chdir(topdir)
        cls.app = UniDecApp()
        cls.testmzml = os.path.join(testdir, "JAW.mzML")
        cls.testhdf5 = os.path.join(testdir, "JAW.hdf5")
        cls.testnew = os.path.join(testdir, "test.hdf5")
        cls.testnewpaths = [os.path.join(testdir, "0.txt")]

    def tearDown(self):
        wx.CallAfter(self.app.quit_application)
        self.app.start()

    # Engine Tests
    def test_unidec_paths(self):
        self.assertTrue(os.path.isfile(self.app.eng.config.UniDecPath))

    def test_new_file(self):
        self.app.new_file(self.testnew)
        self.app.add_files(self.testnewpaths)
        self.app.add_files(self.testnewpaths)
        self.app.on_delete_spectrum([1])
        self.app.view.on_defaults(self.app.view.menu.menuDefault0)
        self.app.eng.config.startz = 25
        self.app.eng.config.endz = 40
        self.app.eng.config.minmz = 5800
        self.app.eng.config.maxmz = 7900
        self.app.eng.config.mzbins = 10
        self.assertEqual(self.app.eng.data.spectra[0].rawdata.shape, (94579L, 2L))
        self.app.on_auto()

    def test_auto_MS(self):
        os.chdir(testdir)
        if os.path.isfile(self.testhdf5):
            os.remove(self.testhdf5)

        paths = [self.testmzml]
        self.app.import_mzml(paths, 2.0)

        self.app.open_file(self.testhdf5)

        self.app.eng.config.startz = 5
        self.app.eng.config.endz = 12
        self.app.eng.config.massub = 24000
        self.app.eng.config.masslb = 20000
        self.app.import_config()

        self.app.on_dataprep_button()
        self.assertEqual(self.app.eng.out, 0)
        # Test Auto Peak Width
        self.app.on_auto_peak_width()
        self.assertGreater(self.app.eng.config.mzsig, 0.5)
        self.assertLess(self.app.eng.config.mzsig, 2)
        # Test Auto Run
        self.app.on_unidec_button()
        self.assertEqual(self.app.eng.out, 0)
        self.app.on_pick_peaks()
        self.assertEqual(self.app.eng.out, 0)

        intensities = [p.height for p in self.app.eng.pks.peaks]
        masses = [p.mass for p in self.app.eng.pks.peaks]
        maxmass = masses[np.argmax(intensities)]

        self.assertAlmostEqual(int(1000 * round(float(maxmass) / 1000)), 22000, places=0)

        # test replot and related
        self.app.on_replot()

        # Test Plot Grids
        self.app.make2dplots()

        # Test Plot Composite
        self.app.on_plot_composite(0)

        # Test make top
        self.app.make_top(0)
        self.assertEqual(self.app.eng.data.spectra[0].data2[0, 1], self.app.eng.data.data2[0, 1])

        # Test ignore and repopulate
        self.app.on_ignore([0])
        self.app.on_isolate([0])
        self.app.on_repopulate()

        # Test Analysis Scripts
        self.app.on_export_params()

        # Test Exit
        self.app.view.on_exit(0)


    def test_presets_resets(self):
        self.app.view.on_defaults(self.app.view.menu.menuDefault2)
        self.assertEqual(self.app.eng.config.endz, 1)
        self.app.view.on_defaults(self.app.view.menu.menuDefault0)
        self.assertEqual(self.app.eng.config.endz, 100)
        self.app.view.on_defaults(self.app.view.menu.menuDefault3)
        self.assertEqual(self.app.eng.config.endz, 30)
        self.app.view.on_defaults(self.app.view.menu.menuDefault1)
        self.assertEqual(self.app.eng.config.endz, 50)
        self.app.on_save_default(0)
        self.app.on_reset(0)
        self.assertEqual(self.app.eng.config.endz, 100)
        self.app.on_load_default(0)
        self.assertEqual(self.app.eng.config.endz, 50)

    def test_flip(self):
        self.app.on_flip_tabbed()


if __name__ is "__main__":
    unittest.main()

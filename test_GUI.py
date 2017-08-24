import unittest
from GUniDec import UniDecApp
import os
import wx
import numpy as np

global topdir, testdir
topdir = os.getcwd()
testdir = os.path.join(topdir, "TestSpectra")


class TestGUI(unittest.TestCase):
    @classmethod
    def setUp(cls):
        os.chdir(topdir)
        cls.app = UniDecApp()

        cls.testms = "0.txt"
        cls.testimms = "imms.dat"
        cls.testmzml = "myo.mzML"

    def tearDown(self):
        wx.CallAfter(self.app.quit_application)
        self.app.start()

    # Engine Tests
    def test_unidec_paths(self):
        self.assertTrue(os.path.isfile(self.app.eng.config.UniDecPath))
        #self.assertTrue(os.path.isfile(self.app.eng.config.UniDecIMPath))
        self.assertTrue(os.path.isfile(self.app.eng.config.rawreaderpath))
        self.assertTrue(os.path.isfile(self.app.eng.config.cdcreaderpath))

    # Open Files
    def test_mzml(self):
        self.app.on_open_file(self.testmzml, testdir, clean=True)
        self.assertEqual(self.app.eng.data.rawdata.shape, (94540L, 2L))

        with self.assertRaises(IOError):
            self.app.on_open_file("junk.junk", testdir)

    def test_auto_MS(self):
        # Test MS open File
        self.app.on_open_file(self.testms, testdir, clean=True)
        self.assertEqual(self.app.eng.data.rawdata.shape, (94579L, 2L))
        self.app.on_reset(0)
        self.app.eng.config.startz = 25
        self.app.eng.config.endz = 40
        self.app.eng.config.minmz = 5800
        self.app.eng.config.maxmz = 7900
        self.app.import_config()
        self.app.on_dataprep_button(0)
        # self.eng.config.masslb = 211000
        # self.eng.config.massub = 222000

        # Test Auto Peak Width
        self.app.on_auto_peak_width(0)
        self.assertGreater(self.app.eng.config.mzsig,5)
        self.assertLess(self.app.eng.config.mzsig, 15)

        # Test Auto Run
        self.app.on_auto(0)
        self.app.on_plot_peaks(0)
        self.app.on_plot_composite(0)
        self.app.on_replot(0)


        intensities = [p.height for p in self.app.eng.pks.peaks]
        masses = [p.mass for p in self.app.eng.pks.peaks]
        maxmass = masses[np.argmax(intensities)]

        self.assertAlmostEqual(int(1000 * round(float(maxmass) / 1000)), 212000, places=0)

        # Test Analysis Scripts
        self.app.on_integrate()
        self.app.on_export_params(0)
        self.app.on_fit_masses(0)
        self.app.on_plot_offsets(0)
        self.app.on_center_of_mass(0)
        self.app.on_charge_plot(0)



        # Test Save and Load File
        if True:
            sname = "test.zip"
            sname = os.path.join(testdir, self.app.eng.config.outfname + ".zip")
            print sname
            self.app.on_save_state(0, sname)
            self.app.on_load_state(0, sname)

        if True:
            self.app.view.on_save_figure_eps(0)
            self.assertTrue(os.path.isfile("0_Figure1.eps"))
            self.app.view.on_save_figure_png(0)
            self.assertTrue(os.path.isfile("0_Figure1.png"))
            self.app.view.on_save_figure_pdf(0)
            self.assertTrue(os.path.isfile("0_Figure1.pdf"))
            self.app.view.on_save_figure_small(0)
            self.assertTrue(os.path.isfile("0_Thumb_Figure1.pdf"))
            self.app.on_pdf_report(0)
            self.assertTrue(os.path.isfile("0_report.pdf"))

        # Test Exit
        self.app.view.on_exit(0)

    def test_auto_IMMS(self):
        self.app.on_open_file(self.testimms, testdir, clean=True)
        self.assertEqual(self.app.eng.data.rawdata.shape, (1872L, 2L))
        self.assertEqual(self.app.eng.data.rawdata3.shape, (52416L, 3L))

        self.app.eng.config.startz = 10
        self.app.eng.config.endz = 18
        self.app.eng.config.mzbins = 4
        # TODO: Fix so that error is raised or mzbins is fixed to not go below the minimum distance for IM.
        self.app.import_config()
        self.app.on_auto(0)

        intensities = [p.height for p in self.app.eng.pks.peaks]
        masses = [p.mass for p in self.app.eng.pks.peaks]
        maxmass = masses[np.argmax(intensities)]

        self.assertAlmostEqual(maxmass / 10000., 98818 / 10000., places=0)

        self.app.on_plot_peaks(0)
        self.app.on_plot_nativeccs(0)
        self.app.on_plot_composite(0)
        self.app.on_replot(0)
        self.app.make_cube_plot(0)

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

    def test_open_raw_ms(self):
        self.rawfile = os.path.join(testdir,"test_ms.raw")
        self.app.on_raw_open(0, self.rawfile)
        print "Removing:",self.app.eng.config.filename
        print os.getcwd()
        os.chdir("..")
        os.remove(self.app.eng.config.filename)

    def test_open_raw_imms(self):
        self.rawfile = os.path.join(testdir, "test_imms.raw")
        self.app.on_raw_open(0, self.rawfile)
        print "Removing:", self.app.eng.config.filename
        os.chdir("..")
        os.remove(self.app.eng.config.filename)




if __name__ is "__main__":
    unittest.main()

# TODO: Add more tests
'''
General:
Check header test on open file

File:
Get from clipboard
Save Figure As

Tools:
Batch processing modes
Wizard
Peak Width Tool
Manual Assignment
Oligomer and Mass Tools

Analysis:
Native Charge Tools
Data Collector
Kendrick
2D extract
Autocorr
Auto Match

Advanced:
File paths stuff
Mode
Intensity Scale

Experimental:
All

IM:
Parameters
Extraction

'''

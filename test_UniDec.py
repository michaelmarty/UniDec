from unittest import TestCase
from unidec import UniDec
from unidec_modules.unidecstructure import DataContainer
import os
import numpy as np


class TestUniDec(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.eng = UniDec()
        cls.testdir = os.path.join(os.getcwd(), "TestSpectra")
        os.chdir(cls.testdir)

        cls.testms = "0.txt"
        cls.testimms = "imms.dat"
        cls.testmzml = "myo.mzML"

    # Check
    def test_unidec_paths(self):
        """
        Test whether the paths specified actually have the correct files there.
        :return:
        """
        self.assertTrue(os.path.isfile(self.eng.config.UniDecPath))
        self.assertTrue(os.path.isfile(self.eng.config.UniDecIMPath))
        self.assertTrue(os.path.isfile(self.eng.config.rawreaderpath))
        self.assertTrue(os.path.isfile(self.eng.config.cdcreaderpath))

    def test_reset_config(self):
        """
        Check wether a specific value is being reset while another is not.
        :return:
        """
        self.eng.config.imflag = 1
        self.eng.config.detectoreffva = 100
        self.eng.reset_config()
        self.assertEqual(self.eng.config.imflag, 1)
        self.assertNotEqual(self.eng.config.detectoreffva, 100)

    def test_export_import_config(self):
        """
        Write the default config to the default config file.
        Switch a value.
        Import the default back.
        Check if the value has switched.
        :return:
        """
        self.eng.initialize()
        self.eng.export_config(self.eng.config.defaultconfig)
        self.eng.config.mzbins = 1000
        self.eng.load_config(self.eng.config.defaultconfig)
        self.assertNotEqual(self.eng.config.mzbins, 1000)

    #Check
    def test_open_file(self):
        """
        Open several test spectra from various data types in self.testdir.
        Check the that dimensions are the same as expected.
        Check that it fails on a file that isn't there.
        :return:
        """
        self.eng.open_file(self.testms, self.testdir)
        self.assertEqual(self.eng.data.rawdata.shape, (94579L, 2L))

        self.eng.open_file(self.testimms, self.testdir)
        self.assertEqual(self.eng.data.rawdata.shape, (1872L, 2L))
        self.assertEqual(self.eng.data.rawdata3.shape, (52416L, 3L))

        self.eng.open_file(self.testmzml, self.testdir)
        self.assertEqual(self.eng.data.rawdata.shape, (177703L, 2L))

        with self.assertRaises(IOError):
            self.eng.open_file("junk.junk", self.testdir)

    # Check
    def test_auto_MS(self):
        """
        Open a file, process the data, and pick the peaks.
        Check that the maximum peak is what is expected.
        :return:
        """
        self.eng.open_file(self.testms, self.testdir)
        self.eng.config.startz = 25
        self.eng.config.endz = 40
        self.eng.config.minmz = 5800
        self.eng.config.maxmz = 7900
        # self.eng.config.masslb = 211000
        # self.eng.config.massub = 222000
        self.eng.autorun()

        intensities = [p.height for p in self.eng.pks.peaks]
        masses = [p.mass for p in self.eng.pks.peaks]
        maxmass = masses[np.argmax(intensities)]

        self.assertAlmostEqual(int(50 * round(float(maxmass)/50)), 211850, places=0)

    #Check
    def test_auto_IMMS(self):
        """
        Open a file, process the data, and pick the peaks.
        Check that the maximum peak is what is expected.
        :return:
        """
        self.eng.open_file(self.testimms, self.testdir)
        self.eng.config.startz = 10
        self.eng.config.endz = 18
        self.eng.autorun()

        intensities = [p.height for p in self.eng.pks.peaks]
        masses = [p.mass for p in self.eng.pks.peaks]
        maxmass = masses[np.argmax(intensities)]

        self.assertAlmostEqual(maxmass / 10000., 98818 / 10000., places=0)

    def test_save_load(self):
        self.eng.open_file(self.testms, self.testdir)
        self.eng.autorun()
        savepath = os.path.join(self.testdir, self.eng.config.outfname + ".zip")
        massdatlength1 = len(self.eng.data.massdat)
        self.eng.save_state(savepath)
        self.eng.data = DataContainer
        self.eng.load_state(savepath)
        massdatlength2 = len(self.eng.data.massdat)
        self.assertEqual(massdatlength1, massdatlength2)





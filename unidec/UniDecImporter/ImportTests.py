import unittest
from unidec.UniDecImporter.Agilent.AgilentImporter import AgilentImporter
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.UniDecImporter.MZML.mzML import *
from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter
from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
from unidec.UniDecImporter.Waters.Waters import WatersDataImporter

#Waters
#Thermo
#mzxml
#Mzml
#Agilent
class ImportTests(unittest.TestCase):
    def test_waters_create(self):
        path = ("Z:\\Group Share\\JGP\\Lipidomics\\EquiSPLASH_msp_construction\\Negative_data\\Kimber_DDA\\"
                "2024-07-09_08_EquiSPLASH_ESIneg_TrapDDA.raw")
        ""
        importer = ImporterFactory.create_importer(path)
        test = WatersDataImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_thermos_create(self):
        path  = "Z:\\Group Share\\JGP\\js8b05641_si_001\\test.raw"
        importer = ImporterFactory.create_importer(path)
        test = ThermoImporter(path)
        self.assertEqual(type(test), type(importer))


    def test_mzxml_create(self):
        path = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control\\23_04_21_PEPPI_1B_A.mzXML"
        importer = ImporterFactory.create_importer(path)
        test = MZXMLImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_mzml_create(self):
        path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.mzML"
        importer = ImporterFactory.create_importer(path)
        test = MZMLImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_agilent_create(self):
        path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d"
        importer = ImporterFactory.create_importer(path)
        test = AgilentImporter(path)
        self.assertEqual(type(test), type(importer))


    def test_populated_waters(self):
        path = ("Z:\\Group Share\\JGP\\Lipidomics\\EquiSPLASH_msp_construction\\Negative_data\\Kimber_DDA\\"
                "2024-07-09_08_EquiSPLASH_ESIneg_TrapDDA.raw")
        ""
        reader= ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        self.assertTrue(len(spectrum) > 0)


    def test_populated_thermo(self):
        path  = "Z:\\Group Share\\JGP\\js8b05641_si_001\\test.raw"
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.grab_centroid_data(1)
        self.assertTrue(len(spectrum) > 0)


    def test_populated_mzxml(self):
        path = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control\\23_04_21_PEPPI_1B_A.mzXML"
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.grab_scan_data(0)
        self.assertTrue(len(spectrum) > 0)

    def test_populated_mzml(self):
        path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.mzML"
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.grab_scan_data(0)
        self.assertTrue(len(spectrum) > 0)


class ImporterTools(unittest.TestCase):
    #Test to see if anything is returned
    def test_merge_spectra(self):
        path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d"
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        merged = merge_spectra([spectrum])
        self.assertTrue(len(merged) > 0)
    def test_merge_spectra_dim(self):
        path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\AgilentData\\2019_05_15_bsa_ccs_02.d"
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        merged = merge_spectra([spectrum])
        self.assertTrue(merged.shape[1]== 2)

#Pull the right scans
#Averaging is working
#Does grabbing scan at 0 index break

#Get scan_times
#min_max_scans
#XIC-> Tics
#

    def test_merge_im_spectra(self):
        pass








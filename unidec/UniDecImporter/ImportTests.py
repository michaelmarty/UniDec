import unittest
try:
    from unidec.UniDecImporter.Agilent.AgilentImporter import AgilentImporter
except Exception as e:
    print(e)
    print("Unable to import AgilentImporter")
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.UniDecImporter.MZML.mzML import *
from unidec.UniDecImporter.MZXML.mzXML import MZXMLImporter
try:
    from unidec.UniDecImporter.Thermo.Thermo import ThermoImporter
except Exception as e:
    print(e)
    print("Unable to import ThermoImporter")
try:
    from unidec.UniDecImporter.Waters.Waters import WatersDataImporter
except Exception as e:
    print(e)
    print("Unable to import WatersDataImporter")
from unidec.tools import traverse_to_unidec3


#Waters
#Thermo
#mzxml
#Mzml
#Agilent
#Grab the root and changedir to it




class ImportTests(unittest.TestCase):
    def test_waters_create(self):
        path = traverse_to_unidec3("UniDec3")
        # Join the returned path with the new file or directory
        path = os.path.join(path, "TestSpectra", "2024-07-09_08_EquiSPLASH_ESIneg_TrapDDA.raw")
        importer = ImporterFactory.create_importer(path)
        test = WatersDataImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_thermos_create(self):

        base_path = traverse_to_unidec3("UniDec3")
        path = os.path.join(base_path, "TestSpectra", "test.raw")
        importer = ImporterFactory.create_importer(path)
        test = ThermoImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_mzxml_create(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "1500_scans_200K_16 fills-qb1.mzXML")
        importer = ImporterFactory.create_importer(path)
        test = MZXMLImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_mzml_create(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzML")
        importer = ImporterFactory.create_importer(path)
        test = MZMLImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_agilent_create(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "2019_05_15_bsa_ccs_02.d")
        importer = ImporterFactory.create_importer(path)
        test = AgilentImporter(path)
        self.assertEqual(type(test), type(importer))

    def test_populated_waters(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "2024-07-09_08_EquiSPLASH_ESIneg_TrapDDA.raw")
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        self.assertTrue(len(spectrum) > 0)

    def test_populated_thermo(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "test.raw")
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.grab_centroid_data(1)
        self.assertTrue(len(spectrum) > 0)

    def test_populated_mzxml(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "1500_scans_200K_16 fills-qb1.mzXML")
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.grab_scan_data(0)
        self.assertTrue(len(spectrum) > 0)

    def test_populated_mzml(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzML")
        reader = ImporterFactory.create_importer(path)
        #print("Path for mzml: ", path)
        spectrum = reader.grab_scan_data(0)
        self.assertTrue(len(spectrum) > 0)




class MZMLTests(unittest.TestCase):

    def test_merge_spectra_functional(self):
        datalist = np.array([[[1]], [[2]]])
        expected_output = np.array([[1,0]])
        merged = merge_spectra(datalist, type = "Interpolate")
        #(2 + 1) / 2 = 1 for idx 0
        #No values left, so appends a 0 to idx 1
        self.assertTrue(np.array_equal(merged, expected_output))

    def test_merge_spectra_populated(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "2019_05_15_bsa_ccs_02.d")
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        merged = merge_spectra([spectrum])
        self.assertTrue(len(merged) > 0)

    def test_merge_spectra_dim(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "2019_05_15_bsa_ccs_02.d")
        reader = ImporterFactory.create_importer(path)
        spectrum = reader.get_data()
        merged = merge_spectra([spectrum])
        self.assertTrue(merged.shape[1] == 2)

    def test_get_scan_times(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzMl")
        importer = ImporterFactory.create_importer(path)
        test_time = importer.get_scan_times(100)
        self.assertTrue(test_time > 0)

    def test_max_scans(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzMl")
        importer = ImporterFactory.create_importer(path)
        test_max = importer.get_max_scans()
        self.assertTrue(test_max > 0)

    def test_get_scans_from_times(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzMl")
        importer = ImporterFactory.create_importer(path)
        test_scans = importer.get_scans_from_times([0, 150])
        self.assertTrue(len(test_scans) > 0)

    def test_tic_ndim(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzMl")
        importer = ImporterFactory.create_importer(path)
        test_tic = importer.get_tic()
        self.assertTrue(test_tic.ndim > 1)

    def test_tic_populated(self):
        path = traverse_to_unidec3("UniDec3")
        path = os.path.join(path, "TestSpectra", "JAW.mzMl")
        importer = ImporterFactory.create_importer(path)
        test_tic = importer.get_tic()
        self.assertTrue(len(test_tic) > 0)




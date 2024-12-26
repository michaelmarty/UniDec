import unittest
import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import numpy as np

topdir = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection"
chrom_dat_examples = ["test_mzml.mzML", "test_mzmlgz.mzML.gz", "test_mzxml.mzXML",
                      "test_thermo.RAW", "test_waters.raw", "test_agilent.d"]
chrom_paths = [os.path.join(topdir, p) for p in chrom_dat_examples]

ss_dir = "SingleScan"
ss_dat_examples = ["test_csv.csv", "test_dat.dat", "test_txt.txt"]

ss_paths = []
for p in ss_dat_examples:
    ss_paths.append(os.path.join(topdir, ss_dir, p))
for p in chrom_dat_examples:
    ss_paths.append(os.path.join(topdir, p))

cdms_dir = "CDMS"
cdms_dat_examples = ["test_csv_cdms.csv", "test_raw_cdms.raw", "test_npz_cdms.npz", "test_dmt_cdms.dmt"]
cdms_paths = [os.path.join(topdir, cdms_dir, p) for p in cdms_dat_examples]

imms_dir = "IMMS"
imms_dat_examples = ["test_waters.raw", "test_agilentimms_mzml.mzML", "test_watersimms_txt.txt"]
imms_paths = [os.path.join(topdir, imms_dir, p) for p in imms_dat_examples]


class ImporterTests(unittest.TestCase):
    def test_tic(self):
        for p in chrom_paths:
            print("Testing TIC for", p)
            importer = ImporterFactory.create_importer(p)
            test_tic = importer.get_tic()
            print("Shape:", np.shape(test_tic))
            self.assertTrue(np.shape(test_tic)[1] == 2)
            self.assertTrue(len(test_tic) > 20)
            self.assertTrue(importer.chrom_support)

    def test_full_avg_scan(self):
        for p in ss_paths:
            print("Testing Avg Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_avg = importer.get_avg_scan()
            print("Shape:", np.shape(test_avg))
            self.assertTrue(np.shape(test_avg)[1] == 2)
            self.assertTrue(len(test_avg) > 500)
            importer.close()

    def test_avg_scan_scanrange(self):
        for p in ss_paths:
            print("Testing Avg Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_avg = importer.get_avg_scan(scan_range=[1, 5])
            print("Shape:", np.shape(test_avg))
            self.assertTrue(np.shape(test_avg)[1] == 2)
            self.assertTrue(len(test_avg) > 500)

    def test_avg_scan_norange(self):
        for p in ss_paths:
            print("Testing Avg Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_avg = importer.get_avg_scan(scan_range=[2,2])
            print("Shape:", np.shape(test_avg))
            self.assertTrue(np.shape(test_avg)[1] == 2)
            self.assertTrue(len(test_avg) > 500)

    def test_avg_time_norange(self):
        for p in ss_paths:
            print("Testing Avg Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_avg = importer.get_avg_scan(time_range=[0.1,0.1])
            print("Shape:", np.shape(test_avg))
            self.assertTrue(np.shape(test_avg)[1] == 2)
            self.assertTrue(len(test_avg) > 500)

    def test_avg_scan_timerange(self):
        for p in ss_paths:
            print("Testing Avg Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_avg = importer.get_avg_scan(time_range=[0, 0.25])
            print("Shape:", np.shape(test_avg))
            self.assertTrue(np.shape(test_avg)[1] == 2)
            self.assertTrue(len(test_avg) > 500)

    def test_single_scan(self):
        for p in chrom_paths:
            print("Testing Single Scan for", p)
            importer = ImporterFactory.create_importer(p)
            test_single = importer.get_single_scan(1)
            print("Shape:", np.shape(test_single))
            self.assertTrue(np.shape(test_single)[1] == 2)
            self.assertTrue(len(test_single) > 20)

    def test_all_scans(self):
        for p in chrom_paths:
            print("Testing All Scans for", p)
            importer = ImporterFactory.create_importer(p)
            test_all = importer.get_all_scans()
            print("Len:", len(test_all))
            self.assertTrue(type(test_all) == list)
            self.assertTrue(len(test_all) == len(importer.scans))

    def test_get_polarity(self):
        for p in chrom_paths:
            print("Testing Polarity for", p)
            importer = ImporterFactory.create_importer(p)
            test_polarity = importer.get_polarity()
            self.assertTrue(test_polarity in ["Positive", "Negative"])

    def test_get_ms_order(self):
        for p in chrom_paths:
            print("Testing MS Order for", p)
            importer = ImporterFactory.create_importer(p)
            test_order = importer.get_ms_order()
            print("MS Order:", test_order)
            self.assertTrue(test_order >= 1)

    def test_check_centroided(self):
        for p in ss_paths:
            print("Testing Centroided for", p)
            importer = ImporterFactory.create_importer(p)
            test_centroided = importer.check_centroided()
            print("Centroided:", test_centroided)
            self.assertTrue(type(test_centroided) == bool)

    def test_get_scans_times(self):
        for p in chrom_paths:
            print("Testing Max Scans and Times for", p)
            importer = ImporterFactory.create_importer(p)
            max_scan = importer.get_max_scan()
            max_time = importer.get_max_time()
            print("Max Scan:", max_scan)
            print("Max Time:", max_time)
            self.assertTrue(max_scan > 20)
            self.assertTrue(max_time > 0.4)

            max_scan_time = importer.get_scan_time(max_scan)
            print("Max Scan Time:", max_scan_time)
            self.assertEqual(max_scan_time, max_time)

            test_scan = max_scan - 2
            test_scan_range = [test_scan, max_scan]
            print("Test Scan Range:", test_scan_range)
            test_times = importer.get_times_from_scans(test_scan_range)
            test_times = [test_times[0], test_times[2]]
            print("Test Times:", test_times)
            test_scan_return = importer.get_scans_from_times(test_times)
            print("Test Scans:", test_scan_return)
            self.assertTrue(test_scan_return[0] == test_scan)
            self.assertTrue(test_scan_return[1] == max_scan)

    def test_repeat_avg(self):
        for p in chrom_paths:
            print("Testing Repeat Avg for", p)
            importer = ImporterFactory.create_importer(p)
            scan_range = [1, 5]
            test_avg = importer.get_avg_scan(scan_range=scan_range)
            test_avg2 = importer.get_avg_scan(scan_range=scan_range)
            self.assertTrue(np.allclose(test_avg, test_avg2))

    def test_get_cdms(self):
        for p in cdms_paths:
            importer = ImporterFactory.create_importer(p)
            curr_dat = importer.get_cdms_data()
            self.assertTrue(curr_dat.ndim == 2)
            self.assertTrue(len(curr_dat[:, 0]) == len(curr_dat[:, 1]) == len(curr_dat[:, 2]) == len(curr_dat[:, 3]))
            self.assertTrue(curr_dat[:, 0].shape == curr_dat[:, 1].shape == curr_dat[:, 2].shape == curr_dat[:, 3].shape)





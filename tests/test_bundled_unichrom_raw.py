"""Conditional UniChrom coverage for the packaged Thermo RAW example."""

import platform
import tempfile
import unittest

import numpy as np

from UniDecImporter.ImporterFactory import VendorReaderUnavailableError
from tests._test_support import copy_unidec_example
from unidec.modules.ChromEng import ChromEngine


WINDOWS_X64 = platform.system() == "Windows" and platform.machine().lower() in {
    "amd64", "x86_64", "x64",
}


@unittest.skipUnless(WINDOWS_X64, "Packaged Thermo RAW reader requires Windows x64")
class TestBundledUniChromRaw(unittest.TestCase):
    def test_time_selection_uses_inclusive_scans_and_vendor_aggregation(self):
        with tempfile.TemporaryDirectory(prefix="unichrom-raw-") as directory:
            path = copy_unidec_example(
                directory, "UniChrom", "SEC_Native_Bispecific_Special.raw"
            )
            engine = ChromEngine()
            try:
                try:
                    loaded_existing = engine.open_chrom(str(path), refresh=True)
                except VendorReaderUnavailableError as error:
                    self.skipTest(f"Thermo RAW runtime is unavailable: {error}")

                self.assertFalse(loaded_existing)
                scans = np.asarray(engine.chromdat.scans)
                times = np.asarray(engine.chromdat.times)
                tic = np.asarray(engine.ticdat)
                self.assertEqual(tic.shape, (len(scans), 2))
                self.assertEqual(len(times), len(scans))
                self.assertTrue(np.isfinite(tic).all())
                self.assertTrue(np.all(np.diff(scans) > 0))
                self.assertTrue(np.all(np.diff(times) > 0))
                np.testing.assert_allclose(tic[:, 0], times)

                first = len(times) // 3
                last = first + 3
                time_range = [float(times[first]), float(times[last])]
                expected_scans = [int(scans[first]), int(scans[last])]
                spectrum = engine.get_data_from_times(*time_range)

                self.assertEqual(engine.scans[:2], expected_scans)
                np.testing.assert_allclose(engine.scans[2:], time_range)
                self.assertEqual(engine.attrs["scanstart"], expected_scans[0])
                self.assertEqual(engine.attrs["scanend"], expected_scans[1])
                self.assertTrue(np.isfinite(spectrum).all())
                self.assertTrue(np.all(np.diff(spectrum[:, 0]) > 0))

                direct_average = engine.chromdat.get_avg_scan(
                    scan_range=expected_scans
                )
                direct_sum = engine.chromdat.get_avg_scan(
                    scan_range=expected_scans, sum_mode=True
                )
                selected_count = np.count_nonzero(
                    (scans >= expected_scans[0]) & (scans <= expected_scans[1])
                )
                self.assertEqual(selected_count, 4)
                np.testing.assert_array_equal(spectrum, direct_average)
                np.testing.assert_array_equal(direct_sum[:, 0], direct_average[:, 0])
                np.testing.assert_allclose(
                    direct_sum[:, 1], direct_average[:, 1] * selected_count
                )
            finally:
                if engine.chromdat is not None:
                    engine.chromdat.close()


if __name__ == "__main__":
    unittest.main()

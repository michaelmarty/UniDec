import os
import tempfile
import unittest

import numpy as np

from unidec.modules.CDEng import read_uccd_binary, write_uccd_binary


class TestUCCDBinary(unittest.TestCase):
    def test_writer_clips_negative_demultiplexing_sidelobes(self):
        chromaxis = np.array([0.0, 1.0])
        mzaxis = np.array([1000.0, 1100.0])
        zaxis = np.array([10.0, 11.0])
        hstack = np.array([
            [[1.0, -0.25], [0.0, 2.0]],
            [[-0.5, 3.0], [4.0, 0.0]],
        ])

        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "test_uccd.bin")
            write_uccd_binary(filename, chromaxis, mzaxis, zaxis, hstack)
            output, axes = read_uccd_binary(
                filename,
                expected_axes=(chromaxis, mzaxis, zaxis),
            )

        np.testing.assert_array_equal(output, np.clip(hstack, 0, None))
        np.testing.assert_array_equal(axes[0], chromaxis)
        np.testing.assert_array_equal(axes[1], mzaxis)
        np.testing.assert_array_equal(axes[2], zaxis)


if __name__ == "__main__":
    unittest.main()

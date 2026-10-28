import os
import struct
import tempfile
import unittest

import numpy as np

from unidec.modules.CDEng import (
    UCCD_BINARY_HEADER,
    UCCD_BINARY_MAGIC,
    UCCD_BINARY_RECORD,
    read_uccd_binary,
    write_uccd_binary,
)


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

    def test_unequal_axes_round_trip_distinguishable_cells(self):
        chromaxis = np.array([0.0, 2.0])
        mzaxis = np.array([100.0, 200.0, 300.0])
        zaxis = np.array([1.0, 2.0, 3.0, 4.0])
        hstack = np.arange(1.0, 25.0).reshape(2, 4, 3)
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "unequal.bin")
            write_uccd_binary(filename, chromaxis, mzaxis, zaxis, hstack)
            output, axes = read_uccd_binary(filename)
        np.testing.assert_array_equal(output, hstack)
        for actual, expected in zip(axes, (chromaxis, mzaxis, zaxis)):
            np.testing.assert_array_equal(actual, expected)

    def test_all_zero_and_singleton_cubes(self):
        cases = (
            (np.array([0.0, 1.0]), np.array([100.0, 200.0]), np.array([1.0]), np.zeros((2, 1, 2))),
            (np.array([0.0]), np.array([100.0]), np.array([1.0]), np.array([[[7.0]]])),
        )
        for chromaxis, mzaxis, zaxis, hstack in cases:
            with self.subTest(shape=hstack.shape), tempfile.TemporaryDirectory() as directory:
                filename = os.path.join(directory, "cube.bin")
                write_uccd_binary(filename, chromaxis, mzaxis, zaxis, hstack)
                output, _ = read_uccd_binary(filename)
                np.testing.assert_array_equal(output, hstack)

    def test_writer_rejects_wrong_shape(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "must have shape"):
                write_uccd_binary(
                    os.path.join(directory, "wrong.bin"),
                    np.array([0.0, 1.0]), np.array([100.0, 200.0, 300.0]), np.array([1.0, 2.0]),
                    np.zeros((2, 3, 2)),
                )

    def test_sum_mz_matches_independent_axis_sum(self):
        chromaxis = np.array([0.0, 1.0])
        mzaxis = np.array([100.0, 200.0, 300.0])
        zaxis = np.array([1.0, 2.0])
        hstack = np.arange(1.0, 13.0).reshape(2, 2, 3)
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "sum.bin")
            write_uccd_binary(filename, chromaxis, mzaxis, zaxis, hstack)
            output, _ = read_uccd_binary(filename, sum_mz=True)
        np.testing.assert_array_equal(output, np.sum(hstack, axis=(0, 1)))

    def test_incomplete_header_axes_and_sparse_records(self):
        with tempfile.TemporaryDirectory() as directory:
            header_file = os.path.join(directory, "header.bin")
            with open(header_file, "wb") as stream:
                stream.write(b"short")
            with self.assertRaisesRegex(ValueError, "Incomplete UCCD header"):
                read_uccd_binary(header_file)

            axes_file = os.path.join(directory, "axes.bin")
            with open(axes_file, "wb") as stream:
                stream.write(UCCD_BINARY_HEADER.pack(UCCD_BINARY_MAGIC, 1, 1, 1, 0))
                stream.write(struct.pack("<f", 0.0))
            with self.assertRaisesRegex(ValueError, "Incomplete UCCD axes"):
                read_uccd_binary(axes_file)

            sparse_file = os.path.join(directory, "sparse.bin")
            with open(sparse_file, "wb") as stream:
                stream.write(UCCD_BINARY_HEADER.pack(UCCD_BINARY_MAGIC, 1, 1, 1, 1))
                stream.write(struct.pack("<fff", 0.0, 100.0, 1.0))
                stream.write(b"short")
            with self.assertRaisesRegex(ValueError, "Incomplete UCCD sparse data"):
                read_uccd_binary(sparse_file)

    def test_invalid_magic_out_of_range_index_and_expected_axis_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            invalid_file = os.path.join(directory, "invalid.bin")
            with open(invalid_file, "wb") as stream:
                stream.write(UCCD_BINARY_HEADER.pack(b"BADMAGIC", 1, 1, 1, 0))
            with self.assertRaisesRegex(ValueError, "Invalid UCCD binary"):
                read_uccd_binary(invalid_file)

            index_file = os.path.join(directory, "index.bin")
            record = np.array([(1, 2.0)], dtype=UCCD_BINARY_RECORD)
            with open(index_file, "wb") as stream:
                stream.write(UCCD_BINARY_HEADER.pack(UCCD_BINARY_MAGIC, 1, 1, 1, 1))
                stream.write(struct.pack("<fff", 0.0, 100.0, 1.0))
                stream.write(record.tobytes())
            with self.assertRaisesRegex(ValueError, "Out-of-range sparse index"):
                read_uccd_binary(index_file)

            valid_file = os.path.join(directory, "valid.bin")
            axes = (np.array([0.0]), np.array([100.0]), np.array([1.0]))
            write_uccd_binary(valid_file, *axes, np.ones((1, 1, 1)))
            with self.assertRaisesRegex(ValueError, "m/z axis does not match"):
                read_uccd_binary(valid_file, expected_axes=(axes[0], np.array([200.0]), axes[2]))


if __name__ == "__main__":
    unittest.main()

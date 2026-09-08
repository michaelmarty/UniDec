"""Native HDF5 regression against the original charge-cube FFT formulation.

Set UNIDEC_TEST_EXECUTABLE to test an isolated native build. Otherwise use the
packaged executable, as in test_native_build. All inputs/outputs are temporary.
"""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import h5py
import numpy as np


def write_input(path, scans=5, charges=4, **settings):
    config = dict(
        metamode=-1, startz=2, endz=charges + 1, numit=-4,
        dtsig=2.5, mzsig=1.2, psfun=0, rawflag=1, datanorm=0,
        zzsig=0., msig=0., psig=0., beta=0.,
        masslb=-2008., massub=-4040., massbins=1., adductmass=0.,
        nativezlb=-100., nativezub=100., intthresh=0.,
        minmz=0., maxmz=2000., mzbins=0., subbuff=0., reductionpercent=0.,
    )
    if charges == 1:
        config["masslb"] = -100.
    config.update(settings)
    mz = np.arange(1000., 1012., dtype=np.float32)
    with h5py.File(path, "w") as hdf:
        attrs = hdf.create_group("config").attrs
        for name, value in config.items():
            attrs[name] = np.int32(value) if isinstance(value, int) else np.float32(value)
        dataset = hdf.create_group("ms_dataset")
        dataset.attrs["num"] = np.int32(scans)
        for scan in range(scans):
            signal = .1 + np.exp(-((mz - 1002. - scan) / 1.5) ** 2)
            signal += .3 * np.exp(-((mz - 1009.) / 1.2) ** 2)
            dataset.create_group(str(scan)).create_dataset(
                "raw_data", data=np.column_stack((mz, signal)).astype(np.float32))
    return config


def peak(center, axis, width, shape):
    if width == 0:
        return (axis == center).astype(float)
    gaussian = np.exp(-(center - axis) ** 2 / (2 * width ** 2))
    lorentzian = (width / 2) ** 2 / ((center - axis) ** 2 + (width / 2) ** 2)
    if shape == 0:
        return gaussian
    if shape == 1:
        return lorentzian
    return np.where(axis < center,
                    np.exp(-(center - axis) ** 2 / (2 * width ** 2 * .180337)),
                    lorentzian)


def cube_reference(initial, mz, config):
    """Keep the old full 3-D forward/sum/broadcast/adjoint operator as oracle."""
    charges = np.arange(config["startz"], config["endz"] + 1)
    mass = ((mz.astype(np.float32)[:, None] - np.float32(config["adductmass"]))
            * charges.astype(np.float32))
    allowed = (mass > abs(config["masslb"])) & (mass < abs(config["massub"]))
    allowed &= initial.sum(axis=0)[:, None] > 0
    # A zero-iteration native run supplies the preprocessed initialization,
    # avoiding a second implementation of the shared cubic HDF5 merger.
    count = allowed.sum(axis=1)
    per_charge = np.divide(initial, count, out=np.zeros_like(initial), where=count > 0)
    cube = per_charge[:, :, None] * allowed
    observed = per_charge * (len(charges) + 2)
    t = np.arange(len(initial), dtype=float)
    mzsig = config["mzsig"] / (2.35482 if config["psfun"] == 0 else 1)
    dtsig = config["dtsig"] / 2.35482
    km = peak(mz[0], mz, mzsig, config["psfun"])
    km += peak(2 * mz[-1] - mz[-2], mz, mzsig, config["psfun"])
    kt = (np.ones(1) if len(t) == 1 else
          peak(0, t, dtsig, config["psfun"]) + peak(len(t), t, dtsig, config["psfun"]))
    kernel = np.zeros_like(cube)
    kernel[:, :, 0] = kt[:, None] * km
    spectrum = np.fft.fftn(kernel)
    for _ in range(abs(config["numit"])):
        predicted_cube = np.fft.ifftn(np.fft.fftn(cube) * spectrum).real
        predicted = predicted_cube.sum(axis=2)
        ratio = np.divide(observed, predicted, out=np.zeros_like(observed), where=predicted > 0)
        broadcast = np.broadcast_to(ratio[:, :, None], cube.shape)
        cube *= np.fft.ifftn(np.fft.fftn(broadcast) * spectrum.conj()).real
        cube = np.maximum(cube, 0) * allowed
    if config["rawflag"] in (0, 2):
        cube = np.fft.ifftn(np.fft.fftn(cube) * spectrum / kernel.sum()).real
        cube = np.maximum(cube, 0) * allowed
    if config["datanorm"] == 1 and cube.max() > 0:
        cube *= observed.max() / cube.max()
    return cube, mass


class TestUniChromNative(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        executable = os.environ.get("UNIDEC_TEST_EXECUTABLE")
        if not executable:
            from unidec.modules.unidecstructure import UniDecConfig
            config = UniDecConfig()
            config.initialize_system_paths()
            executable = config.UniDecPath
        cls.executable = str(Path(executable).resolve())
        if not Path(cls.executable).is_file():
            raise unittest.SkipTest("Build the native executable or set UNIDEC_TEST_EXECUTABLE")

    def run_native(self, path):
        result = subprocess.run([self.executable, str(path), "-nthreads", "2"],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def test_matches_charge_cube_reference(self):
        cases = [dict(psfun=shape, rawflag=rawflag, datanorm=norm)
                 for shape in range(3) for rawflag in range(4) for norm in range(2)]
        cases += [dict(scans=1), dict(charges=1), dict(mzsig=0.),
                  dict(scans=4, charges=3), dict(dtsig=.05)]
        for case in cases:
            with self.subTest(**case), tempfile.TemporaryDirectory() as directory:
                path = Path(directory, "chrom.hdf5")
                config = write_input(path, **case)
                with h5py.File(path, "a") as hdf:
                    hdf["config"].attrs.modify("numit", np.int32(0))
                    hdf["config"].attrs.modify("rawflag", np.int32(1))
                    # Initialization must precede cube-wide output normalization.
                    hdf["config"].attrs.modify("datanorm", np.int32(0))
                self.run_native(path)
                with h5py.File(path, "r") as hdf:
                    mz = hdf["ms_dataset/mz_axis"][:].astype(float)
                    initial = hdf["ms_dataset/mz_grid"][:].reshape(-1, len(mz)).astype(float)
                # datanorm=1 normalizes each observed scan before iteration.
                if config["datanorm"] == 1:
                    with h5py.File(path, "r") as hdf:
                        maxima = [hdf[f"ms_dataset/{i}/raw_data"][:, 1].max()
                                  for i in range(len(initial))]
                    initial /= np.asarray(maxima)[:, None]
                cube, mass = cube_reference(initial, mz, config)
                write_input(path, **case)
                self.run_native(path)
                with h5py.File(path, "r") as hdf:
                    data = hdf["ms_dataset"]
                    expected_mz = cube.sum(axis=2)
                    actual_mz = data["mz_grid"][:].reshape(expected_mz.shape)
                    np.testing.assert_allclose(actual_mz, expected_mz, rtol=1e-4, atol=2e-6)
                    axis = data["mass_axis"][:]
                    mass_grid = np.zeros((len(initial), len(axis)))
                    position = (mass - axis[0]) / config["massbins"]
                    lower = np.floor(position).astype(int)
                    fraction = position - lower
                    for scan in range(len(initial)):
                        for offset, weight in ((0, 1 - fraction), (1, fraction)):
                            indexes = lower + offset
                            valid = (indexes >= 0) & (indexes < len(axis))
                            np.add.at(mass_grid[scan], indexes[valid], (cube[scan] * weight)[valid])
                    np.testing.assert_allclose(data["mass_grid"][:].reshape(mass_grid.shape),
                                               mass_grid, rtol=1e-4, atol=2e-5)
                    for name, grid in (("mz", expected_mz), ("mass", mass_grid)):
                        total = grid.sum(axis=0)
                        total *= (1. if config["datanorm"] == 1 else grid.max()) / total.max()
                        np.testing.assert_allclose(data[f"{name}_sum"][:], total,
                                                   rtol=1e-4, atol=2e-5)
                    for scan in range(len(initial)):
                        np.testing.assert_array_equal(data[f"{scan}/mass_data"][:, 0], axis)
                        np.testing.assert_allclose(data[f"{scan}/mass_data"][:, 1], mass_grid[scan],
                                                   rtol=1e-4, atol=2e-5)
                        self.assertEqual(data[str(scan)].attrs["length_mz"], len(mz))
                        self.assertEqual(data[str(scan)].attrs["length_mass"], len(axis))
                    self.assertEqual(hdf["config"].attrs["gridsflag"], 1)


if __name__ == "__main__":
    unittest.main()

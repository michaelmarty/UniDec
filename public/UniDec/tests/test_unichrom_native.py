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


def write_input(path, scans=5, charges=4, wide_mz=False, edge_only=False, **settings):
    config = dict(
        metamode=-1, startz=2, endz=charges + 1, numit=-4,
        dtsig=2.5, unichromzeropad=0, mzsig=1.2, psfun=0, rawflag=1, datanorm=0,
        zzsig=0., msig=0., psig=0., beta=0.,
        masslb=-2008., massub=-4040., massbins=1., adductmass=0.,
        nativezlb=-100., nativezub=100., intthresh=0.,
        minmz=0., maxmz=2000., mzbins=0., subbuff=0., reductionpercent=0.,
    )
    if charges == 1:
        config["masslb"] = -100.
    config.update(settings)
    mz = (np.linspace(500., 2000., 64, dtype=np.float32) if wide_mz else
          np.arange(1000., 1012., dtype=np.float32))
    with h5py.File(path, "w") as hdf:
        attrs = hdf.create_group("config").attrs
        for name, value in config.items():
            attrs[name] = np.int32(value) if isinstance(value, int) else np.float32(value)
        dataset = hdf.create_group("ms_dataset")
        dataset.attrs["num"] = np.int32(scans)
        for scan in range(scans):
            if wide_mz:
                signal = .1 + np.exp(-((mz - 900. - 10 * scan) / 100.) ** 2)
                signal += .3 * np.exp(-((mz - 1600.) / 80.) ** 2)
            else:
                signal = .1 + np.exp(-((mz - 1002. - scan) / 1.5) ** 2)
                signal += .3 * np.exp(-((mz - 1009.) / 1.2) ** 2)
            if edge_only and scan > 0:
                signal = np.full_like(mz, 1e-8)
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


def charge_smooth_reference(cube, mz, charges, config, floor):
    """Apply one charge-smoothing pass from an unchanged input cube."""
    mz_count, charge_count = cube.shape[1:]
    up = np.empty(mz_count * charge_count, dtype=int)
    down = np.empty_like(up)
    for mz_index, mz_value in enumerate(mz):
        for charge_index, charge in enumerate(charges):
            index = mz_index * charge_count + charge_index
            mass = (mz_value - config["adductmass"]) * charge
            upper_mz = mass / (charge + 1) + config["adductmass"]
            lower_mz = (mass / (charge - 1) + config["adductmass"]
                        if charge != 1 else 0.)

            def mapped_index(value):
                position = ((value - mz[0]) / (mz[-1] - mz[0]) * (mz_count - 1))
                return int(np.floor(position + .5))

            upper_index = mapped_index(upper_mz)
            if charge_index == charge_count - 1 or not 0 <= upper_index < mz_count:
                up[index] = index
            else:
                up[index] = upper_index * charge_count + charge_index + 1
            lower_index = mapped_index(lower_mz)
            if charge_index == 0 or not 0 <= lower_index < mz_count:
                down[index] = index
            else:
                down[index] = lower_index * charge_count + charge_index - 1

    flat = cube.reshape(len(cube), -1)
    if floor <= 0:
        return ((flat + flat[:, down] * abs(floor) +
                 flat[:, up] * abs(floor)) / 3).reshape(cube.shape)
    logs = np.log(flat + floor)
    logs[~np.isfinite(logs)] = 0
    return np.maximum(np.exp((logs + logs[:, down] + logs[:, up]) / 3) - floor,
                      0).reshape(cube.shape)


def scan_source_indexes(scan_count, padding, zero_padding):
    indexes = np.arange(scan_count + 2 * padding) - padding
    if zero_padding:
        return np.where((indexes >= 0) & (indexes < scan_count), indexes, -1)
    if scan_count == 1:
        return np.zeros_like(indexes)
    indexes %= 2 * scan_count
    return np.where(indexes < scan_count, indexes, 2 * scan_count - indexes - 1)


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
    mzsig = config["mzsig"] / (2.35482 if config["psfun"] == 0 else 1)
    dtsig = config["dtsig"] / 2.35482
    padding = int(np.ceil(3 * dtsig))
    padded_scans = len(initial) + 2 * padding
    t = np.arange(padded_scans, dtype=float)
    km = peak(mz[0], mz, mzsig, config["psfun"])
    km += peak(2 * mz[-1] - mz[-2], mz, mzsig, config["psfun"])
    kt = peak(0, t, dtsig, config["psfun"]) + peak(len(t), t, dtsig, config["psfun"])
    kernel = np.zeros((padded_scans, len(mz), len(charges)))
    kernel[:, :, 0] = kt[:, None] * km
    spectrum = np.fft.fftn(kernel)
    sources = scan_source_indexes(len(initial), padding, config["unichromzeropad"])

    def forward(values):
        extended = np.zeros_like(kernel)
        valid = sources >= 0
        extended[valid] = values[sources[valid]]
        result = np.fft.ifftn(np.fft.fftn(extended) * spectrum).real
        return result[p:c]

    def adjoint(values):
        extended = np.zeros_like(kernel)
        extended[p:c] = values
        padded = np.fft.ifftn(np.fft.fftn(extended) * spectrum.conj()).real
        result = np.zeros_like(values)
        for scan, source in enumerate(sources):
            if source >= 0:
                result[source] += padded[scan]
        return result

    p, c = padding, padding + len(initial)
    sensitivity = adjoint(np.ones_like(cube))
    for iteration in range(abs(config["numit"])):
        if config["psig"] >= 1 and iteration > 0:
            width = int(config["psig"])
            previous = cube.copy()
            for mz_index in range(len(mz)):
                average = previous[:, max(0, mz_index - width):
                                   mz_index + width + 1].sum(axis=1) / (1 + 2 * width)
                cube[:, mz_index] = np.where(allowed[mz_index], average,
                                             previous[:, mz_index])
        if config["zzsig"] != 0:
            cube = charge_smooth_reference(
                cube, mz, charges, config, config["zzsig"] * observed.max())
            cube *= allowed
        predicted_cube = forward(cube)
        predicted = predicted_cube.sum(axis=2)
        ratio = np.divide(observed, predicted, out=np.zeros_like(observed), where=predicted > 0)
        broadcast = np.broadcast_to(ratio[:, :, None], cube.shape)
        cube *= adjoint(broadcast) * kernel.sum() / sensitivity
        cube = np.maximum(cube, 0) * allowed
    if config["rawflag"] in (0, 2):
        cube = forward(cube) / kernel.sum()
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

    def run_native(self, path, *args):
        result = subprocess.run([self.executable, str(path), *(args or ("-nthreads", "2"))],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def test_grid_refresh_preserves_coupled_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "refresh.hdf5")
            write_input(path)
            # Missing grids are generated by the coupled solver once.
            self.assertIn("iterating", self.run_native(path, "-grids"))
            with h5py.File(path, "r") as hdf:
                before = {name: hdf["ms_dataset/" + name][:]
                          for name in ("mz_grid", "mass_grid", "mz_axis", "mass_axis")}
            output = self.run_native(path, "-grids")
            self.assertNotIn("iterating", output)
            self.assertIn("Grids Already Made", output)
            with h5py.File(path, "r") as hdf:
                self.assertIn("peaks/peakdata", hdf)
                for name, values in before.items():
                    np.testing.assert_array_equal(hdf["ms_dataset/" + name][:], values)

    def test_coupled_grid_refresh_rejects_unavailable_charge_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "charges.hdf5")
            for mode in (6, 7):
                for command in ("-grids", "-all"):
                    write_input(path, exchoice=mode)
                    result = subprocess.run([self.executable, str(path), command],
                                            capture_output=True, text=True, timeout=60)
                    self.assertEqual(result.returncode, 12)
                    self.assertIn("charge-resolved outputs", result.stderr)
                    self.assertNotIn("iterating", result.stdout)

    def test_positive_width_command_contracts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "commands.hdf5")
            write_input(path)
            output = self.run_native(path, "-proc")
            self.assertNotIn("iterating", output)
            with h5py.File(path, "r") as hdf:
                self.assertIn("ms_dataset/0/processed_data", hdf)
                self.assertNotIn("ms_dataset/mz_grid", hdf)

            write_input(path)
            output = self.run_native(path, "-all")
            self.assertEqual(output.count("iterating"), 1)
            with h5py.File(path, "r") as hdf:
                self.assertIn("ms_dataset/mz_grid", hdf)
                self.assertIn("ms_dataset/mass_grid", hdf)
                self.assertIn("peaks/peakdata", hdf)

            for command in ("-extract", "-peaks"):
                self.assertNotIn("iterating", self.run_native(path, command))

            write_input(path)
            output = self.run_native(path, "-newgrids")
            self.assertEqual(output.count("iterating"), 1)
            with h5py.File(path, "r") as hdf:
                self.assertIn("peaks/peakdata", hdf)

            for command in ("-ultraextract", "-charges", "-scanpeaks"):
                result = subprocess.run([self.executable, str(path), command],
                                        capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 12)
                self.assertIn("charge-resolved per-scan outputs", result.stderr)
                self.assertNotIn("iterating", result.stdout)

    def test_scan_padding_prevents_first_to_last_wrap(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "edge.hdf5")
            write_input(path, edge_only=True, numit=0, rawflag=0)
            output = self.run_native(path)
            self.assertIn("mirrored padding", output)
            with h5py.File(path, "r") as hdf:
                mz_count = len(hdf["ms_dataset/mz_axis"])
                grid = hdf["ms_dataset/mz_grid"][:].reshape(-1, mz_count)
            self.assertLess(grid[-1].max(), grid[0].max() * .01)

    def test_matches_charge_cube_reference(self):
        cases = [dict(psfun=shape, rawflag=rawflag, datanorm=norm)
                 for shape in range(3) for rawflag in range(4) for norm in range(2)]
        cases += [dict(scans=1), dict(charges=1), dict(mzsig=0.),
                  dict(scans=4, charges=3), dict(dtsig=.05),
                  dict(unichromzeropad=1),
                  dict(charges=7, rawflag=0),
                  dict(psig=1.), dict(psig=3., charges=7),
                  dict(psig=20., scans=1, charges=1),
                  dict(wide_mz=True, charges=10, psig=1., zzsig=1., numit=-8,
                       masslb=-100., massub=-25000., massbins=10.),
                  dict(wide_mz=True, charges=10, zzsig=-1., numit=-8,
                       masslb=-100., massub=-25000., massbins=10.),
                  dict(wide_mz=True, charges=10, zzsig=1., numit=-8,
                       masslb=-100., massub=-25000., massbins=10.)]
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

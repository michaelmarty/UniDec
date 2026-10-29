import tempfile
import unittest
import os
from pathlib import Path

import h5py
import numpy as np

from unidec.engine import UniDec
from unidec.modules.unidecstructure import UniDecConfig


class TestUniDecFileValidation(unittest.TestCase):
    def test_empty_and_malformed_spectra_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="unidec-invalid-") as directory:
            for name, contents in (("empty.txt", ""), ("malformed.txt", "not numeric\n")):
                with self.subTest(name=name):
                    path = Path(directory, name)
                    path.write_text(contents, encoding="utf-8")
                    engine = UniDec(silent=True)

                    with self.assertRaisesRegex(ValueError, "at least two columns"):
                        engine.open_file(path.name, directory, clean=True, silent=True)

                    self.assertEqual(engine.data.rawdata.size, 0)


class TestUniDecConfigSerialization(unittest.TestCase):
    values = {
        "startz": 7,
        "endz": 42,
        "mzbins": 0.25,
        "masslb": 1234.5,
        "massub": 98765.5,
        "rawflag": 1,
        "rescale_newblur": 0,
        "suppression_topn": 4,
        "suppression_topx": 0.35,
        "suppression_harmonic": 1,
        "cmap": "viridis",
        "smooth": 1.5,
        "subbuff": 3.0,
        "subtype": 1,
        "intthresh": 0.02,
        "datanorm": 0,
        "mzsig": 1.2,
        "psfun": 2,
        "psfunz": 1,
        "zzsig": 2.2,
        "psig": 3.3,
        "beta": 5.5,
        "pusher": 1.25,
        "mindt": 1.0,
        "maxdt": 8.0,
        "twaveflag": 2,
        "temp": 30.0,
        "pressure": 3.0,
        "dtsig": 2.5,
        "chromatogram_type": 1,
        "zerolog": -9.5,
    }
    hdf5_only_values = {
        "exchoice": 4,
        "exchoicez": 2,
        "exthresh": 0.4,
        "exnorm": 3,
        "exnormz": 4,
        "exwindow": 12.5,
        "UClineardecon": 0,
        "UCtype": 1,
        "unichromzeropad": 3,
    }
    text_only_values = {
        "CDslope": 0.42,
        "CDzbins": 2.0,
        "CDres": 1.5,
        "CDScanStart": 3,
        "CDScanEnd": 30,
        "HTksmooth": 4.0,
        "HTcycleindex": 17,
        "demultiplexmode": "MRS",
    }

    def test_newblur_rescaling_defaults_on(self):
        self.assertEqual(UniDecConfig().rescale_newblur, 1)

    def test_managed_output_names_support_short_and_legacy_modes(self):
        with tempfile.TemporaryDirectory(prefix="unidec-names-") as directory:
            config = UniDecConfig()
            config.outfname = str(Path(directory, "sample"))
            config.default_file_names()
            self.assertEqual(config.output_file("mass.txt"), str(Path(directory, "sample_mass.txt")))

            config.outfname = os.path.join(directory, "")
            config.default_file_names(s="")
            self.assertEqual(config.output_file("mass.txt"), str(Path(directory, "mass.txt")))

    def test_text_config_uses_relative_managed_paths_and_preserves_provenance(self):
        with tempfile.TemporaryDirectory(prefix="unidec-relative-config-") as directory:
            run = Path(directory, "sample_unidecfiles")
            run.mkdir()
            source = UniDecConfig()
            source.udir = str(run)
            source.outfname = os.path.join(str(run), "")
            source.default_file_names(s="")
            source.original_file = "sample with spaces.raw"
            source.config_export(source.confname)

            text = Path(source.confname).read_text(encoding="utf-8")
            self.assertIn("input input.dat\n", text)
            self.assertIn("output ./\n", text)
            self.assertIn("original_file sample with spaces.raw\n", text)

            restored = UniDecConfig()
            restored.config_import(source.confname)
            self.assertEqual(restored.original_file, "sample with spaces.raw")

            restored.original_file = "destination.raw"
            restored.config_import(source.confname)
            self.assertEqual(restored.original_file, "destination.raw")

    def test_text_config_uses_exact_parameter_keys(self):
        with tempfile.TemporaryDirectory(prefix="unidec-exact-config-") as directory:
            path = Path(directory, "config.dat")
            path.write_text(
                "numit_extra 1\n"
                "startz_suffix 20\n"
                "aggressive_mode 2\n"
                "peakthresh 0.25\n",
                encoding="utf-8",
            )

            restored = UniDecConfig()
            restored.config_import(str(path))

            self.assertEqual(restored.numit, 100)
            self.assertEqual(restored.startz, 1)
            self.assertEqual(restored.aggressiveflag, 0)
            self.assertEqual(restored.peakthresh, 0.25)


    def test_new_run_uses_short_names_and_versioned_npz_cache(self):
        with tempfile.TemporaryDirectory(prefix="unidec-short-run-") as directory:
            spectrum = Path(directory, "sample.txt")
            np.savetxt(spectrum, [[500.0, 1.0], [501.0, 2.0], [502.0, 1.0]])

            engine = UniDec(ignore_args=True, silent=True)
            engine.open_file(spectrum.name, directory, silent=True)
            output = Path(directory, "sample_unidecfiles")

            self.assertEqual(Path(engine.config.confname), output / "conf.dat")
            self.assertEqual(Path(engine.config.infname), output / "input.dat")
            self.assertFalse((output / "sample_conf.dat").exists())
            with np.load(output / "rawdata.npz", allow_pickle=False) as cache:
                np.testing.assert_allclose(cache["data"], engine.data.rawdata)
                self.assertEqual(int(cache["unidec_rawdata_version"]), 1)
                self.assertEqual(str(cache["mode"]), "ms")
                self.assertEqual(str(cache["original_file"]), spectrum.name)

            restored = UniDec(ignore_args=True, silent=True)
            restored.open_file("rawdata.npz", str(output), silent=True)
            self.assertEqual(restored.config.original_file, spectrum.name)
            np.testing.assert_allclose(restored.data.rawdata, engine.data.rawdata)

            (output / "rawdata.npz").write_bytes(b"truncated")
            repaired = UniDec(ignore_args=True, silent=True)
            repaired.open_file(spectrum.name, directory, silent=True)
            with np.load(output / "rawdata.npz", allow_pickle=False) as cache:
                np.testing.assert_allclose(cache["data"], engine.data.rawdata)

    def test_text_and_hdf5_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="unidec-config-") as directory:
            for extension in ("dat", "hdf5"):
                with self.subTest(extension=extension):
                    path = Path(directory, f"config.{extension}")
                    source = UniDecConfig()
                    source.outfname = str(Path(directory, "source"))
                    source.default_file_names()
                    expected_values = dict(self.values)
                    if extension == "dat":
                        expected_values.update(self.text_only_values)
                    else:
                        expected_values.update(self.hdf5_only_values)
                    for name, value in expected_values.items():
                        setattr(source, name, value)

                    if extension == "dat":
                        source.config_export(str(path))
                    else:
                        source.write_hdf5(str(path))

                    restored = UniDecConfig()
                    restored.outfname = str(Path(directory, "restored"))
                    restored.default_file_names()
                    if extension == "dat":
                        restored.config_import(str(path))
                    else:
                        restored.read_hdf5(str(path))

                    for name, value in expected_values.items():
                        if isinstance(value, float):
                            self.assertAlmostEqual(getattr(restored, name), value, places=6)
                        else:
                            self.assertEqual(getattr(restored, name), value)

    def test_hdf5_names_types_datasets_and_missing_optional_fallbacks(self):
        with tempfile.TemporaryDirectory(prefix="unidec-hdf-contract-") as directory:
            path = Path(directory, "config.hdf5")
            source = UniDecConfig()
            source.startz = 7
            source.dtsig = 2.5
            source.zerolog = -9.5
            source.masslist = np.array([100.0, 200.0])
            source.manuallist = np.array([[100.0, 5.0, 1.0]])
            source.write_hdf5(str(path))

            with h5py.File(path, "r") as hdf:
                attrs = hdf["config"].attrs
                self.assertEqual(attrs["startz"], 7)
                self.assertEqual(attrs["dtsig"], 2.5)
                self.assertEqual(attrs["zerolog"], -9.5)
                self.assertTrue(np.issubdtype(type(attrs["startz"]), np.integer))
                np.testing.assert_array_equal(hdf["config/masslist"], [100.0, 200.0])
                self.assertEqual(hdf["config/manuallist"].shape, (1, 3))

            legacy = Path(directory, "legacy.hdf5")
            with h5py.File(legacy, "w") as hdf:
                hdf.create_group("config").attrs["startz"] = 9
            restored = UniDecConfig()
            default_dtsig = restored.dtsig
            default_suppression = restored.suppression_topn
            default_rescale_newblur = restored.rescale_newblur
            restored.read_hdf5(str(legacy))
            self.assertEqual(restored.startz, 9)
            self.assertEqual(restored.dtsig, default_dtsig)
            self.assertEqual(restored.suppression_topn, default_suppression)
            self.assertEqual(restored.rescale_newblur, default_rescale_newblur)

    def test_check_badness_repairs_malformed_limits_and_rejects_reversed_range(self):
        config = UniDecConfig()
        config.minmz = "bad"
        config.maxmz = "bad"
        badness, warning = config.check_badness()
        self.assertEqual((config.minmz, config.maxmz), (0, 100000))
        self.assertEqual(badness, 0)
        self.assertEqual(warning, "")

        config.minmz = 200
        config.maxmz = 100
        badness, warning = config.check_badness()
        self.assertEqual(badness, 1)
        self.assertIn("Max m/z", warning)


if __name__ == "__main__":
    unittest.main()

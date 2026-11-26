import tempfile
import unittest
import os
from pathlib import Path

import h5py
import numpy as np

from unidec.engine import UniDec
from unidec.modules.IMEng import UniDecIM
from unidec.modules.unidecstructure import DataContainer, UniDecConfig


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
    def test_legacy_peak_shape_inflation_is_dropped_on_export(self):
        with tempfile.TemporaryDirectory() as directory:
            text_path = Path(directory, "conf.dat")
            text_path.write_text("mzsig 2.5\npeakshapeinflate 3\n", encoding="utf-8")
            hdf_path = Path(directory, "config.hdf5")
            with h5py.File(hdf_path, "w") as hdf:
                hdf.create_group("config").attrs.update(mzsig=2.5, peakshapeinflate=3)

            config = UniDecConfig()
            config.config_import(str(text_path))
            self.assertEqual(config.mzsig, 2.5)
            self.assertFalse(hasattr(config, "inflate"))
            config.config_export(str(text_path))
            self.assertNotIn("peakshapeinflate", text_path.read_text(encoding="utf-8"))

            config = UniDecConfig()
            config.read_hdf5(str(hdf_path))
            self.assertEqual(config.mzsig, 2.5)
            self.assertFalse(hasattr(config, "inflate"))
            config.write_hdf5(str(hdf_path))
            with h5py.File(hdf_path, "r") as hdf:
                self.assertNotIn("peakshapeinflate", hdf["config"].attrs)

    conventional_text_scalar_keys = {
        "version", "imflag", "cdmsflag", "input", "output", "original_file",
        "numit", "numz", "endz", "startz", "zzsig", "psig", "beta",
        "suppression_topn", "suppression_topx", "suppression_percent",
        "suppression_startit", "suppression_harmonic", "suppression_satellite",
        "mzsig", "psfun", "zpsfn", "discreteplot", "smashflag", "massub",
        "masslb", "msig", "molig", "massbins", "mtabsig", "isomode", "minmz",
        "maxmz", "subbuff", "subtype", "smooth", "mzbins", "peakwindow",
        "peakthresh", "normthresh", "peakplotthresh", "plotsep", "intthresh",
        "reductionpercent", "aggressive", "rawflag", "rescale_newblur",
        "adductmass", "nativezub", "nativezlb", "poolflag", "accvol",
        "noiseflag", "linflag", "cmap", "peakcmap",
        "spectracmap", "publicationmode", "isotopemode", "peaknorm", "datanorm",
        "baselineflag", "orbimode", "integratelb", "integrateub", "filterwidth",
        "zerolog", "doubledec", "centroided", "verbose", "variablepw", "minratio",
    }
    conventional_hdf5_scalar_keys = {
        "version", "imflag", "cdmsflag", "original_file", "numit", "endz",
        "startz", "zzsig", "psig", "beta", "suppression_topn",
        "suppression_topx", "suppression_percent", "suppression_startit",
        "suppression_harmonic", "suppression_satellite", "mzsig", "psfun",
        "psfunz", "variablepw", "minratio", "discreteplot", "massub", "masslb",
        "msig", "molig", "massbins", "mtabsig", "minmz", "maxmz", "subbuff",
        "subtype", "smooth", "mzbins", "peakwindow", "peakthresh", "normthresh",
        "peakplotthresh", "plotsep", "intthresh", "reductionpercent", "aggressive",
        "rawflag", "rescale_newblur", "adductmass", "nativezub", "nativezlb",
        "poolflag", "accvol", "noiseflag", "linflag", "cmap",
        "peakcmap", "spectracmap", "publicationmode", "isotopemode", "peaknorm",
        "datanorm", "baselineflag", "orbimode", "integratelb", "integrateub",
        "filterwidth", "zerolog", "manualfileflag", "mfileflag", "doubledec",
        "smashflag",
    }
    values = {
        "startz": 7,
        "endz": 42,
        "mzbins": 0.25,
        "masslb": 1234.5,
        "massub": 98765.5,
        "rawflag": 1,
        "aggressiveflag": 2,
        "separation": 0.125,
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
        "variablepw": 1,
        "minratio": 0.2,
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
        "CDzbins": 2,
        "CDres": 1.5,
        "CDScanStart": 3,
        "CDScanEnd": 30,
        "CDiitflag": True,
        "HTksmooth": 4.0,
        "HTcycleindex": 17,
        "demultiplexmode": "MRS",
    }

    im_values = {
        "zout": 1,
        "pusher": 1.25,
        "mindt": 0.5,
        "maxdt": 9.5,
        "ccsub": 2500.0,
        "ccslb": 100.0,
        "dtsig": 2.5,
        "csig": 1.5,
        "ccsbins": 2.0,
        "subbufdt": 1.5,
        "smoothdt": 2.5,
        "nativeccsub": 150.0,
        "nativeccslb": -150.0,
        "twaveflag": 2,
        "temp": 30.0,
        "pressure": 3.0,
        "volt": 150.0,
        "gasmass": 28.0,
        "to": 0.5,
        "driftlength": 1.2,
        "tcal1": 2.1,
        "tcal2": 3.2,
        "tcal3": 4.3,
        "tcal4": 5.4,
        "edc": 0.7,
    }
    im_serialized_keys = {
        "zout", "pusher", "mindt", "maxdt", "ccsub", "ccslb", "dtsig", "csig",
        "ccsbins", "subbufdt", "smoothdt", "ubnativeccs", "lbnativeccs", "twaveflag",
        "temp", "pressure", "volt", "gasmass", "tnaught", "driftlength", "tcal1",
        "tcal2", "tcal3", "tcal4", "edc",
    }
    cd_ht_text_values = {
        "CDslope": 0.42,
        "CDzbins": 2,
        "CDres": 1.5,
        "CDScanStart": 3,
        "CDScanEnd": 30,
        "CDiitflag": True,
        "CDprethresh": 0.2,
        "HTksmooth": 4.0,
        "CDScanCompress": 5,
        "chromatogram_type": 1,
        "chromatogram_norm": 1,
        "HTmaxscans": 50,
        "HToutputlb": 2.2,
        "HToutputub": 9.9,
        "HTtimepad": 1.1,
        "HTanalysistime": 37.0,
        "HTxaxis": "Scan",
        "HTcycleindex": 17,
        "HTcycletime": 3.0,
        "HTtimeshift": 6.0,
        "htbit": "101",
        "FTstart": 1.0,
        "FTend": 9.0,
        "FTflatten": True,
        "FTapodize": 2,
        "demultiplexmode": "MRS",
        "FTsmooth": 3.0,
        "HTmaskn": 99,
        "HTwin": 7,
    }
    cd_ht_text_keys = {
        "CDslope", "CDzbins", "CDres", "CDScanStart", "CDScanEnd", "CDiitflag",
        "CDprethresh", "HTksmooth", "CDScanCompress", "chromatogram_type",
        "chromatogram_norm", "HTmaxscans", "HToutlb", "HToutub", "HTtimepad",
        "HTanalysistime", "HTxaxis", "HTcycleindex", "HTcylctime", "HTtimeshift",
        "htbit", "FTstart", "FTend", "FTflatten", "FTapodize", "demultiplexmode",
        "FTsmooth", "HTmaskn", "HTwin",
    }
    hdf5_mode_values = {
        "metamode": -1,
        "exwindow": 12.5,
        "exchoice": 4,
        "exchoicez": 2,
        "exthresh": 0.4,
        "exnorm": 3,
        "exnormz": 4,
        "UClineardecon": 0,
        "UCtype": 1,
        "unichromzeropad": 3,
        "calc_unichrom_scores": True,
        "time_window": 2.5,
        "chrom_peak_width": 0.8,
        "chromatogram_type": 1,
        "chromatogram_norm": 1,
        "sw_time_window": 3.0,
        "sw_scan_offset": 5,
        "time_start": 4.0,
        "time_end": 14.0,
        "scan_window": 6.0,
    }
    hdf5_mode_keys = {
        "metamode", "exwindow", "exchoice", "exchoicez", "exthresh", "exnorm", "exnormz",
        "UClineardecon", "UCtype", "unichromzeropad", "calc_unichrom_scores",
        "chrom_time_window", "chrom_peak_width", "chromatogram_type",
        "chromatogram_norm", "sw_time_window", "sw_scan_offset", "time_start",
        "time_end", "scan_window",
    }
    text_only_mode_keys = cd_ht_text_keys - {"chromatogram_type", "chromatogram_norm"}
    hdf5_only_mode_keys = hdf5_mode_keys - {"chromatogram_type", "chromatogram_norm"}

    def test_newblur_rescaling_defaults_on(self):
        self.assertEqual(UniDecConfig().rescale_newblur, 1)

    def test_conventional_text_and_hdf5_scalar_key_inventories(self):
        """Literal expected sets protect the reviewed inventory from serializer omissions."""
        with tempfile.TemporaryDirectory(prefix="unidec-key-inventory-") as directory:
            root = Path(directory)
            source = UniDecConfig()
            source.outfname = str(root / "source")
            source.default_file_names()
            source.original_file = "sample.raw"
            source.minmz = 100.0
            source.maxmz = 1000.0
            source.integratelb = 10000.0
            source.integrateub = 20000.0

            text_path = root / "config.dat"
            source.config_export(str(text_path))
            text_keys = {
                line.split(maxsplit=1)[0]
                for line in text_path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            }
            self.assertEqual(self.conventional_text_scalar_keys - text_keys, set())

            hdf_path = root / "config.hdf5"
            source.write_hdf5(str(hdf_path))
            with h5py.File(hdf_path, "r") as hdf:
                hdf5_keys = set(hdf["config"].attrs)
            self.assertEqual(self.conventional_hdf5_scalar_keys - hdf5_keys, set())

    def test_scalar_registry_iterators_shadow_production_serializers(self):
        with tempfile.TemporaryDirectory(prefix="unidec-registry-shadow-") as directory:
            root = Path(directory)
            source = UniDecConfig()
            for name, value in self.values.items():
                setattr(source, name, value)
            source.minmz = ""
            source.maxmz = ""
            source.integratelb = None
            source.integrateub = None

            before = {
                spec.attribute: getattr(source, spec.attribute)
                for spec in source.scalar_parameter_specs()
            }
            text_values = {
                key: (value, dtype)
                for key, value, dtype in source.iter_serialized_parameters("text")
            }
            hdf5_values = {
                key: (value, dtype)
                for key, value, dtype in source.iter_serialized_parameters("hdf5")
            }
            expected_text_keys = self.conventional_text_scalar_keys - {
                "version", "input", "output", "original_file", "numz", "isomode",
                "centroided", "verbose",
            }
            expected_hdf5_keys = self.conventional_hdf5_scalar_keys - {
                "version", "original_file", "integratelb", "integrateub",
            }
            self.assertEqual(set(text_values), expected_text_keys)
            self.assertEqual(set(hdf5_values), expected_hdf5_keys)

            text_path = root / "config.dat"
            source.config_export(str(text_path))
            production_text = {
                line.split(maxsplit=1)[0]: line.split(maxsplit=1)[1]
                for line in text_path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            }
            for key, (value, _) in text_values.items():
                self.assertIn(key, production_text)
                self.assertEqual(production_text[key], str(value), key)

            hdf5_path = root / "config.hdf5"
            source.write_hdf5(str(hdf5_path))
            with h5py.File(hdf5_path, "r") as hdf:
                production_hdf5 = hdf["config"].attrs
                for key, (value, dtype) in hdf5_values.items():
                    self.assertIn(key, production_hdf5)
                    if dtype is float:
                        self.assertAlmostEqual(production_hdf5[key], value, places=12, msg=key)
                    else:
                        self.assertEqual(production_hdf5[key], value, key)

            after = {
                spec.attribute: getattr(source, spec.attribute)
                for spec in source.scalar_parameter_specs()
            }
            self.assertEqual(after, before)
            self.assertNotIn("aggressiveflag", hdf5_values)
            self.assertNotIn("separation", hdf5_values)
            self.assertNotIn("zpsfn", hdf5_values)

    def test_mode_specific_text_inventories_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="unidec-mode-text-") as directory:
            path = Path(directory, "config.dat")
            source = UniDecConfig()
            expected = self.im_values | self.cd_ht_text_values
            for name, value in expected.items():
                setattr(source, name, value)

            source.config_export(str(path))
            keys = {
                line.split(maxsplit=1)[0]
                for line in path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            }
            self.assertEqual(self.im_serialized_keys | self.cd_ht_text_keys, keys &
                             (self.im_serialized_keys | self.cd_ht_text_keys))
            self.assertEqual(keys & self.hdf5_only_mode_keys, set())

            restored = UniDecConfig()
            restored.config_import(str(path))
            for name, value in expected.items():
                self.assertEqual(getattr(restored, name), value, name)
            for name in ("zout", "CDzbins", "CDScanStart", "CDScanEnd", "CDScanCompress",
                         "HTmaxscans", "HTcycleindex", "FTapodize", "HTmaskn", "HTwin"):
                self.assertIsInstance(getattr(restored, name), int, name)
            for name in ("CDiitflag", "FTflatten"):
                self.assertIsInstance(getattr(restored, name), bool, name)

    def test_mode_specific_hdf5_inventories_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="unidec-mode-hdf5-") as directory:
            path = Path(directory, "config.hdf5")
            source = UniDecConfig()
            expected = self.im_values | self.hdf5_mode_values
            for name, value in expected.items():
                setattr(source, name, value)
            source.write_hdf5(str(path))

            expected_keys = self.im_serialized_keys | self.hdf5_mode_keys
            with h5py.File(path, "a") as hdf:
                attrs = hdf["config"].attrs
                self.assertEqual(expected_keys, set(attrs) & expected_keys)
                self.assertEqual(set(attrs) & self.text_only_mode_keys, set())
                attrs["nativeccsub"] = 999.0
                attrs["nativeccslb"] = -999.0
                attrs["future_mode_setting"] = 123

            restored = UniDecConfig()
            restored.read_hdf5(str(path))
            for name, value in expected.items():
                self.assertEqual(getattr(restored, name), value, name)
            self.assertFalse(hasattr(restored, "future_mode_setting"))

            with h5py.File(path, "a") as hdf:
                del hdf["config"].attrs["ubnativeccs"]
                del hdf["config"].attrs["lbnativeccs"]
            legacy = UniDecConfig()
            legacy.read_hdf5(str(path))
            self.assertEqual(legacy.nativeccsub, 999.0)
            self.assertEqual(legacy.nativeccslb, -999.0)

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

    def test_config_paths_recognize_absolute_paths_from_supported_platforms(self):
        with tempfile.TemporaryDirectory(prefix="unidec-platform-paths-") as directory:
            config_dir = os.path.abspath(directory)
            absolute_paths = (
                r"C:\data\sample.dat",
                r"\\server\share\sample.dat",
                "/var/data/sample.dat",
            )
            for path in absolute_paths:
                with self.subTest(path=path):
                    self.assertEqual(UniDecConfig._resolve_config_path(path, config_dir), path)
                    self.assertEqual(UniDecConfig._path_for_config(path, config_dir), path)

            self.assertEqual(
                UniDecConfig._resolve_config_path("sidecars/sample.dat", config_dir),
                os.path.normpath(os.path.join(config_dir, "sidecars/sample.dat")),
            )

    def test_variable_peak_width_text_sidecar_path_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="unidec-fwhm-config-") as directory:
            path = Path(directory, "config.dat")
            source = UniDecConfig()
            source.variablepw = 1
            source.fwhmfile = str(Path(directory, "fwhm.dat"))
            source.config_export(str(path))

            self.assertIn("fwhmfile fwhm.dat\n", path.read_text(encoding="utf-8"))
            restored = UniDecConfig()
            restored.config_import(str(path))
            self.assertEqual(Path(restored.fwhmfile), Path(directory, "fwhm.dat"))

    def test_data_hdf5_round_trips_fwhm_dataset(self):
        with tempfile.TemporaryDirectory(prefix="unidec-fwhm-hdf-") as directory:
            path = Path(directory, "data.hdf5")
            source = DataContainer()
            source.data2 = np.array([[100.0, 1.0], [101.0, 0.5]])
            source.fwhmlist = np.array([0.2, 0.3])
            source.write_hdf5(str(path))

            with h5py.File(path, "r") as hdf:
                np.testing.assert_array_equal(hdf["ms_data/fwhm"], source.fwhmlist)

            restored = DataContainer()
            restored.read_hdf5(str(path))
            np.testing.assert_array_equal(restored.fwhmlist, source.fwhmlist)

    def test_hdf5_preserves_doubledec_smash_and_source_provenance(self):
        with tempfile.TemporaryDirectory(prefix="unidec-hdf-sidecars-") as directory:
            root = Path(directory)
            path = root / "config.hdf5"
            kernel_path = root / "kernel.txt"
            kernel = np.array([[-1.0, 0.25], [0.0, 1.0], [1.0, 0.25]])
            smashlist = np.array([[100.0, 110.0], [200.0, 205.0]])
            np.savetxt(kernel_path, kernel)

            source = UniDecConfig()
            source.outfname = str(root / "source")
            source.default_file_names()
            source.doubledec = True
            source.kernel = str(kernel_path)
            source.smashflag = 1
            source.smashlist = smashlist
            source.original_file = "sample.raw"
            source.write_hdf5(str(path))

            with h5py.File(path, "r") as hdf:
                attrs = hdf["config"].attrs
                self.assertTrue(bool(attrs["doubledec"]))
                self.assertEqual(attrs["smashflag"], 1)
                self.assertEqual(attrs["original_file"], "sample.raw")
                np.testing.assert_array_equal(hdf["config/doubledec_kernel"], kernel)
                np.testing.assert_array_equal(hdf["config/smashlist"], smashlist)

            restored = UniDecConfig()
            restored.outfname = str(root / "restored")
            restored.default_file_names()
            restored.read_hdf5(str(path))
            self.assertTrue(restored.doubledec)
            self.assertEqual(restored.smashflag, 1)
            self.assertEqual(restored.original_file, "sample.raw")
            np.testing.assert_array_equal(restored.kerneldata, kernel)
            np.testing.assert_array_equal(restored.smashlist, smashlist)

            materialized = Path(restored.materialize_doubledec_kernel())
            self.assertEqual(materialized, root / "restored_doubledec_kernel.txt")
            np.testing.assert_array_equal(np.loadtxt(materialized), kernel)

            source.doubledec = False
            source.kernel = ""
            source.kerneldata = np.array([])
            source.smashflag = 0
            source.smashlist = np.array([])
            source.write_hdf5(str(path))
            with h5py.File(path, "r") as hdf:
                self.assertNotIn("doubledec_kernel", hdf["config"])
                self.assertNotIn("smashlist", hdf["config"])
                self.assertFalse(bool(hdf["config"].attrs["doubledec"]))
                self.assertEqual(hdf["config"].attrs["smashflag"], 0)

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

    def test_text_charge_range_prefers_endz_and_accepts_legacy_numz(self):
        with tempfile.TemporaryDirectory(prefix="unidec-charge-config-") as directory:
            path = Path(directory, "config.dat")
            path.write_text("startz 5\nnumz 100\nendz 12\n", encoding="utf-8")
            restored = UniDecConfig()
            restored.config_import(str(path))
            self.assertEqual(restored.startz, 5)
            self.assertEqual(restored.endz, 12)
            self.assertEqual(restored.numz, 8)

            path.write_text("startz 5\nnumz 8\n", encoding="utf-8")
            legacy = UniDecConfig()
            legacy.config_import(str(path))
            self.assertEqual(legacy.endz, 12)
            self.assertEqual(legacy.numz, 8)

    def test_text_config_restores_explicit_relative_sidecar_paths(self):
        with tempfile.TemporaryDirectory(prefix="unidec-sidecar-config-") as directory:
            root = Path(directory)
            sidecars = root / "sidecars"
            sidecars.mkdir()
            masses = sidecars / "custom_masses.dat"
            assignments = sidecars / "custom_manual.dat"
            kernel = sidecars / "custom_kernel.dat"
            np.savetxt(masses, [100.0, 200.0])
            np.savetxt(assignments, [[100.0, 5.0, 1.0]])
            np.savetxt(kernel, [[-1.0, 0.5], [0.0, 1.0], [1.0, 0.5]])
            path = root / "config.dat"
            path.write_text(
                "mfile sidecars/custom_masses.dat\n"
                "manualfile sidecars/custom_manual.dat\n"
                "kernel sidecars/custom_kernel.dat\n",
                encoding="utf-8",
            )

            restored = UniDecConfig()
            restored.config_import(str(path))
            self.assertEqual(Path(restored.mfile), masses)
            self.assertEqual(Path(restored.manualfile), assignments)
            self.assertEqual(Path(restored.kernel), kernel)
            np.testing.assert_array_equal(restored.masslist, [100.0, 200.0])
            np.testing.assert_array_equal(restored.manuallist, [[100.0, 5.0, 1.0]])


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

    def test_mismatched_raw_cache_is_regenerated_from_active_source(self):
        with tempfile.TemporaryDirectory(prefix="unidec-cache-mismatch-") as directory:
            spectrum = Path(directory, "sample.txt")
            source_data = np.array([[500.0, 1.0], [501.0, 2.0], [502.0, 1.0]])
            np.savetxt(spectrum, source_data)

            engine = UniDec(ignore_args=True, silent=True)
            engine.open_file(spectrum.name, directory, silent=True)
            cache_path = Path(directory, "sample_unidecfiles", "rawdata.npz")

            for mode, original_file in (("imms", spectrum.name), ("ms", "other.raw")):
                with self.subTest(mode=mode, original_file=original_file):
                    np.savez_compressed(
                        cache_path,
                        data=np.array([[900.0, 9.0], [901.0, 8.0]]),
                        unidec_rawdata_version=np.array(1),
                        mode=np.array(mode),
                        original_file=np.array(original_file),
                    )
                    restored = UniDec(ignore_args=True, silent=True)
                    restored.open_file(spectrum.name, directory, silent=True)
                    with np.load(cache_path, allow_pickle=False) as cache:
                        self.assertEqual(str(cache["mode"]), "ms")
                        self.assertEqual(str(cache["original_file"]), spectrum.name)
                        np.testing.assert_allclose(cache["data"], source_data)

    def test_imms_managed_cache_round_trip_preserves_mode_and_provenance(self):
        with tempfile.TemporaryDirectory(prefix="unidec-imms-cache-") as directory:
            cache_path = Path(directory, "rawdata.npz")
            sparse_data = np.array([
                [100.0, 1.0, 1.0],
                [100.0, 2.0, 2.0],
                [101.0, 1.0, 3.0],
                [101.0, 2.0, 4.0],
            ])
            np.savez_compressed(
                cache_path,
                data=sparse_data,
                unidec_rawdata_version=np.array(1),
                mode=np.array("imms"),
                original_file=np.array("sample.raw"),
            )

            engine = UniDecIM(ignore_args=True, silent=True)
            engine.open_file(cache_path.name, directory, silent=True)

            self.assertEqual(engine.config.original_file, "sample.raw")
            np.testing.assert_allclose(engine.data.rawdata3, sparse_data)

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
            source.aggressiveflag = 2
            source.separation = 0.125
            source.variablepw = 1
            source.minratio = 0.2
            source.masslist = np.array([100.0, 200.0])
            source.manuallist = np.array([[100.0, 5.0, 1.0]])
            source.write_hdf5(str(path))

            with h5py.File(path, "r") as hdf:
                attrs = hdf["config"].attrs
                self.assertEqual(attrs["startz"], 7)
                self.assertEqual(attrs["dtsig"], 2.5)
                self.assertEqual(attrs["zerolog"], -9.5)
                self.assertEqual(attrs["aggressive"], 2)
                self.assertEqual(attrs["plotsep"], 0.125)
                self.assertEqual(attrs["variablepw"], 1)
                self.assertEqual(attrs["minratio"], 0.2)
                self.assertTrue(np.issubdtype(type(attrs["startz"]), np.integer))
                np.testing.assert_array_equal(hdf["config/masslist"], [100.0, 200.0])
                self.assertEqual(hdf["config/manuallist"].shape, (1, 3))

            legacy = Path(directory, "legacy.hdf5")
            with h5py.File(legacy, "w") as hdf:
                attrs = hdf.create_group("config").attrs
                attrs["startz"] = 9
                attrs["aggressiveflag"] = 1
                attrs["separation"] = 0.5
                attrs["zpsfn"] = 2
            restored = UniDecConfig()
            default_dtsig = restored.dtsig
            default_suppression = restored.suppression_topn
            default_rescale_newblur = restored.rescale_newblur
            restored.read_hdf5(str(legacy))
            self.assertEqual(restored.startz, 9)
            self.assertEqual(restored.aggressiveflag, 1)
            self.assertEqual(restored.separation, 0.5)
            self.assertEqual(restored.psfunz, 2)
            self.assertEqual(restored.dtsig, default_dtsig)
            self.assertEqual(restored.suppression_topn, default_suppression)
            self.assertEqual(restored.rescale_newblur, default_rescale_newblur)

            mixed = Path(directory, "mixed-aliases.hdf5")
            with h5py.File(mixed, "w") as hdf:
                attrs = hdf.create_group("config").attrs
                attrs["aggressive"] = 2
                attrs["aggressiveflag"] = 1
                attrs["plotsep"] = 0.25
                attrs["separation"] = 0.5
                attrs["psfunz"] = 1
                attrs["zpsfn"] = 2
            restored.read_hdf5(str(mixed))
            self.assertEqual(restored.aggressiveflag, 2)
            self.assertEqual(restored.separation, 0.25)
            self.assertEqual(restored.psfunz, 1)

    def test_hdf5_optional_defaults_do_not_mutate_config(self):
        with tempfile.TemporaryDirectory(prefix="unidec-hdf-optionals-") as directory:
            path = Path(directory, "config.hdf5")
            source = UniDecConfig()
            source.minmz = ""
            source.maxmz = ""
            source.tcal3 = ""
            source.tcal4 = ""
            source.integratelb = None
            source.integrateub = None
            source.massbins = 10
            source.peakwindow = 500
            source.write_hdf5(str(path))

            self.assertEqual(source.minmz, "")
            self.assertEqual(source.maxmz, "")
            self.assertEqual(source.tcal3, "")
            self.assertEqual(source.tcal4, "")
            self.assertIsNone(source.integratelb)
            self.assertIsNone(source.integrateub)
            with h5py.File(path, "r") as hdf:
                attrs = hdf["config"].attrs
                self.assertEqual(attrs["minmz"], 0.0)
                self.assertEqual(attrs["maxmz"], 1000000.0)
                self.assertEqual(attrs["tcal3"], 0.0)
                self.assertEqual(attrs["tcal4"], 0.0)
                self.assertNotIn("integratelb", attrs)
                self.assertNotIn("integrateub", attrs)
                self.assertTrue(np.issubdtype(type(attrs["massbins"]), np.floating))
                self.assertTrue(np.issubdtype(type(attrs["peakwindow"]), np.floating))

    def test_automatic_integration_bounds_round_trip_and_accept_legacy_pair(self):
        with tempfile.TemporaryDirectory(prefix="unidec-integration-bounds-") as directory:
            root = Path(directory)
            text_path = root / "config.dat"
            source = UniDecConfig()
            source.integratelb = None
            source.integrateub = None
            source.config_export(str(text_path))

            text = text_path.read_text(encoding="utf-8")
            self.assertIn("integratelb auto\n", text)
            self.assertIn("integrateub auto\n", text)
            restored = UniDecConfig()
            restored.integratelb = 10.0
            restored.integrateub = 20.0
            restored.config_import(str(text_path))
            self.assertIsNone(restored.integratelb)
            self.assertIsNone(restored.integrateub)

            legacy_path = root / "legacy.dat"
            legacy_path.write_text(
                text.replace("integratelb auto", "integratelb -1")
                .replace("integrateub auto", "integrateub -1"),
                encoding="utf-8",
            )
            restored.config_import(str(legacy_path))
            self.assertIsNone(restored.integratelb)
            self.assertIsNone(restored.integrateub)

    def test_hdf5_automatic_bounds_remove_stale_values_and_accept_legacy_pair(self):
        with tempfile.TemporaryDirectory(prefix="unidec-integration-hdf5-") as directory:
            path = Path(directory, "config.hdf5")
            source = UniDecConfig()
            source.integratelb = -2.0
            source.integrateub = 3.0
            source.write_hdf5(str(path))
            source.integratelb = None
            source.integrateub = None
            source.write_hdf5(str(path))

            with h5py.File(path, "r+") as hdf:
                attrs = hdf["config"].attrs
                self.assertNotIn("integratelb", attrs)
                self.assertNotIn("integrateub", attrs)
                attrs["integratelb"] = -1.0
                attrs["integrateub"] = -1.0

            restored = UniDecConfig()
            restored.integratelb = 10.0
            restored.integrateub = 20.0
            restored.read_hdf5(str(path))
            self.assertIsNone(restored.integratelb)
            self.assertIsNone(restored.integrateub)

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

        config.minmz = 100
        config.maxmz = 200
        config.integratelb = 2.0
        config.integrateub = 1.0
        badness, warning = config.check_badness()
        self.assertEqual(badness, 1)
        self.assertIn("Integration lower bound", warning)

        config.integratelb = 0.0
        config.integrateub = 0.0
        badness, warning = config.check_badness()
        self.assertEqual((badness, warning), (0, ""))


if __name__ == "__main__":
    unittest.main()

import tempfile
import unittest
from pathlib import Path

import numpy as np

from unidec.modules.unidecstructure import UniDecConfig

from tests._test_support import copy_unidec_example


class TestBundledConfigCompatibility(unittest.TestCase):
    cases = (
        (
            ("ADH_unidecfiles", "ADH_conf.dat"),
            {
                "inputversion": "9.0.0b1",
                "imflag": 0,
                "cdmsflag": 0,
                "startz": 1,
                "endz": 70,
                "mzbins": 1.0,
                "mzsig": 16.73188,
                "psfun": 2,
                "subtype": 2,
            },
        ),
        (
            ("IMMS", "aqpz_unidecfiles", "aqpz_conf.dat"),
            {
                "inputversion": "9.0.0b1",
                "imflag": 1,
                "cdmsflag": 0,
                "startz": 10,
                "endz": 18,
                "mzbins": 4.0,
                "mindt": 10.08,
                "maxdt": 15.12,
                "dtsig": 0.2,
            },
        ),
        (
            ("CDMS", "BgalGroEL_unidecfiles", "BgalGroEL_conf.dat"),
            {
                "inputversion": "8.1.3",
                "imflag": 0,
                "cdmsflag": 1,
                "startz": 30,
                "endz": 80,
                "mzbins": 10.0,
                "CDslope": 110750.0,
                "CDScanStart": -1,
                "CDScanEnd": -1,
            },
        ),
    )

    def test_known_bundled_values_load_with_expected_types(self):
        with tempfile.TemporaryDirectory(prefix="unidec-bundled-config-") as directory:
            for parts, expected in self.cases:
                with self.subTest(config=parts[-1], directory=parts[:-1]):
                    path = copy_unidec_example(directory, *parts)
                    config = UniDecConfig()
                    config.config_import(str(path))

                    for name, value in expected.items():
                        actual = getattr(config, name)
                        if isinstance(value, float):
                            self.assertIsInstance(actual, (float, np.floating))
                            self.assertAlmostEqual(actual, value, places=6)
                        else:
                            self.assertEqual(actual, value)

    def test_older_cd_config_preserves_defaults_for_absent_suppression_fields(self):
        with tempfile.TemporaryDirectory(prefix="unidec-legacy-config-") as directory:
            path = copy_unidec_example(
                directory, "CDMS", "BgalGroEL_unidecfiles", "BgalGroEL_conf.dat"
            )
            config = UniDecConfig()
            defaults = {
                "suppression_topn": config.suppression_topn,
                "suppression_topx": config.suppression_topx,
                "suppression_harmonic": config.suppression_harmonic,
            }

            config.config_import(str(path))

            for name, value in defaults.items():
                self.assertEqual(getattr(config, name), value)

    def test_import_then_export_rebases_files_and_round_trips_values(self):
        with tempfile.TemporaryDirectory(prefix="unidec-config-roundtrip-") as directory:
            root = Path(directory)
            source_path = copy_unidec_example(root, "ADH_unidecfiles", "ADH_conf.dat")
            source = UniDecConfig()
            source.config_import(str(source_path))
            source.outfname = str(root / "rebased" / "sample")
            Path(source.outfname).parent.mkdir()
            source.default_file_names()
            source.startz = 3
            source.endz = 27
            source.mzbins = 0.5
            source.suppression_topn = 4
            source.cmap = "viridis"
            source.config_export(source.confname)

            restored = UniDecConfig()
            restored.outfname = str(root / "restored" / "sample")
            Path(restored.outfname).parent.mkdir()
            restored.default_file_names()
            restored.config_import(source.confname)

            self.assertEqual(restored.startz, 3)
            self.assertEqual(restored.endz, 27)
            self.assertEqual(restored.numz, 25)
            self.assertEqual(restored.mzbins, 0.5)
            self.assertEqual(restored.suppression_topn, 4)
            self.assertEqual(restored.cmap, "viridis")
            self.assertEqual(Path(restored.confname).parent, root / "restored")
            self.assertTrue(Path(source.confname).is_file())

    def test_bundled_companion_lists_load_and_round_trip(self):
        """Literal bundled oligomer/manual rows survive the real config loaders."""
        cases = (
            ("ADH", "ofile.dat", "oligomerlist",
             np.array([["0.0", "36985.0", "0", "4", "ADH"]])),
            ("GroEL HCD UniDec", "manualfile.dat", "manuallist",
             np.array([[11956.807435, 94.469166, 67.0]])),
        )
        with tempfile.TemporaryDirectory(prefix="unidec-companion-") as directory:
            for index, (stem, suffix, field, expected) in enumerate(cases):
                with self.subTest(field=field):
                    root = Path(directory, str(index))
                    root.mkdir()
                    parts = (stem + "_unidecfiles",)
                    conf = copy_unidec_example(root, *parts, stem + "_conf.dat")
                    companion = copy_unidec_example(root, *parts, stem + "_" + suffix)
                    source = UniDecConfig()
                    source.outfname = str(root / stem)
                    source.default_file_names()
                    source.config_import(str(conf))
                    actual = getattr(source, field)
                    if field == "oligomerlist":
                        np.testing.assert_array_equal(actual, expected)
                    else:
                        np.testing.assert_allclose(actual, expected, atol=1e-6, rtol=0)
                    self.assertEqual(Path(source.ofile if field == "oligomerlist"
                                          else source.manualfile), companion)

                    # Export to a second directory so loading cannot reuse the
                    # original companion; singleton rows must remain 2D.
                    output = root / "exported"
                    output.mkdir()
                    source.outfname = str(output / stem)
                    source.default_file_names()
                    source.config_export(source.confname)
                    restored = UniDecConfig()
                    restored.outfname = source.outfname
                    restored.default_file_names()
                    restored.config_import(source.confname)
                    if field == "oligomerlist":
                        np.testing.assert_array_equal(restored.oligomerlist, expected)
                    else:
                        np.testing.assert_allclose(restored.manuallist, expected, atol=1e-6, rtol=0)


if __name__ == "__main__":
    unittest.main()

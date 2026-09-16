import tempfile
import unittest
from pathlib import Path

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
        "suppression_topn": 4,
        "suppression_topx": 0.35,
        "suppression_harmonic": 1,
        "cmap": "viridis",
    }

    def test_text_and_hdf5_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="unidec-config-") as directory:
            for extension in ("dat", "hdf5"):
                with self.subTest(extension=extension):
                    path = Path(directory, f"config.{extension}")
                    source = UniDecConfig()
                    source.outfname = str(Path(directory, "source"))
                    source.default_file_names()
                    for name, value in self.values.items():
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

                    for name, value in self.values.items():
                        if isinstance(value, float):
                            self.assertAlmostEqual(getattr(restored, name), value, places=6)
                        else:
                            self.assertEqual(getattr(restored, name), value)


if __name__ == "__main__":
    unittest.main()

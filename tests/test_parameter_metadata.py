import unittest

import numpy as np

from unidec.modules.unidecstructure import UniDecConfig


class TestConventionalParameterMetadata(unittest.TestCase):
    def test_schema_is_serializable_and_records_aliases_and_invalidation(self):
        config = UniDecConfig()
        schema = config.parameter_schema()

        self.assertEqual(schema["psfunz"]["text_key"], "zpsfn")
        self.assertEqual(schema["psfunz"]["hdf5_key"], "psfunz")
        self.assertEqual(schema["startz"]["dtype"], "int")
        self.assertEqual(schema["minmz"]["affects"], "processing")
        self.assertEqual(schema["peakthresh"]["affects"], "peaks")
        self.assertNotIn("kernel", schema)
        self.assertNotIn("smashfile", schema)

    def test_parameter_values_expose_unset_optionals_as_none(self):
        config = UniDecConfig()
        self.assertEqual(config.parameter_values(("minmz", "maxmz")), {
            "minmz": None,
            "maxmz": None,
        })

    def test_validation_rejects_unknown_types_bounds_choices_and_nonfinite(self):
        config = UniDecConfig()
        cases = (
            ({"kernel": "unsafe.txt"}, "Unsupported parameters"),
            ({"startz": 2.5}, "must be an integer"),
            ({"startz": True}, "not a boolean"),
            ({"peakthresh": 1.1}, "at most 1"),
            ({"psfunz": 4}, "must be one of"),
            ({"mzsig": np.inf}, "must be finite"),
        )
        for values, message in cases:
            with self.subTest(values=values):
                with self.assertRaisesRegex(ValueError, message):
                    config.validate_parameters(values)

    def test_cross_field_validation_uses_proposed_and_existing_values(self):
        config = UniDecConfig()
        with self.assertRaisesRegex(ValueError, "startz must be less"):
            config.validate_parameters({"startz": 20, "endz": 10})
        with self.assertRaisesRegex(ValueError, "minmz must be less"):
            config.validate_parameters({"minmz": 1000, "maxmz": 100})
        with self.assertRaisesRegex(ValueError, "masslb must be less"):
            config.validate_parameters({"masslb": config.massub})

    def test_update_is_atomic_and_keeps_charge_count_derived(self):
        config = UniDecConfig()
        original_startz = config.startz
        original_endz = config.endz
        with self.assertRaises(ValueError):
            config.update_parameters({"startz": 8, "peakthresh": 2})
        self.assertEqual(config.startz, original_startz)
        self.assertEqual(config.endz, original_endz)

        updated = config.update_parameters({
            "startz": 8,
            "endz": 24,
            "massbins": 5,
            "minmz": None,
        })
        self.assertEqual(updated["massbins"], 5.0)
        self.assertEqual(config.startz, 8)
        self.assertEqual(config.endz, 24)
        self.assertEqual(config.numz, 17)
        self.assertEqual(config.minmz, "")


if __name__ == "__main__":
    unittest.main()

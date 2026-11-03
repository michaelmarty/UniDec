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


class TestIonMobilityParameterMetadata(unittest.TestCase):
    def setUp(self):
        self.config = UniDecConfig()
        self.config.imflag = 1

    def test_schema_adds_only_reviewed_arrival_time_and_ccs_parameters(self):
        schema = self.config.parameter_schema()
        self.assertEqual(schema["mindt"]["modes"], ("imms",))
        self.assertEqual(schema["mindt"]["affects"], "processing")
        self.assertEqual(schema["ccsbins"]["units"], "CCS")
        self.assertNotIn("pusher", schema)
        self.assertNotIn("twaveflag", schema)
        self.assertNotIn("tcal1", schema)

        conventional = UniDecConfig().parameter_schema()
        self.assertNotIn("mindt", conventional)
        self.assertNotIn("ccsbins", conventional)

    def test_im_ranges_validate_atomically(self):
        original_mindt = self.config.mindt
        with self.assertRaisesRegex(ValueError, "mindt must be less"):
            self.config.update_parameters({"mindt": 8, "maxdt": 2})
        self.assertEqual(self.config.mindt, original_mindt)

        with self.assertRaisesRegex(ValueError, "ccslb must be less"):
            self.config.update_parameters({"ccslb": 500, "ccsub": 100})
        self.assertEqual(self.config.ccslb, 100)
        self.assertEqual(self.config.ccsub, 25000)

        updated = self.config.update_parameters({
            "mindt": 1,
            "maxdt": 8,
            "smoothdt": 2,
            "dtsig": 0.5,
            "ccslb": 200,
            "ccsub": 2000,
            "ccsbins": 10,
        })
        self.assertEqual(updated["mindt"], 1.0)
        self.assertEqual(self.config.maxdt, 8.0)
        self.assertEqual(self.config.ccsbins, 10.0)


class TestCDMSParameterMetadata(unittest.TestCase):
    def setUp(self):
        self.config = UniDecConfig()
        self.config.cdmsflag = 1

    def test_schema_adds_reviewed_cdms_parameters_only_in_cdms_mode(self):
        schema = self.config.parameter_schema()
        self.assertIn("CDslope", schema)
        self.assertIn("CDzbins", schema)
        self.assertIn("csig", schema)
        self.assertNotIn("CDslope", UniDecConfig().parameter_schema())

    def test_zero_slope_is_rejected_atomically(self):
        original = self.config.CDslope
        with self.assertRaisesRegex(ValueError, "CDslope must be nonzero"):
            self.config.update_parameters({"CDslope": 0, "CDzbins": 2})
        self.assertEqual(self.config.CDslope, original)
        self.assertEqual(self.config.CDzbins, 1)


class TestUniChromParameterMetadata(unittest.TestCase):
    def setUp(self):
        self.config = UniDecConfig()
        self.config.agent_mode = "unichrom"

    def test_schema_adds_reviewed_chromatogram_parameters(self):
        schema = self.config.parameter_schema()
        self.assertEqual(schema["chromatogram_type"]["choices"], (0, 1))
        self.assertEqual(schema["time_window"]["affects"], "windowing")
        self.assertIn("dtsig", schema)
        self.assertNotIn("time_window", UniDecConfig().parameter_schema())

    def test_time_range_is_validated_atomically(self):
        with self.assertRaisesRegex(ValueError, "both be set"):
            self.config.update_parameters({"time_start": 1, "time_end": None})
        self.assertEqual(self.config.time_start, "")
        self.assertEqual(self.config.time_end, "")

        with self.assertRaisesRegex(ValueError, "less than"):
            self.config.update_parameters({"time_start": 5, "time_end": 2})
        updated = self.config.update_parameters({
            "chromatogram_type": 1,
            "time_window": 0.5,
            "time_start": 1,
            "time_end": 5,
            "dtsig": 0,
        })
        self.assertEqual(updated["time_window"], 0.5)
        self.assertEqual(self.config.chromatogram_type, 1)
        self.assertEqual(self.config.time_end, 5.0)


if __name__ == "__main__":
    unittest.main()

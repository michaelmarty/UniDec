import tempfile
import unittest
from pathlib import Path

from unidec.modules.IMEng import UniDecIM
from unidec.session import UniDecSession

from tests._test_support import copy_unidec_example


class TestUniDecIMSession(unittest.TestCase):
    def test_bundled_imms_complete_workflow_uses_mode_specific_contract(self):
        with tempfile.TemporaryDirectory(prefix="unidec-im-session-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "IMMS", "aqpz.dat")
            session = UniDecSession(
                engine=UniDecIM(ignore_args=True, silent=True),
                working_directory=root / "session",
                allowed_input_roots=[root],
            )

            opened = session.open_spectrum(spectrum)
            self.assertTrue(opened["ok"], opened)
            self.assertEqual(session.get_state()["mode"], "imms")
            schema = session.get_parameters()["schema"]
            self.assertIn("mindt", schema)
            self.assertIn("ccsbins", schema)
            self.assertNotIn("twaveflag", schema)

            updated = session.set_parameters({
                "minmz": 4504,
                "maxmz": 11988,
                "mzbins": 4,
                "startz": 10,
                "endz": 18,
                "masslb": 5000,
                "massub": 500000,
                "massbins": 10,
                "mzsig": 5,
                "mindt": 10.08,
                "maxdt": 15.12,
                "subbufdt": 2,
                "smoothdt": 0,
                "dtsig": 0.2,
                "ccslb": 100,
                "ccsub": 25000,
                "ccsbins": 20,
            })
            self.assertTrue(updated["ok"], updated)

            run = session.run("all")
            self.assertTrue(run["ok"], run)
            self.assertTrue(session.get_state()["available"]["pick_peaks"])
            results = session.get_results(max_peaks=5)
            self.assertTrue(results["ok"], results)
            self.assertGreater(results["peak_count"], 0)

            processed_plot = session.get_plot("processed_spectrum")
            self.assertTrue(processed_plot["ok"], processed_plot)
            self.assertTrue(Path(processed_plot["path"]).is_file())
            ccs_plot = session.get_plot("ccs_distribution")
            self.assertTrue(ccs_plot["ok"], ccs_plot)
            self.assertTrue(Path(ccs_plot["path"]).is_file())
            self.assertEqual(ccs_plot["xlabel"], "CCS")

            manifest = session.get_run_manifest(run["operation_id"])
            self.assertTrue(manifest["ok"], manifest)
            self.assertEqual(manifest["manifest"]["mode"], "imms")

            checkpoint = session.save_state("imms-checkpoint.zip")
            self.assertTrue(checkpoint["ok"], checkpoint)
            restored = UniDecSession(
                engine=UniDecIM(ignore_args=True, silent=True),
                working_directory=root / "restored",
                allowed_input_roots=[root],
            )
            loaded = restored.load_state(checkpoint["path"])
            self.assertTrue(loaded["ok"], loaded)
            self.assertEqual(restored.get_state()["mode"], "imms")
            self.assertEqual(restored.get_state()["status"], "ready")
            self.assertEqual(restored.get_results()["error"]["code"], "stale_results")
            self.assertAlmostEqual(restored.get_parameters()["values"]["mindt"], 10.08)
            self.assertEqual(
                restored.save_state(root / "outside-imms.zip")["error"]["code"],
                "path_not_allowed",
            )
            rerun = restored.run("all")
            self.assertTrue(rerun["ok"], rerun)
            self.assertGreater(restored.get_results()["peak_count"], 0)


if __name__ == "__main__":
    unittest.main()

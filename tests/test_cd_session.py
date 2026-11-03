import tempfile
import unittest
from pathlib import Path

from unidec.modules.CDEng import UniDecCD
from unidec.session import UniDecCDSessionEngineAdapter, UniDecSession

from tests._test_support import copy_unidec_example


class TestUniDecCDSession(unittest.TestCase):
    def test_bundled_cdms_complete_workflow_uses_managed_outputs(self):
        with tempfile.TemporaryDirectory(prefix="unidec-cd-session-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "CDMS", "BgalGroEL.dmt")
            engine = UniDecCD()
            session = UniDecSession(
                engine=engine,
                engine_adapter=UniDecCDSessionEngineAdapter(engine),
                working_directory=root / "session",
                allowed_input_roots=[root],
            )
            try:
                opened = session.open_spectrum(spectrum)
                self.assertTrue(opened["ok"], opened)
                self.assertEqual(session.get_state()["mode"], "cdms")
                self.assertTrue(Path(engine.config.udir).is_relative_to(root / "session"))
                self.assertFalse((root / "BgalGroEL_unidecfiles").exists())

                updated = session.set_parameters({
                    "minmz": 5000.103515625,
                    "maxmz": 19937.166015625,
                    "mzbins": 10,
                    "startz": 30,
                    "endz": 80,
                    "masslb": 5000,
                    "massub": 1000000,
                    "massbins": 100,
                    "numit": 100,
                    "mzsig": 5,
                    "psig": 1,
                    "psfun": 0,
                    "psfunz": 0,
                    "CDslope": 110750,
                    "CDzbins": 1,
                    "CDres": 0,
                    "csig": 5,
                    "zzsig": 1e-6,
                    "subtype": 1,
                })
                self.assertTrue(updated["ok"], updated)

                run = session.run("all")
                self.assertTrue(run["ok"], run)
                results = session.get_results(max_peaks=5)
                self.assertTrue(results["ok"], results)
                self.assertGreater(results["peak_count"], 0)
                self.assertGreaterEqual(results["r_squared"], 0)
                self.assertTrue(Path(session.get_plot("mass_distribution")["path"]).is_file())
                checkpoint = session.save_state("cdms.zip")
                self.assertTrue(checkpoint["ok"], checkpoint)
                restored_engine = UniDecCD()
                restored = UniDecSession(
                    engine=restored_engine,
                    engine_adapter=UniDecCDSessionEngineAdapter(restored_engine),
                    working_directory=root / "restored",
                    allowed_input_roots=[root],
                )
                loaded = restored.load_state(checkpoint["path"])
                self.assertTrue(loaded["ok"], loaded)
                self.assertEqual(restored.get_state()["mode"], "cdms")
                self.assertEqual(restored.get_results()["error"]["code"], "stale_results")
                rerun = restored.run("all")
                self.assertTrue(rerun["ok"], rerun)
                self.assertGreater(restored.get_results()["peak_count"], 0)
                restored_engine.importer.close()
            finally:
                importer = getattr(engine, "importer", None)
                if importer is not None:
                    importer.close()


if __name__ == "__main__":
    unittest.main()

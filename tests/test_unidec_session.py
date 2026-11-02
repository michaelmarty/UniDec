import tempfile
import unittest
from pathlib import Path

from unidec.session import UniDecSession

from tests._test_support import copy_unidec_example


class TestUniDecSession(unittest.TestCase):
    def test_empty_stale_invalid_and_busy_responses_are_structured(self):
        with tempfile.TemporaryDirectory(prefix="unidec-session-unit-") as directory:
            session = UniDecSession(working_directory=directory)
            self.assertEqual(session.get_state()["status"], "empty")
            self.assertEqual(session.run("all")["error"]["code"], "empty")
            self.assertEqual(session.run("unknown")["error"]["code"], "invalid_stage")
            self.assertEqual(session.get_results()["error"]["code"], "stale_results")

            session._lock.acquire()
            try:
                self.assertEqual(session.set_parameters({"startz": 8})["error"]["code"], "busy")
            finally:
                session._lock.release()

    def test_parameter_updates_are_atomic_and_invalidate_results(self):
        with tempfile.TemporaryDirectory(prefix="unidec-session-params-") as directory:
            session = UniDecSession(working_directory=directory)
            original = session.engine.config.startz
            failed = session.set_parameters({"startz": 8, "peakthresh": 2})
            self.assertFalse(failed["ok"])
            self.assertEqual(session.engine.config.startz, original)
            self.assertEqual(session.revision, 0)

            updated = session.set_parameters({"startz": 8, "endz": 24, "massbins": 5})
            self.assertTrue(updated["ok"])
            self.assertEqual(session.revision, 1)
            self.assertEqual(session.engine.config.numz, 17)
            self.assertFalse(any(session.get_state()["available"].values()))

    def test_bundled_adh_headless_rerun_produces_fresh_bounded_results(self):
        with tempfile.TemporaryDirectory(prefix="unidec-session-adh-") as directory:
            root = Path(directory)
            spectrum = copy_unidec_example(root, "ADH.txt")
            session = UniDecSession(working_directory=root / "session")

            opened = session.open_spectrum(spectrum)
            self.assertTrue(opened["ok"], opened)
            first_revision = session.revision
            first_run = session.run("all")
            self.assertTrue(first_run["ok"], first_run)
            first_results = session.get_results(max_peaks=2)
            self.assertTrue(first_results["ok"])
            self.assertLessEqual(len(first_results["peaks"]), 2)
            plot = session.get_plot("mass_distribution")
            self.assertTrue(plot["ok"], plot)
            self.assertTrue(Path(plot["path"]).is_file())
            self.assertGreater(plot["source_points"], 0)
            self.assertEqual(plot["xlabel"], "Mass (Da)")

            updated = session.set_parameters({"startz": 8, "endz": 24})
            self.assertTrue(updated["ok"], updated)
            self.assertGreater(session.revision, first_revision)
            self.assertEqual(session.get_results()["error"]["code"], "stale_results")
            self.assertEqual(session.get_plot("mass_distribution")["error"]["code"], "stale_plot")
            second_run = session.run("all")
            self.assertTrue(second_run["ok"], second_run)
            second_results = session.get_results()
            self.assertTrue(second_results["ok"])
            self.assertEqual(second_results["revision"], session.revision)
            self.assertTrue(session.get_state()["available"]["pick_peaks"])


if __name__ == "__main__":
    unittest.main()

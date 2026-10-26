import json
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from unidec import gui_self_test


class TestGuiSelfTest(unittest.TestCase):
    def test_worker_registry_is_unique_and_fixtures_exist(self):
        names = [worker.name for worker in gui_self_test.WORKERS]
        self.assertEqual(len(names), len(set(names)))
        for worker in gui_self_test.WORKERS:
            self.assertTrue(gui_self_test.resolve_fixture(worker.fixture).exists())

    def test_copy_fixture_supports_files_and_directories(self):
        with tempfile.TemporaryDirectory() as source_dir, tempfile.TemporaryDirectory() as output_dir:
            source_dir = Path(source_dir)
            source_file = source_dir / "fixture.txt"
            source_file.write_text("data", encoding="utf-8")
            nested_dir = source_dir / "fixture.raw"
            nested_dir.mkdir()
            (nested_dir / "data.bin").write_bytes(b"data")

            copied_file = gui_self_test.copy_fixture(source_file, output_dir)
            copied_dir = gui_self_test.copy_fixture(nested_dir, output_dir)

            self.assertEqual(copied_file.read_text(encoding="utf-8"), "data")
            self.assertEqual((copied_dir / "data.bin").read_bytes(), b"data")

    def test_resolve_fixture_rejects_lfs_pointer(self):
        with tempfile.TemporaryDirectory() as directory:
            pointer = Path(directory) / "fixture.raw"
            pointer.write_text(
                "version https://git-lfs.github.com/spec/v1\n"
                "oid sha256:0000000000000000000000000000000000000000000000000000000000000000\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "Git LFS pointer"):
                gui_self_test.resolve_fixture(pointer.name, directory)

    def test_process_timeout_is_reported(self):
        command = [sys.executable, "-c", "import time; time.sleep(2)"]
        started = time.monotonic()
        return_code, stdout, stderr, status = gui_self_test._run_process(command, 0.05)
        self.assertIsNone(return_code)
        self.assertEqual(status, gui_self_test.TIMEOUT)
        self.assertLess(time.monotonic() - started, 1)

    def test_controller_writes_combined_report(self):
        worker = gui_self_test.WORKERS[0]

        def run_process(command, timeout, cancel_event=None):
            report_path = Path(command[command.index("--report") + 1])
            report_path.write_text(json.dumps({
                "name": worker.name,
                "label": worker.label,
                "fixture": worker.fixture,
                "status": gui_self_test.PASS,
                "message": "Opened and validated",
            }), encoding="utf-8")
            return 0, "worker output", "", None

        with tempfile.TemporaryDirectory() as directory:
            report_path = Path(directory) / "report.json"
            with patch.object(gui_self_test, "WORKERS", (worker,)), \
                    patch.object(gui_self_test, "_run_process", side_effect=run_process):
                report = gui_self_test.run_controller(report_path)

            self.assertTrue(report["success"])
            self.assertEqual(report["summary"][gui_self_test.PASS], 1)
            self.assertEqual(report["results"][0]["stdout"], "worker output")
            self.assertEqual(json.loads(report_path.read_text(encoding="utf-8")), report)

    def test_controller_reports_cancellation_before_first_worker(self):
        cancel_event = threading.Event()
        cancel_event.set()
        with tempfile.TemporaryDirectory() as directory:
            report = gui_self_test.run_controller(
                Path(directory) / "report.json",
                cancel_event=cancel_event,
            )
        self.assertFalse(report["success"])
        self.assertEqual(report["summary"][gui_self_test.CANCELLED], 1)
        self.assertEqual(report["results"][0]["name"], "controller")

    def test_unknown_worker_writes_failure_report(self):
        with tempfile.TemporaryDirectory() as directory:
            report_path = Path(directory) / "worker.json"
            self.assertEqual(gui_self_test.run_worker("missing", report_path), 1)
            result = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(result["status"], gui_self_test.FAIL)


if __name__ == "__main__":
    unittest.main()

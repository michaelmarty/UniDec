import unittest

from unidec.run_manifest import RunManifestHistory


class TestRunManifestHistory(unittest.TestCase):
    def test_records_peaks_and_history_are_bounded_and_comparable(self):
        history = RunManifestHistory()
        for operation_id in range(1, 102):
            history.append(
                unidec_version="test",
                session_id="session",
                mode="metaunidec",
                operation_id=operation_id,
                revision=operation_id,
                input_identity={"kind": "test"},
                parameters={"startz": operation_id},
                workload={"spectrum_count": 2},
                stage="all",
                started_at="2026-01-01T00:00:00+00:00",
                duration=0.1,
                response={
                    "ok": True,
                    "results": {
                        "spectrum_count": 2,
                        "peak_count": 30,
                        "mass_points": operation_id,
                        "peaks": [
                            {"mass": float(index), "height": float(30 - index)}
                            for index in range(30)
                        ],
                        "peaks_truncated": False,
                    },
                },
            )

        self.assertIsNone(history.get(1))
        manifest = history.get(101)
        self.assertEqual(len(manifest["results"]["peaks"]), 25)
        self.assertTrue(manifest["results"]["peaks_truncated"])
        comparison = history.compare(100, 101)
        self.assertEqual(comparison["metrics"]["mass_points"]["delta"], 1)
        self.assertEqual(
            comparison["parameter_changes"]["startz"], {"before": 100, "after": 101})

    def test_invalid_operation_identifier_is_rejected(self):
        history = RunManifestHistory()
        with self.assertRaisesRegex(ValueError, "integer or null"):
            history.get(True)


if __name__ == "__main__":
    unittest.main()

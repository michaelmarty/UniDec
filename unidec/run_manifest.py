"""Bounded run provenance shared by headless UniDec workflow sessions."""

from copy import deepcopy


class RunManifestHistory:
    """Retain bounded JSON-compatible run records without owning a session lifecycle."""

    MAX_RUN_RECORDS = 100
    MAX_RECORDED_PEAKS = 25
    MAX_FIXTURE_ERROR_MESSAGE = 2000

    def __init__(self):
        self._records = []

    def append(self, *, unidec_version, session_id, mode, operation_id, revision,
               input_identity, parameters, workload, stage, started_at, duration, response):
        previous_parameters = self._records[-1]["parameters"] if self._records else {}
        parameter_changes = {
            name: {"before": previous_parameters.get(name), "after": value}
            for name, value in parameters.items()
            if previous_parameters.get(name) != value
        }
        results = deepcopy(response.get("results")) if response.get("ok") else None
        if results is not None and "peaks" in results:
            peaks = results["peaks"]
            results["peaks"] = peaks[:self.MAX_RECORDED_PEAKS]
            results["peaks_truncated"] = (
                bool(results.get("peaks_truncated"))
                or len(peaks) > self.MAX_RECORDED_PEAKS
            )
        record = {
            "manifest_version": 1,
            "unidec_version": unidec_version,
            "session_id": session_id,
            "mode": mode,
            "operation_id": operation_id,
            "revision": revision,
            "input": deepcopy(input_identity),
            "parameters": deepcopy(parameters),
            "parameter_changes": parameter_changes,
            "workload": deepcopy(workload),
            "requested_stage": stage,
            "started_at": started_at,
            "duration": duration,
            "status": "success" if response.get("ok") else "error",
            "error": deepcopy(response.get("error")),
            "warnings": [],
            "results": results,
        }
        self._records.append(record)
        del self._records[:-self.MAX_RUN_RECORDS]

    def get(self, operation_id=None):
        if operation_id is None:
            record = self._records[-1] if self._records else None
        elif isinstance(operation_id, int) and not isinstance(operation_id, bool):
            record = next(
                (item for item in self._records if item["operation_id"] == operation_id), None)
        else:
            raise ValueError("operation_id must be an integer or null")
        return None if record is None else deepcopy(record)

    def compare(self, before_operation_id, after_operation_id):
        before = self.get(before_operation_id)
        after = self.get(after_operation_id)
        if before is None or after is None:
            return None
        before_results = before.get("results") or {}
        after_results = after.get("results") or {}
        before_workload = before.get("workload") or {}
        after_workload = after.get("workload") or {}

        def comparison(first, second):
            return {
                "before": first,
                "after": second,
                "delta": None if first is None or second is None else second - first,
            }

        metric_names = (
            "spectrum_count", "window_count", "scan_count", "peak_count", "mass_points",
            "uniscore",
        )
        metrics = {}
        for name in metric_names:
            first = before_results.get(name, before_workload.get(name))
            second = after_results.get(name, after_workload.get(name))
            if first is not None or second is not None:
                metrics[name] = comparison(first, second)
        before_peaks = before_results.get("peaks") or []
        after_peaks = after_results.get("peaks") or []
        first_mass = before_peaks[0]["mass"] if before_peaks else None
        second_mass = after_peaks[0]["mass"] if after_peaks else None
        metrics["dominant_peak_mass"] = comparison(first_mass, second_mass)
        parameter_names = before["parameters"].keys() & after["parameters"].keys()
        return {
            "before_operation_id": before_operation_id,
            "after_operation_id": after_operation_id,
            "parameter_changes": {
                name: {
                    "before": before["parameters"][name],
                    "after": after["parameters"][name],
                }
                for name in parameter_names
                if before["parameters"][name] != after["parameters"][name]
            },
            "metrics": metrics,
        }

    def regression_fixture(self, manifest):
        source_error = manifest["error"] or {}
        message = str(source_error.get("message", ""))
        return {
            "fixture_schema": "unidec-regression-fixture/v1",
            "source": {
                "unidec_version": manifest["unidec_version"],
                "mode": manifest["mode"],
                "run_operation_id": manifest["operation_id"],
                "revision": manifest["revision"],
                "input": manifest["input"],
                "parameters": manifest["parameters"],
                "workload": manifest["workload"],
                "requested_stage": manifest["requested_stage"],
            },
            "expected_error": {
                "code": source_error.get("code"),
                "message": message[:self.MAX_FIXTURE_ERROR_MESSAGE],
                "message_truncated": len(message) > self.MAX_FIXTURE_ERROR_MESSAGE,
            },
        }

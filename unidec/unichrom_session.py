"""Managed, headless open boundary for UniChrom chromatograms."""

import json
from pathlib import Path
import tempfile
import threading
import time
import uuid
from datetime import datetime, timezone

import numpy as np

from unidec.modules.ChromEng import ChromEngine
from unidec.run_manifest import RunManifestHistory


class UniChromSession:
    """Open one raw chromatogram while keeping its derived HDF5 project managed."""

    MAX_RESULT_PEAKS = 100

    def __init__(self, engine=None, working_directory=None, allowed_input_roots=None):
        self.engine = engine or ChromEngine()
        if working_directory is None:
            working_directory = tempfile.mkdtemp(prefix="unichrom-session-")
        self.working_directory = Path(working_directory).resolve()
        self.working_directory.mkdir(parents=True, exist_ok=True)
        if allowed_input_roots is None:
            allowed_input_roots = (self.working_directory,)
        self.allowed_input_roots = tuple(Path(root).resolve() for root in allowed_input_roots)
        if not self.allowed_input_roots:
            raise ValueError("At least one allowed input root is required")
        invalid_roots = [root for root in self.allowed_input_roots if not root.is_dir()]
        if invalid_roots:
            raise ValueError(f"Input roots must be existing directories: {invalid_roots}")

        self.session_id = uuid.uuid4().hex
        self.status = "empty"
        self.revision = 0
        self.operation_id = 0
        self.last_error = None
        self._source = None
        self._managed_project = None
        self._lock = threading.Lock()
        self._stage_revision = {"process": None, "deconvolve": None, "pick_peaks": None}
        self._run_history = RunManifestHistory()

    @staticmethod
    def _is_within(path, root):
        try:
            path.relative_to(root)
            return True
        except ValueError:
            return False

    def _input_path(self, path):
        resolved = Path(path).resolve()
        if not any(self._is_within(resolved, root) for root in self.allowed_input_roots):
            raise ValueError(f"Path is outside the permitted input roots: {resolved}")
        return resolved

    def _output_path(self, path):
        candidate = Path(path)
        if not candidate.is_absolute():
            candidate = self.working_directory / candidate
        resolved = candidate.resolve()
        if not self._is_within(resolved, self.working_directory):
            raise ValueError(f"Path is outside the working directory: {resolved}")
        return resolved

    @staticmethod
    def _axis_summary(values):
        data = np.asarray(values)
        if data.ndim != 2 or data.shape[1] < 2 or len(data) == 0:
            return {"points": 0, "time_range": None, "maximum_intensity": None}
        return {
            "points": int(len(data)),
            "time_range": [float(np.min(data[:, 0])), float(np.max(data[:, 0]))],
            "maximum_intensity": float(np.max(data[:, 1])),
        }

    def _chromatogram_summary(self):
        scans = np.asarray(self.engine.fullscans)
        return {
            "polarity": self.engine.polarity,
            "scan_count": int(scans.size),
            "scan_range": None if scans.size == 0 else [
                int(np.min(scans)), int(np.max(scans))],
            "tic": self._axis_summary(self.engine.ticdat),
            "bpc": self._axis_summary(self.engine.bpcdat),
            "window_count": len(self.engine.data.spectra),
        }

    def _state(self):
        return {
            "session_id": self.session_id,
            "mode": "unichrom",
            "status": self.status,
            "revision": self.revision,
            "operation_id": self.operation_id,
            "last_error": self.last_error,
            "source": self._source,
            "chromatogram": None if self._managed_project is None else self._chromatogram_summary(),
            "available": {
                stage: stage_revision == self.revision
                for stage, stage_revision in self._stage_revision.items()
            },
        }

    def get_state(self):
        return self._state()

    def _manifest_workload(self):
        return {
            "scan_count": int(np.asarray(self.engine.fullscans).size),
            "window_count": len(self.engine.data.spectra),
        }

    def _record_run(self, operation_id, stage, started_at, started, response):
        self._run_history.append(
            unidec_version=self.engine.config.version,
            session_id=self.session_id,
            mode="unichrom",
            operation_id=operation_id,
            revision=self.revision,
            input_identity=self._source,
            parameters=self.engine.config.parameter_values(),
            workload=self._manifest_workload(),
            stage=stage,
            started_at=started_at,
            duration=time.perf_counter() - started,
            response=response,
        )

    def get_run_manifest(self, operation_id=None):
        try:
            manifest = self._run_history.get(operation_id)
        except ValueError as error:
            return self._failure("invalid_operation_id", error, set_error=False)
        if manifest is None:
            return self._failure("run_not_found", "No matching run manifest is available",
                                 set_error=False)
        return {"ok": True, "manifest": manifest}

    def compare_runs(self, before_operation_id, after_operation_id):
        try:
            comparison = self._run_history.compare(before_operation_id, after_operation_id)
        except ValueError as error:
            return self._failure("invalid_operation_id", error, set_error=False)
        if comparison is None:
            return self._failure("run_not_found", "No matching run manifest is available",
                                 set_error=False)
        return {"ok": True, **comparison}

    def get_parameters(self, names=None):
        if self._managed_project is None:
            return self._failure("no_chromatogram", "Open a chromatogram before reading parameters",
                                 set_error=False)
        try:
            return {
                "ok": True,
                "session_id": self.session_id,
                "revision": self.revision,
                "values": self.engine.config.parameter_values(names),
                "schema": self.engine.config.parameter_schema(),
            }
        except ValueError as error:
            return self._failure("invalid_parameters", error, set_error=False)

    def _failure(self, code, message, operation_id=None, set_error=True):
        error = {"code": code, "message": str(message)}
        if operation_id is not None:
            error["operation_id"] = operation_id
        if set_error:
            self.last_error = error
            self.status = "error"
        return {"ok": False, "error": error, "state": self._state()}

    def _invalidate_results(self):
        for stage in self._stage_revision:
            self._stage_revision[stage] = None

    def open_chromatogram(self, path):
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        operation_id = self.operation_id
        self.status = "running"
        self.last_error = None
        try:
            try:
                source = self._input_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, operation_id)
            if not source.exists() or not (source.is_file() or source.is_dir()):
                return self._failure("not_found", f"Chromatogram does not exist: {source}", operation_id)

            operation_directory = self.working_directory / f"operation-{operation_id}"
            operation_directory.mkdir(parents=True, exist_ok=False)
            managed_project = operation_directory / f"{source.stem}.hdf5"
            self.engine.open_chrom(
                str(source), load_hdf5=False, refresh=True, output_path=managed_project)

            self._source = {
                "kind": "directory" if source.is_dir() else "file",
                "path": str(source),
                "size": None if source.is_dir() else source.stat().st_size,
            }
            self._managed_project = managed_project
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return {"ok": True, "operation_id": operation_id, "state": self._state()}
        except Exception as error:
            return self._failure("open_failed", error, operation_id)
        finally:
            self._lock.release()

    def create_time_windows(self, width, start=None, end=None):
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        operation_id = self.operation_id
        self.status = "running"
        self.last_error = None
        try:
            if self._managed_project is None:
                return self._failure("no_chromatogram", "Open a chromatogram before creating windows",
                                     operation_id)
            if (isinstance(width, bool) or not isinstance(width, (int, float))
                    or not np.isfinite(width) or width <= 0):
                return self._failure("invalid_windows", "width must be a positive finite number",
                                     operation_id)
            for name, value in (("start", start), ("end", end)):
                if (value is not None and (isinstance(value, bool)
                        or not isinstance(value, (int, float)) or not np.isfinite(value))):
                    return self._failure(
                        "invalid_windows", f"{name} must be a finite number or null", operation_id)
            if (start is None) != (end is None):
                return self._failure(
                    "invalid_windows", "start and end must either both be set or both be null", operation_id)
            if start is not None and start >= end:
                return self._failure("invalid_windows", "start must be less than end", operation_id)

            self.engine.config.time_window = float(width)
            self.engine.config.time_start = None if start is None else float(start)
            self.engine.config.time_end = None if end is None else float(end)
            self.engine.add_regular_times()
            self.engine.reset_vars(export=True)
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return {"ok": True, "operation_id": operation_id, "state": self._state()}
        except Exception as error:
            return self._failure("window_creation_failed", error, operation_id)
        finally:
            self._lock.release()

    def set_parameters(self, values):
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        operation_id = self.operation_id
        self.status = "running"
        self.last_error = None
        try:
            if self._managed_project is None:
                return self._failure("no_chromatogram", "Open a chromatogram before setting parameters",
                                     operation_id)
            try:
                updated = self.engine.config.update_parameters(values)
            except ValueError as error:
                return self._failure("invalid_parameters", error, operation_id)
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return {
                "ok": True,
                "operation_id": operation_id,
                "updated": updated,
                "state": self._state(),
            }
        except Exception as error:
            return self._failure("parameter_update_failed", error, operation_id)
        finally:
            self._lock.release()

    def run(self, stage="all"):
        if stage != "all":
            return self._failure("invalid_stage", "UniChrom currently supports only stage='all'",
                                 set_error=False)
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        operation_id = self.operation_id
        self.status = "running"
        self.last_error = None
        started = time.perf_counter()
        started_at = datetime.now(timezone.utc).isoformat()

        def finish(response):
            self._record_run(operation_id, stage, started_at, started, response)
            return response

        try:
            if self._managed_project is None:
                return finish(self._failure(
                    "no_chromatogram", "Open a chromatogram before running", operation_id))
            if not self.engine.data.spectra:
                return finish(self._failure(
                    "no_windows", "Create chromatogram windows before running", operation_id))

            self.engine.process_data()
            process_result = getattr(self.engine, "out", None)
            if process_result not in (None, 0):
                raise RuntimeError(f"UniChrom processing returned {process_result}")
            self._stage_revision["process"] = self.revision

            result = self.engine.run_unidec()
            if result not in (None, 0):
                raise RuntimeError(f"UniChrom deconvolution returned {result}")
            self._stage_revision["deconvolve"] = self.revision

            self.engine.pick_peaks()
            self._stage_revision["pick_peaks"] = self.revision
            self.status = "ready"
            return finish({
                "ok": True,
                "operation_id": operation_id,
                "state": self._state(),
                "results": self._results(self.MAX_RESULT_PEAKS),
            })
        except Exception as error:
            return finish(self._failure("run_failed", error, operation_id))
        finally:
            self._lock.release()

    def _results(self, max_peaks):
        peaks = sorted(
            self.engine.pks.peaks,
            key=lambda peak: float(getattr(peak, "height", 0)),
            reverse=True,
        )
        bounded = [{
            "mass": float(peak.mass),
            "height": float(peak.height),
            "dscore": float(getattr(peak, "dscore", 0)),
        } for peak in peaks[:max_peaks]]
        mass_data = np.asarray(self.engine.data.massdat)
        return {
            "revision": self.revision,
            "window_count": len(self.engine.data.spectra),
            "peak_count": len(peaks),
            "peaks": bounded,
            "peaks_truncated": len(peaks) > max_peaks,
            "mass_points": int(len(mass_data)),
            "mass_range": None if len(mass_data) == 0 else [
                float(np.min(mass_data[:, 0])), float(np.max(mass_data[:, 0]))],
            "uniscore": float(getattr(self.engine.pks, "uniscore", 0)),
        }

    def get_results(self, max_peaks=100):
        if not isinstance(max_peaks, int) or isinstance(max_peaks, bool) or max_peaks < 1:
            return self._failure("invalid_limit", "max_peaks must be a positive integer",
                                 set_error=False)
        if self._stage_revision["pick_peaks"] != self.revision:
            return self._failure("stale_results", "Run UniChrom before reading results",
                                 set_error=False)
        return {
            "ok": True,
            "session_id": self.session_id,
            "results": self._results(min(max_peaks, self.MAX_RESULT_PEAKS)),
        }

    def export_regression_fixture(self, operation_id, path):
        previous_status = self.status
        previous_error = self.last_error
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        export_operation_id = self.operation_id
        try:
            manifest_response = self.get_run_manifest(operation_id)
            if not manifest_response.get("ok"):
                error = manifest_response["error"]
                return self._failure(error["code"], error["message"], export_operation_id)
            manifest = manifest_response["manifest"]
            if manifest["status"] != "error":
                return self._failure(
                    "run_succeeded", "Regression fixtures can be exported only for failed runs",
                    export_operation_id)
            try:
                target = self._output_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, export_operation_id)
            if target.suffix.lower() != ".json":
                return self._failure("invalid_path", "Regression fixtures must use a .json suffix",
                                     export_operation_id)
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                with target.open("x", encoding="utf-8") as stream:
                    json.dump(self._run_history.regression_fixture(manifest), stream,
                              indent=2, sort_keys=True)
                    stream.write("\n")
            except FileExistsError:
                return self._failure("already_exists", f"Refusing to overwrite {target}",
                                     export_operation_id)
            self.status = previous_status
            self.last_error = previous_error
            return {"ok": True, "operation_id": export_operation_id, "path": str(target),
                    "source_operation_id": manifest["operation_id"], "state": self._state()}
        except Exception as error:
            return self._failure("export_failed", error, export_operation_id)
        finally:
            self._lock.release()

"""Managed, headless session boundary for MetaUniDec collections."""

import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import threading
import time
import uuid
from datetime import datetime, timezone

import numpy as np

from unidec.metaunidec.mudeng import MetaUniDec
from unidec.run_manifest import RunManifestHistory


class MetaUniDecCollectionSession:
    """Open and summarize one MetaUniDec HDF5 project without mutating its source."""

    STAGES = ("all",)
    MAX_SPECTRUM_SUMMARIES = 100
    MAX_RESULT_PEAKS = 100
    MAX_HASH_BYTES = 256 * 1024 * 1024

    def __init__(self, engine=None, working_directory=None, allowed_input_roots=None):
        self.engine = engine or MetaUniDec()
        if working_directory is None:
            working_directory = tempfile.mkdtemp(prefix="metaunidec-session-")
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
        self._lock = threading.Lock()
        self._source_identity = None
        self._active_project = None
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
        roots = self.allowed_input_roots + (self.working_directory,)
        if not any(self._is_within(resolved, root) for root in roots):
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

    @classmethod
    def _file_identity(cls, path):
        size = path.stat().st_size
        identity = {
            "kind": "metaunidec_collection",
            "path": str(path),
            "size": size,
            "sha256": None,
            "hash_status": "omitted_large_file" if size > cls.MAX_HASH_BYTES else "unavailable",
        }
        if size <= cls.MAX_HASH_BYTES:
            try:
                digest = hashlib.sha256()
                with path.open("rb") as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                        digest.update(chunk)
                identity["sha256"] = digest.hexdigest()
                identity["hash_status"] = "complete"
            except OSError:
                pass
        return identity

    @staticmethod
    def _json_value(value):
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        return str(value)

    def _collection_summary(self):
        data = self.engine.data
        spectra = []
        for spectrum in data.spectra[:self.MAX_SPECTRUM_SUMMARIES]:
            spectra.append({
                "index": self._json_value(spectrum.index),
                "name": self._json_value(spectrum.name),
                "variable_1": self._json_value(spectrum.var1),
                "variable_2": self._json_value(spectrum.var2),
                "ignored": bool(spectrum.ignore),
                "points": {
                    "raw": int(len(spectrum.rawdata)),
                    "processed": int(len(spectrum.data2)),
                    "mass": int(len(spectrum.massdat)),
                    "charge": int(len(spectrum.zdata)),
                },
            })
        return {
            "spectrum_count": len(data.spectra),
            "variable_names": [str(data.v1name), str(data.v2name)],
            "spectra": spectra,
            "spectra_truncated": len(data.spectra) > self.MAX_SPECTRUM_SUMMARIES,
            "aggregate": {
                "mass_points": int(len(data.massdat)),
                "mz_points": int(len(data.mzdat)),
                "peak_count": int(len(data.peaks)),
                "extract_shape": list(np.shape(data.exgrid)),
            },
        }

    def _state(self):
        return {
            "session_id": self.session_id,
            "mode": "metaunidec",
            "status": self.status,
            "revision": self.revision,
            "operation_id": self.operation_id,
            "last_error": self.last_error,
            "source": self._source_identity,
            "collection": None if self._active_project is None else self._collection_summary(),
            "available": {
                stage: stage_revision == self.revision
                for stage, stage_revision in self._stage_revision.items()
            },
        }

    def get_state(self):
        return self._state()

    def _manifest_workload(self):
        return {"spectrum_count": 0 if self._active_project is None else len(self.engine.data.spectra)}

    def _record_run(self, operation_id, stage, started_at, started, response):
        self._run_history.append(
            unidec_version=self.engine.config.version,
            session_id=self.session_id,
            mode="metaunidec",
            operation_id=operation_id,
            revision=self.revision,
            input_identity=self._source_identity,
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
        if self._active_project is None:
            return self._failure("no_collection", "Open a collection before reading parameters",
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

    def open_collection(self, path):
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
            if not source.is_file():
                return self._failure("not_found", f"Collection does not exist: {source}", operation_id)

            operation_directory = self.working_directory / f"operation-{operation_id}"
            operation_directory.mkdir(parents=True, exist_ok=False)
            managed_project = operation_directory / source.name
            source_identity = self._file_identity(source)
            shutil.copy2(source, managed_project)
            self.engine.open(str(managed_project))

            self._source_identity = source_identity
            self._active_project = managed_project
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return {
                "ok": True,
                "operation_id": operation_id,
                "state": self._state(),
            }
        except Exception as error:
            return self._failure("open_failed", error, operation_id)
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
            if self._active_project is None:
                return self._failure("no_collection", "Open a collection before setting parameters",
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
        if stage not in self.STAGES:
            return self._failure("invalid_stage", "MetaUniDec currently supports only stage='all'",
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
            if self._active_project is None:
                return finish(self._failure(
                    "no_collection", "Open a collection before running", operation_id))

            self.engine.process_data()
            process_result = getattr(self.engine, "out", None)
            if process_result not in (None, 0):
                raise RuntimeError(f"MetaUniDec processing returned {process_result}")
            self._stage_revision["process"] = self.revision

            self.engine.run_unidec()
            deconvolution_result = getattr(self.engine, "out", None)
            if deconvolution_result not in (None, 0):
                raise RuntimeError(f"MetaUniDec deconvolution returned {deconvolution_result}")
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
        bounded = []
        for peak in peaks[:max_peaks]:
            bounded.append({
                "mass": float(peak.mass),
                "height": float(peak.height),
                "dscore": float(getattr(peak, "dscore", 0)),
            })
        mass_data = np.asarray(self.engine.data.massdat)
        return {
            "revision": self.revision,
            "spectrum_count": len(self.engine.data.spectra),
            "peak_count": len(peaks),
            "peaks": bounded,
            "peaks_truncated": len(peaks) > max_peaks,
            "mass_points": int(len(mass_data)),
            "mass_range": None if len(mass_data) == 0 else [
                float(np.min(mass_data[:, 0])),
                float(np.max(mass_data[:, 0])),
            ],
        }

    def get_results(self, max_peaks=100):
        if not isinstance(max_peaks, int) or isinstance(max_peaks, bool) or max_peaks < 1:
            return self._failure("invalid_limit", "max_peaks must be a positive integer",
                                 set_error=False)
        if self._stage_revision["pick_peaks"] != self.revision:
            return self._failure("stale_results", "Run the collection before reading results",
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

    def save_checkpoint(self, path):
        if not self._lock.acquire(blocking=False):
            return self._failure("busy", "Another session mutation is running", set_error=False)
        self.operation_id += 1
        operation_id = self.operation_id
        self.status = "running"
        self.last_error = None
        try:
            if self._active_project is None:
                return self._failure("no_collection", "Open a collection before saving a checkpoint",
                                     operation_id)
            try:
                target = self._output_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, operation_id)
            if target.suffix.lower() not in (".h5", ".hdf5"):
                return self._failure("invalid_path", "Collection checkpoints must use .h5 or .hdf5",
                                     operation_id)
            if target.exists():
                return self._failure("already_exists", f"Checkpoint already exists: {target}", operation_id)
            target.parent.mkdir(parents=True, exist_ok=True)
            self.engine.config.write_hdf5(str(self._active_project))
            shutil.copy2(self._active_project, target)
            self.status = "ready"
            return {
                "ok": True,
                "operation_id": operation_id,
                "path": str(target),
                "state": self._state(),
            }
        except Exception as error:
            return self._failure("checkpoint_save_failed", error, operation_id)
        finally:
            self._lock.release()

    def load_checkpoint(self, path):
        return self.open_collection(path)

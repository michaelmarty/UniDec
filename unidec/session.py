"""Persistent, GUI-independent control surface for conventional UniDec."""

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import threading
import time
import uuid
import zipfile

import numpy as np

from unidec.engine import UniDec


class UniDecSessionEngineAdapter:
    """Normalize the engine calls and plot sources used by a session."""

    def __init__(self, engine):
        self.engine = engine

    @property
    def mode(self):
        return "imms" if self.engine.config.imflag == 1 else "ms"

    def open_spectrum(self, spectrum, output_prefix):
        self.engine.outfile = str(output_prefix)
        self.engine.open_file(
            spectrum.name,
            str(spectrum.parent),
            simple_output=True,
            silent=True,
        )

    def process(self):
        result = self.engine.process_data(silent=True)
        if result not in (None, 0):
            raise RuntimeError(f"UniDec processing returned {result}")

    def deconvolve(self):
        return self.engine.run_unidec(silent=True)

    def pick_peaks(self):
        self.engine.pick_peaks()

    def save_state(self, target):
        self.engine.save_state(str(target))

    def load_state(self, source):
        return self.engine.load_state(str(source))

    def plot_definitions(self, process_fresh):
        processed_data = self.engine.data.data2
        if self.mode == "imms" and process_fresh:
            data3 = np.asarray(self.engine.data.data3)
            if data3.ndim == 2 and data3.shape[1] >= 3 and len(data3) > 0:
                mz, inverse = np.unique(data3[:, 0], return_inverse=True)
                processed_data = np.column_stack((mz, np.bincount(inverse, weights=data3[:, 2])))
        plots = {
            "processed_spectrum": ("process", processed_data, "m/z", "Intensity"),
            "mass_distribution": (
                "deconvolve", self.engine.data.massdat, "Mass (Da)", "Intensity"),
        }
        if self.mode == "imms":
            plots["ccs_distribution"] = (
                "deconvolve", self.engine.data.ccsdata, "CCS", "Intensity")
        return plots


class UniDecCDSessionEngineAdapter(UniDecSessionEngineAdapter):
    """Adapt the distinct CD-MS engine workflow to the session contract."""

    @property
    def mode(self):
        return "cdms"

    def open_spectrum(self, spectrum, output_prefix):
        importer = getattr(self.engine, "importer", None)
        if importer is not None:
            importer.close()
        self.engine.outfile = str(output_prefix)
        self.engine.open_file(str(spectrum), refresh=True)

    def process(self):
        if self.engine.process_data(transform=True) != 1:
            raise RuntimeError("CD-MS processing did not produce a histogram")

    def deconvolve(self):
        self.engine.run_deconvolution(process_data=False)
        return 0

    def pick_peaks(self):
        self.engine.pick_peaks(calc_dscore=False)
        self.engine.dscore()

    def save_state(self, target):
        self.engine.save_state(str(target))

    def load_state(self, source):
        extract_directory = source.parent / "cdms-checkpoint"
        extract_directory.mkdir()
        with zipfile.ZipFile(source) as archive:
            archive.extractall(extract_directory)
        raw_cache = extract_directory / "rawdata.npz"
        if not raw_cache.is_file():
            return False
        self.engine.outfile = str(extract_directory / source.stem)
        self.engine.open_file(str(raw_cache), refresh=False)
        return True


class UniDecSession:
    """Own one conventional-MS engine and serialize its mutable operations."""

    STAGES = ("process", "deconvolve", "pick_peaks", "all")
    MAX_RUN_RECORDS = 100
    MAX_RECORDED_PEAKS = 25
    MAX_HASH_BYTES = 256 * 1024 * 1024
    MAX_FIXTURE_ERROR_MESSAGE = 2000

    def __init__(self, engine=None, working_directory=None, allowed_input_roots=None,
                 engine_adapter=None):
        self.engine = engine or UniDec(ignore_args=True, silent=True)
        self.engine_adapter = engine_adapter or UniDecSessionEngineAdapter(self.engine)
        if self.engine_adapter.engine is not self.engine:
            raise ValueError("The session and engine adapter must own the same engine")
        if working_directory is None:
            working_directory = tempfile.mkdtemp(prefix="unidec-session-")
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
        self._stage_revision = {"process": None, "deconvolve": None, "pick_peaks": None}
        self._last_duration = None
        self._input_identity = None
        self._run_records = []

    def _state(self):
        config = self.engine.config
        return {
            "session_id": self.session_id,
            "mode": self._mode(),
            "file": config.filename or None,
            "working_directory": str(self.working_directory),
            "allowed_input_roots": [str(root) for root in self.allowed_input_roots],
            "status": self.status,
            "revision": self.revision,
            "operation_id": self.operation_id,
            "last_error": self.last_error,
            "available": {
                stage: self._stage_revision[stage] == self.revision
                for stage in self._stage_revision
            },
        }

    def get_state(self):
        return self._state()

    def get_parameters(self, names=None):
        return {
            "session_id": self.session_id,
            "revision": self.revision,
            "values": self.engine.config.parameter_values(names),
            "schema": self.engine.config.parameter_schema(),
        }

    def _failure(self, code, message, operation_id=None, set_error=True):
        error = {"code": code, "message": str(message)}
        if operation_id is not None:
            error["operation_id"] = operation_id
        if set_error:
            self.last_error = error
            self.status = "error"
        return {"ok": False, "error": error, "state": self._state()}

    def _begin_mutation(self):
        if not self._lock.acquire(blocking=False):
            return None
        self.operation_id += 1
        self.status = "running"
        self.last_error = None
        return self.operation_id

    def _finish_success(self, operation_id, **values):
        response = {"ok": True, "operation_id": operation_id, "state": self._state()}
        response.update(values)
        return response

    def _invalidate_results(self):
        for stage in self._stage_revision:
            self._stage_revision[stage] = None

    def _mode(self):
        return self.engine_adapter.mode

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

    def open_spectrum(self, path):
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            try:
                spectrum = self._input_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, operation_id)
            if not spectrum.is_file():
                return self._failure("not_found", f"Spectrum does not exist: {spectrum}", operation_id)
            session_directory = self.working_directory / f"operation-{operation_id}"
            session_directory.mkdir(parents=True, exist_ok=True)
            self.engine_adapter.open_spectrum(spectrum, session_directory / spectrum.stem)
            self._input_identity = self._file_identity(spectrum, "spectrum")
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return self._finish_success(operation_id)
        except Exception as error:
            return self._failure("open_failed", error, operation_id)
        finally:
            self._lock.release()

    def set_parameters(self, values):
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            updated = self.engine.config.update_parameters(values)
            self.revision += 1
            self._invalidate_results()
            self.status = "ready" if self.engine.config.filename else "empty"
            return self._finish_success(operation_id, values=updated)
        except ValueError as error:
            return self._failure("invalid_parameters", error, operation_id)
        finally:
            self._lock.release()

    def save_state(self, path):
        previous_status = self.status
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            if not self.engine.config.filename:
                return self._failure("empty", "Open a spectrum before saving state", operation_id)
            try:
                target = self._output_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, operation_id)
            if target.suffix.lower() != ".zip":
                return self._failure("invalid_path", "State checkpoints must use a .zip suffix", operation_id)
            if target.exists():
                return self._failure("already_exists", f"Refusing to overwrite {target}", operation_id)
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                self.engine_adapter.save_state(target)
            except NotImplementedError as error:
                return self._failure("unsupported_mode", error, operation_id)
            self.status = previous_status
            return self._finish_success(operation_id, path=str(target))
        except Exception as error:
            return self._failure("save_failed", error, operation_id)
        finally:
            self._lock.release()

    def load_state(self, path):
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            try:
                source = self._input_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, operation_id)
            if not source.is_file():
                return self._failure("not_found", f"State checkpoint does not exist: {source}", operation_id)
            if source.suffix.lower() != ".zip":
                return self._failure("invalid_path", "State checkpoints must use a .zip suffix", operation_id)
            with zipfile.ZipFile(source) as archive:
                for member in archive.infolist():
                    member_path = Path(member.filename.replace("\\", "/"))
                    if member_path.is_absolute() or ".." in member_path.parts:
                        return self._failure(
                            "invalid_checkpoint", f"Unsafe checkpoint member: {member.filename}", operation_id)
            operation_directory = self.working_directory / f"operation-{operation_id}"
            operation_directory.mkdir(parents=True, exist_ok=True)
            managed_copy = operation_directory / source.name
            shutil.copy2(source, managed_copy)
            try:
                loaded = self.engine_adapter.load_state(managed_copy)
            except NotImplementedError as error:
                return self._failure("unsupported_mode", error, operation_id)
            if not loaded:
                return self._failure("load_failed", f"Could not load {source}", operation_id)
            self._input_identity = self._file_identity(source, "checkpoint")
            self.revision += 1
            self._invalidate_results()
            self.status = "ready"
            return self._finish_success(operation_id, source=str(source))
        except Exception as error:
            return self._failure("load_failed", error, operation_id)
        finally:
            self._lock.release()

    def _require_fresh(self, stage):
        return self._stage_revision[stage] == self.revision

    @classmethod
    def _file_identity(cls, path, kind):
        size = path.stat().st_size
        identity = {
            "kind": kind,
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

    def _record_run(self, operation_id, stage, started_at, duration, response):
        parameters = self.engine.config.parameter_values()
        previous = self._run_records[-1] if self._run_records else None
        previous_parameters = {} if previous is None else previous["parameters"]
        parameter_changes = {
            name: {"before": previous_parameters.get(name), "after": value}
            for name, value in parameters.items()
            if previous_parameters.get(name) != value
        }
        record = {
            "manifest_version": 1,
            "unidec_version": self.engine.config.version,
            "session_id": self.session_id,
            "mode": self._mode(),
            "operation_id": operation_id,
            "revision": self.revision,
            "input": None if self._input_identity is None else dict(self._input_identity),
            "parameters": parameters,
            "parameter_changes": parameter_changes,
            "requested_stage": stage,
            "started_at": started_at,
            "duration": duration,
            "status": "success" if response.get("ok") else "error",
            "error": response.get("error"),
            "warnings": [],
            "output_prefix": self.engine.outfile or None,
            "results": None,
        }
        results = {}
        if self._require_fresh("deconvolve"):
            results["r_squared"] = float(self.engine.config.error)
        if self._require_fresh("pick_peaks"):
            ranked = sorted(self.engine.pks.peaks, key=lambda peak: peak.height, reverse=True)
            results["peak_count"] = len(ranked)
            results["peaks"] = [
                {
                    "mass": float(peak.mass),
                    "height": float(peak.height),
                    "integral": float(peak.integral),
                    "dscore": float(peak.dscore),
                }
                for peak in ranked[:self.MAX_RECORDED_PEAKS]
            ]
            results["peaks_truncated"] = len(ranked) > self.MAX_RECORDED_PEAKS
        if results:
            record["results"] = results
        self._run_records.append(record)
        del self._run_records[:-self.MAX_RUN_RECORDS]

    def get_run_manifest(self, operation_id=None):
        if operation_id is None:
            record = self._run_records[-1] if self._run_records else None
        elif isinstance(operation_id, int) and not isinstance(operation_id, bool):
            record = next(
                (item for item in self._run_records if item["operation_id"] == operation_id),
                None,
            )
        else:
            return self._failure(
                "invalid_operation_id", "operation_id must be an integer or null", set_error=False)
        if record is None:
            return self._failure("run_not_found", "No matching run manifest is available", set_error=False)
        return {"ok": True, "manifest": deepcopy(record)}

    def compare_runs(self, before_operation_id, after_operation_id):
        before_response = self.get_run_manifest(before_operation_id)
        after_response = self.get_run_manifest(after_operation_id)
        if not before_response.get("ok"):
            return before_response
        if not after_response.get("ok"):
            return after_response
        before = before_response["manifest"]
        after = after_response["manifest"]
        before_results = before.get("results") or {}
        after_results = after.get("results") or {}

        def comparison(name):
            first = before_results.get(name)
            second = after_results.get(name)
            return {
                "before": first,
                "after": second,
                "delta": None if first is None or second is None else second - first,
            }

        before_peaks = before_results.get("peaks") or []
        after_peaks = after_results.get("peaks") or []
        before_mass = before_peaks[0]["mass"] if before_peaks else None
        after_mass = after_peaks[0]["mass"] if after_peaks else None
        parameter_changes = {
            name: {"before": before["parameters"][name], "after": after["parameters"][name]}
            for name in before["parameters"]
            if before["parameters"][name] != after["parameters"][name]
        }
        return {
            "ok": True,
            "before_operation_id": before_operation_id,
            "after_operation_id": after_operation_id,
            "parameter_changes": parameter_changes,
            "r_squared": comparison("r_squared"),
            "peak_count": comparison("peak_count"),
            "dominant_peak_mass": {
                "before": before_mass,
                "after": after_mass,
                "delta": None if before_mass is None or after_mass is None else after_mass - before_mass,
            },
        }

    def export_regression_fixture(self, operation_id, path):
        previous_status = self.status
        previous_error = self.last_error
        export_operation_id = self._begin_mutation()
        if export_operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            manifest_response = self.get_run_manifest(operation_id)
            if not manifest_response.get("ok"):
                error = manifest_response["error"]
                return self._failure(error["code"], error["message"], export_operation_id)
            manifest = manifest_response["manifest"]
            if manifest["status"] != "error":
                return self._failure(
                    "run_succeeded",
                    "Regression fixture descriptions can be exported only for failed runs",
                    export_operation_id,
                )
            try:
                target = self._output_path(path)
            except ValueError as error:
                return self._failure("path_not_allowed", error, export_operation_id)
            if target.suffix.lower() != ".json":
                return self._failure(
                    "invalid_path",
                    "Regression fixture descriptions must use a .json suffix",
                    export_operation_id,
                )

            source_error = manifest["error"] or {}
            message = str(source_error.get("message", ""))
            fixture = {
                "fixture_schema": "unidec-regression-fixture/v1",
                "source": {
                    "unidec_version": manifest["unidec_version"],
                    "mode": manifest["mode"],
                    "run_operation_id": manifest["operation_id"],
                    "revision": manifest["revision"],
                    "input": manifest["input"],
                    "parameters": manifest["parameters"],
                    "requested_stage": manifest["requested_stage"],
                },
                "expected_error": {
                    "code": source_error.get("code"),
                    "message": message[:self.MAX_FIXTURE_ERROR_MESSAGE],
                    "message_truncated": len(message) > self.MAX_FIXTURE_ERROR_MESSAGE,
                },
            }
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                with target.open("x", encoding="utf-8") as stream:
                    json.dump(fixture, stream, indent=2, sort_keys=True)
                    stream.write("\n")
            except FileExistsError:
                return self._failure(
                    "already_exists", f"Refusing to overwrite {target}", export_operation_id)
            self.status = previous_status
            self.last_error = previous_error
            return self._finish_success(
                export_operation_id,
                path=str(target),
                source_operation_id=manifest["operation_id"],
            )
        except Exception as error:
            return self._failure("export_failed", error, export_operation_id)
        finally:
            self._lock.release()

    def run(self, stage="all"):
        if stage not in self.STAGES:
            return self._failure("invalid_stage", f"Unsupported run stage: {stage}", set_error=False)
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        started = time.perf_counter()
        started_at = datetime.now(timezone.utc).isoformat()

        def finish(response):
            duration = time.perf_counter() - started
            self._record_run(operation_id, stage, started_at, duration, response)
            return response

        try:
            if not self.engine.config.filename:
                return finish(self._failure("empty", "Open a spectrum before running UniDec", operation_id))

            stages = ("process", "deconvolve", "pick_peaks") if stage == "all" else (stage,)
            for current in stages:
                if current == "process":
                    self._invalidate_results()
                    try:
                        self.engine_adapter.process()
                    except RuntimeError as error:
                        return finish(self._failure("processing_failed", error, operation_id))
                    self._stage_revision["process"] = self.revision
                elif current == "deconvolve":
                    if not self._require_fresh("process"):
                        return finish(self._failure(
                            "stale_prerequisite", "Process data at the current revision first", operation_id)
                        )
                    self._stage_revision["deconvolve"] = None
                    self._stage_revision["pick_peaks"] = None
                    result = self.engine_adapter.deconvolve()
                    if result != 0:
                        diagnostic = getattr(self.engine, "last_run_diagnostic", None)
                        return finish(self._failure(
                            "native_failed",
                            f"UniDec returned {result}" + (f": {diagnostic}" if diagnostic else ""),
                            operation_id,
                        ))
                    self._stage_revision["deconvolve"] = self.revision
                    self._stage_revision["pick_peaks"] = None
                else:
                    if not self._require_fresh("deconvolve"):
                        return finish(self._failure(
                            "stale_prerequisite", "Deconvolve at the current revision first", operation_id)
                        )
                    self._stage_revision["pick_peaks"] = None
                    self.engine_adapter.pick_peaks()
                    self._stage_revision["pick_peaks"] = self.revision

            self._last_duration = time.perf_counter() - started
            self.status = "complete" if self._require_fresh("pick_peaks") else "ready"
            return finish(self._finish_success(operation_id, duration=self._last_duration))
        except Exception as error:
            return finish(self._failure("run_failed", error, operation_id))
        finally:
            self._lock.release()

    def get_results(self, max_peaks=100):
        if not isinstance(max_peaks, int) or isinstance(max_peaks, bool) or max_peaks < 0:
            return self._failure("invalid_limit", "max_peaks must be a nonnegative integer", set_error=False)
        if not self._require_fresh("pick_peaks"):
            return self._failure("stale_results", "Run all stages at the current revision", set_error=False)
        peaks = []
        for peak in self.engine.pks.peaks[:max_peaks]:
            peaks.append({
                "mass": float(peak.mass),
                "height": float(peak.height),
                "integral": float(peak.integral),
                "dscore": float(peak.dscore),
            })
        return {
            "ok": True,
            "session_id": self.session_id,
            "revision": self.revision,
            "operation_id": self.operation_id,
            "runtime": self._last_duration,
            "r_squared": float(self.engine.config.error),
            "peak_count": len(self.engine.pks.peaks),
            "peaks": peaks,
            "truncated": len(self.engine.pks.peaks) > max_peaks,
        }

    def get_plot(self, name):
        """Render a named plot to a managed PNG artifact without constructing wxPython."""
        plots = self.engine_adapter.plot_definitions(self._require_fresh("process"))
        if name not in plots:
            return self._failure("unknown_plot", f"Unsupported plot: {name}", set_error=False)
        required_stage, data, xlabel, ylabel = plots[name]
        if not self._require_fresh(required_stage):
            return self._failure(
                "stale_plot", f"Run {required_stage} at the current revision first", set_error=False)
        data = np.asarray(data)
        if data.ndim != 2 or data.shape[1] < 2 or len(data) == 0:
            return self._failure("plot_unavailable", f"No data are available for {name}", set_error=False)

        from matplotlib.figure import Figure

        artifacts = self.working_directory / "artifacts"
        artifacts.mkdir(parents=True, exist_ok=True)
        path = artifacts / f"revision-{self.revision}-{name}.png"
        figure = Figure(figsize=(8, 5))
        axis = figure.subplots()
        axis.plot(data[:, 0], data[:, 1])
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_title(name.replace("_", " ").title())
        figure.savefig(path, format="png", dpi=120, bbox_inches="tight")
        return {
            "ok": True,
            "session_id": self.session_id,
            "revision": self.revision,
            "operation_id": self.operation_id,
            "name": name,
            "path": str(path),
            "mime_type": "image/png",
            "source_points": len(data),
            "xlabel": xlabel,
            "ylabel": ylabel,
        }

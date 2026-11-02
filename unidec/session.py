"""Persistent, GUI-independent control surface for conventional UniDec."""

from pathlib import Path
import tempfile
import threading
import time
import uuid

import numpy as np

from unidec.engine import UniDec


class UniDecSession:
    """Own one conventional-MS engine and serialize its mutable operations."""

    STAGES = ("process", "deconvolve", "pick_peaks", "all")

    def __init__(self, engine=None, working_directory=None):
        self.engine = engine or UniDec(ignore_args=True, silent=True)
        if working_directory is None:
            working_directory = tempfile.mkdtemp(prefix="unidec-session-")
        self.working_directory = Path(working_directory).resolve()
        self.working_directory.mkdir(parents=True, exist_ok=True)
        self.session_id = uuid.uuid4().hex
        self.status = "empty"
        self.revision = 0
        self.operation_id = 0
        self.last_error = None
        self._lock = threading.Lock()
        self._stage_revision = {"process": None, "deconvolve": None, "pick_peaks": None}
        self._last_duration = None

    def _state(self):
        config = self.engine.config
        return {
            "session_id": self.session_id,
            "mode": "ms",
            "file": config.filename or None,
            "working_directory": str(self.working_directory),
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

    def open_spectrum(self, path):
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        try:
            spectrum = Path(path).resolve()
            if not spectrum.exists():
                return self._failure("not_found", f"Spectrum does not exist: {spectrum}", operation_id)
            session_directory = self.working_directory / f"operation-{operation_id}"
            session_directory.mkdir(parents=True, exist_ok=True)
            self.engine.outfile = str(session_directory / spectrum.stem)
            self.engine.open_file(
                spectrum.name,
                str(spectrum.parent),
                simple_output=True,
                silent=True,
            )
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

    def _require_fresh(self, stage):
        return self._stage_revision[stage] == self.revision

    def run(self, stage="all"):
        if stage not in self.STAGES:
            return self._failure("invalid_stage", f"Unsupported run stage: {stage}", set_error=False)
        operation_id = self._begin_mutation()
        if operation_id is None:
            return self._failure("busy", "Another session mutation is running", set_error=False)
        started = time.perf_counter()
        try:
            if not self.engine.config.filename:
                return self._failure("empty", "Open a spectrum before running UniDec", operation_id)

            stages = ("process", "deconvolve", "pick_peaks") if stage == "all" else (stage,)
            for current in stages:
                if current == "process":
                    self.engine.process_data(silent=True)
                    self._stage_revision["process"] = self.revision
                    self._stage_revision["deconvolve"] = None
                    self._stage_revision["pick_peaks"] = None
                elif current == "deconvolve":
                    if not self._require_fresh("process"):
                        return self._failure(
                            "stale_prerequisite", "Process data at the current revision first", operation_id)
                    result = self.engine.run_unidec(silent=True)
                    if result != 0:
                        return self._failure("native_failed", f"UniDec returned {result}", operation_id)
                    self._stage_revision["deconvolve"] = self.revision
                    self._stage_revision["pick_peaks"] = None
                else:
                    if not self._require_fresh("deconvolve"):
                        return self._failure(
                            "stale_prerequisite", "Deconvolve at the current revision first", operation_id)
                    self.engine.pick_peaks()
                    self._stage_revision["pick_peaks"] = self.revision

            self._last_duration = time.perf_counter() - started
            self.status = "complete" if self._require_fresh("pick_peaks") else "ready"
            return self._finish_success(operation_id, duration=self._last_duration)
        except Exception as error:
            return self._failure("run_failed", error, operation_id)
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
        plots = {
            "processed_spectrum": ("process", self.engine.data.data2, "m/z", "Intensity"),
            "mass_distribution": ("deconvolve", self.engine.data.massdat, "Mass (Da)", "Intensity"),
        }
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

"""Managed, headless session boundary for MetaUniDec collections."""

import hashlib
from pathlib import Path
import shutil
import tempfile
import threading
import uuid

import numpy as np

from unidec.metaunidec.mudeng import MetaUniDec


class MetaUniDecCollectionSession:
    """Open and summarize one MetaUniDec HDF5 project without mutating its source."""

    MAX_SPECTRUM_SUMMARIES = 100
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
        }

    def get_state(self):
        return self._state()

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
            self.status = "ready"
            return {
                "ok": True,
                "operation_id": operation_id,
                "updated": updated,
                "state": self._state(),
            }
        finally:
            self._lock.release()

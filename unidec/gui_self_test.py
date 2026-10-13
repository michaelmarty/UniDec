"""Packaged GUI smoke tests for the UniDec launcher."""

from dataclasses import dataclass
import gc
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np


PASS = "pass"
SKIP = "skip"
FAIL = "failure"
TIMEOUT = "timeout"
CANCELLED = "cancelled"


@dataclass(frozen=True)
class WorkerSpec:
    name: str
    label: str
    fixture: str
    application: str
    optional: bool = False
    companions: tuple = ()


WORKERS = (
    WorkerSpec("unidec", "UniDec", "ADH.txt", "unidec"),
    WorkerSpec("unidec-im", "UniDec IM", "IMMS/aqpz.dat", "unidec-im"),
    WorkerSpec("unidec-cd", "UniDecCD", "CDMS/BgalGroEL.dmt", "unidec-cd"),
    WorkerSpec("unidec-cd-thermo", "UniDecCD Thermo", "CDMS/GroEL_CDMS_1.RAW",
               "unidec-cd", optional=True),
    WorkerSpec("isodec", "IsoDec", "IsoDec/test2.txt", "isodec"),
    WorkerSpec("unichrom", "UniChrom native", "UniChrom/SEC_Native_Bispecific_Special.raw",
               "unichrom"),
    WorkerSpec("unichrom-state", "UniChrom saved state",
               "UniChrom/SEC_Native_Bispecific_Special.hdf5", "unichrom",
               companions=("UniChrom/SEC_Native_Bispecific_Special.raw",)),
    WorkerSpec("metaunidec", "MetaUniDec", "UniChrom/SEC_Native_Bispecific_Special.hdf5",
               "metaunidec"),
    WorkerSpec("unichrom-cd", "UniChromCD", "CDMS/BgalGroEL.dmt", "unichrom-cd"),
)
WORKER_BY_NAME = {worker.name: worker for worker in WORKERS}


def example_data_dir():
    return Path(__file__).resolve().parent / "bin" / "Example Data"


def resolve_fixture(relative_path, root=None):
    path = (Path(root) if root else example_data_dir()) / Path(relative_path)
    if not path.exists():
        raise FileNotFoundError(f"Packaged fixture is missing: {relative_path}")
    if path.is_file():
        with path.open("rb") as stream:
            header = stream.read(128)
        if header.startswith(b"version https://git-lfs.github.com/spec"):
            raise RuntimeError(f"Packaged fixture is an unresolved Git LFS pointer: {relative_path}")
    return path


def copy_fixture(source, workspace):
    destination = Path(workspace) / source.name
    if source.is_dir():
        shutil.copytree(source, destination)
    else:
        shutil.copy2(source, destination)
    return destination


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def _finite_nonempty(value, label):
    array = np.asarray(value)
    if array.size == 0:
        raise AssertionError(f"{label} is empty")
    try:
        if not np.all(np.isfinite(array.astype(float))):
            raise AssertionError(f"{label} contains non-finite values")
    except (TypeError, ValueError):
        raise AssertionError(f"{label} is not numeric") from None


def _construct_application(application):
    if application == "unidec":
        from unidec.GUniDec import UniDecApp
        return UniDecApp(ignore_args=True)
    if application == "unidec-im":
        from unidec.UniDecIM import UniDecIMApp
        return UniDecIMApp(ignore_args=True)
    if application == "unidec-cd":
        from unidec.UniDecCD import UniDecCDApp
        return UniDecCDApp(ignore_args=True)
    if application == "isodec":
        from unidec.IsoDecGUI import IsoDecPres
        return IsoDecPres()
    if application == "unichrom":
        from unidec.UniChrom import ChromApp
        return ChromApp(ignore_args=True)
    if application == "metaunidec":
        from unidec.MetaUniDec import UniDecApp
        return UniDecApp(ignore_args=True)
    if application == "unichrom-cd":
        from unidec.UniChromCD import UniChromCDApp
        return UniChromCDApp(ignore_args=True)
    raise KeyError(f"Unknown self-test application: {application}")


def _open_fixture(app, spec, fixture):
    if spec.application == "unidec-cd":
        # Opening CDMS data normally starts processing immediately. This worker
        # covers the file reader; processing has separate native tests.
        app.on_dataprep_button = lambda *args, **kwargs: None
        app.on_open_file(fixture.name, str(fixture.parent))
    elif spec.application in ("unidec", "unidec-im", "isodec", "unichrom-cd"):
        app.on_open_file(fixture.name, str(fixture.parent))
    elif spec.application == "unichrom":
        app.open_file(str(fixture))
    elif spec.application == "metaunidec":
        app.open_file(str(fixture))


def _validate_application(app, spec, fixture):
    if spec.application in ("unidec", "isodec"):
        _finite_nonempty(app.eng.data.data2, "spectrum")
    elif spec.application == "unidec-im":
        _finite_nonempty(app.eng.data.data3, "ion-mobility spectrum")
    elif spec.application in ("unidec-cd", "unichrom-cd"):
        _finite_nonempty(app.eng.farray, "charge-detection data")
    elif spec.application == "unichrom":
        _finite_nonempty(app.eng.ticdat, "chromatogram")
    elif spec.application == "metaunidec":
        spectra = app.eng.data.get_spectra()
        if not spectra:
            raise AssertionError("MetaUniDec dataset is empty")
        _finite_nonempty(spectra[0].data2, "MetaUniDec spectrum")

    recorded = getattr(app, "top_path", None)
    if spec.application == "metaunidec":
        recorded = app.eng.config.hdf_file
    elif spec.application == "unichrom":
        recorded = app.eng.path
    if not recorded or Path(recorded).resolve().parent != fixture.resolve().parent:
        raise AssertionError("application did not record the copied fixture path")


def run_worker(name, report_path):
    started = time.time()
    spec = WORKER_BY_NAME.get(name)
    if spec is None:
        result = {"name": name, "status": FAIL, "message": "Unknown worker"}
        _write_json(report_path, result)
        return 1

    workspace = Path(tempfile.mkdtemp(prefix=f"unidec-self-test-{name}-"))
    app = None
    result = {
        "name": spec.name,
        "label": spec.label,
        "fixture": spec.fixture,
        "status": FAIL,
        "workspace": str(workspace),
    }
    return_code = 1
    try:
        try:
            fixture = copy_fixture(resolve_fixture(spec.fixture), workspace)
            for companion in spec.companions:
                copy_fixture(resolve_fixture(companion), workspace)
        except FileNotFoundError as error:
            if spec.optional:
                result.update(status=SKIP, message=str(error))
                return_code = 2
            else:
                raise
        else:
            from unittest.mock import patch
            from unidec.modules.unidec_presbase import UniDecPres

            with patch.object(UniDecPres, "read_recent", return_value=[]), \
                    patch.object(UniDecPres, "write_to_recent", return_value=None), \
                    patch.object(UniDecPres, "cleanup_recent_file", return_value=None), \
                    patch.object(UniDecPres, "on_load_default", return_value=None):
                try:
                    app = _construct_application(spec.application)
                    _open_fixture(app, spec, fixture)
                    app.wx_app.Yield()
                    _validate_application(app, spec, fixture)
                except (ImportError, ModuleNotFoundError) as error:
                    if spec.optional:
                        result.update(status=SKIP, message=f"Optional reader unavailable: {error}")
                        return_code = 2
                    else:
                        raise
                else:
                    result.update(status=PASS, message="Opened and validated")
                    return_code = 0
    except Exception as error:
        result.update(status=FAIL, message=str(error), traceback=traceback.format_exc())
        return_code = 1
    finally:
        if app is not None and app.view is not None:
            try:
                for reader_name in ("importer", "chromdat"):
                    reader = getattr(app.eng, reader_name, None)
                    if reader is not None and hasattr(reader, "close"):
                        reader.close()
                        setattr(app.eng, reader_name, None)
                app.view.Destroy()
                app.wx_app.Yield()
            except Exception:
                pass
        app = None
        gc.collect()
        if result["status"] in (PASS, SKIP):
            try:
                shutil.rmtree(workspace)
                result["workspace"] = None
            except OSError as error:
                result.update(
                    status=FAIL,
                    message=f"Opened successfully but could not remove temporary workspace: {error}",
                )
                return_code = 1
        result["duration_seconds"] = round(time.time() - started, 3)
        _write_json(report_path, result)
    return return_code


def _worker_command(name, report_path):
    arguments = ["--gui-self-test-worker", name, "--report", str(report_path)]
    if getattr(sys, "frozen", False):
        return [sys.executable, *arguments]
    return [sys.executable, "-m", "unidec.gui_self_test", *arguments]


def _run_process(command, timeout, cancel_event=None):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    deadline = time.monotonic() + timeout
    while True:
        if cancel_event is not None and cancel_event.is_set():
            process.terminate()
            try:
                stdout, stderr = process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                stdout, stderr = process.communicate()
            return None, stdout, stderr, CANCELLED
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            process.kill()
            stdout, stderr = process.communicate()
            return None, stdout, stderr, TIMEOUT
        try:
            stdout, stderr = process.communicate(timeout=min(0.2, remaining))
            return process.returncode, stdout, stderr, None
        except subprocess.TimeoutExpired:
            pass


def default_report_path():
    stamp = time.strftime("%Y%m%d-%H%M%S")
    return Path(tempfile.gettempdir()) / f"unidec-gui-self-test-{stamp}.json"


def run_controller(report_path=None, timeout=120, callback=None, cancel_event=None):
    report_path = Path(report_path) if report_path else default_report_path()
    results = []
    with tempfile.TemporaryDirectory(prefix="unidec-self-test-reports-") as report_directory:
        for index, spec in enumerate(WORKERS, start=1):
            if cancel_event is not None and cancel_event.is_set():
                break
            worker_report = Path(report_directory) / f"{spec.name}.json"
            return_code, stdout, stderr, process_status = _run_process(
                _worker_command(spec.name, worker_report), timeout, cancel_event)
            if process_status:
                result = {
                    "name": spec.name,
                    "label": spec.label,
                    "fixture": spec.fixture,
                    "status": process_status,
                    "message": f"Worker {process_status}",
                }
            elif worker_report.is_file():
                result = json.loads(worker_report.read_text(encoding="utf-8"))
            else:
                result = {
                    "name": spec.name,
                    "label": spec.label,
                    "fixture": spec.fixture,
                    "status": FAIL,
                    "message": f"Worker exited {return_code} without a report",
                }
            if stdout:
                result["stdout"] = stdout
            if stderr:
                result["stderr"] = stderr
            results.append(result)
            if callback:
                callback(index, result)
            if process_status == CANCELLED:
                break

    if (cancel_event is not None and cancel_event.is_set()
            and (not results or results[-1]["status"] != CANCELLED)):
        result = {
            "name": "controller",
            "label": "GUI test suite",
            "fixture": None,
            "status": CANCELLED,
            "message": "Cancelled before the next worker started",
        }
        results.append(result)
        if callback:
            callback(len(results), result)

    summary = {
        status: sum(result["status"] == status for result in results)
        for status in (PASS, SKIP, FAIL, TIMEOUT, CANCELLED)
    }
    report = {
        "report_path": str(report_path),
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "results": results,
        "summary": summary,
        "success": summary[FAIL] == 0 and summary[TIMEOUT] == 0 and summary[CANCELLED] == 0,
    }
    _write_json(report_path, report)
    return report


def _argument_value(arguments, option, default=None):
    if option not in arguments:
        return default
    index = arguments.index(option)
    if index + 1 >= len(arguments):
        raise ValueError(f"{option} requires a value")
    return arguments[index + 1]


def main(arguments=None):
    arguments = list(sys.argv[1:] if arguments is None else arguments)
    report_path = _argument_value(arguments, "--report")
    if "--gui-self-test-worker" in arguments:
        name = _argument_value(arguments, "--gui-self-test-worker")
        if not report_path:
            raise ValueError("--gui-self-test-worker requires --report")
        return run_worker(name, report_path)
    if "--gui-self-test" in arguments:
        timeout = float(_argument_value(arguments, "--timeout", 120))
        report = run_controller(report_path, timeout=timeout)
        print(json.dumps(report["summary"], sort_keys=True))
        print(f"Report: {report['report_path']}")
        return 0 if report["success"] else 1
    raise ValueError("Expected --gui-self-test or --gui-self-test-worker")


if __name__ == "__main__":
    raise SystemExit(main())

# Conventional MS from Python

Use `unidec.engine.UniDec` for direct access to NumPy arrays and `unidec.session.UniDecSession` for managed files and structured results. Both examples use the packaged ADH spectrum and the native UniDec executable. Install UniDec in the Python environment you intend to use. The bundled text input requires no vendor runtime or GUI. Run the code blocks in order.

## Run the direct engine

For your own text input, supply one row per point with finite m/z and intensity values in ascending m/z order. Copy the input into an isolated directory. By default, `open_file` creates a sibling directory named after the spectrum with a `_unidecfiles` suffix. When reusing that directory it can load cached processed data and an existing configuration. Pass `refresh=True` to refresh the input instead. Set parameters after opening so an existing configuration cannot replace your choices.

```python
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from unidec.engine import UniDec
from unidec.session import UniDecSession


def copy_demo_spectrum(directory):
    source = files("unidec").joinpath("bin", "Example Data", "ADH.txt")
    path = directory / "ADH.txt"
    path.write_bytes(source.read_bytes())
    return path


with TemporaryDirectory() as directory:
    root = Path(directory)
    spectrum = copy_demo_spectrum(root)
    engine = UniDec(ignore_args=True, silent=True)
    engine.open_file(str(spectrum), silent=True)

    schema = engine.config.parameter_schema()
    assert "startz" in schema and "endz" in schema
    engine.config.update_parameters({"startz": 8, "endz": 24})
    if engine.process_data(silent=True) == 1:
        raise RuntimeError("UniDec rejected the processing settings")
    status = engine.run_unidec(silent=True)
    if status:
        raise RuntimeError(engine.last_run_diagnostic or f"UniDec exited with {status}")

    assert engine.data.data2.ndim == 2 and engine.data.data2.shape[1] == 2
    assert engine.data.massdat.ndim == 2 and engine.data.massdat.shape[1] == 2
    peaks = engine.pick_peaks()
    assert peaks.size and peaks.shape[1] == 2
    print(peaks[:3])  # Mass in Da and normalized height
    print("Fit R-squared:", engine.config.error)

    restore_directory = root / "restore"
    restore_directory.mkdir()
    archive = restore_directory / "checkpoint.zip"
    engine.save_state(str(archive))
    restored_engine = UniDec(ignore_args=True, silent=True)
    assert restored_engine.load_state(str(archive))
    np.testing.assert_array_equal(restored_engine.data.ztab, engine.data.ztab)
```

`open_file` resets the engine's data and peak containers, reads the input, and sets up working paths. Read `engine.data.rawdata` for the imported spectrum. `process_data` exports configuration and writes the processed spectrum for the native executable. `run_unidec` imports the fit and mass distribution on exit status zero. A failed result import can raise an exception even after the executable succeeds. `pick_peaks` writes the peak list when peaks are found and can rescale the mass distribution's intensity column.

Direct configuration updates change the config object in memory. Explicitly reprocess and rerun after changing settings to get results for the new configuration. The direct engine can replace an existing archive when saving. Use a trusted archive in a dedicated directory when restoring. Restoration extracts the working files and writes a reconstructed text spectrum beside the ZIP.

## Run a managed session

A session tracks which stages are available at the current revision. The example runs the three stages together, inspects a bounded peak list, writes a plot, and saves and restores a checkpoint. Use individual stages `process`, `deconvolve`, and `pick_peaks` when you need to inspect intermediate results.

```python
def require_ok(response):
    if not response["ok"]:
        raise RuntimeError(response["error"])
    return response


with TemporaryDirectory() as directory:
    root = Path(directory)
    spectrum = copy_demo_spectrum(root)
    session = UniDecSession(
        working_directory=root / "analysis",
        allowed_input_roots=[root],
    )
    require_ok(session.open_spectrum(spectrum))
    settings = session.get_parameters(["startz", "endz"])
    print("Initial charge range:", settings["values"])
    require_ok(session.set_parameters({"startz": 8, "endz": 24}))
    require_ok(session.run("all"))
    assert session.get_state()["available"]["pick_peaks"]

    results = require_ok(session.get_results(max_peaks=3))
    assert 0 < len(results["peaks"]) <= 3
    print(results["r_squared"], results["peaks"])
    plot = require_ok(session.get_plot("mass_distribution"))
    assert Path(plot["path"]).is_file()

    checkpoint = require_ok(session.save_state("checkpoint.zip"))
    restored_session = UniDecSession(
        working_directory=root / "restored",
        allowed_input_roots=[root],
    )
    require_ok(restored_session.load_state(checkpoint["path"]))
    assert restored_session.get_results()["error"]["code"] == "stale_results"
    require_ok(restored_session.run("all"))
    restored_results = require_ok(restored_session.get_results(max_peaks=3))
    np.testing.assert_allclose(
        [peak["mass"] for peak in restored_results["peaks"]],
        [peak["mass"] for peak in results["peaks"]],
    )
```

The session writes operation files and plots below its working directory. `processed_spectrum` needs a fresh processing stage and `mass_distribution` needs a fresh deconvolution stage. `get_results` requires fresh peak picking. A successful open, restore, or parameter update clears stage availability. A failed parameter update reports `invalid_parameters` and leaves the revision and parameter values as they were. Check `error.code` for failures such as `stale_prerequisite`, `native_failed`, or `busy`. Read-only `get_state` and `get_parameters` return their payloads directly. An unknown parameter name passed to `get_parameters` raises `ValueError`.

Checkpoint saves require a new ZIP path under the workspace. Restore copies the archive into a managed operation directory. The restored configuration and data are available to rerun, and `get_results` reports stale results until you do. Both examples remove their temporary files when the directory context exits.

See [configuration](../api/configuration.md) for validation rules and the generated parameter reference, and [data and peaks](../api/data-and-peaks.md) for the interpretation of the returned arrays and summaries.

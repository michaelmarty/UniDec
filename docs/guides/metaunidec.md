# MetaUniDec collections from Python

Use `unidec.meta_session.MetaUniDecCollectionSession` to deconvolve multiple related spectra on a common mass axis and extract peak signals across the collection. The headless workflow needs the native UniDec executable. For MCP, start the [agent interface](agent-interface.md) with `--mode metaunidec`.

## Inputs and project ownership

`open_collection(path)` expects a MetaUniDec HDF5 project in an allowed input directory. The project stores its configuration under `/config` and spectra under numbered groups in `/ms_dataset`. Each spectrum has a `raw_data` array shaped `(N, 2)` with m/z and intensity columns. Spectra can have different m/z axes and lengths. Spectrum names and experimental variables can be stored as attributes on the numbered groups. Use the MetaUniDec importer or the engine's collection writer to build this structure.

The session copies the HDF5 project to a new operation directory below its workspace. Native processing and analysis update that copy, and auxiliary exports go in its sibling `_unidecfiles` directory. If you use the direct engine, `unidec.metaunidec.mudeng.MetaUniDec.open(path)` instead works on the given project. The direct stages are `process_data()`, `run_unidec()`, and `pick_peaks()`.

## Run a two-spectrum collection

This example builds the input project with the bundled ADH and BSA spectra. All files are created inside a temporary directory. The session's only run stage, `"all"`, preprocesses spectra, runs native deconvolution, and builds aggregate grids and peak extracts. Apply common public settings through `set_parameters` before running.

```python
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from unidec.meta_session import MetaUniDecCollectionSession
from unidec.metaunidec.mudeng import MetaUniDec


def require_ok(response):
    if not response["ok"]:
        raise RuntimeError(response["error"])
    return response


with TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "collection.hdf5"
    examples = files("unidec").joinpath("bin", "Example Data")

    creator = MetaUniDec()
    creator.data.new_file(str(source))
    for name in ("ADH.txt", "BSA.txt"):
        spectrum = root / name
        spectrum.write_bytes(examples.joinpath(name).read_bytes())
        creator.data.add_file(path=str(spectrum))
    creator.config.write_hdf5(str(source))
    original = source.read_bytes()

    workspace = root / "analysis"
    session = MetaUniDecCollectionSession(
        working_directory=workspace, allowed_input_roots=[root],
    )
    opened = require_ok(session.open_collection(source))
    assert opened["state"]["collection"]["spectrum_count"] == 2
    require_ok(session.set_parameters({
        "startz": 5, "endz": 50, "masslb": 10000, "massub": 250000, "mzbins": 1,
    }))
    require_ok(session.run("all"))
    results = require_ok(session.get_results(max_peaks=5))["results"]
    assert 0 < len(results["peaks"]) <= 5
    print(results["mass_range"], results["peaks"])

    project = Path(session.engine.config.hdf_file)
    assert project.is_relative_to(workspace)
    with h5py.File(project, "r") as hdf:
        axis = hdf["/ms_dataset/mass_axis"]
        assert axis.shape == hdf["/ms_dataset/mass_sum"].shape
        assert axis.size == results["mass_points"]
        assert hdf["/ms_dataset/mass_grid"].size == 2 * axis.size
        assert hdf["/peaks/extracts"].shape == (2, results["peak_count"])
        assert np.isfinite(hdf["/peaks/extracts"][:]).all()

    checkpoint = require_ok(session.save_checkpoint("collection.hdf5"))
    restored = MetaUniDecCollectionSession(
        working_directory=root / "restored", allowed_input_roots=[workspace],
    )
    require_ok(restored.load_checkpoint(checkpoint["path"]))
    assert restored.get_results()["error"]["code"] == "stale_results"
    require_ok(restored.run())
    restored_results = require_ok(restored.get_results(max_peaks=5))["results"]
    assert restored_results["peak_count"] == results["peak_count"]
    np.testing.assert_allclose(
        [peak["mass"] for peak in restored_results["peaks"]],
        [peak["mass"] for peak in results["peaks"]],
    )
    assert source.read_bytes() == original
```

## Read the outputs

- `get_state()["collection"]` gives spectrum counts, variable names, and up to 100 spectrum summaries. `get_results(max_peaks=...)` gives aggregate mass bounds (Da), mass point count, total peak count, and up to 100 peaks ordered by height. Peak summaries contain mass, height, and DScore when available.
- In the managed HDF5 project, `/ms_dataset/mass_axis` and `mass_sum` hold the shared mass axis (Da) and summed intensity. The flat `mass_grid` stores one consecutive block of mass intensities per spectrum. On the engine, `data.massdat` is the two-column mass/intensity sum.
- `/peaks/peakdata` contains aggregate mass and height, with a score column when native scoring supplies it. `/peaks/extracts` has rows for spectra and columns for aggregate peaks. The engine exposes the transposed array as `data.exgrid`, so each `Peak.extracts` is one peak's signal across the collection.
- A successful open or settings update makes session results stale. `save_checkpoint` creates a new `.h5` or `.hdf5` snapshot inside the workspace. `load_checkpoint` opens it through the same managed-copy workflow. Rerun to make results current in the new session.

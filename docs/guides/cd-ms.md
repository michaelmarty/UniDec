# Charge detection MS from Python

Use `unidec.modules.CDEng.UniDecCD` for single-ion events. Its default workflow uses the native UniDec executable and runs without a GUI. The bundled `BgalGroEL.dmt` input can be read without a vendor runtime.

## Inputs and charge calibration

Use a supported event file such as DMT or a text export with one row per ion. For text files, put **m/z and measured intensity** first. The optional remaining columns are scan number, injection-time field, and acquisition time, in that order. The importer fills missing metadata for single-scan files. It interprets metadata by format, so preserve the original units and scaling when exporting. The engine keeps only events with positive intensity. The intensity supplies the charge estimate and needs the calibration appropriate to the input.

With a positive `CDslope` and `subtype=1`, charge is intensity divided by the slope. The noise-scaled branch (a negative slope or `subtype=0`) divides by the absolute slope times the estimated noise. Configure the m/z and charge ranges and both bin widths before processing. A nonzero `zzsig` requires unit charge bins (`CDzbins=1`). See the [parameter reference](../api/configuration.md) for validated settings.

## Run the bundled event file

Construct the session with the CD-MS adapter and pass it the same engine instance. The adapter uses the dedicated CD-MS stages and computes CD peak scores from the two-dimensional fit. For MCP, choose `--mode cdms` at startup. See the [agent interface](agent-interface.md).

```python
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from unidec.modules.CDEng import UniDecCD
from unidec.session import UniDecCDSessionEngineAdapter, UniDecSession


def require_ok(response):
    if not response["ok"]:
        raise RuntimeError(response["error"])
    return response


with TemporaryDirectory() as directory:
    root = Path(directory)
    source = files("unidec").joinpath("bin", "Example Data", "CDMS", "BgalGroEL.dmt")
    spectrum = root / "BgalGroEL.dmt"
    spectrum.write_bytes(source.read_bytes())

    engine = UniDecCD()
    session = UniDecSession(
        engine=engine,
        engine_adapter=UniDecCDSessionEngineAdapter(engine),
        working_directory=root / "analysis",
        allowed_input_roots=[root],
    )
    try:
        require_ok(session.open_spectrum(spectrum))
        assert session.get_state()["mode"] == "cdms"
        assert engine.darray.shape[1] == 5
        require_ok(session.set_parameters({
            "minmz": 5000.103515625, "maxmz": 19937.166015625, "mzbins": 10,
            "startz": 30, "endz": 80, "masslb": 5000, "massub": 1000000,
            "massbins": 100, "numit": 100, "mzsig": 5, "psig": 1,
            "psfun": 0, "psfunz": 0, "CDslope": 110750, "CDzbins": 1,
            "CDres": 0, "csig": 5, "zzsig": 1e-6, "subtype": 1,
        }))

        require_ok(session.run("process"))
        histogram_shape = (len(engine.ztab), len(engine.mz))
        assert engine.harray.shape == histogram_shape
        np.testing.assert_allclose(engine.data.data2[:, 1], engine.harray.sum(axis=0))

        require_ok(session.run("deconvolve"))
        assert engine.observed.shape == engine.fit2d.shape == histogram_shape
        assert np.isfinite(engine.fit2d).all()
        assert engine.data.massdat.shape[1] == 2
        assert np.isfinite(engine.data.massdat).all()

        require_ok(session.run("pick_peaks"))
        results = require_ok(session.get_results(max_peaks=3))
        assert 0 < len(results["peaks"]) <= 3
        assert all(0 <= peak["dscore"] <= 1 for peak in results["peaks"])
        print(results["r_squared"], results["peaks"])

        plot = require_ok(session.get_plot("mass_distribution"))
        assert Path(plot["path"]).is_file()
    finally:
        if getattr(engine, "importer", None) is not None:
            engine.importer.close()
```

For direct engine use, call `open_cdms_file(str(path))`, then apply calibration and settings. `process_data(transform=True)` returns `1` for a populated histogram or `0` when no histogram is produced. Then call `run_deconvolution(process_data=False)` on that histogram. Calling `run_deconvolution()` directly also performs preprocessing. Check that deconvolution produced a fit before calling `pick_peaks(calc_dscore=False)` and `dscore()`. The native route raises `RuntimeError` for a nonzero executable exit code.

## Interpret the results

- `engine.darray` holds the imported events after the positive-intensity filter. `engine.farray` is the working event array after preprocessing filters. `engine.zarray` contains the corresponding estimated charges.
- `engine.harray` has shape `(number of charge bins, number of m/z bins)`, with axes `engine.ztab` and `engine.mz`. The working histogram is replaced by the deconvolved result. The aligned `engine.observed` and `engine.fit2d` arrays retain the histogram used for fitting and its fitted intensities.
- `engine.data.data2` is a two-column m/z/intensity projection of the current histogram. `engine.data.massdat` contains mass in Da and intensity. After deconvolution, `engine.data.massgrid` can be reshaped to `(len(engine.data.massdat), len(engine.ztab))` for mass by charge.
- Session results summarize mass peaks and the clipped fit R-squared from the 2D fit. DScore describes peak quality. It ranges from zero to one. The underlying event arrays and histograms are available on the engine.

The direct engine normally writes to a sibling `_unidecfiles` directory. The session writes its versioned `rawdata.npz` event cache, configuration, native files, and PNGs below its own workspace. A session ZIP restores through that event cache. Restore with a CD-MS adapter and rerun before reading fresh results.

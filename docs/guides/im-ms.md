# Ion mobility MS from Python

Use `unidec.modules.IMEng.UniDecIM` for ion mobility spectra. The engine can run headlessly. Deconvolution requires the native UniDec executable. The example uses the bundled `aqpz.dat` text spectrum and needs no vendor runtime.

## Inputs and calibration

For text input, supply three numeric columns in this order: **m/z, drift or arrival time (ms), intensity**. A sparse list can omit zero-intensity bins. UniDec expands the input onto a dense rectangular grid and builds the raw m/z projection by summing intensities over time. Use finite values and coordinates from the instrument's sampled axes. Support for other formats, such as combined-ion-mobility mzML, depends on the importer.

Set instrument calibration before running real samples. The native engine uses this calibration to convert arrival time to collision cross section (CCS). Drift-tube and traveling-wave calibrations have different settings. The fixture below matches the engine's built-in drift-tube calibration defaults. The public [parameter reference](../api/configuration.md) covers time selection, smoothing and background subtraction, time peak width, and CCS bounds and bin size. Instrument calibration fields such as `twaveflag` belong to the engine configuration and are outside the validated session schema.

## Run the bundled spectrum

This example opens the file with an IM-MS engine inside `UniDecSession`. The temporary workspace contains the raw cache, configuration, native input and output, and PNG. With the direct engine, use `open_file`, `process_data`, `run_unidec`, and `pick_peaks` in that order. Direct opening normally uses a sibling `_unidecfiles` directory. For MCP, choose `--mode imms` at startup. See the [agent interface](agent-interface.md).

```python
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from unidec.modules.IMEng import UniDecIM
from unidec.session import UniDecSession


def require_ok(response):
    if not response["ok"]:
        raise RuntimeError(response["error"])
    return response


with TemporaryDirectory() as directory:
    root = Path(directory)
    source = files("unidec").joinpath("bin", "Example Data", "IMMS", "aqpz.dat")
    spectrum = root / "aqpz.dat"
    spectrum.write_bytes(source.read_bytes())

    engine = UniDecIM(ignore_args=True, silent=True)
    session = UniDecSession(
        engine=engine,
        working_directory=root / "analysis",
        allowed_input_roots=[root],
    )
    require_ok(session.open_spectrum(spectrum))
    assert session.get_state()["mode"] == "imms"
    require_ok(session.set_parameters({
        "minmz": 4504, "maxmz": 11988, "mzbins": 4,
        "startz": 10, "endz": 18, "mzsig": 5,
        "mindt": 10.08, "maxdt": 15.12, "subbufdt": 2, "dtsig": 0.2,
        "ccslb": 100, "ccsub": 25000, "ccsbins": 20,
    }))

    require_ok(session.run("process"))
    grid = engine.data.data3
    mz_count = len(np.unique(grid[:, 0]))
    time_count = len(np.unique(grid[:, 1]))
    assert grid.shape == (mz_count * time_count, 3)

    require_ok(session.run("deconvolve"))
    data = engine.data
    assert data.massdat.shape[1] == data.ccsdata.shape[1] == 2
    assert np.isfinite(data.ccsdata).all()
    assert data.massccs.shape == (len(data.massdat), len(data.ccsdata))
    assert data.ccsz.shape == (len(data.ztab), len(data.ccsdata))
    assert data.mztgrid.shape == (mz_count, time_count, len(data.ztab))

    require_ok(session.run("pick_peaks"))
    results = require_ok(session.get_results(max_peaks=3))
    assert 0 < len(results["peaks"]) <= 3
    print(results["peaks"])

    plot = require_ok(session.get_plot("ccs_distribution"))
    assert Path(plot["path"]).is_file()
```

## Interpret the results

- `data.rawdata3` and `data.data3` hold raw and processed triplets, respectively. Each row is one m/z and time bin. Time varies fastest within each m/z block. `data.rawdata` is the raw two-column m/z projection. For a processed projection, sum `data.data3[:, 2]` across time for each m/z or request the session's `processed_spectrum` plot.
- `data.massdat` and `data.ccsdata` are floating-point axis/intensity pairs. The axes are mass in Da and CCS in Å². The session's `mass_distribution` and `ccs_distribution` plots are available after deconvolution.
- In a full import, `data.massccs` is indexed by mass then CCS, `data.ccsz` by charge then CCS, and `data.mztgrid` by m/z then time then charge. `data.fitdat2d` has the same coordinates as the processed triplets with fitted intensities in column three. `run_unidec(efficiency=True)` skips the larger grids and fits.
- `get_results` summarizes mass peaks and fit R-squared. The CCS arrays are available on the engine. Peak picking can rescale the mass distribution's intensities.

Session ZIP checkpoints preserve files needed to rerun. Restore into a session constructed with `UniDecIM`. Run the stages again before requesting fresh results or plots.

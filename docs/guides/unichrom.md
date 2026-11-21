# UniChrom time windows from Python

Use `unidec.unichrom_session.UniChromSession` to turn a raw LC-MS chromatogram into a collection of averaged spectra and deconvolve it. This headless workflow requires the native UniDec executable and an importer for the source format. For MCP, start the [agent interface](agent-interface.md) with `--mode unichrom`.

## Inputs and managed projects

Pass a raw source such as a multi-scan mzML file or an importer-supported vendor directory under an allowed input root. Opening the source reads scan metadata and builds total ion chromatogram (TIC) and base peak chromatogram (BPC) summaries. It creates an initially empty HDF5 project under a new operation directory in the workspace. The original source stays in place and the importer reads from it when generating windows.

`create_time_windows(width, start=None, end=None)` averages scans into regular windows in the importer's time units, minutes for mzML. Use finite positive widths. Supply both optional bounds or omit both to use the chromatogram's range. The bounds select window start times. Each window requests a full interval of length `width`, so the last interval can extend past the specified `end`. Creating windows replaces the spectra in the managed project and makes previous results stale. Each numbered spectrum in `/ms_dataset` stores an `(N, 2)` m/z and intensity array in `raw_data` and attributes such as `timemid` and `scanstart`.

Use `run("all")` after window creation. It preprocesses the window spectra, runs native deconvolution, and extracts aggregate mass peaks. After changing window settings, call `create_time_windows` again to rebuild the spectra. Set `dtsig=0` for independent window deconvolution. A positive `dtsig` selects coupled chromatographic deconvolution. Its aggregate grids and extracts use the same project structure, while charge-resolved per-window outputs are unavailable in that mode.

## Run a small mzML chromatogram

The example creates a real six-scan mzML from the packaged ADH spectrum with fixed scan times and intensity multipliers. Install the optional fixture writer with `python -m pip install psims` to run it. The regular analysis workflow can open an existing mzML without that writer.

```python
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from psims.mzml.writer import MzMLWriter

from unidec.unichrom_session import UniChromSession


def require_ok(response):
    if not response["ok"]:
        raise RuntimeError(response["error"])
    return response


with TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "chromatogram.mzML"
    reference = files("unidec").joinpath("bin", "Example Data", "ADH.txt")
    with reference.open("rb") as stream:
        mz, intensity = np.loadtxt(stream, unpack=True)

    scales = np.array([0.25, 0.5, 1.0, 0.75, 0.5, 0.25])
    times = np.arange(len(scales)) * 0.5
    with source.open("wb") as stream, MzMLWriter(stream) as writer:
        writer.controlled_vocabularies()
        writer.file_description(["MS1 spectrum"], [])
        writer.software_list([{
            "id": "fixture", "version": "1",
            "params": ["custom unreleased software tool"],
        }])
        writer.instrument_configuration_list([])
        writer.data_processing_list([{
            "id": "fixture_processing",
            "processing_methods": [{"order": 0, "software_reference": "fixture",
                                    "params": ["Conversion to mzML"]}],
        }])
        with writer.run(id="fixture_run"):
            with writer.spectrum_list(count=len(scales)):
                for index, (time, scale) in enumerate(zip(times, scales), start=1):
                    writer.write_spectrum(
                        mz, intensity * scale, id=f"scan={index}",
                        scan_start_time=time, params=[{"ms level": 1}],
                        centroided=False,
                    )
            with writer.chromatogram_list(count=1):
                writer.write_chromatogram(
                    times, intensity.sum() * scales, id="TIC",
                    chromatogram_type="total ion current chromatogram",
                )

    original = source.read_bytes()
    workspace = root / "analysis"
    session = UniChromSession(working_directory=workspace, allowed_input_roots=[root])
    try:
        opened = require_ok(session.open_chromatogram(source))
        assert opened["state"]["chromatogram"]["scan_count"] == len(scales)
        np.testing.assert_allclose(session.engine.ticdat[:, 0], times)
        windows = require_ok(session.create_time_windows(width=1.0, start=0.0, end=1.0))
        assert windows["state"]["chromatogram"]["window_count"] == 2
        require_ok(session.set_parameters({
            "startz": 5, "endz": 50, "masslb": 10000, "massub": 250000,
            "mzbins": 1, "dtsig": 0,
        }))
        require_ok(session.run())
        results = require_ok(session.get_results(max_peaks=5))["results"]
        assert results["window_count"] == 2
        assert 0 < len(results["peaks"]) <= 5
        print(results["mass_range"], results["peaks"])

        project = Path(session.engine.config.hdf_file)
        assert project.is_relative_to(workspace)
        with h5py.File(project, "r") as hdf:
            assert hdf["/ms_dataset/0/raw_data"].shape[1] == 2
            assert hdf["/ms_dataset/0"].attrs["timemid"] < hdf["/ms_dataset/1"].attrs["timemid"]
            assert hdf["/ms_dataset/mass_axis"].size == results["mass_points"]
            assert hdf["/peaks/extracts"].shape == (2, results["peak_count"])
            assert np.isfinite(hdf["/peaks/extracts"][:]).all()
        assert source.read_bytes() == original
        assert not source.with_suffix(".hdf5").exists()
    finally:
        importer = getattr(session.engine, "chromdat", None)
        if importer is not None:
            importer.close()
```

## Read the outputs

`get_state()["chromatogram"]` gives bounded TIC/BPC metadata, scan range and count, and window count. For full chromatogram arrays, read `engine.ticdat` and `engine.bpcdat`. They contain time and intensity pairs. `get_results(max_peaks=...)` returns the window count, aggregate mass bounds (Da), a mass point count, UniScore, and up to 100 peaks ordered by height. The HDF5 project stores the aggregate mass axis and intensities in `/ms_dataset` and one extraction row per window in `/peaks/extracts`.

The managed HDF5 project is a derived collection. Keep the original raw source to regenerate windows. UniChrom's session has no portable checkpoint operation. For direct engine use, `unidec.modules.ChromEng.ChromEngine.open_chrom` accepts an explicit `output_path`. The default derives the HDF5 path from the source name.

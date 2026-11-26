# Python-to-native boundary

The normal Python engines launch the UniDec executable using an argument list. Conventional MS, IM-MS, and native CD-MS use `unidec.tools.unidec_call`. It passes an absolute text configuration path, sets the subprocess working directory to the configuration directory, and captures a diagnostic and numeric exit status. MetaUniDec and UniChrom call the same executable with an HDF5 project path and, where needed, a stage argument.

## Configuration and inputs

| Boundary | Python writer and native reader |
| --- | --- |
| Text configuration | `UniDecConfig.config_export` writes space-separated key/value lines. `LoadConfig` in `udstruct.c` reads the keys. Managed paths such as `input` and `output` are written relative to the configuration file when possible and resolved from its directory. Other settings can refer to text sidecars. The configuration object supplies paths for the current short or legacy-prefixed output names. |
| HDF5 configuration | `UniDecConfig.write_hdf5` writes typed attributes under `/config`. `mh5LoadConfig` in `h5io.c` reads the native subset and retains initialized values for missing attributes. Array settings such as a mass list use datasets. The Python reader also restores settings used only by Python or the GUI. |

The reviewed scalar persistence registry supplies types, serialized keys, aliases, and unset values for a subset of configuration fields. It is distinct from the public parameter schema, which filters automation settings by acquisition mode. Check both serializers and both C readers before changing a field. For example, the Python setting `psfunz` writes `zpsfn` in text and `psfunz` in HDF5. Conventional MS derives charge count from `endz - startz + 1`. For variable peak widths, supply one finite positive FWHM per processed m/z point through a text sidecar or the HDF5 spectrum group's `fwhm` dataset.

The conventional text route reads m/z and intensity from the prepared `input.dat` or its legacy-prefixed equivalent. IM-MS writes flattened binary float32 triplets for m/z, arrival time, and intensity. The single-histogram CD-MS path writes a text table for m/z, charge, and histogram intensity from the calibrated event histogram. In HDF5, conventional single-spectrum deconvolution reads `/ms_data/processed_data` with two columns. A collection stores per-spectrum arrays under `/ms_dataset/<index>`, including `raw_data` for preprocessing and `processed_data` for deconvolution.

## Native dispatch

The executable's `ImportConfig` currently selects HDF5 when the path contains `.hdf5`. Use lowercase `.hdf5` for project files passed to the executable, including projects passed through a Python engine. After loading configuration, `main` in `UniDec.c` dispatches as follows.

| Configuration | Command path |
| --- | --- |
| Single-spectrum mode (`metamode=-2`), conventional settings | Conventional `run_unidec`, using the supplied text or HDF5 inputs. |
| Single-spectrum text mode with `imflag=1` | `run_unidec_IM`. |
| Single-spectrum text mode with `cdmsflag=1` or `2` | `run_unidec_CD` for a single histogram or `run_unidec_UCCD` for a chromatographic CD-MS stack. |
| HDF5 collection, `dtsig=0` or absent | `run_metaunidec`. The normal no-flag call preprocesses and deconvolves per spectrum. `-proc` preprocesses, `-decon` uses prepared spectra, and `-grids` builds or reads aggregate grids and extracts peaks. |
| HDF5 collection, positive `dtsig` | Coupled `run_chromatogram`. Preprocessing and extraction commands retain the collection path. A grid refresh can reuse the coupled mass grids. Charge-resolved per-window extraction commands are unavailable in the coupled route. |

## Result imports

- **Conventional text.** `UniDec.unidec_imports` requires a nonempty mass/intensity table. In a full import, the native float32 fit has one value per prepared m/z point. The flat mass grid has one value per mass/charge pair, with charge varying fastest. Python combines the m/z grid values with their m/z and charge coordinates. The conventional reader uses the error file and processed spectrum to calculate fit R-squared. Optional outputs depend on the mode and settings, including `rawflag` and Python's `efficiency` import option.
- **IM-MS and CD-MS text.** IM-MS imports mass and CCS distributions and reshapes mobility grids using the processed m/z and arrival-time axes and charge count. CD-MS reshapes the native fit and deconvolved output back to the charge-by-m/z histogram layout and transforms the output to the mass axis. These routines also interpret run-specific diagnostics or fit evidence. See the [IM-MS guide](../guides/im-ms.md) and [CD-MS guide](../guides/cd-ms.md) for Python result shapes.
- **HDF5.** Conventional native output includes `/ms_data/mass_data` and, for full outputs, flat mass and m/z grids. The native HDF5 route and the text importer have different output sets. Collection Python code reloads the per-spectrum arrays and reads the common mass axis, mass sum, and extracts. `MetaDataSet` turns the mass axis and sum into `data.massdat` and transposes the extraction matrix from spectrum-by-peak on disk to peak-by-spectrum in memory.

## Example: compare conventional text and HDF5 runs

This diagnostic example runs the same prepared spectrum through both configuration formats. It uses the bundled BSA input and executable, with all working files in a temporary directory. The HDF5 path calls the executable directly and reads its datasets explicitly.

```python
from importlib.resources import files
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from unidec.engine import UniDec

with TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "BSA.txt"
    sample = files("unidec").joinpath("bin", "Example Data", "BSA.txt")
    source.write_bytes(sample.read_bytes())

    engine = UniDec(ignore_args=True, silent=True)
    engine.open_file(str(source), silent=True)
    config = engine.config
    config.update_parameters({"startz": 8, "endz": 24, "mzbins": 1})
    if engine.process_data(silent=True) == 1:
        raise RuntimeError("Preprocessing failed")
    status = engine.run_unidec(silent=True)
    if status:
        raise RuntimeError(engine.last_run_diagnostic)

    prepared = np.loadtxt(config.infname, ndmin=2)
    charge_count = config.endz - config.startz + 1
    assert engine.data.fitdat.shape == (len(prepared),)
    assert engine.data.mzgrid.shape == (len(prepared) * charge_count, 3)

    project = root / "same-input.hdf5"
    config.write_hdf5(str(project))
    with h5py.File(project, "a") as hdf:
        assert hdf["config"].attrs["metamode"] == -2
        assert hdf["config"].attrs["endz"] == config.endz
        hdf.require_group("ms_data").create_dataset("processed_data", data=prepared)

    result = subprocess.run(
        [config.UniDecPath, str(project)], cwd=root,
        capture_output=True, text=True, timeout=60,
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)

    with h5py.File(project, "r") as hdf:
        group = hdf["ms_data"]
        mass = group["mass_data"][:]
        assert group["mass_grid"].size == len(mass) * charge_count
        assert group["mz_grid"].size == len(prepared) * charge_count
        assert np.isfinite(mass).all()
        np.testing.assert_allclose(
            group["mass_grid"][:], engine.data.massgrid, rtol=1e-5, atol=1e-6,
        )
        np.testing.assert_allclose(
            group["mz_grid"][:], engine.data.mzgrid[:, 2], rtol=1e-5, atol=1e-6,
        )
    np.testing.assert_allclose(mass, engine.data.massdat, rtol=1e-5, atol=1e-6)
```

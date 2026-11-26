# Data and peaks

**Tier 1.** Read these containers to inspect results from an engine. The containers are mutable Python objects that can be inspected without wxPython or the native executable. The stages that populate them have their own [requirements](modes.md).

| Import path | Where to look |
| --- | --- |
| `unidec.modules.unidecstructure.DataContainer` | On a single-spectrum engine, `engine.data` holds arrays such as `rawdata`, `data2` (processed m/z spectrum), and `massdat` (deconvolved mass spectrum). IM-MS adds mobility arrays such as `rawdata3` and `data3`. |
| `unidec.modules.peakstructure.Peaks` and `unidec.modules.peakstructure.Peak` | After peak picking, `engine.pks.peaks` gives individual peaks with fields such as `mass` and `height`. |
| Session `get_state()` and `get_results(max_peaks=...)` | Use these to inspect freshness and retrieve bounded result dictionaries for single spectra, collections, and chromatograms. |

Arrays can be empty before their stage runs. Read the results after a successful run and copy arrays you need to retain independently of later engine operations. Reading a container or a session result summary does not write an export.

See the [IM-MS guide](../guides/im-ms.md) for mobility and CCS arrays and the [CD-MS guide](../guides/cd-ms.md) for event, histogram, and fit arrays.

See the [MetaUniDec guide](../guides/metaunidec.md) for aggregate mass grids and peak extracts and the [UniChrom guide](../guides/unichrom.md) for chromatogram and window outputs.

## Conventional-MS result fields

- `rawdata` and `data2` are floating-point arrays with m/z in the first column and intensity in the second. The processed spectrum can have a different point count or grid after selection and preprocessing.
- `massdat` has mass in daltons in column zero and signal intensity in column one. Its intensity scale can change during peak picking. `ztab` contains the integer charge states from `startz` through `endz`, inclusive.
- In a full native import, `fitdat` is a one-dimensional floating-point array of fitted intensities aligned with the processed m/z axis in `data2[:, 0]`. `engine.config.error` is the fit R-squared.
- `pick_peaks()` returns a floating-point array of detected mass and height pairs when peaks are found. It returns an empty array when there are no detections. `Peak.mass` is in daltons, `Peak.height` follows the configured normalization, and `Peak.integral` records the integrated area when calculated separately. `Peak.dscore` holds the peak quality score when calculated.
- Session `get_results(max_peaks=N)` returns the first `N` peaks in engine order, the total peak count, a truncation flag, and fit R-squared. Use a nonnegative integer. A limit of zero returns summary data with an empty peak list.

To inspect the peak container independently of a deconvolution, populate it with synthetic mass and height pairs.

```{code-block} python
:name: inspect-peak-container

from unidec.modules.peakstructure import Peaks

peaks = Peaks()
peaks.add_peaks([[100.0, 4.0], [110.0, 2.0]])
assert [(peak.mass, peak.height) for peak in peaks.peaks] == [(100.0, 4.0), (110.0, 2.0)]
```

## Supported Python reference

| Entry point | Supported operations |
| --- | --- |
| {any}`DataContainer <unidec.modules.unidecstructure.DataContainer>` | Holds mutable spectrum, grid, and fit arrays. |
| {any}`Peak <unidec.modules.peakstructure.Peak>` | Holds measurements for a single detected peak. |
| {any}`Peaks <unidec.modules.peakstructure.Peaks>` | {any}`add_peaks <unidec.modules.peakstructure.Peaks.add_peaks>` |

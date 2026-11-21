# Scientific analysis from Python

These Tier 2 APIs operate on spectra, peak lists, and extracted traces. The examples use synthetic data and run without wxPython, a vendor runtime, or the native executable. Engine scoring methods can also use data produced by native deconvolution. Start with the [mode guides](../api/modes.md) to produce those inputs.

## Process and quantify a spectrum

Import these helpers from `unidec.tools`. Supply a floating-point array with shape `(N, 2)`, an ascending x axis, and intensity in the second column. The x axis can be m/z or mass in Da.

| Callable | Contract |
| --- | --- |
| `datachop(data, newmin, newmax)` | Selects rows within the inclusive x bounds and keeps input order. |
| `normalize(data)` | Divides intensity by its maximum in place and returns the same array. Zero-maximum data are left alone. Use a floating-point copy to retain the original scale. |
| `gsmooth(data, sig)` | Applies Gaussian smoothing to the intensity column in place. The standard deviation is in samples, so use a regularly spaced axis for a constant physical width. |
| `integrate(data, start, end)` | Returns the trapezoidal area and the selected rows with `start <= x < end`. The area has intensity times x-axis units. An empty selection gives zero. A single selected point returns its intensity. |
| `peakdetect(data, window=..., threshold=..., norm=...)` | With no configuration override or ppm window, use an integer half-window in samples. For `norm=True`, a threshold between zero and one is relative to maximum intensity. For `norm=False`, use an absolute intensity threshold. The cutoff is strict. The result is in input-axis order, with columns for x and original height. No detections produce an empty one-dimensional array. |

```{code-block} python
:name: preprocess-spectrum

import numpy as np

from unidec import tools

spectrum = np.column_stack((np.arange(5.0), [0.0, 2.0, 4.0, 2.0, 0.0]))
cropped = tools.datachop(spectrum, 1.0, 3.0)
area, selected = tools.integrate(spectrum, 1.0, 4.0)
np.testing.assert_array_equal(cropped, selected)
assert area == 6.0

scaled = tools.normalize(spectrum.copy())
np.testing.assert_allclose(scaled[:, 1], [0.0, 0.5, 1.0, 0.5, 0.0])
np.testing.assert_array_equal(
    tools.peakdetect(spectrum, window=1, threshold=0.5, norm=True),
    [[2.0, 4.0]],
)
```

## Inspect fit and peak scores

For conventional MS, use `unidec.engine.UniDec.dscore()` after deconvolution with the full output import and after picking peaks. It updates peak fields in place. `Peak.dscore` multiplies the mass-shape, uniqueness, FWHM, and charge-distribution components. `pks.uniscore` is the mean of peak DScores weighted by peak height raised to `power`, multiplied by fit R-squared (`config.error`). These quantities describe fit and peak evidence. Consider how acquisition mode, scoring settings, and signal quality affect comparisons across runs.

For CD-MS, use `unidec.modules.CDEng.UniDecCD.dscore()` after native or Python single-histogram deconvolution and peak picking. It requires aligned `observed` and `fit2d` histograms when peaks are present. Without that evidence it raises `ValueError`. The optional `uscore_method` selects convolution-based fit and purity scoring (`1`, default), local area agreement (`2`), or a neutral UScore of one (`0`). The score uses a positive finite FWHM window multiplier and nonnegative finite `power`. The [CD-MS guide](cd-ms.md) demonstrates fit validation and peak scoring.

For an extracted chromatographic peak, use `unidec.modules.ChromEng.chromatographic_peak_scores(time, intensity)` directly. Supply equal-length finite vectors with at least three samples and strictly increasing times. The first score is the ratio of the shorter to the longer half width at half maximum. The second penalizes signal rises when moving away from the apex within that region. Both are between zero and one. If either half-height crossing is missing, or the trace is invalid or has no positive peak, the function returns `(0.0, 0.0)`. It neither changes the inputs nor writes files.

```{code-block} python
:name: score-chromatographic-peak

import numpy as np

from unidec.modules.ChromEng import chromatographic_peak_scores

scores = chromatographic_peak_scores(
    [0, 1, 2, 3, 4, 5], [0.0, 0.5, 1.0, 0.75, 0.5, 0.0],
)
np.testing.assert_allclose(scores, [0.5, 1.0])
assert chromatographic_peak_scores([0, 1, 2], [0, 0.5, 1]) == (0.0, 0.0)
```

## Fit single peaks and mixtures

Use `unidec.modules.fitting.fit_peak` for one isolated peak plus a constant background. Pass equal-length one-dimensional x and intensity arrays and initial guesses for the center, positive FWHM, area, and background. Shape codes are Gaussian (`0`), Lorentzian (`1`), and split Gaussian/Lorentzian (`2`). The returned tuple contains the parameters in **FWHM, center, area, background** order, their estimated standard errors in the same order, and fitted intensities on the original x axis. Width and center have the units of x. Area has intensity times x units. Use `isolated_peak_fit` to estimate starting values from an isolated slice. SciPy can raise if optimization fails or warn if it cannot estimate parameter uncertainty.

For multiple peaks with a shared shape code and a baseline already removed, use `unidec.modules.MassFitter.MassFitter`. Supply an `(N, 2)` mass distribution. A two-column guess array contains center and an intensity seed per peak. The `"smallguess"` constructor option seeds each width at 0.5 x units. `perform_fit("sort")` stores and returns the fitted intensities and a center-sorted `(P, 3)` parameter array in **center, FWHM, area** order. Both fitting examples operate in memory.

```{code-block} python
:name: fit-peak-models

import numpy as np

from unidec import tools
from unidec.modules.fitting import fit_peak, psfit
from unidec.modules.MassFitter import MassFitter

mass = np.linspace(92.0, 108.0, 321)
single = psfit(mass, s=1.2, m=100.0, a=3.0, b=0.2, psfun=0)
params, errors, fitted = fit_peak(
    mass, single, psfun=0, midguess=99.8, fwhmguess=1.0, aguess=2.5, bguess=0.1,
)
np.testing.assert_allclose(params, [1.2, 100.0, 3.0, 0.2], atol=1e-6)
assert errors.shape == (4,)
np.testing.assert_allclose(fitted, single, atol=1e-6)

mixture = (
    3.0 * tools.make_peak_shape(mass, 0, 1.0, 98.0, norm_area=True)
    + 2.0 * tools.make_peak_shape(mass, 0, 1.5, 102.0, norm_area=True)
)
guesses = np.array([[98.2, 3.0], [101.8, 2.0]])
fitter = MassFitter(np.column_stack((mass, mixture)), guesses, 0, "smallguess")
fit_data, peak_params = fitter.perform_fit("sort")
np.testing.assert_allclose(peak_params, [[98.0, 1.0, 3.0], [102.0, 1.5, 2.0]], atol=1e-5)
np.testing.assert_allclose(fit_data, mixture, atol=1e-5)
```

## Match picked peaks to candidate masses

Use `unidec.modules.matchtools.site_match` to annotate a `unidec.modules.peakstructure.Peaks` object. Pass a nonempty candidate mass vector, an integer matrix with one candidate per row and one index per site, a list of option names for each site, and the site labels. Use consistent mass units for observed and candidate masses, typically Da. The closest candidate sets `Peak.match` and the signed observed-minus-candidate error sets `Peak.matcherror`. A label is assigned only when the absolute error is strictly below `tolerance`, or when no tolerance is given. The return array has four rows (observed masses, matched masses, errors, labels). NumPy can coerce mixed values to strings. Read numerical match values from the peak fields.

`get_sitematch_list` can build that candidate matrix from a pandas table of option names, masses, and site probabilities. It enumerates combinations, so the candidate count grows with the product of option counts. `sitematch_to_target` filters candidate rows within a strict mass tolerance without mutating peaks. For UniDec Processing Pipeline tables, `calc_pairs` produces candidate masses and labels and `UPP_check_peaks` annotates the peaks and adds quantitation fields to the supplied row. It uses peak heights by default or already-populated integrals with `integrate=True`. These matching calls operate in memory. `sitematch_to_excel` writes a spreadsheet.

```{code-block} python
:name: match-candidate-masses

import numpy as np

from unidec.modules.matchtools import site_match
from unidec.modules.peakstructure import Peaks

peaks = Peaks()
peaks.add_peaks(np.array([[100.25, 3.0], [109.0, 2.0]]))
matches = site_match(
    peaks, np.array([100.0, 110.0]), np.array([[0], [1]]),
    onames=[["A", "B"]], sitenames=["Site1"], tolerance=1.0,
)
assert matches.shape == (4, 2)
assert [peak.label for peak in peaks.peaks] == ["[Site1:A]", ""]
np.testing.assert_allclose([peak.match for peak in peaks.peaks], [100.0, 110.0])
np.testing.assert_allclose([peak.matcherror for peak in peaks.peaks], [0.25, -1.0])
```

## Calculate mass defects and maps

In `unidec.tools.simple_mass_defect`, divide mass by a positive reference repeat mass and subtract its nominal repeat count. Use `centermode=1` for flooring (fractional defect from zero up to one) or `centermode=0` for rounding (centered around zero). With `normtype=1` the residual is converted back to Da when masses are in Da. With `normtype=0` it stays fractional.

Use `unidec.tools.kendrick_analysis` for a mass distribution. Provide a finite, ascending `(N, 2)` array of mass and nonnegative intensity with positive signal in the selected range, plus a positive reference repeat mass. Optional `massrange` crops the data inclusively. `nbins` sets the number of defect-axis bins. `transformmode=0` sums intensities into nearest defect bins and `transformmode=1` interpolates intensities onto the map. The intensity grid is normalized by its maximum. Use a distribution with variation in the defects and signal on the mapped grid to obtain a finite, resolved map.

The return tuple contains the defect projection `(defect, summed normalized intensity)`, the flattened map `(nominal mass, defect, normalized intensity)`, and three aligned grids for nominal mass, defect, and intensity. The grids have one row per distinct nominal repeat count and `nbins` columns. Nominal mass is always repeat count times reference mass. `xaxistype=0` selects fractional **defect** and `xaxistype=1` selects defect in mass units.

```{code-block} python
:name: map-mass-defects

import numpy as np

from unidec.tools import kendrick_analysis, simple_mass_defect

assert simple_mass_defect(87.0, 8.0, centermode=0, normtype=1) == -1.0
spectrum = np.array([[81.0, 1.0], [82.0, 2.0], [89.0, 3.0], [90.0, 4.0]])
projection, points, mass_grid, defect_grid, intensity_grid = kendrick_analysis(
    spectrum, 8.0, centermode=1, nbins=2, transformmode=0, xaxistype=1,
)
np.testing.assert_allclose(mass_grid[:, 0], [80.0, 88.0])
np.testing.assert_allclose(defect_grid[0], [1.0, 2.0])
np.testing.assert_allclose(intensity_grid, [[0.25, 0.5], [0.75, 1.0]])
np.testing.assert_allclose(projection, [[1.0, 1.0], [2.0, 1.5]])
assert points.shape == (4, 3)
```

Both helpers operate in memory. On an engine, `kendrick_peaks` updates the reference setting and each picked peak's defect fields. `kendrick_continuous` reads `data.massdat` and writes `1D_Mass_Defects.txt` and `2D_Mass_Defects.txt` in the engine's output directory.

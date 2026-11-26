# UniChrom and UniDec Harmonization Plan

This plan covers deconvolution and post-fit output in coupled UniChrom
(`dtsig > 0`). With `dtsig = 0`, MetaUniDec runs the standard UniDec
deconvolution core independently on each window. Preprocessing and peak
picking are outside the scope. [Native dispatch](UniDec.c#L70)

## 1. Apply charge scaling and fit constraints

- **`orbimode`**. Reuse UniDec's charge scaling by `1/z` after fitting and
  before reconvolution and conversion to mass. The iterative fit stays the
  same. [UniDec output handling](udmain.c#L356)
- **`manualfileflag`, `mfileflag`, and `mtabsig`**. Reuse manual charge
  assignments and mass-list rules when constructing the allowed m/z–charge
  mask. Build it on the common m/z axis for the FFT path and on each window's
  processed axis for the direct path. Preserve UniDec's rule precedence and
  the zero-tolerance behavior that selects the nearest m/z bin for each
  target mass and charge. [Manual assignments](udio.c#L586) ·
  [Mass-list rules](udcore.c#L1117) ·
  [FFT mask](UniChrom_Main.c#L283) ·
  [Direct mask](UniChrom_Nonlinear.c#L665)

## 2. Define output scaling and mass deconvolution

- **`rescale_newblur`**. Use the existing forward operator to calculate the
  fitted signal and derive the scale factor for reconvolved output. Decide
  whether that factor applies to the whole run or to each window. A run-wide
  factor preserves relative chromatographic intensities. Per-window factors
  follow the standalone normalization more closely. [UniDec fit scaling](udcore.c#L1061)
- **`doubledec` and `kernel`**. Reuse the 1-D mass deconvolution as a
  post-fit operation. The proposed behavior is to apply it per window,
  align the resulting mass axes, and rebuild the aggregate using the
  configured normalization. This keeps the aggregate and mass traces based
  on the same transformed spectra. [UniDec mass deconvolution](udcore.c#L2145) ·
  [UniChrom mass-grid output](UniChrom_Main.c#L495)

## 3. Extend the coupled solver

- **`variablepw` and `fwhmfile`**. Start with the direct path, which already
  stores individual response weights. Define how widths map to each window's
  processed m/z points and use them to set the operator weights and
  neighborhoods. A single FFT kernel assumes a shift-invariant m/z response,
  so supporting variable widths there would require a different operator.
  [UniDec width setup](udmain.c#L88) ·
  [Direct operator](UniChrom_Nonlinear.c#L180) ·
  [FFT kernel](UniChrom_Main.c#L226)
- **`baselineflag` and `aggressiveflag`**. Mode 2 is the simpler port because
  it estimates and subtracts a fixed baseline before iteration. Mode 1
  requires adding a baseline component to the forward prediction and updating
  it during the coupled fit. [Fixed baseline](udmain.c#L196) ·
  [Joint baseline update](udcore.c#L523)
- **Signed `mzsig`**. Coupled UniChrom uses the magnitude as the peak width and
  always applies the adjoint. Standalone UniDec uses the sign to select its
  m/z back-convolution. Matching that behavior requires defining how the
  correction still propagates across windows when the m/z back-convolution
  is disabled. [UniDec update](udcore.c#L557) ·
  [Coupled update](UniChrom_Main.c#L923)

## Verification

Check both coupled paths with each supported setting enabled and disabled.
Preserve the default results and the independent-window dispatch. Use focused
cases for overlapping assignments, zero-tolerance mass lists, charge scaling,
relative window intensities, and agreement between mass traces and aggregate
outputs under the configured normalization. For operator changes, check the
forward and adjoint together on small spectra with known widths or baselines.
[Native workflow tests](../../tests/test_unichrom_native.py)

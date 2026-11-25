# Retired UCCD-style blur in UniChrom

UniChrom previously offered UCCD-style charge and mass blurring for comparisons
in its coupled HDF5 workflow. Both coupled solvers now use UniDec blurring.
Existing HDF5 projects can still be opened. UniChrom ignores the retired
attribute when reading them. The standalone UCCD command still uses its own blur.

## Former guide text

> Coupled deconvolution uses the UniDec charge and mass blur settings (`zzsig`
> and `msig`) by default on both common and per-scan m/z axes. To compare with
> legacy UCCD blurring, pass `{"unichromuccdblur": 1}` to `set_parameters` or set
> the `unichromuccdblur` attribute on the HDF5 `/config` group. Omit it or set it
> to `0` to use UniDec blurring.

## Floor comparison

The UCCD option interpreted `zzsig` and `msig` as fractions of the observed
peak for a positive floor. It added the floor before taking the logarithm and
subtracted it after exponentiation. UniDec uses widths for the Gaussian
neighbor weights and substitutes `zerolog` for zeros and missing neighbors.
An absolute UCCD floor of `1e-12` therefore differs from UniDec's
`blur_it_mean` at `zzsig=1`. For the values `[1, 0, 1]`, the smoothed center
is approximately `1e-4` with the UCCD floor and `7.73e-3` with UniDec at
`zerolog=-12`.

The [historical performance plan](UNICHROM_PERFORMANCE_PLAN.md) records earlier
work on the UCCD-style smoother.

## Restore the comparison mode

[UNICHROM_UCCD_BLUR.patch](UNICHROM_UCCD_BLUR.patch) preserves both native blur
paths, the Python flag and HDF5 configuration, and the regression tests. To
check and apply the patch from the repository root, run:

```shell
git apply --check public/UniDec/unidec/src/Archive/UNICHROM_UCCD_BLUR.patch
git apply public/UniDec/unidec/src/Archive/UNICHROM_UCCD_BLUR.patch
```

Rebuild the native executable after applying the patch. The restored flag
defaults to `0`.

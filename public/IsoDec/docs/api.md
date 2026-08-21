# Python API

## `IsoDecRuntime`

`IsoDecRuntime(phaseres=8, verbose=False)` owns a mutable `IsoDecConfig`, a
native `IsoDecWrapper`, and the accumulated `MatchedCollection`.

- `batch_process_spectrum(data, centroided=False, refresh=False)` processes an
  `N x 2` array and returns a `MatchedCollection`.
- `process_file(path, assume_centroided=False)` processes each scan exposed by
  a built-in or UniDec reader.
- `pks_to_mass(binsize=0.1)` creates a binned zero-charge mass spectrum.
- `export_peaks(type="tsv", filename=...)` exports assignments.

## `IsoDecWrapper`

`IsoDecWrapper` is the direct ctypes interface to `isodeclib`. Its
`predict_charge`, `encode`, and `process_spectrum` methods accept NumPy arrays.

## `IsoDecConfig`

`IsoDecConfig` contains peak detection, isotope matching, charge model, and
scan metadata parameters. It is independent of UniDec and can be modified
before processing.

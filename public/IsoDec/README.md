# IsoDec

IsoDec assigns charge states and deconvolves isotopically resolved mass
spectrometry data. Its native core uses [IsoGen](https://github.com/michaelmarty/IsoGen)
for peptide and RNA isotope distributions.

## Installation

Install a published wheel:

```shell
python -m pip install isodec
```

IsoDec supports Python 3.9 and newer on 64-bit Windows, Linux, and macOS,
including Apple silicon. Native builds are also tested on ARM64 Windows and
Linux. Compatible wheels include the native IsoDec and IsoGen libraries. A
source build requires CMake 3.22.1 or newer, a C/C++ compiler, FFTW, and an
initialized IsoGen submodule:

```shell
git clone --recurse-submodules https://github.com/michaelmarty/IsoDec.git
python -m pip install ./IsoDec
```

## Python usage

```python
import numpy as np
from isodec import IsoDecRuntime

spectrum = np.loadtxt("spectrum.txt")
engine = IsoDecRuntime()
peaks = engine.batch_process_spectrum(spectrum, centroided=True)

for peak in peaks:
    print(peak.mz, peak.z, peak.monoiso, peak.matchedintensity)
```

TXT, DAT, CSV, and NPZ single-scan spectra are supported without UniDec. Raw
vendor files, mzML/mzXML, and I2MS import use UniDec's extensive importer stack:

```shell
python -m pip install "isodec[unidec]"
```

## Command line

```shell
isodec spectrum.txt --centroided -o assigned_peaks.tsv
python -m isodec --help
```

## Development

```shell
git submodule update --init --recursive
python -m pip install -e ".[test]"
python -m pytest
```

Documentation sources live in `docs/` and are built with MkDocs. Release and
wheel-building details are in [PUBLISHING.md](PUBLISHING.md).

## License and citation

IsoDec is released under the BSD 3-Clause License. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for IsoGen, FFTW, and Intel
runtime notices. When using IsoDec, please cite Pavek et al.,
“[A Fast Neural Network for Isotopic Charge State Assignment](https://doi.org/10.1021/jacs.5c03162),”
*J. Am. Chem. Soc.* **2025**, *147*, 21610–21620. Machine-readable metadata is
provided in [CITATION.cff](CITATION.cff).

## Changelog

### 2.0.1 (2026-08-18)

- Improved native processing performance with optimized AVX matrix operations,
  linear-time peak detection, faster isotope matching, and reusable work
  buffers.
- Moved installed native libraries and executables into `isodec/bin`; Windows
  binaries are included with the package.
- Added processing-time and matched-peak reporting for native spectrum runs.
- Added automated MSVC-versus-IntelLLVM performance benchmarking.

### 2.0.0 (2026-08-17)

- Split IsoDec from UniDec into the standalone `isodec` package, with IsoGen as
  a submodule and an independent version and release cycle.
- Added experimental charge-state assignment directly from profile-mode data.
- Added cross-platform tests—including native ARM64 Windows and Linux
  coverage—documentation, and automated wheel and source distribution
  publishing.

### 1.0.0

- Initial IsoDec release with IsoGen integration for isotopic charge-state
  assignment and deconvolution of mass spectrometry data, distributed as part
  of UniDec.

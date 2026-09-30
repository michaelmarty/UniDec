# CA ETD regression input

These are the user-supplied centroided `tests/ca_etd.dat` and modified-sequence FASTA
used for the native brute-force evaluation in `docs/native-brute-force-plan.md`.
The spectrum is copied without resampling or intensity normalization; it has
50,121 rows (m/z, intensity). The FASTA is `seq.fasta` in this output directory,
following IsoDecGUI's sequence-file layout. It includes the N-terminal acetylation.

Spectrum SHA-256:
`4c681c374eb7949a3d689e0ae156a5c353a2e237cd58b2550b6c24c5aa14ec5b`.

Reproduce the grouping evaluation from the IsoDec directory:

```shell
python -m benchmarks.evaluate_mass_grouping tests/ca_etd.dat tests/ca_etd_unidecfiles/seq.fasta
python -m pytest tests/test_native_mass_grouping.py
```

Fragment labels are algorithmic assignments, not independent experimental
ground truth. This one input is insufficient to select a global default order.

`reference_original.json.gz` freezes the corrected Python-only matcher and
grouping output from 2026-09-29: 794 hits and 298 original-order groups, including
hit records, membership, and every mass-group field. It is compressed UTF-8 JSON
(no pickle). The test compares native output both to this saved reference and to
a separate Python run. Regenerate deliberately after reviewing behavior changes:

```python
import gzip
import json
import runpy
from pathlib import Path
import numpy as np
from isodec.brute_force_seq_match import brute_force_pep_match
from isodec.config import IsoDecConfig

data = Path("tests/ca_etd_unidecfiles")
sequence = "".join(s.strip() for s in (data / "seq.fasta").read_text().splitlines()
                   if not s.startswith(">"))
config = IsoDecConfig()
config.activescan, config.activescanrt, config.activescanorder = 1, 0, 2
result = brute_force_pep_match(sequence, np.loadtxt("tests/ca_etd.dat"),
    config=config, native=False, centroided=True, fragmentation_type="ETD")
snapshot = runpy.run_path("tests/test_native_mass_grouping.py")["_snapshot_arrays"]
values = {key: value.tolist() for key, value in snapshot(result).items()}
(data / "reference_original.json.gz").write_bytes(gzip.compress(
    json.dumps(values, separators=(",", ":")).encode("utf-8"), mtime=0))
```

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

## `match_fragments`

```python
match_fragments(
    pks,
    sequence,
    ion_types=None,
    monoisotopic=True,
    ppm_tolerance=5,
    allow_duplicate_assignments=False,
    **isogen_kwargs,
)
```

`match_fragments` compares the neutral masses in a `MatchedCollection` with
the theoretical backbone fragments of a protein or peptide. It returns the
same collection after adding annotations and summary results.

```python
from isodec import match_fragments

match_fragments(
    pks,
    "S[Acetylation]PEPTIDE",
    fragmentation_type="ETD",
    ppm_tolerance=20,
)

print(pks.fragment_match_percent)
print(pks.sequence_coverage)
print(pks.fragment_matches)
```

### Ion types

Pass a compact string such as `"by"` or `"cz'"`, or an iterable such as
`("a", "x+1", "y-1")`, to `ion_types`.

| Ion type | Terminus | Neutral terminal shift from the residue sum |
| --- | --- | --- |
| `a` | N | -CO |
| `a+1` | N | -CO+H |
| `b` | N | none |
| `c` | N | +NH3 |
| `x` | C | +CO2 |
| `x+1` | C | +CO2+H |
| `y` | C | +H2O |
| `y-1` | C | +H2O-H |
| `z` | C | +H2O-NH3 |
| `z'` | C | +H2O-NH2, one neutral hydrogen above `z` |

`z+1`, `z•`, `z·`, and `z.` are accepted aliases for `z'`. When
`ion_types` is omitted, the default is `b` and `y` unless a named
fragmentation method is selected.

### Fragmentation methods

Pass `fragmentation_type` through the IsoGen options to select conventional
ion series:

| `fragmentation_type` | Ion series |
| --- | --- |
| `CID`, `HCD`, `SID`, `IRMPD` | `b`, `y` |
| `ETD`, `ECD` | `c`, `z'` |
| `EThcD`, `BYCZ*` | `b`, `y`, `c`, `z'` |
| `UVPD` | `a`, `b`, `c`, `x`, `y`, `z'` |
| `UVPD4` | `a`, `a+1`, `x+1`, `y-1` |
| `UVPD6` | `a`, `a+1`, `x+1`, `x`, `y-1`, `z'` |
| `UVPD9` | `a`, `a+1`, `b`, `c`, `x`, `x+1`, `y`, `y-1`, `z'` |

Explicit `ion_types` override `fragmentation_type`.

### Options

| Argument | Meaning |
| --- | --- |
| `pks` | A peak collection whose peaks provide `monoisos` or `monoiso` neutral masses. |
| `sequence` | An amino acid sequence, optionally containing supported ProForma annotations. |
| `monoisotopic` | Use monoisotopic fragment masses when true and average masses when false. |
| `ppm_tolerance` | Maximum absolute mass error in ppm. Must be nonnegative. |
| `allow_duplicate_assignments` | Store every valid label on each peak instead of only its smallest-error label. |
| `fragmentation_type` | Select a named ion-series set when `ion_types` is omitted. |
| `ambiguous_rule` | `"reject"` omits ambiguous fragment masses; `"both"` emits numbered possibilities such as `b2#1` and `b2#2`. |

Localized, terminal, global fixed, region-localized, labile, and unlocalized
ProForma modifications are handled by IsoGen 1.1.1. A modification contributes
to a fragment only when its site is retained. Ambiguous modifications and the
ambiguous residues `B` and `Z` follow `ambiguous_rule`.

By default, each peak receives its valid annotation with the smallest absolute
ppm error as `peak.sequence_match`, or `None` when it has no match. The same
theoretical ion may annotate peaks from multiple scans. With
`allow_duplicate_assignments=True`, `sequence_match` is a list of labels.

The collection receives three summary attributes:

| Attribute | Value |
| --- | --- |
| `fragment_match_percent` | Percentage of peaks with at least one fragment annotation. |
| `sequence_coverage` | Fraction of the sequence's backbone cleavage sites covered by at least one annotation. |
| `fragment_matches` | DataFrame indexed by cleavage site, with theoretical `*_mass`, observed `*_match`, and `match_count` columns. |

N-terminal ions are placed at their fragment length in `fragment_matches`.
C-terminal ions are placed at the corresponding N-terminal cleavage position,
so complementary ions contribute to the same coverage site.

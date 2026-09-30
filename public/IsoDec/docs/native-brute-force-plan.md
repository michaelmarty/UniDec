# Native brute-force fragment matching plan

## Current closeout assessment (2026-09-29)

The native fragment matcher and initial-batch C mass merge are implemented and
integrated into the runtime and IsoDecGUI. No further major algorithm migration
is required. Python appends later scans by design; native scan appending would
be a separate feature. Original insertion order
remains the default. The sections below preserve the design and experiment
history; their earlier counts, proposed ports, and release status are superseded
by this assessment and the dated follow-ups.

Current acceptance evidence:

- Sequence-composition fragment batches include supported ion termini and
  known-formula modifications; unknown-composition modifications are rejected.
- Corrected ETD results agree with the frozen Python reference: 794 hits,
  298 original-order groups, or 297 matched-intensity-order groups. The normal
  fixture has 277 hits and 197 original-order groups.
- Membership, mass distributions, intensities, scans, exports, and multi-scan
  accumulation have native/Python parity coverage. C ownership and error cases
  have CTest and sanitizer coverage for the mass-group ABI.
- This review reran all 80 IsoDec tests on Windows/Python 3.12: 80 passed,
  no skips, in 25.96 seconds. This is a source-checkout test with the local
  batch-capable IsoGen and DLL, not a fresh installed-wheel acceptance run.
  A separate probe confirmed the IsoGen batch API and all three native
  fragment-match, grouping, and centroid-rematch entry points are available.
- Linux x86_64 Release build, CTest, and an isolated installed wheel previously
  passed all 80 tests with the pinned IsoGen wheel. That validation used system
  FFTW under WSL; it does not prove repaired-wheel portability or ARM64 support.
- Native grouping avoids the Python JIT merge routines, but the corrected warm
  ETD measurements do not establish a speed gain. Historical GUI timings are
  not current performance claims after the subsequent GUI changes.

### Remaining major release work

1. **Release and require the batch-capable PyIsoGen dependency.** IsoDec's
   `pyisogen>=1.1.2` requirement and the pinned batch-capable IsoGen source both
   identify version 1.1.2. The published 1.1.2 used during Linux validation lacks
   `calc_pep_fragment_isodists`; this review also downloaded the current
   published Windows x64 wheel and confirmed that API is absent. The current
   PyPI release still identifies as 1.1.2. Publish the required API under a
   distinguishable version, update the minimum dependency and submodule pointer
   deliberately, then prove a clean
   pip install works without an editable source checkout or replacement wheel.
   See [PyIsoGen on PyPI](https://pypi.org/project/pyisogen/).
2. **Finish installed-wheel acceptance on the supported platform matrix.**
   Validate Windows x64/ARM64, macOS Intel/Apple silicon, and repaired Linux
   x86_64/ARM64 wheels using the released dependency. Include native loading,
   frozen ETD output, normal processing, multi-scan/export parity, and CTest.
   Require the fragment-match and grouping entry points in release tests:
   `native=True` currently permits a Python fallback, and grouping tests skip
   when their ABI is absent, so a green parity test alone is insufficient.
   The macOS publish jobs currently repair and upload wheels without installing
   and testing them. The IsoDec workflows live under `public/IsoDec/.github`
   in this combined repository; ensure they actually run for the standalone
   package revision, since root workflows currently cover UniDec.

### Bounded robustness work before release

The six brute-force tests cover successful matching, charge limits, malformed
public inputs, modified-fragment parity, nearby centroids, and coverage summaries.
They do not yet directly exercise the fragment C ABI's full edge contract.
Add focused native/Python tests for empty batches, zero or missing isotopes,
threshold equality, spectrum-boundary windows, invalid envelope buffers,
older-library fallback, and native error propagation. Keep the fallback tests
separate from release tests that require the new native symbols.

Audit the fragment binding's dimensions/counts before passing pointers, and
guard `fragment_match.c`'s signed `capacity * 2` growth against overflow. The
grouping ABI's stronger validation does not automatically cover this separate
entry point. These are targeted hardening and acceptance tasks, not reasons
to redesign the matcher or mass merge.

### Deferred, optional research

Further C scan appending, a new grouping default, fragment-identity grouping,
more chemistry support, and additional speed optimization are separate work.
More independently labeled spectra would help evaluate identification accuracy
and insertion-order choices; the current ETD labels are algorithmic assignments,
not experimental ground truth. Reprofile the current full GUI only if another
performance change or a new speed claim is proposed.

## Original goal and design

Move the repeated isotope-envelope search for top-down fragments into IsoDec's
C library and predict each fragment's isotope envelope from its sequence, not
from its monoisotopic mass. Give IsoGen the protein sequence once and have it
return an ordered batch of fragment labels, neutral monoisotopic masses, and
isotope envelopes. Pass that batch and IsoDec's prepared centroid spectrum to
the C matcher once; return matches identified by fragment index and charge.
Keep `MatchedPeak` construction and sequence coverage reporting in Python.

At the start of this work, the Python reference path in
`isodec/brute_force_seq_match.py` used `isogen.calc_pep_fragments()` for
fragment masses and labels, `calc_isotope_dist_dual()` for mass-based peptide
isotope envelopes, and IsoDec's `find_matches()` and score/area rules. The
implemented Python reference now consumes the same sequence-composition batch
as the native matcher. IsoDec's
native library already links to IsoGen and calls `fft_pep_mass_to_dist()`;
the pinned IsoGen revision also exports `fft_pep_seq_to_dist()` for sequence
input. Sequence mode intentionally changes predicted intensities and may
change which fragments pass the matching thresholds. The current mass mode is
therefore a comparison baseline, not an exact match-set oracle.

Before implementation, record cold and warm times for the Python matcher and
the full GUI action on the supplied centroided input and sequence. Record the
set of `(fragment label, charge)` matches, score, projected intensities, and
sequence coverage. Use that saved result as a comparison and performance
baseline. In a recent local run, the matcher took about 5.9 s cold and 0.82 s
warm, while subsequent fragment annotation took about 0.49 s; repeat these
measurements on the same revision and input before setting a speed target.

## Batching decision

Use one public IsoGen batch call per protein/fragmentation setting. Its Python
entry point should parse the sequence and ProForma annotations once, derive
fragment compositions with cumulative N- and C-terminal residue counts, and
invoke a native batch isotope calculator on packed elemental compositions.
This keeps the existing Python ProForma parser as the source of truth while
avoiding a separate native prediction call for every fragment. IsoGen should
return a contiguous intensity matrix plus aligned labels and neutral masses;
IsoDec can construct each neutral mass axis and convert it for each charge.

The previous proposed IsoDec C loop already crossed Python-to-C only once;
its repeated IsoGen calls would have stayed inside native code. Therefore a
batch IsoGen API is primarily a way to reuse sequence parsing/compositions,
centralize fragment chemistry, and amortize IsoGen setup. Passing the result
back through Python adds one bulk transfer before matching. Benchmark this
two-call path against a direct native IsoDec-to-IsoGen batch call using the
same packed compositions before claiming a speedup. Keep the public IsoGen
batch result useful independently of IsoDec.

## IsoGen changes

1. **Verify sequence-mode FFT.** At the pinned IsoGen submodule revision,
   compare native `fft_pep_seq_to_dist()` with Python
   `isogen.isodist(fragment_sequence, type="PEPTIDE", method="FFT",
   dist_only=True)`. Establish accepted residues, isotope ordering and
   normalization, offsets, supported lengths, invalid-input behavior, and
   thread safety. Use sequence mode for the intensity vector; retain the
   fragment's monoisotopic mass for its neutral mass-axis origin.
2. **Add a batch fragment API.** Expose a Python function accepting one protein
   sequence and the existing fragmentation/ion-type options. Return stable,
   aligned labels, monoisotopic neutral masses, and an `(n_fragments, isolen)`
   contiguous intensity array. Parse once and build prefix/suffix elemental
   compositions; pass those compositions to one native IsoGen batch call.
   Specify isotope-length, threshold/normalization, array ownership, ordering,
   empty-input, and error behavior. Retain existing single-fragment APIs.
3. **Make fragment chemistry explicit.** The current IsoGen documentation says
   `ion_type` adjusts the mass axis while sequence-mode intensities retain
   standard terminal composition; ProForma modifications also change the mass
   axis but are stripped before intensities are calculated. For an exact
   sequence-composition envelope, include a/b/c/x/y/z terminal composition and
   supported modification formulas in each packed fragment composition before
   the native batch calculation. A mass-only modification cannot determine
   its isotope composition; define whether it is rejected, approximated with
   a visible warning, or handled by a documented fallback. Do not silently
   label such results as exact sequence predictions.
4. **Test sequence and chemistry parity in IsoGen.** Compare batched output with
   individual sequence-mode output where their chemistry is equivalent. Compare
   native and Python sequence mode for short and long fragments, each supported ion series,
   residue composition differences at similar mass, terminal modifications,
   known-formula modifications, and unsupported mass-only changes. Check both
   monoisotopic origin and isotope intensities against independently assembled
   elemental formulas where possible. Test batch ordering, duplicate masses,
   large proteins, and output lifetime. Publish the IsoGen batch API first,
   then update IsoDec's pinned submodule revision and build docs.

**IsoGen deliverable:** one tested protein-to-fragment-envelope batch API,
including exact supported ion/composition chemistry and an explicit policy for
modifications lacking elemental formulas.

## IsoDec changes

1. **Define the batch ABI.** Add a C entry point in `isodeclib.h` and
   `isodeclib.c` accepting sorted centroid m/z and intensity arrays, IsoGen's
   aligned neutral monoisotopic masses and contiguous isotope-intensity matrix,
   spectrum bounds, maximum charge, and the existing match parameters
   (`adductmass`, isotope threshold, ppm tolerance, minimum matched peaks,
   cosine threshold, area threshold, and `minusoneaszero`). Use `double` for
   masses and m/z. Return a status, required/output hit count, and a caller
   sized array of records containing fragment index, charge, score, isotope
   indexes, centroid indexes, and the values needed to reconstruct projected
   intensities. Define ownership, capacity/overflow behavior, and an ABI
   version or feature check. Keep Python strings out of the C ABI.
2. **Port the current matching rules.** Consume IsoGen's batch envelopes in
   fragment order, apply the same strict isotope threshold, build each neutral
   mass axis from the fragment's monoisotopic origin and validated isotope
   spacing, and convert it to m/z
   for every charge whose envelope intersects the observed range.
   Preserve `find_matches()` peak selection and tolerance semantics, distinct
   matched-peak count, cosine score, intensity scaling, area coverage, and
   top-intensity exception. Reuse scratch buffers within a call. Avoid a
   second centroiding step; the GUI/runtime already supplies prepared
   centroids.
3. **Bind and integrate.** Add the new structure and function signature in
   `isodec/c_interface.py`. Call IsoGen's protein-to-fragment batch API once
   with the GUI's fragmentation settings. Pass its aligned masses and
   intensities to IsoDec C, then reconstruct `MatchedPeak` and
   `MatchedCollection` with the original labels, matched centroids, isotope and
   mass distributions, measured intensity scale, scan metadata, and mass-list
   entries. Route `IsoDecRuntime.brute_force_pep_match()` and the existing GUI
   button through this path. Keep a Python matcher fallback when an older
   native library lacks the entry point; it must consume the same IsoGen batch
   rather than revert silently to mass mode. Surface native failures and
   unsupported modification chemistry explicitly.
4. **Prove parity before optimizing further.** Add focused native and pytest
   cases for empty inputs, unsorted or invalid inputs at the Python boundary,
   duplicate/nearby centroid peaks, missing isotopes, threshold equality,
   zero-intensity peaks, charge limits, modified sequences, and envelope edges
   at the spectrum bounds. Compare native results with a Python sequence-mode
   reference for match identities, scores, observed and projected intensities,
   mass entries, and sequence coverage on synthetic cases and the supplied
   real spectrum. Compare mass mode separately to explain changes in coverage
   and score; require correct sequence chemistry and peak alignment rather
   than identical acceptance decisions across the two isotope models.
5. **Benchmark and package.** Measure IsoGen batch generation, transfer to
   IsoDec C, matching, and the full GUI action separately on the same input,
   centroid array, and build type. Compare cold and warm times with the
   per-fragment native-call prototype. Record candidate count, hit count, and
   coverage alongside time. Verify the CMake build and wheel include the
   native entry point and its IsoGen
   dependency on supported platforms. Do not broaden the C port to GUI
   annotation unless the measured end-to-end result shows that it is the next
   bottleneck.

**IsoDec deliverable:** a compatible batch matcher with parity tests, a Python
binding and fallback, runtime/GUI integration, and a before/after benchmark.

## Implementation order and acceptance

1. Save mass-mode outputs and timings; settle the IsoGen sequence and fragment
   chemistry contract, then implement and test its protein-to-fragment batch
   API.
2. Build a Python matcher using IsoGen's batch as the sequence-mode reference;
   implement and test the IsoDec C batch matcher independently of the GUI.
3. Bind it in Python and compare complete `MatchedPeak` output with the Python
   sequence-mode reference.
4. Switch the runtime and GUI dispatch, then verify coverage, projected
   intensity, timing printout, and shutdown behavior in the GUI.
5. Benchmark the release build and document the measured gain and remaining
   time spent in annotation.

Accept the migration when native and Python sequence-mode results agree on
the supplied data, projected intensities use the observed centroid scale,
fragment mass axes and chemistry are correct, focused edge cases pass, and
measured end-to-end time improves. Report coverage under both isotope models;
do not treat a difference from mass mode alone as a regression. Keep the Python
sequence-mode implementation available as a regression oracle until those
checks pass across supported builds.

## Local implementation result

The IsoGen source change is recorded in submodule commit `bb848ea` and the
IsoDec source/binaries in this working tree use that revision. The supplied
50,121-centroid ETD spectrum produced 780 hits in both the native and Python
sequence matchers. Their fragment/charge pairs and projected intensities were
identical; the largest score difference was below `6e-16`. Sequence coverage
was 81.8%. On this Windows x64 build, IsoGen generated 516 fragment envelopes
in about 0.20 s, native matching took about 0.02 s, and a warm complete call
including `MatchedPeak`/mass-list construction took about 0.58 s. First calls
were slower because existing Numba-based mass-list helpers compile on demand.

The IsoGen submodule commit must be published before the new pointer can be
used from a fresh clone. The local IsoGen checkout is installed in editable
mode into the Python 3.12 interpreter used by the GUI; no `PYTHONPATH` override
is needed locally. Other installed GUI environments still need a released
IsoGen package with this API. The native IsoDec binary has been rebuilt from
source in this working tree. Both commits remain local pending publication.

## GUI acceptance profile (2026-09-27)

On Windows x64 with Python 3.12, loaded the supplied `ca_etd.dat` (50,121
prepared centroids) into the GUI's `data2` field and used the `seq.fasta` in
its associated UniDec files directory, ETD fragmentation, and 5 ppm tolerance.
Invoked the actual `on_brute_force_match()` handler twice in one GUI process,
with the window shown. This measures the GUI matching and display path
without altering the supplied data files; opening and data preparation were
not part of these timings.

| Run | Console matching timer | Whole GUI handler | IsoGen batch | Native match | Sequence annotation | Spectrum plot | Coverage plot |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| First | 12.46 s | 14.90 s | 0.261 s | 0.020 s | 0.747 s | 1.373 s | 0.872 s |
| Second | 3.30 s | 7.52 s | 0.606 s | 0.052 s | 1.400 s | 2.060 s | 1.688 s |

The console timer ends after sequence annotation, before the plots and peak
display are updated. Both runs displayed **780 peaks and 81.8% sequence
coverage**. Every matched centroid intensity was identical to the value in
`ca_etd.dat` at the same m/z. For every matched peak, the projected envelope's
maximum intensity equaled its largest matched observed intensity exactly.
The coverage plot was populated, and closing the window through `on_exit()`
ended the subprocess with exit code 0.

A separate run on the same input compared native and Python sequence-mode
matchers: both returned the same 780 ordered `(fragment label, charge)` pairs
and 81.7829% coverage, with zero projected-intensity difference and a maximum
score difference of `5.6e-16`. Its native call included first-use Numba
compilation and its Python call ran warm, so those two call times are not a
fair speed comparison. The GUI stage timings above are single observations;
plotting and startup costs varied between the two runs.

## GUI speed pass (2026-09-27)

The brute-force result now retains IsoGen's aligned fragment labels and masses.
The GUI builds its sequence-coverage table from the accepted labels instead of
calling the generic mass-based `match_fragments()` pass again. It also preserves
the spectrum plot when it already displays the same prepared `data2` array
without extra peak markers, labels, or other plot artists;
the mass and coverage plots replace their own contents without clearing every
plot first.

On the same supplied ETD input, a fresh GUI process took 17.63 s for the first
whole handler call and 1.72 s for the second. On the second call, IsoGen batch
generation took 0.233 s, native matching 0.022 s, coverage summarization
0.250 s, mass plotting 0.207 s, and coverage plotting 0.632 s. The spectrum
plot was reused. The previous two-call profile measured 7.52 s for its second
whole handler call, including 1.40 s in generic sequence annotation and
2.06 s in spectrum plotting. These are individual observations, not a
controlled multi-run benchmark.

A follow-up real-file run after the overlay guard measured 1.96 s for its
second whole handler call, with no spectrum redraw. It again returned 780
hits, 81.8% coverage, and a clean GUI shutdown.

Both new GUI calls still returned 780 hits and 81.8% sequence coverage. The
matched centroid intensities remained identical to the input, and every
projected envelope retained the observed intensity scale. GUI shutdown exited
with code 0. Focused IsoDec and UniDec GUI workflow tests passed.

On a fresh interpreter, the real spectrum started with no compiled Numba
signatures for `fastnearest`, `fastwithin_abstol_withnearest`,
`within_ppm_plus_mm`, `merge_massdist`, `merge_decon_centroids`,
`isodist_match`, or `calc_css_from_data`. After constructing the 780 matched
peaks and 294 mass-list entries, each had compiled at least one signature.
These are examples of first-use compilation in the mass-list path; the
signature check does not attribute a separate elapsed time to each function.

## Mass-list startup profile (2026-09-27)

Instrumented `MatchedCollection.add_pk_to_masses()` around each of the 780
insertions on the same supplied centroid spectrum. Before caching, mass-list
insertion took 9.11 s in a fresh interpreter and 0.11 s on a second call in
that interpreter. The full matcher calls took 9.48 s and 0.42 s respectively.
The initial first-use delay is concentrated in mass-list construction.

Enabled Numba disk caching for the seven compiled helpers listed above. The
first call after this code change still compiled and populated the cache:
10.20 s in mass-list insertion. In a new interpreter after cache population,
the first call took 1.37 s in mass-list insertion and 1.73 s for the full
matcher; its second call took 0.10 s and 0.44 s respectively. All runs
returned 780 peaks and 294 mass-list entries. Focused brute-force and
fragment-matching tests passed (12 tests). These are single-run timings on
one Windows x64 machine; a different Python/Numba version or argument
signature must populate its own cache.

The GUI sets its 5 ppm tolerance as a `float`, while the direct matcher uses
the default integer `5`. Numba created distinct signatures for
`within_ppm_plus_mm`, `merge_massdist`, and `isodist_match` for these two
types. Accordingly, the first GUI call after only the integer-tolerance cache
had been populated still took 15.46 s; its warm call took 1.97 s. After that
GUI call populated the float-tolerance cache, a new GUI process took 5.33 s
for its first whole handler call, with 1.62 s spent in mass-list insertion.
The remaining time includes GUI display work. All these GUI calls returned
780 peaks and 81.8% coverage. Cache gains therefore depend on the argument
types used by the actual GUI path.

This removes most of the measured startup penalty without a C mass-list port.
If first-use time still matters after deployment testing, profile the
remaining roughly 1.4 s of cached first-call mass-list work before designing
a single batched C mass-list builder.

## C mass-list port proposal (historical; completed)

**Scope and payoff.** This would be an IsoDec change; IsoGen's fragment batch
API needs no change. Target the one-spectrum brute-force path first and keep
`MatchedCollection.add_pk_to_masses()` for other workflows. With the Numba
cache populated, the measured mass-list cost is about 1.4-1.6 s on a first
call in a new process and 0.1 s on a warm call. Benchmark again on target
machines before starting: this is the maximum direct saving available from
the port, and the GUI still spends time on plot updates.

**Existing C code audit.** `fragment_match.c` already has binary-bound searches,
isotope-envelope scoring, and the `FragmentHit` batch output. `isodeclib.c`
has `nearfast()` for finding a nearby spectrum m/z and native peak/shift
matching. These are useful starting points for search and score operations,
but they do not group accepted peaks into `MatchedMass` objects. In both
paths, Python currently calls `pks.add_pk_to_masses()` for every native hit:
`brute_force_seq_match.py` does so in `_add_fragment_peak()`, and
`c_interface.py` does so after `process_spectrum()`. The existing `nearfast()`
accepts a `float` query and searches spectrum m/z, so a mass-list search
still needs a `double`-precision implementation checked against Python's
nearest/tie behavior. The C port can reuse the batch-hit format and search
patterns; grouping and merged-mass state are new native work.

1. **Freeze the Python contract.** Save, in insertion order, the 780 hit
   records and all fields of the resulting 294 `MatchedMass` objects for the
   supplied spectrum. Include `monoiso`, `monoisos`, `massdist`, merged
   deconvolved centroids, `mzs`, `zs`, isotope distributions, intensities,
   scan fields, and the `clusters` membership. Add small fixtures that force
   duplicate charges, close masses, isotope shifts, unequal intensities, and
   ambiguous group choices. Use these as the reference for the C path.
2. **Add a versioned native batch result.** Extend the current C fragment
   matcher with a new entry point that groups its hits before returning them;
   keep `match_fragment_batch()` unchanged for older callers. Supply the
   spectrum, fragment envelopes, matching and grouping parameters, and scan
   metadata once. Return the existing hit records plus a group ID for each hit,
   final mass scalars, and flattened variable-length arrays with offsets for
   merged centroids, mass distributions, and monoisotopic candidates. Define
   counts, an ABI version/structure size, precision (`double` for masses and
   intensities), ownership, a matching free function, and error/overflow
   behavior in `isodeclib.h`.
3. **Port the order-dependent merge rules.** In `isodeclib.c`, process hits in
   the same order as Python. Preserve nearest-mass candidate selection,
   isotope-shift ppm and cosine checks, tie handling, first-seen charge-state
   distributions, weighted monoisotopic updates, centroid merging, mass-axis
   fitting, and apex/total/scan intensity accounting. Reuse temporary buffers
   and preserve Python's group insertion and update order, using the corrected
   arithmetic validated below. Add a C test for empty input and bounds, then fixtures
   for each merge branch.
4. **Reconstruct Python objects once.** Bind the new ABI in `c_interface.py`.
   Continue making `MatchedPeak` objects from the hit records, then create
   `MatchedMass` objects from the native summaries and hit-to-group IDs. Set
   `clusters` to the actual Python peak objects; assemble simple per-group
   charge, scan, and intensity arrays from those members without running the
   Python merge rules. Preserve the public `pks.masses`/`pks.monoisos` layout
   used by plotting and exports. Do not call `add_pk_to_masses()` for each hit
   on this path. Feature-detect the new entry
   point and use the current Python path with older native libraries.
5. **Prove parity and measure the GUI.** Compare group membership and every
   mass-list field against the Python reference, with floating values checked
   at a defined tolerance. Check the mass plot and TSV/MSAlign exports as well
   as 780 hits, 294 masses, 81.8% coverage, projected intensities, and GUI
   shutdown on the supplied data. Build and test Windows, Linux, and macOS
   libraries and wheels. Time first-ever, cache-ready fresh-process, and warm
   calls separately, including the whole GUI action. Adopt the C path only if
   it improves the measured user-visible time without changing mass-list
   behavior.

The hardest part is step 3: mass grouping changes as each peak is inserted,
so a simple C grouping pass over fixed monoisotopic masses would not reproduce
the current output. Porting only `fastnearest()` would leave Python's
per-peak merge loop and most of the interface cost in place.

### How the current mass grouping works

`brute_force_seq_match.py` turns each accepted fragment/charge hit into a
`MatchedPeak` and immediately calls `MatchedCollection.add_pk_to_masses()`.
That method searches the collection's ordered `monoisos` lookup array for the
nearest mass, then tests every group within `maxshift * 1.25` Da. Each
candidate must pass `MatchedMass.check_if_match()`: a coarse mass bound,
scan and optional retention-time bounds, an isotope-shift-aware ppm test,
and isotope-envelope cosine similarity against the group's merged
deconvolved centroids. With no passing group, it inserts a new group at the
nearest-mass position. With one, it merges there. With several, it chooses
the group with the closest apex retention time; the first candidate wins
an equal-distance tie.

On a merge, `MatchedMass.merge_in_pk()` adds the hit to `clusters`, updates
apex and scan/intensity totals, and retains the first isotope distribution
seen for each charge. It combines the new peak's neutral-mass centroids with
the group's centroids, merging close points with intensity-weighted masses.
It also updates monoisotopic candidates and the representative `monoiso`:
near-equal representatives are intensity averaged; otherwise the theoretical
distribution with the better cosine fit to the merged centroids supplies the
representative mass. Finally, it fits the chosen mass distribution's axis
and intensity to those centroids. The next hit therefore sees changed group
mass, centroids, apex, and fit; grouping is order dependent.

Two details belong in the parity fixtures. The collection's
`monoisos` lookup entry is set when a group is inserted and is not refreshed
when that group's representative `monoiso` changes. Candidate lookup can
therefore use the original insertion mass while `check_if_match()` uses the
updated group mass. Before the 2026-09-28 fix, `MatchedMass.monoisos` evaluated
`old * total + incoming * intensity / (total + intensity)`, whereas the
representative `monoiso` update divides the entire weighted sum. Both now
divide the entire weighted sum; use this corrected behavior for C parity.

### Batch grouping experiment on the supplied ETD spectrum (2026-09-27)

Loaded `ca_etd.dat` as the GUI's prepared centroid spectrum (50,121 points),
used the associated `seq.fasta`, and collected the same 780 native fragment
hits before mass insertion. Replayed those hit objects through the existing
Python grouping method in native hit order and in stable monoisotopic-mass
order, using independent copies of the hits for each replay. Both orders
produced **294 groups with exactly the same membership**. All compared fields,
including `monoiso`, `monoisos`, merged centroids, fitted `massdist`, charge
arrays, and intensity/scan fields, agreed within `1e-10` relative and
absolute tolerance on this spectrum. An initial comparison falsely reported
215 different `monoisos` arrays because `MatchedMass` shares the first peak's
candidate list by reference; reusing the same peak objects across replays
contaminated the second result. The independent-copy result above supersedes
that comparison.

Grouping by exact fragment mass alone would produce 293 groups: 292 of the
294 current groups match exact-mass groups, but four hits at mass
`4221.229560615` (`z'36`, charges 3-6) are currently split into groups of
three and one. Exact-mass grouping would combine them without applying the
isotope-envelope check.

The original order is **not intensity based**. IsoGen emits ETD fragments by
ion series (`c` then `z'`) and increasing fragment length. For each fragment,
the C matcher appends accepted charges in ascending order. Python retains
this hit order when inserting masses. The first eight hits in this run had
matched-intensity ranks 382, 220, 45, 290, 77, 204, 724, and 265 among
all 780 hits.

To test a strong-peak anchor, replayed the hits in descending
`matchedintensity` order (the sum of matched observed isotope intensities),
and separately in descending `peakint` order (the strongest matched centroid).
Each ordering produced **293 groups**: the same four `z'36` hits at charges
3-6 became one group rather than the current groups of three and one. The
charge-3 hit has matched intensity 0.062597 and cosine score 0.7470, whereas
charges 4 and 5 have matched intensities 0.276047 and 0.253138 with scores
0.9969 and 0.9790. Starting with charge 4 changes the evolving group enough
to accept charge 3. The 780 individual hits and their sequence annotations
remain the same, so sequence coverage does not change.

Sequence coverage counts cleavage rows with at least one assigned fragment
label, using `summarize_assigned_fragments()` on the individual `MatchedPeak`
objects. The supplied run has 211 covered rows out of 258, or **81.7829%**.
Reordering or regrouping the same 780 accepted hits does not alter any
label, so native hit order, mass order, descending matched intensity,
descending apex intensity, and exact-mass grouping all give **81.7829%**.
The GUI rounds each to 81.8%. Only a change that accepts, rejects, or
relabels individual fragment hits would change this coverage measure.

Among the 292 groups with unchanged membership, descending matched intensity
changed `monoisos` and the first-seen charge arrays in 153 groups, fitted
`massdist` and merged centroids in 106, and scan-intensity sums in 31 at the
stated tolerance. Descending `peakint` changed those fields in 146, 104, and
28 groups respectively. Representative `monoiso`, total intensity, and apex
fields remained equal in the shared groups. Thus a strong-peak anchor is a
plausible alternative grouping policy, but its group summaries and one group
assignment differ from current output; it should be evaluated against more
spectra and export expectations before becoming the default.

**Practical option.** Return the complete hit batch and bucket it by mass to
narrow candidate searches. Replaying hits in original order within those
buckets preserves current behavior; stable mass order also matched every
checked output field on this dataset, though other spectra may differ. Pure
exact-mass grouping skips the isotope compatibility test and merges one
additional group here. Intensity-first grouping can be an explicit mode if
better anchors prove useful on a broader validation set. Sorting does not
itself remove the per-hit merge work; a native batch builder is still needed
for the performance goal.

### Accuracy before the arithmetic fixes (2026-09-28)

Evaluated all five policies on the same 780 accepted hits, with independent
deep copies for every replay. `benchmarks/evaluate_mass_grouping.py` reproduces
the evaluation from a prepared centroid spectrum and FASTA file:

```shell
python -m benchmarks.evaluate_mass_grouping tests/ca_etd.dat tests/ca_etd_unidecfiles/seq.fasta --fragmentation ETD
```

Each group's reference is the theoretical c or z' mass associated with its
member labels. No policy combined different theoretical fragment masses on
this dataset. All representative `monoiso` values agree with their reference
within numerical rounding (maximum error below `5e-10` ppm). This is a
consistency check: brute-force matching initializes `monoiso` from the
theoretical mass, so that field cannot independently establish accuracy.

The fitted `massdist` axis can move during merging. Recovering its implied
monoisotopic origin by subtracting the retained first isotope's mass offset
gives the following errors relative to the source fragment mass. The exact
mass policy forces groups by identical theoretical mass, then uses the
existing merge code in original hit order to construct their summaries.

| Policy | Groups | Fitted origins within 5 ppm | Median absolute error, ppm | 95th percentile, ppm | Maximum, ppm |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original order | 294 | 243 (82.7%) | 1.239 | 8.453 | 16.228 |
| Mass order | 294 | 243 (82.7%) | 1.239 | 8.453 | 16.228 |
| Descending matched intensity | 293 | 283 (96.6%) | 1.099 | 4.406 | 9.286 |
| Descending apex intensity | 293 | 280 (95.6%) | 1.137 | 4.645 | 9.286 |
| Exact mass, original merge order | 293 | 241 (82.3%) | 1.258 | 8.467 | 16.228 |

As a check against the measured data, estimated each group's monoisotopic
mass from its matched centroid m/z values: convert to neutral mass, subtract
each assigned isotope's offset, then average using observed intensity.
All estimates lie within 5 ppm. Original and mass order have median absolute
error 0.992 ppm and maximum 4.562 ppm; the other three policies have median
0.991 ppm and maximum 3.854 ppm. That small change comes from combining the
four z'36 hits. This estimate still uses the accepted isotope assignments and
therefore does not establish whether a fragment identification is correct.

The alternative `monoisos` candidate arrays fail badly under every policy:

| Policy | Candidate entries within 5 ppm | Groups retaining any valid candidate |
| --- | ---: | ---: |
| Original or mass order | 172 / 476 | 169 / 294 |
| Descending matched or apex intensity | 172 / 476 | 168 / 293 |
| Exact mass | 170 / 475 | 167 / 293 |

Two current arithmetic defects must be addressed before choosing a grouping
policy or preserving the code in C. The candidate update divides only the
second term of its weighted sum, corrupting masses even when both inputs
are identical. Also, `merge_massdist()` adds the absolute distance to a
centroid, so a centroid below the predicted mass moves the predicted axis
upward; repeated merges can accumulate positive drift. These findings revise
the earlier parity-first recommendation: establish corrected Python
arithmetic as the reference, then repeat this comparison before the C port.
Descending matched intensity is the best tested order for the current fitted
axes, but its advantage may partly compensate for the fitting defect. No
production grouping or mass arithmetic was changed by this evaluation.

### Corrected arithmetic and accuracy rerun (2026-09-28)

Fixed the candidate update to divide the entire weighted sum. `MatchedMass`
now owns a copy of the first peak's candidate array, so merging no longer
changes that peak's candidates. Fixed `merge_massdist()` to shift toward the
observed centroid using a signed displacement, with its magnitude checked
against the existing tolerance. Original hit order remains the production
default. The preceding accuracy table records the behavior before these fixes.

Reran all five policies on `ca_etd.dat` with the same sequence. Every policy
now has exactly one valid candidate per group: 294/294 for original and mass
order, 293/293 for the other policies. Representative and candidate masses
agree with their source c/z' masses within numerical rounding. Fitted-axis
accuracy is now:

| Policy | Fitted origins within 5 ppm | Median absolute error, ppm | 95th percentile, ppm | Maximum, ppm |
| --- | ---: | ---: | ---: | ---: |
| Original order | 293 / 294 | 0.983 | 3.293 | 5.168 |
| Mass order | 293 / 294 | 0.983 | 3.293 | 5.168 |
| Descending matched intensity | 292 / 293 | 0.998 | 3.294 | 5.168 |
| Descending apex intensity | 292 / 293 | 0.998 | 3.294 | 5.168 |
| Exact mass, original merge order | 292 / 293 | 0.998 | 3.294 | 5.168 |

The large apparent accuracy advantage of intensity-first insertion disappears
after correcting the fitting arithmetic. The remaining fitted-axis outlier
is z'128 (charges 10-13). Its individually matched centroids remain within
5 ppm of the theoretical neutral mass after isotope offsets are removed.
The group fitter also uses local unmatched centroids and checks displacement
against the evolving fitted axis; it does not enforce a final 5 ppm bound
relative to the original theoretical monoisotopic mass. The observed-centroid
mass estimates remain within 5 ppm for every group under every policy.

The normal brute-force run still returns 780 hits and 211/258 covered
cleavages (81.7829%). Validation passed: 19 focused grouping, brute-force,
and fragment tests; three spectrum/output tests; and two GUI workflow tests.
New regression cases cover identical and unequal candidate masses with unequal
weights, candidate ownership, positive/negative/zero axis shifts, repeated
fitting without drift, and out-of-tolerance centroids. The corrected Python
implementation is the reference for the proposed native port; these results
provide no accuracy reason to switch the default insertion order.

### Grouping runtime after the fixes (2026-09-28)

Benchmarked the corrected Python grouping on the same 780 accepted hits using
`--timing-repeats 15` with `benchmarks.evaluate_mass_grouping`. Warmed each
policy first, then randomized policy order within each of 15 rounds. Each
trial used independent copies of the input hits. Timings include policy
sorting and mass-group construction; they exclude input copying, IsoGen,
native fragment matching, accuracy analysis, and GUI work. These are warm
Windows/Python 3.12 measurements, not startup or end-to-end GUI timings.

| Policy | Median, ms | Minimum, ms | Maximum, ms |
| --- | ---: | ---: | ---: |
| Original order | 48.52 | 44.57 | 52.93 |
| Mass order | 45.35 | 42.92 | 49.11 |
| Descending matched intensity | 51.62 | 50.39 | 122.45 |
| Descending apex intensity | 51.87 | 49.03 | 54.84 |
| Exact mass | 40.61 | 37.01 | 45.52 |

Exact-mass grouping is fastest in this run, saving about 7.9 ms (16%) of
grouping time, but skips isotope compatibility checks and combines the
previously split z'36 group. Mass order is fastest among policies retaining
the checks, saving about 3.2 ms (6.5%); its timing range overlaps original
order. Intensity-first ordering offers no speed advantage here. These small
warm-call savings do not justify changing the default on their own; any
native port should target measured startup costs and the complete merge loop.

The arithmetic corrections are documented in the README changelog for
IsoDec 2.0.5, with the package and citation versions updated together.

### Similarity and stability checks before the C port (2026-09-28)

`benchmarks.compare_mass_grouping` replays independent copies of accepted
hits, compares exact group membership and charge sets, and aligns projected
isotope intensities by fragment label and isotope number. It reports global
cosine similarity and normalized absolute intensity difference. Two inputs
from the supplied carbonic anhydrase experiment were tested: `ca_etd.dat`
(50,121 prepared centroids, 780 hits), and the associated profile TXT after
IsoDec centroiding (4,209 centroids, 294 hits). These are different data
preparations of the same experiment, not independent proteins or methods.

| Policy | ETD groups / changed hits | ETD projection cosine / relative L1 | Profile groups / changed hits | Profile projection cosine |
| --- | ---: | ---: | ---: | ---: |
| Original order | 294 / 0 | 1 / 0 | 146 / 0 | 1 |
| Mass order | 294 / 0 | 1 / 0 | 146 / 0 | 1 |
| Descending matched intensity | 293 / 4 | 0.999972 / 0.119% | 145 / 3 | 1 |
| Descending apex intensity | 293 / 4 | 0.999972 / 0.119% | 145 / 3 | 1 |
| Exact mass | 293 / 4 | 1 / approximately 0 | 145 / 3 | 1 |

On both inputs, the changed hits are charges of the same z'36 fragment:
charges 3-6 in `ca_etd.dat` and charges 3-5 in the profile-derived spectrum.
Mass order reproduces all original groups, isotope projections, charge order,
and group fields on both inputs. Exact mass preserves the summed projection
despite merging z'36 groups, but that merge changes group membership and
therefore charge grouping in mass-list output. The intensity policies preserve
the charge *sets* and total intensity of each shared group, but reorder their
first-seen charge arrays in 146-153 ETD groups and 55-57 profile groups.
The maximum fitted-axis difference among shared ETD groups is 1.794 ppm for
the intensity policies; it is zero for mass and exact-mass order. A global
cosine near one therefore does not establish group-level parity.

With a fixed random seed, 30 shuffled ETD insertions produced 293 groups in
26 runs and 294 in four; the only changed membership involved the same four
z'36 hits. Relative projection L1 ranged from 0.030% to 0.434%, and global
cosine from 0.999792 to 0.999998. Ten shuffled insertions of the profile
centroids produced 145 groups in six runs and 146 in four, changing the same
three z'36 hits. This is a specific order-sensitive case, not evidence that
the full grouping is independent of insertion order.

`benchmarks.synthetic_grouping_comparison` generated two known ETD fragments
(`c4` and `z'4`) of `PEPTIDEK`, each at charges 2 and 3. Clean spectra, one
missing isotope per envelope, weak random noise with 0.35 ppm mass jitter,
and nearby low-intensity shoulders each returned the four intended hits in
two groups. All five policies gave identical membership and projected
intensities on these four cases.

Two harder synthetic decoys expose the grouping rule's limit. With two known
fragment sources 1.0033 Da apart and broad five-isotope envelopes, original,
mass, and both intensity orders merged all four charge hits into one
mixed-identity group. Exact-mass grouping kept the two known identities
separate. With two distinct sources of exactly the same mass, all five
policies—including exact mass—merged them. A brute-force-specific group key
based on the IsoGen fragment index would preserve known source identity in
both cases. This does not prove that label grouping will improve every
experimental result, but it is a more meaningful alternative to evaluate
than sorting the accepted hits. Both decoys are reproducible with
`benchmarks.synthetic_grouping_comparison`.

The remaining z'128 fitted-axis outlier was traced through charges 10-13:
its axis moved 0, 4.135, 5.599, and 5.168 ppm from the theoretical origin
as each peak merged. The charge-11 window contained 70 local centroids but
only 16 accepted isotope matches; the charge-12 window had 77 and 17.
`calc_mass_dists()` currently contributes *all* local centroids to grouping.
In a diagnostic replay using only matched centroids, all 294 fitted axes
were within 5 ppm, with a maximum error of 4.425 ppm and unchanged group
count. This is a separate candidate behavior change: test membership,
intensities, and other spectra before adopting it as the C reference.

The actual GUI handler was also timed in a fresh Python 3.12 process with
the prepared ETD spectrum already loaded: 5.57 s for its first call and
0.737 s for its second. Both calls returned 780 hits and 81.7829% sequence
coverage; the window closed cleanly. These are single observations including
matching and plots, excluding file load and data preparation. An approximately
3-8 ms warm grouping-order gain is a small part of this user-visible time.

**Port decision.** Preserve corrected Python insertion order for an exact C
parity implementation. For a behavior-changing native brute-force path,
prototype grouping by IsoGen fragment index and measure it against the
current groups and exports. Test matched-centroid-only fitting as a separate
change. Mass sorting offers little measured end-to-end benefit, and
exact-mass grouping cannot distinguish isobaric source identities. Independent
proteins and fragmentation methods are still needed before promoting either
behavior change; the supplied local alternate data come from the same
experiment.

### Matched-centroid grouping and final C recommendation (2026-09-28)

Implemented the matched-centroid rule in shared Python `MatchedPeak.calc_mass_dists()`.
Brute-force peaks already carry accepted centroids and use those directly.
Normal native IsoDec peaks do not populate that Python field, so the shared
method now finds their accepted centroids in the local window against the
returned isotope distribution before converting to neutral mass. A peak
without an isotope distribution can still use its full local window. The
original `centroids` field remains available for plotting and inspection.
This is a behavior change for normal IsoDec as well as the brute-force path.

| Corrected grouping policy | ETD groups | Fitted origins within 5 ppm | Maximum fitted-origin error | Warm grouping median (15 runs) |
| --- | ---: | ---: | ---: | ---: |
| Original order | 294 | 294 / 294 | 4.425 ppm | 41.92 ms |
| Mass order | 294 | 294 / 294 | 4.425 ppm | 40.38 ms |
| Descending matched intensity | 293 | 293 / 293 | 4.425 ppm | 46.90 ms |
| Descending apex intensity | 293 | 293 / 293 | 4.425 ppm | 46.64 ms |
| Exact mass | 293 | 293 / 293 | 4.425 ppm | 35.60 ms |

All candidate masses still match their known fragment masses. On the prepared
ETD input, mass order retains all 294 original groups and has identical
projected intensities. The intensity policies still combine the four z'36
hits: projection cosine is 0.999964 for matched-intensity order and 0.999975
for apex-intensity order. Exact-mass grouping combines those hits but its
summed projected intensity matches original order. The profile-derived
input still yields 294 hits and 146 original groups; its policy membership
and projection comparisons are unchanged. Randomized insertion retains the
same 26/30 and 6/10 split/merge counts on the two inputs. The simple
synthetic cases and both close-mass decoys retain their earlier outcomes.

Normal IsoDec's supplied 3,537-centroid regression spectrum still returns
277 matched peaks. Matched-centroid grouping yields 197 mass groups instead
of 195: two old pairs split, one separated by roughly three isotope masses
and the other by roughly two. Both are plausible improvements in specificity,
though this spectrum lacks independent sequence labels to prove their
identities. In 16 alternating warm runs, complete direct processing took a
median 46.04 ms with matched-centroid merging versus 42.85 ms with the old
full-window merging. These timings include normal native matching and Python
mass grouping; the roughly 3 ms increase reflects Python rematching work.
The GUI brute-force handler took 5.35 s on a fresh-process first call and
0.791 s warm, returning 780 hits and 81.7829% coverage both times. These
are single GUI observations with the prepared spectrum already loaded.

All 44 IsoDec tests and the two focused UniDec GUI workflow tests passed.
The suite now includes a regression showing that an unrelated local
centroid cannot enter a merged mass group, plus a normal-style peak with
no prepopulated `matchedcentroids` field. The normal mass-group count
regression was updated from 195 to 197 after inspecting the two split pairs.

**Final recommendation:** Port the corrected, order-dependent merge to C as
one source-neutral grouping implementation. Feed it each accepted peak's
matched centroid masses, theoretical isotope distribution, charge,
monoisotopic candidates, intensity, scan, and retention time. The brute-force
`FragmentHit` records already include matched centroid indexes and counts;
the normal native matcher already computes indexes internally. Preserve those
indexes or an explicit matched-centroid list through the native boundary so
normal IsoDec does not repeat peak matching in Python. Invoke the same native
grouping implementation from both the brute-force and normal IsoDec paths;
support original-order parity and a selectable descending-intensity order,
and verify all group fields and exports against the corresponding Python
reference. A fragment-identity grouping option
can be evaluated later for sequence-driven workflows, but it cannot serve
as the shared default because normal IsoDec has no known fragment identity.
Mass sorting saves only about 1.5 ms of warm grouping here, and exact-mass
grouping changes grouping semantics. Target batched native grouping and
startup costs rather than a standalone nearest-mass lookup.

### Which grouping is preferable for z'36? (2026-09-28)

The four accepted z'36 hits at charges 3-6 have the *same* IsoGen
monoisotopic mass, and each passed the fragment-envelope threshold. The
charge-3 hit is weaker (matched intensity 0.062597, cosine 0.747) than
charges 4 and 5 (0.276047/0.9969 and 0.253138/0.979). The original order
seeds the group with charge 3 and leaves it separate; starting with a
stronger hit groups all four. **One group is the better interpretation for
this labeled fragment**, subject to the weaker charge-3 hit being a genuine
identification. This is a reason to consider strongest-first insertion and
not merely a cosmetic difference in output order.

The same order change does not have the same evidential support in normal
IsoDec. Replaying its 277 unlabeled native peaks after matched-centroid
selection gives 197 groups in original order, 198 in mass order, and 196 in
either descending-intensity order. Intensity-first changes only two peaks:
it combines a charge-6 candidate at 6035.0874 Da with a charge-5 candidate
at 6033.1060 Da. Their approximately 1.9814-Da separation is near a
two-isotope shift, but sequence identity is unavailable to decide whether
this merge is correct. The broad-envelope synthetic decoy also shows that
intensity-first insertion can merge two *different* known sources one
isotope mass apart. Thus the z'36 improvement alone does not justify
switching the global default for normal IsoDec.

For the shared C implementation, make insertion order an explicit policy:
original order for parity and intensity order for evaluation. Both policies
must use the same matched-centroid merge and scoring code. Use descending
**matched isotope intensity** as the proposed stronger-hit priority; compare
its group identities and exports against known-source data before making it
the default for either workflow. A separate optional fragment-identity guard
could protect sequence-driven matching from the synthetic mixed-identity
case, while the core grouping API continues to work without labels for
normal IsoDec. The C port should not hard-code original order as a statement
of correctness.

### `fastnearest()` C feasibility trial (2026-09-27)

Built a temporary `double`-precision C version of `fastnearest()` and called
it through `ctypes` with the same conversion required by the Python mass-list
path. It returned the same indexes as the Numba function for all 486 direct
lookups from the supplied ETD run, plus empty, singleton, exact-hit, and tie
cases. Of those 486 real lookups, 358 received Python lists and needed a
contiguous NumPy conversion for the C call.

On repeated warm lookups, the best of four blocks of 486 calls took 2.78 ms
with Numba and 2.71 ms with C plus `ctypes` and conversion. Alternating
full-matcher trials varied by much more than this 0.07 ms difference, so they
did not establish a user-visible speedup. The other nearest searches occur
inside Numba-compiled `fastwithin_abstol_withnearest()` and `merge_massdist()`;
substituting a `ctypes` function for the Python-level call does not replace
those compiled calls. No standalone C lookup was retained in the source tree.
The existing Numba cache remains the faster path to a material startup gain.

### Normal native match indexes are not yet a grouping input (2026-09-28)

Before using normal IsoDec's `MatchedPeak.matchedindsexp` array for the C mass
grouping port, decoded its increasing isotope-index prefix and treated the
paired centroid indexes as offsets into the exported `startindex:endindex`
window. On the supplied 3,537-centroid regression spectrum, this changed the
matched centroid set in **71 of 277 peaks**. The group count stayed at 197 and
group membership stayed the same, but fitted mass distributions differed.
The direct substitution was therefore reverted. The Python rematch of the
exported window remains the grouping reference.

For example, one charge-1 peak has theoretical m/z values near 649.2872 and
650.2905. Its stored native centroid offsets are 5 and 7, which index
649.2629 and 650.2237 in the exported window. The Python rematch selects
649.2872 and 650.2906. This shows that the native offsets cannot simply be
interpreted relative to the exported window; it does not establish whether
the mismatch comes from the original search-window base, result bookkeeping,
or another cause. A regression test now protects this peak's mass-grouping
centroids.

The shared C port needs either absolute accepted-centroid indexes with an
explicit count and documented base, or a native rematch against the same
exported window using Python `find_matches()` semantics. Compare the selected
centroid sets and all fitted group fields with the current Python path before
removing its rematch. Native acceptance indexes alone are not a parity oracle
for mass grouping.

### Shared native grouping implementation and order decision (2026-09-28)

Added a version-1 C batch grouping ABI in `mass_group.c`, used by normal
IsoDec and sequence brute-force matching. Normal peaks first go through one
native batch rematch of the exported centroid windows; fragment peaks retain
their accepted centroid sets. Both then use the same order-dependent C merge.
Python reconstructs `MatchedMass` objects and their peak memberships from
the native summaries. An older library without these entry points continues
to use Python grouping. `IsoDecConfig.mass_group_order` accepts `"original"`
and `"matched_intensity"`; original order is the default. The UniDec GUI's
separate config class also defaults to original order when this field is
absent.

With a Windows x64 release build, the 3,537-centroid normal regression input
returned the same 277 peaks and 197 original-order groups as Python. Selected
centroid sets, group membership, fitted mass distributions, merged centroids,
charge and scan arrays, tabular output, and binned mass spectrum agreed. The
normal native peak intensities are float32, and representative/candidate
masses can differ by less than 0.001 Da because the C grouping arithmetic
uses double for masses. Six native grouping tests, all 51 IsoDec tests, and
the two focused UniDec GUI workflow tests passed with the rebuilt local DLL.
The rebuilt `isodeclib.dll` is in the package's `bin` directory.

In 15 alternating warm measurements of complete normal spectrum processing
after three warmups, the Python path had a 35.08 ms median and the native path
a 32.21 ms median. The native path had one 91 ms outlier; the result suggests
only a small warm gain on this input. Linux, macOS, and wheel builds remain
to be validated.

The supplied 50,121-centroid ETD input returned the same 780 labeled hits and
294 original-order groups through both the native and Python grouping paths.
Group membership, lookup masses, representative and candidate masses, fitted
mass distributions, merged centroids, and total intensities agreed at
`rtol=1e-6, atol=1e-3`. In eight alternating warm complete matcher calls per
path, after two warmups, the median was 184.71 ms with C grouping and
206.25 ms with Python grouping. With the existing Numba disk cache, three
fresh-interpreter native calls took 211-213 ms each; three Python-grouping
calls took 633-841 ms each. These
calls exclude file loading and GUI plotting.

The real GUI handler returned 780 hits, 294 groups, and 81.7829% sequence
coverage, with a clean shutdown. One fresh GUI process took 1.20 s for its
first handler call and 0.65 s for its second. In a separate alternating
warm comparison of eight handler calls per path, the median was 690.0 ms
with native grouping and 700.7 ms with Python grouping. Plot-time variation
was larger than that 10.7 ms median difference, so the direct matcher
benchmark is stronger evidence of a speed gain than the whole-GUI timing.

Descending matched-intensity order reproduces the Python policy and yields
196 rather than 197 groups on the normal regression spectrum. It combines a
charge-6 candidate near 6035.0874 Da and a charge-5 candidate near
6033.1060 Da; there is no sequence identity for this spectrum to judge that
merge. The earlier labeled ETD case supports combining all four z'36 hits,
and the new native intensity-order run reproduced the Python policy's 293
groups without mixing theoretical masses on that spectrum. The known-source
synthetic decoy demonstrates that isotope-shift grouping can still combine
distinct fragments. Keep original order as the shared
default until independent labeled inputs establish that intensity-first
improves identification across both workflows. Intensity order remains
available explicitly for evaluation.

### Review fixes and corrected reference (2026-09-29)

Fixed the review findings in the Python reference, native implementations, and
binding rather than merely loosening parity checks:

- Grouping recovers missing accepted centroids using `calc_mass_dists`, refreshes
  cached neutral centroids, and supports missing isotope distributions.
- Coincident zero-weight centroids retain a finite mass in both implementations.
- The binding validates buffer shapes, finite values, lengths, and scan/charge
  integer ranges before exposing memory to C. Native pair/count overflow checks
  reject oversized input; valid output pointers reset on errors, and scan
  subtraction cannot overflow a signed integer.
- Native result slots clear the unused isotope-array tails, and Python reads
  only `realisolength`. The normal fixture's formerly unsorted envelope is now
  monotonic; a synthetic stale-tail regression also checks intensity sums.
- Cosine scoring no longer wraps a missing preceding isotope to the final
  observed isotope. This fixes Python scoring, C grouping, and C fragment
  matching. An actually present preceding isotope still contributes its penalty.
  Perfect one-isotope matches now score 1, not approximately 0.7071.
- Group masses, distributions, total intensities, and scan intensities now use
  double precision consistently, including the Python incremental path. Scores
  and acceptance thresholds were not retuned to preserve erroneous old counts.

The normal fixture retains 277 hits, 197 original-order groups, and 196
intensity-order groups. The repository ETD fixture now yields 794 accepted hits,
298 original-order groups, and 297 intensity-order groups through independently
executed Python and native paths. Both orders have zero groups mixing theoretical
fragment masses on this input; every fitted-axis origin remains within 5 ppm
(maximum 4.425 ppm). These labels are not independent experimental ground truth,
so original order remains the default.

`tests/ca_etd.dat` and `tests/ca_etd_unidecfiles/` contain the supplied centroid spectrum and FASTA, their
provenance, and a compressed JSON reference generated by the corrected Python-only
path. The reference freezes hit records and every group field. Expanded tests
compare separate native and Python objects, all scan/RT fields, memberships,
missing envelopes, zero intensities, malformed buffers, tolerance boundaries,
shuffled insertions, duplicate charges/scans, and multiple scans. All-group
comparisons use absolute mass tolerance 1e-9 Da and intensity relative tolerance
1e-12; integer fields and membership are exact. Actual TSV/MSAlign exports and
mass spectra are compared, in addition to the frozen original-order reference.

Appending scans intentionally remains a Python incremental operation after the
first native batch, tested against an all-Python grouping/rematching fallback.
Intensity ordering is per incoming spectrum, not a global replay of prior scans.
Both evaluation scripts now suppress the batch finalizer when collecting pristine
hits; patching `add_pk_to_masses` alone no longer suppresses native grouping.

Windows x64 release builds, all 80 IsoDec tests (with importlib collection), and
both focused UniDec GUI workflow tests pass. The 33 focused native-grouping and
integration tests also pass with deprecation warnings treated as errors. Linux GCC
standalone grouping ABI checks pass with `-Wall -Wextra -Werror` and AddressSanitizer
plus UndefinedBehaviorSanitizer, including repeated allocation/free, oversized
counts, zero weights, and extreme scan numbers. The fragment matcher also compiles
with GCC warnings treated as errors. The WSL environment lacks CMake and FFTW
development dependencies, so this is not a full Linux package validation. A macOS
host is unavailable. Windows wheel validation needs the missing packaging tools
(`build`, `scikit-build-core`, and `wheel`) in an isolated environment.

After two warmups and five alternating measurements per path, the full corrected
ETD matcher took median 217.30 ms with native grouping and 216.40 ms with Python
grouping. Input validation removes the earlier small warm advantage on this run;
these results do not substantiate a warm speedup or update the historical GUI
timings. Native grouping still avoids invoking Python's JIT merge routines.

Deprecated trapezoidal integration calls were replaced throughout package code,
scripts, and both teaching notebooks. IsoDec uses SciPy's `trapezoid` to retain
its NumPy 1.23/Python 3.9 dependency floor; UniDec and scripts use `np.trapezoid`.

### Linux native build follow-up (2026-09-29)

The WSL Ubuntu x86_64 environment now provides CMake, Ninja, GCC 11.4.0,
FFTW development files, and Python venv support. A fresh Release CMake build
completed all 14 targets with the existing default options, followed by an
install into an isolated staging prefix outside the repository. No source
changes or tracked native artifact updates were needed.

The mass-group ABI checks now live in `isodec_test.cpp` and run through CTest.
They cover ownership, repeated allocation/free, empty input, oversized counts,
and extreme scan values. The installed library has `$ORIGIN` RUNPATH and
resolves its sibling `isogen.so` and system FFTW without missing dependencies.
The installed `isodec_test`, run with the staged libraries on its loader path,
processed a copy of `tests/test_spectrum.txt` (3,537 rows) and wrote 266 matches.
This executable smoke test uses its own defaults and does not establish parity
with the 277-hit Python-configured workflow.

With pip authorized, built `isodec-2.0.5-py3-none-linux_x86_64.whl` through
`python -m build --wheel` and installed it into an isolated Python 3.10 WSL
virtual environment. `pip check` passed. Tests ran from a separate copy of
the fixtures and test files, using the wheel-installed package and native
library rather than the Windows source package binaries.

The first full run with pip's published `pyisogen==1.1.2` dependency had
71 passes and 9 failures: that release lacks `calc_pep_fragment_isodists`.
Building and installing a PyIsoGen wheel from the existing pinned submodule
(without changing its sources) resolved those failures. **All 80 IsoDec tests
passed in 19.56 seconds**, including normal-spectrum grouping, multi-scan and
export parity, and the frozen ETD reference: 794 hits, 298 original-order
groups, and 297 intensity-order groups.

The pinned PyIsoGen source currently also identifies itself as version 1.1.2,
so the declared `pyisogen>=1.1.2` dependency does not distinguish the required
batch-capable build from the published release. A normal fresh pip install
therefore still lacks fragment batch support; this remains a release dependency
issue, despite passing validation with the pinned source build. No dependency
version or publishing configuration was changed during this validation.

This validates a local Linux x86_64 wheel on WSL Ubuntu with system FFTW;
it does not establish manylinux portability or ARM64 compatibility. Windows
wheel and macOS validation remain outstanding as previously recorded.

### IsoDecGUI names and ETD fixture layout (2026-09-29)

IsoDecGUI now writes `conf.dat`, `input.dat`, `rawdata.txt`, and other plain
output names inside each `<sample>_unidecfiles` directory, following the
UniChromCD pattern. Existing `<sample>_conf.dat` files and their companion
lists still load when a short config is absent. Interactive and batch exports
use `results.tsv` and `results_ms1.msalign`/`results_ms2.msalign` inside that
directory. The sequence remains `seq.fasta` there.

The ETD example is `tests/ca_etd.dat` beside `tests/test_spectrum.txt`;
`tests/ca_etd_unidecfiles` holds `seq.fasta`, the frozen reference, and its
README. The moved spectrum retains its SHA-256. The former standalone C ABI
harness was merged into `isodec_test.cpp` as `--test-mass-group` and registered
with CTest. After these changes, the native CTest passed, all 80 IsoDec tests
passed on Linux, and the focused UniDec naming and GUI workflow tests passed
on Windows.

### C mass merge closeout (2026-09-29)

The C mass merge migration is complete for initial batches in normal IsoDec
and sequence brute-force matching. Python continues to append later scans by
design; the multi-scan parity test covers that boundary. Original insertion
order remains the default, and matched-intensity order remains an evaluation
option. The corrected ETD comparison showed no warm speed gain, so there is
no pending performance claim for this migration. Publishing the required
PyIsoGen API and completing platform wheel validation remain separate release
work described above.

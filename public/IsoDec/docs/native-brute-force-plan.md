# Native brute-force fragment matching plan

## Goal and current boundary

Move the repeated isotope-envelope search for top-down fragments into IsoDec's
C library and predict each fragment's isotope envelope from its sequence, not
from its monoisotopic mass. Give IsoGen the protein sequence once and have it
return an ordered batch of fragment labels, neutral monoisotopic masses, and
isotope envelopes. Pass that batch and IsoDec's prepared centroid spectrum to
the C matcher once; return matches identified by fragment index and charge.
Keep `MatchedPeak` construction and sequence coverage reporting in Python.

The current Python path is the reference implementation in
`isodec/brute_force_seq_match.py`. It uses `isogen.calc_pep_fragments()` for
fragment masses and labels, `calc_isotope_dist_dual()` for mass-based peptide
isotope envelopes, and IsoDec's `find_matches()` and score/area rules. IsoDec's
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
used from a fresh clone. A released IsoGen package with this API is also needed
for ordinary installed-GUI use; until then, local tests use the submodule on
`PYTHONPATH`. The native IsoDec binary has been rebuilt from source in this
working tree.

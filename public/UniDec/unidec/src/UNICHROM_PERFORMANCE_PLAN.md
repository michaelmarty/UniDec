# UniChrom performance investigation and implementation handoff

Updated 2026-09-09 after implementing and validating charge-summed 2-D FFTs.
This is the remaining-work handoff; the completed reduction is no longer an
implementation task. The measurements below come from that implementation
session. This document update runs no application code, builds or benchmarks.

Updated again after fusing the post-adjoint correction with cleanup in
`UniChrom_Main.c`. That exact cube-traffic optimization is also complete;
remaining tasks below begin with phase profiling and larger opportunities.

Updated once more after making point smoothing scan-parallel with private
per-scan rolling sums. This removes nested OpenMP regions from the UniChrom
path while preserving the existing smoothing calculation.

Updated after caching the mass-bin lower index and interpolation fraction per
m/z/charge cell. This removes repeated `calcmass`/`floorf` work across scans;
the change is exact and its isolated timing effect is currently below the
run-to-run noise floor.

A strided FFTW `plan_many` output prototype was also measured and rejected:
it raised the current median from about 0.147 s to 0.163 s without regularizers
and from 0.374 s to 0.402 s with regularizers. The existing per-charge 2-D
transforms are retained; future output work should test contiguous bounded
batches or separable direct passes instead of assuming strided batching wins.

The follow-up contiguous batch uses temporary charge-major planes and one
FFTW `plan_many` forward/inverse pair. It is retained because it was numerically
equivalent and measured 0.144 s without regularizers and 0.375 s with them in
the latest repeated workload. Allocation or plan failure falls back to the
established per-charge path; the temporary buffers are released before output
masking.

Changing the reduced FFT plans from `FFTW_ESTIMATE` to `FFTW_MEASURE` was also
tested and rejected: the representative unregularized median rose to about
0.146 s from 0.144 s, with added plan-time work. Keep `FFTW_ESTIMATE` unless a
workload with many repeated invocations demonstrates that planning amortizes.

The next execution pass added three bounded changes. An optional
`UNICHROM_PROFILE` build now reports preprocessing, merge/grid I/O, solve, and
output-FFT wall-clock phases without changing the default binary. UniChrom now
configures the available MKL thread count for the reduced transforms and can
initialize FFTW threads when that backend exposes its thread API; both settings
are restored safely on cleanup. Finally, masks below 75% occupancy use cached
per-m/z allowed-charge lists for projection, correction, and m/z output sums,
while dense masks retain the contiguous loops. The sparse list allocation is
best-effort and falls back to the dense implementation if it cannot allocate.
The isolated Release build and native 3-D equivalence test pass after these
changes. A direct/separable convolution implementation was not promoted:
without phase data and kernel-support measurements it would add a second
numerical path prematurely; the reduced FFT remains the measured default.
Enable the phase diagnostics in a configured build with
`-DUNIDEC_PROFILE_UNICHROM=ON`; the option defines `UNICHROM_PROFILE` only for
the native executable. The sparse list path still needs a dedicated fixture
with valid positive mass bounds before it can be claimed by durable regression
coverage; the existing broad-window suite exercises the dense fallback.

The profiled executable was rebuilt and exercised successfully. On a small
64-scan/32-charge smoke workload, preprocessing and solve were each only a few
milliseconds; the solve counters attributed nearly all measured solve CPU time
to regularizers, with projection, FFT, update, and convergence rounded below
the six-decimal display precision. This workload is too small to select an
optimization, so the next benchmark must use the larger merged grid from the
timing harness or a representative real chromatogram before direct convolution
is attempted.

The larger 64-scan/1,024-bin/32-charge profile selected the next exact change:
positive-floor `zsig` regularization consumed about 0.209 s of a 0.262 s solve,
while both reduced FFTs consumed about 0.040 s. Adding `#pragma omp simd` to
the independent `blur_it_UCCD` cell loop reduced the same regularizer phase to
about 0.136 s and the solve to about 0.192 s in the isolated Release build
(roughly 35% lower regularizer time; single-run directional evidence). The
native 29-case equivalence regression still passes. This keeps the existing
log/finite handling and is now the leading optimization to remeasure with
repeated runs before considering direct convolution.

Six fresh-input runs (one warmup, five measured) on the same zsig workload
confirmed the direction: the pre-SIMD optimized executable had a 0.380 s median
end-to-end time, while the SIMD build measured 0.285 s (about 25% lower,
0.283–0.292 s after warmup). The comparison uses identical settings and
executables except for the pragma; it is still synthetic and single-machine
evidence. The next regularizer experiment should therefore focus on preserving
this path while testing compiler vector-math options or reducing repeated
positive-floor logarithm work, with numerical regression after every change.

Adding `restrict` qualifiers to the same blur loop was also tested. The new
build passed the native regression and measured a 0.283 s median over five
post-warmup runs versus 0.285 s for the SIMD-only build; that difference is
within run-to-run noise, but the qualifiers are retained because they are
semantically valid and cost-free. No further exact gain has yet been shown.

The attempted two-pass positive-floor optimization in `blur_it_UCCD` was
incorrect and has been reverted. It wrote logarithms into the output scratch
buffer, then overwrote that buffer with exponentiated results while later cells
still gathered neighboring values from it. The result mixed log-domain and
intensity-domain values. The existing native oracle had `zzsig=0`, so it did
not exercise this path. The initial repair restored the SIMD single-pass
implementation. A `zzsig=1` native oracle case now covers this dependency.

On 2026-09-10, safe log caching replaced that repair: the first SIMD pass fills
the existing scratch buffer with sanitized logarithms, and the second gathers
from that immutable buffer while writing results directly into the data array.
This removes two of three logarithm evaluations per cell and the positive-floor
copy-back, without allocating another cube. Both UniChrom and UCCD callers now
use this in-place-data contract. The nonpositive-floor branch retains its
original arithmetic and copies its completed scratch result back internally.

An alternating old/new Release comparison on fresh copies of
`SEC_Native_Bispecific_Special.hdf5`, using four requested workers, one warmup
and three measured `-decon` runs per executable, reduced median process time
from 1.493 s to 1.105 s (26%). The final m/z and mass grids were bit-identical.
Synthetic negative-charge-smoothing output was also identical; consecutive
charge/mass smoothing and zero-padding comparisons had relative L2 errors
below 9.1e-7. These are single-machine measurements. The native oracle now also
covers negative charge smoothing on a wide m/z grid.

An attempted durable sparse-mask case using a narrow positive mass window did
not match the existing oracle, so that fixture was removed pending diagnosis.
The sparse implementation remains unclaimed by regression coverage; broad
window cases continue to exercise the dense fallback. Sparse correctness is
now a follow-up before further sparse tuning.

The mirrored direct-convolution experiment was fast on the synthetic timing
case but produced unacceptable `dtsig=1` GUI results on representative data.
It changed kernel truncation, m/z and scan boundary handling, iteration
normalization, and adjoint edge weighting simultaneously. Attempts to correct
the reflected-boundary sensitivity did not restore the expected result. The
experiment, its CMake option, and its special regression were removed; the
source and installed executable are back on the committed reduced-FFT solver.

On the synthetic GUI command path, the restored FFT build gives a median
normalized m/z-grid correlation of 0.994 between `dtsig=0` and `dtsig=1`
(relative L2 about 10.7%). Treat this as a smoke comparison rather than a
scientific acceptance result.

The scan-boundary fix is now implemented around the validated FFT operator.
Each side is padded by `ceil(3 * dtsig)` after `dtsig` is converted to internal
sigma units. Mirrored padding is the default. The internal configuration/HDF5
attribute `unichromzeropad=1` selects zero-filled padding instead. The forward
operation
extends and crops the scan axis; the adjoint center-pads and folds mirrored
samples back onto their source scans. An `H^T 1` sensitivity correction keeps
the Richardson-Lucy update normalized at the boundaries. Output reconvolution
uses the same extension and crop. A first-scan-only regression verifies that
the last scan remains below 1% of the first, and the numerical oracle covers
both boundary modes.

The exact reported `SEC_Native_Bispecific_Special.hdf5` file identified the
charge-smoothing regression above. Its configuration uses `zzsig=1` for 50
iterations. With the incorrect two-pass smoother, `dtsig=1` ended at a reported
convergence metric of 8517.46 and its mass grid was effectively unrelated to
the `dtsig=0` result (median normalized correlation -0.0025). With the restored
single-pass smoother, convergence is 0.04596; the median normalized correlation
between `dtsig=0` and `dtsig=1` is 0.993 for the m/z grid and 0.946 for the mass
grid. These comparisons use copies of the same source file and the GUI's
`-decon` then `-grids` command sequence.

On the same file after boundary padding, `dtsig=1` uses two padded scans per
side. Single runs measured 1.522 s for mirrored deconvolution and 1.527 s for
zero-padded deconvolution, versus 1.101 s for the separate `dtsig=0` command.
Mirrored and zero-padded outputs are nearly identical away from the boundaries:
their median normalized correlations exceed 0.9999999 for both m/z and mass;
the two edge scans have relative differences of about 0.18% in m/z and 1.28%
in mass. Against `dtsig=0`, mirrored padding retains median correlations of
0.993 for m/z and 0.944 for mass on this file.

Direct convolution was reconsidered after the padding fix using a profiled
Release build on `SEC_Native_Bispecific_Special.hdf5`. Three runs put the
iteration FFT phase at 0.140–0.150 s and output FFT at 0.075–0.077 s, while
all regularizers together took 0.959–0.964 s of the 1.232–1.245 s solve.
That timer includes point smoothing and does not isolate charge smoothing.
Thus, even removing every FFT operation would save at most about 0.22 s on this
workload, and a real replacement would save less. At three internal sigma the
Gaussian kernel would use five scan taps and seven m/z taps, but that compact
stencil would also omit nonzero tails; a full-length exact direct convolution
would not be competitive. No direct backend was added. Reconsider a separable
direct-scan/batched-m/z-FFT hybrid only for a representative workload where
profiling shows convolution, rather than regularization, dominates total time.

## Scope and objective

### 2026-09-10 command-path audit and post-cache profiling

The Python grid/peak importer calls `MetaUniDec.make_grids()`, which invokes
`-grids`. Positive `dtsig` previously bypassed command dispatch and reran
`run_chromatogram` for that request, doing another full solve and skipping the
normal peak refresh. Native `-grids` dispatch now reuses valid coupled grids
and executes MetaUniDec's grid/peak path. Missing grids are generated by
UniChrom once first. Charge extraction modes 6/7 report unsupported coupled
charge outputs instead of silently generating independent per-scan results.
Regression coverage checks missing/existing grids, preserved grid values,
peak generation, and the unsupported charge-extraction path.

An alternating Release A/B on fresh SEC file copies, one warmup and three
measured command sequences per executable, measured median `-decon` plus
`-grids` latency of 2.264 s before and 1.326 s after (41% lower).
Median individual deconvolution times were 1.114/1.111 s and grid-refresh times
1.103/0.215 s. Medians of components need not sum to the median total.
Both final grids were bit-identical. These timings use default worker settings
and include process startup/I/O, but not Python GUI import or rendering.

The profile build now reports point, charge, and mass smoothing separately.
Three post-cache SEC runs measured charge smoothing at 0.354–0.402 s, point
smoothing at 0.307–0.354 s, iteration convolution at 0.143–0.160 s, and output
FFT at 0.082–0.086 s. Mass smoothing was disabled. The aggregate regularizer
timer still includes softmax and suppression when enabled.

A cell-scheduled two-pass charge smoother was prototyped with an explicit
barrier between log-cache construction and neighbor reads. Four alternating
measured runs after warmup changed median process time from 1.135 to 1.122 s
(about 1%) and charge time from 0.346 to 0.322 s. It passed the small oracle,
but the SEC comparison exceeded the existing elementwise tolerance (43 of
162,300 m/z cells; overall relative L2 2.03e-5). This prototype was removed:
the small timing gain did not justify accepting changed numerical behavior.
The established scan-parallel log cache remains installed.

Next measured numerical candidate is point smoothing, now comparable in cost
to charge smoothing. Keep the FFT backend while these two phases dominate.
Other positive-width commands such as `-proc` and `-all` still need a separate
command-contract audit; this fix is deliberately scoped to grid/peak refresh.

Bring positive-`dtsig` processing close to the corresponding `dtsig=0` workflow,
while retaining chromatographic coupling and current numerical/output behavior.
HDF5 UniChrom is the primary track. UCCD is a separate optional track below;
it has not received this optimization and its charge dimension is measured.

Proposed performance acceptance target, to confirm against representative data:
median end-to-end time at most 1.25 times the zero-width baseline on the primary
workload, with 1.5 times as an interim milestone. These are planning targets,
not promises. Also report the full workload range, memory, iteration counts,
and numerical differences. Positive width solves a different problem, so the
zero-width result is a speed reference, not a numerical oracle.

The earlier analysis-only stage was followed by an authorized implementation
and validation stage. This update is documentation-only. Future implementation
tasks should include the focused builds and experiments below; do not infer
a new requirement to obtain permission before routine authorized validation.

Paths below are relative to `public/UniDec` unless explicitly stated otherwise.
Read repository `AGENTS.md` and `SKILLS.md`, preserve user changes, and keep
native binaries/build trees out of source patches. Use existing FFTW/MKL and
OpenMP support; no new sparse-matrix dependency is proposed.

## Current baseline and evidence

`UniChrom_Main.c` now sums latent charge before prediction, computes the
adjoint on the same chromatography/m/z grid, and applies one correction to
each charge row. It reuses that 2-D workspace for per-charge output
reconvolution. The iteration kernel remains unnormalized; only output
reconvolution uses unit-sum normalization. Preserve the mask order, periodic
boundaries, asymmetric peak shapes and cube-wide output normalization.

Validation already performed:

- Fresh isolated IntelLLVM 2025.1/MKL Release builds of old and new sources.
- 55 old/new cases covering peak shapes, output/normalization modes,
  regularizers individually and together, convergence, and zero-width routing.
  Worst relative L2 difference across compared datasets: `7.15e-7`.
- `tests/test_unichrom_native.py`: 29 configurations checked against the
  original 3-D formulation, including singleton scan/charge dimensions,
  zero m/z width, narrow chromatography and masked mass outputs.
- Existing UniChrom routing and UCCD binary tests passed. Packaged binaries
  were not replaced; select the actual newly built executable for future work.
- The post-adjoint correction now checks the allowed mask, finiteness and
  nonnegative result while multiplying, removing the separate full-cube cleanup
  pass. The native regression, UniChrom routing test and UCCD binary test pass
  after this change.
- Point smoothing now allocates `scan_count * numz` scratch sums, runs scans in
  one OpenMP region, and uses a serial charge-block loop per scan. This avoids
  sharing `smoothing_sums` and avoids nested parallel regions in UniChrom;
  other callers retain the original parallel wrapper.
- Mass output mapping now precomputes one lower-bin index and fraction for each
  m/z/charge cell, then reuses them for every scan. Native output equivalence
  remains within the existing tolerance. A repeat of the synthetic timing
  harness measured 0.142 s without regularizers and 0.389 s with them; the
  latter overlaps the earlier 0.374 s result, so this is not a claimed speedup.
- The strided `plan_many` output prototype was numerically correct but slower
  (0.163 s / 0.402 s in the same two synthetic rows), so it was removed. This
  rejects that layout, not batched output in general.

Synthetic native-command timings, seconds (median [minimum, maximum]):

| Settings | Old 3-D implementation | Current 2-D implementation | dtsig=0 command | Old/current speedup |
| --- | --- | --- | --- | --- |
| Regularizers disabled | 0.932 [0.907, 0.936] | 0.144 [0.144, 0.157] | 0.0805 [0.0765, 0.0854] | 6.5x |
| beta=0.5, zzsig=0.1, psig=1 | 1.213 [1.201, 1.247] | 0.375 [0.371, 0.398] | 0.0801 [0.0742, 0.0911] | 3.2x |

Protocol for the latest numbers: 64 scans, 1,024 input m/z samples (about 1,025 merged bins), 32
charges, Gaussian peaks, imported dtsig=2.5 and mzsig=1.2, numit=-20,
rawflag=2, datanorm=0, four requested OpenMP workers, one warmup and five
measured runs with alternating execution order. Each run used a fresh temporary
HDF5 input. Timings include native process startup and I/O, exclude fixture
creation and Python import, and are not phase profiles. Effective MKL thread
count was not separately recorded. No real-data or peak-memory benchmark has
been completed.

Compared with the pre-fusion 2-D run, the cleanup fusion and scan-parallel
smoothing together reduced the latest medians by about 6% without regularizers
and 12% with regularizers. Timing spread and machine load make these directional
results; the numerical result remained unchanged within the existing tolerance.

The contiguous output batch was compared with the retained per-charge output
path on the same synthetic workload. Its regularized median improved by roughly
4%; the unregularized difference was within timing noise. Treat this as a
bounded output optimization, not evidence that all batched FFT layouts help.

Important comparison limit: the default MetaUniDec command deconvolves scans
but does not enter the mode-3/mode-4 branch that generates merged grids. UniChrom
writes merged grids during its command. The observed current/zero ratios
(about 2.05x and 5.73x) therefore compare command latency with different output
work; they do not establish the gap for the same user-visible workflow.
The old/current speedups compare the same positive-width output contract.

The temporary comparison harness and raw timing files were not added to the
repository; their results are summarized above. Make representative benchmarks
reproducible before building on these numbers. Reuse the regression source at
`tests/test_unichrom_native.py`, currently present as a new working-tree file,
and include it with the implementation when preparing the changes for review.

## What the source establishes

| Location | Observation and implication |
| --- | --- |
| `unidec/src/UniDec.c`, `main` | Meta mode routes positive `dtsig` to `run_chromatogram`; zero goes to `run_metaunidec`. CDMS mode 2 routes separately to `run_unidec_UCCD`. |
| `unidec/src/UniChrom_Main.c`, `make_kernel_fft_UniChrom` | Builds the equivalent charge-independent 2-D kernel directly using `PeakDist` and the periodic scan peak formula. |
| Same file, `run_chromatogram` | Each iteration does four Q-cell FFT executions, plus dense charge projection and latent updates. There is no cube-sized ratio broadcast. |
| Same file, FFT initialization | Plans and kernel spectrum are already reused. Real transforms already exist. Recommending either as a new optimization would miss the actual opportunity. Plans use `FFTW_ESTIMATE`; the execution path now applies the available MKL/FFTW thread configuration around the reduced plans. |
| Same file, allocation/update loops | `blur`, `scratch`, `oldblur` still scale with the cube; FFT buffers scale with the measured grid. Pre-forward cleanup remains a full cube pass; post-adjoint correction now fuses cleanup into multiplication. Copies, smoothing and reductions still traverse the cube. |
| Same file, output reconvolution | For rawflag 0/2, charge planes are copied into contiguous temporary batches and transformed with one `plan_many` pair, then scattered back and masked. Allocation/plan failure retains the per-charge fallback. |
| Same file, preprocessing/output | Every call processes all spectra, merges through HDF5, reads the merged arrays back, and eventually writes merged and per-scan mass outputs. These costs can dominate after iteration acceleration. |
| `unidec/src/MetaUniDec_Main.c` | The zero-width baseline has independent-scan parallel fast deconvolution, conditional on `rawflag > 1`, no manual assignments, no double deconvolution, and available workers. Benchmark the actual branch used. |
| `unidec/src/UCCD_Main.c` | Sparse binary input is expanded to dense data. Sparse correction builds ratios only at observed nonzero indexes, then clears/scatters into a dense workspace and executes dense FFTs. Sparse I/O is not sparse convolution. |
| Same file, `batched_2d` | Zero width uses batched m/z-charge 2-D FFTs and compacts empty scans. Positive width uses a full 3-D FFT and retains the original scan geometry. Thread setup and `FFTW_MEASURE` already exist here. |
| `unidec/src/udtools.c`, `MakeKernel2D`, `PeakDist`; `UCCD_Main.c`, `make_kernel3D_UCCD` | The kernel factors into one-dimensional axis kernels. Corner-image sums and the chromatography two-image construction define the exact sampled periodic kernel. |
| `unidec/src/udstruct.c`, `PostImport`; `h5io.c` | Imported widths are converted to internal units; chromatography width is divided by 2.35482. Do not apply the conversion again when building direct kernels. |

For `T` scans, `M` m/z bins, `Z` charges, let `N=T*M*Z` and `Q=T*M`.
The HDF5 iteration convolution work is O(Q log Q), with O(N) projection and
latent updates. Main persistent solver buffers now consume
`12N + 4Q + 16*T*(floor(M/2)+1)` bytes: three float cubes, one float FFT grid,
and two complex spectra. This excludes observed data, regularizer tables,
outputs and FFT-library internal allocations; measure peak resident memory
rather than presenting this subtotal as process memory. Output reconvolution
still costs O(Z*Q log Q), plus charge-plane packing/scattering. The retained
batch reduces FFT plan executions while paying one contiguous temporary copy.

## Priority 1: profile the remaining work, then fix the largest phase

Do this before choosing between direct convolution and sparse storage.

1. Establish a reproducible current-build baseline and the zero-width workflow
   through the point where both supply the grids needed by the same Python
   consumer. Report command latency separately. Trace the actual caller's
   follow-up grid/import calls instead of switching both to `-all`, which may
   also add peak extraction or other work. Record iteration counts, per-scan
   overrides, effective worker counts and merged versus per-scan grid lengths.
2. Time preprocessing/merge, FFT planning, each regularizer, projection,
   forward/adjoint FFTs, update/masks, convergence, output gathering, output
   FFTs, scattering, mass transform and HDF5 writes separately. Use monotonic
   wall time. Keep instrumentation temporary or behind an existing diagnostic
   mechanism; do not add a performance-settings UI.
3. Prioritize regularization and cube traffic on regularized workloads.
   Enabling the tested regularizers now increases current command time by about
   0.228 s, roughly 61% of the 0.374 s total. This is a configuration comparison,
   not a measured phase share or a guaranteed removable cost. Isolate beta,
   charge smoothing and point smoothing in turn. If that increment remains,
   even eliminating all work in the unregularized run would not reach the
   measured zero-width command time. FFT replacement alone is unlikely to
   close this case's gap.
4. Measure output reconvolution explicitly. For I iterations there are 4I
   iteration FFT executions and another 2 output FFT executions for rawflag
   0/2 in the retained batched path. The earlier per-charge path required 2Z;
   counts alone are not a timing
   breakdown. The retained path uses bounded contiguous packing; compare it
   against separable direct passes if this phase remains significant. Keep the
   per-charge fallback for allocation/plan failure. Do not allocate a full
   cube FFT again by default or parallelize the fallback charge loop using its
   shared mutable FFT workspace.
5. Two exact traffic/parallel candidates are complete: post-adjoint correction
   and cleanup are fused, and point smoothing uses one scan-level parallel
   region with private sums. Reprofile before attempting another pass fusion.
   Next candidates are precomputing reused regularizer values or reducing
   measured output FFT/gather costs. Concrete correctness constraints are
   listed below. Reprofile after each one.

Retain every output charge through reconvolution and masking before mass
mapping. A charge-summed output cannot recover mass, and per-plane maximum
normalization would change the result. The current inverse FFT scale is 1/Q;
do not introduce an extra charge factor or RL sensitivity normalization.

## Priority 2: separable direct convolution and a small hybrid

If profiles justify it, evaluate these implementations of H, the current
scan/m/z convolution. The reduced FFT implementation already exists:

| Candidate | HDF5 per-application work | Suitable case |
| --- | --- | --- |
| Reduced 2-D FFT | O(Q log Q) | Broad kernels, long tails; exact sampled reference |
| Direct scan pass + direct m/z pass | O(Q*(Kt+Km)) | Short effective support in both axes |
| Direct scan pass + batched m/z FFT | O(Q*(Kt+log M)) | Narrow chromatography and broad m/z kernel |

`Kt` and `Km` count retained taps, including wrapped offsets. Use two
one-dimensional passes, not a product stencil with `Kt*Km` work. Width zero
is an identity pass. Very small positive width is not automatically zero:
omitting nonzero taps is an approximation even if their effect is tiny.

First construct full-length axis arrays that reproduce `MakeKernel2D` and
`periodic_scan_peak_UCCD`, including origin/end image contributions. Preserve
the existing periodic boundaries. Do not substitute a minimum-distance
Gaussian, zero padding, reflection, or per-edge renormalization. The split
Gaussian/Lorentzian shape is asymmetric; the adjoint uses reversed periodic
offsets. Derive each axis from its actual helper: charge split-shape orientation
differs from m/z. Keep normalization explicit and apply it once.

Full-length direct convolution is an exact arithmetic reference, not necessarily
fast. For a fast compact stencil, choose support from omitted sampled-kernel
L1 mass relative to the full discrete kernel sum. Record tail mass and output
error separately. Gaussian kernels can be short; Lorentzian and split kernels
can have significant long tails. Keep an FFT fallback when support is broad.
Do not assume a Gaussian radius works for all peak shapes or silently truncate
tails. Any nonzero-tail omission needs an explicitly accepted accuracy budget
before becoming default. Retaining all numerically nonzero taps is another
exact option when float underflow already makes support compact.

Reuse kernel taps across iterations. Traverse contiguous m/z rows for scan
passes, use disjoint output tiles, and avoid atomics. Start with one measured
choice between direct and FFT; only add a dimension-specific hybrid if it
wins materially. No public backend selector or tuning framework is needed.

Benchmark iteration and output operators separately: direct convolution on
Q cells might save little iteration time but help the Z-plane output pass,
or vice versa. A kernel width in physical m/z units can span many bins after
merging, so record actual tap counts rather than assuming a small numeric
mzsig implies a cheap direct operator. For retained FFTs, measure any thread
or planning change on the smaller Q-cell domain; carrying over UCCD's thread
cap or `FFTW_MEASURE` policy without measuring startup can make short runs slower.

## Priority 3: exploit sparsity where it survives the algorithm

Distinguish three independent quantities: observed nonzero fraction, fraction
of allowed m/z-charge cells, and latent occupancy after regularization. Sparse
input does not establish sparse intermediate arrays.

1. **Static allowed-cell lists for HDF5.** The existing allowed mask is shared
   by all scans. Initially use per-m/z lists of allowed charge indexes for
   projection and multiplication, retaining dense latent storage. This avoids
   CSR metadata and a full regularizer rewrite. Ensure disallowed values are
   zero at the same stage as before. Keep dense iteration for high occupancy.
2. **Sparse observed-ratio adjoint.** For finite valid data, the ratio is zero
   where observations are zero. A direct compact adjoint can scatter from
   those observations, or gather only relevant contributions. Compare total
   work, including output clearing, indexing, thread reductions and the support
   expansion between separable passes. In 2-D a naive sparse scatter costs
   O(nnz(observed)*Kt*Km), which can lose to dense separable passes.
3. **Compressed latent storage only after evidence.** A static compressed
   m/z-charge pattern repeated over scans could reduce memory from N to
   `T*allowed_count`, but neighbor tables and regularizers need exact remapping.
   Softargmax includes the full charge count and zero entries in its formula;
   deleting those entries changes it. Smoothing can create intermediate values
   before the final mask. Missing neighbors cannot simply be dropped or
   renormalized. Budget the mapping and output expansion costs.

The large cube FFT buffers have already gone, so recalculate the sparse
memory/performance benefit using the current buffer sizes. Profile density
after each regularizer; support growth can erase the advantage. Prefer static
allowed-cell lists and exact loop pruning to maintaining dynamic sparse
matrices. Sparse gather/scatter should earn its place against the current
small dense measured grid, not the removed full-cube FFT implementation.

Do not construct a giant explicit convolution matrix: it duplicates translation-
invariant coefficients and multiplies storage by stencil size. Native compact
taps or a few axis operators are the initial sparse representation. Generic
CSR/CSC is a later option only if the operator becomes irregular and measurements
justify it. Do not prune latent values by an arbitrary intensity cutoff or
delete empty scans when positive `dtsig` couples their original positions.

## Separate UCCD track for the linked source file

UCCD observes `[T,M,Z]`; it has no charge projection/broadcast identity to exploit.
Its kernel still factors as `kt*km*kz`. Prioritize:

1. Retain the current dtsig-zero branch as the speed and regression baseline.
2. Prototype direct chromatography convolution plus the existing batched
   m/z-charge FFT approach for positive width. Preserve all original scans.
   Apply the scan pass and m/z-charge pass in the forward operator and their
   transposes in the adjoint. Keep the current unnormalized iteration kernel,
   `predicted != 0` ratio guard, and output rules. Unlike HDF5 UniChrom, UCCD
   reconvolves only for `rawflag == 0`, with charge width zero at output.
3. Compare fully separable direct passes: O(N*(Kt+Km+Kz)) per operator versus
   O(N log N) for the current cube FFT. Skip identity axes (`csig == 0`, etc.).
4. Reuse observed nonzero indexes for a sparse direct adjoint only if support
   expansion and reduction costs are favorable. Existing sparse correction
   already saves ratio arithmetic; repeating that change will not save FFTs.
5. Quantify how much of the zero-width advantage comes from empty-scan
   compaction. With positive width, scan deletion changes distances and the
   circular domain. Any future tiling must retain original indexes and provide
   halos for both forward and adjoint applications each iteration; independently
   converging cropped windows is a different algorithm.

UCCD already initializes FFTW threads or sets MKL thread limits, caps default
workers, and uses measured plans. Measure planning separately from execution;
do not describe these as missing optimizations. Do not share stateful FFT
contexts across simultaneous workers.

## Candidate changes for the measured bottleneck

- Fuse HDF5 charge projection with required pre-forward cleanup only if a phase
  profile shows the remaining cleanup pass matters. Correction multiplication
  and post-update cleanup are already fused; preserve its allowed-mask,
  finiteness and nonnegative checks, plus the checkpoint schedule and meaning
  of `oldblur` (previous checkpoint, not necessarily previous iteration).
- `blur_it_UCCD` repeats log evaluations for each neighbor in positive-floor
  regularization. Do not reuse its output buffer for a separate log pass: its
  neighbor gathers make that in-place strategy order-dependent. Any future log
  cache requires distinct immutable storage and must pass the `zzsig` oracle.
- `softargmax` has its own OpenMP loop inside HDF5's scan-parallel call; avoid
  introducing another nested level there. Point smoothing's UniChrom path is
  now scan-parallel with private sums; keep its other callers unchanged and
  never reintroduce shared scratch counters under an outer pragma.
- Allocate inactive z/m regularizer tables only when used. Mass-bin lower
  indexes and interpolation weights are now precomputed after mass-axis
  selection; preserve output summation order and non-finite handling under
  review. Further output work should target measured FFT/gather costs.
- Measure preprocessing, merger I/O and output writes before changing them.
  A small in-memory merger handoff may remove a write/read cycle, but preserve
  dataset creation and downstream expectations. Do not skip preprocessing
  without a complete freshness rule for raw data and processing settings.
  Keep HDF5 access on the owning thread, following existing MetaUniDec practice.
- Inspect merged-axis inflation: `make_grid` uses global range and `mzres`.
  Report the resulting M versus per-scan lengths. Coarsening bins or narrowing
  scientific ranges is a quality tradeoff, not a transparent optimization.
- Do not reduce iteration caps, relax convergence, drop regularizers, or use
  independent deconvolution followed by smoothing to claim equivalent speed.
  GPU support, new dependencies and solver acceleration are deferred until
  simpler changes miss the target and the remaining cost is measured.

## Work packages for the next agent team

Use separate branches/worktrees or exclusive file ownership. Numerical agents
must not concurrently edit `UniChrom_Main.c`. The coordinator integrates one
candidate at a time; benchmark and review work can run independently.

| Agent / package | Deliverable | Dependency / gate |
| --- | --- | --- |
| A: baseline and profiling | Reproducible real/synthetic workload runner, matched-output zero-width comparison, phase timings and memory for current sources | Start here; reuse existing numerical fixtures |
| B: regularizers and cube traffic | Small exact changes to the dominant regularizer, parallel structure or remaining memory passes, with per-change measurements; cleanup fusion and scan-parallel point smoothing are complete | A phase profile; exclusive HDF5 source ownership |
| C: convolution and output | Direct/separable iteration candidate; bounded contiguous output batch is complete; assess kernel alternatives and tail/crossover measurements; mass-bin mapping cache is complete | A profile; prototype on a separate branch and integrate sequentially with B |
| D: correctness and sparsity review | Extend existing native regression for changed branches, audit contracts and sparse feasibility; review B/C equivalence | A; compressed storage only if current density/memory measurements justify it |

UCCD's hybrid is a separate package if that workflow is requested; no HDF5
implementation task remains for charge collapse. Each package reports
changed files, semantic assumptions, numerical error, timings, memory and
remaining limitations. Review exact transformations before approximate ones.

## Validation and performance protocol for future execution

Use Release builds with the same compiler, FFT backend and requested worker
count. Follow CMake 3.25/HDF5/FFTW requirements and existing OpenMP options.
Use fresh out-of-tree builds and identify the actual executable launched;
normal CMake builds copy tracked artifacts to `unidec/bin`, so review those
separately. Never modify original HDF5 benchmark inputs: both paths write them.

1. Capture dimensions, peak-shape selectors, imported/internal widths, m/z
   spacing, all regularizers, `rawflag`, `datanorm`, fixed-mass-axis setting,
   allowed fraction, observed density, build identity and worker settings.
   Choose the primary real workload plus small/large, narrow/broad kernel,
   sparse/dense and few/many charge cases. Include empty scan gaps and signals
   at the periodic boundaries.
2. Benchmark current 2-D positive width, the next candidate, and zero width
   on separate copies with otherwise identical settings and matched required
   outputs. Retain an old 3-D build only as an optional historical reference.
   Use a fixed iteration count (negative `numit` suppresses early stopping in these
   loops), then also test normal positive-`numit` convergence. Report both;
   equal caps do not establish equal executed iterations.
3. Measure wall time externally and by phases: preprocessing/merge, allocation
   and planning, regularization, projection, forward, ratio, adjoint, update,
   convergence, output reconvolution, mass transform, I/O and Python import.
   Do not treat HDF5 UniChrom's `clock()` print as portable wall time. Report
   first-run latency and warmed median/range from at least five comparable
   runs, peak resident memory, and iterations. Avoid mixing Python import
   overhead into a native-only result without labeling it.
4. Test operators independently: impulse at center and each edge, asymmetric
   peak shape, constant input, axis sums, full-support direct versus FFT,
   and the adjoint identity `<H x,y> ~= <x,H^T y>`. Preserve existing
   charge-collapse regression and extend it for the changed branch with
   unequal charge intensities and a nontrivial mask. Check both raw output
   and normalized per-charge reconvolution.
5. Test one update and a fixed sequence against the current validated 2-D
   implementation, then convergence/output behavior with beta, psig, zsig,
   msig and suppressions individually and in combination. Initial suggested
   exact-transform gates for normalized fixtures: relative L2 <= 1e-5 for
   operators and <= 1e-4 for fixed-iteration outputs, with an explicit
   scale-aware absolute tolerance near zero. These are provisional gates;
   investigate failures rather than silently loosening them. Inspect weak
   peaks, mass locations, areas and chromatographic widths, not only norms.
   Approximate truncation requires its own accepted scientific error budget.
6. Exercise T=1, Z=1, smallest supported M, odd/even sizes, empty input,
   all-disallowed constraints, non-finite input, zero/invalid widths, numit=0,
   zero denominators and large dimension products. Current UniChrom handles
   singleton charge in its new kernel builder and validates cube size in a
   wide type before allocating its latent buffers. Preserve these checks.
   Upstream `make_grid` still needs a size-arithmetic audit for oversized
   inputs; the shared `MakeKernel2D` used by other paths still reads the
   second-last axis element for singleton dimensions. Do not reintroduce that
   helper's singleton access in HDF5, or use undefined behavior as an oracle.
7. Verify HDF5 root datasets `mz_grid`, `mz_axis`, `mz_sum`, `mass_grid`,
   `mass_axis`, `mass_sum` below `/ms_dataset`, flattened scan-major shapes,
   per-scan `mass_data`, length/width attributes and got-grids marker.
   Read through `unidec/metaunidec/mudstruct.py::import_grids` and spectrum
   import. Preserve global normalization and fixed/dynamic mass-axis rules.
   For UCCD verify sparse binary headers, original axes, record indexes and
   expanded scan ordering. Test zero-width dispatch separately from numerical
   equivalence; `test_unichrom_workflow.py` mocks the native call whereas
   `test_unichrom_native.py` exercises native HDF5 outputs numerically.
8. Put regression coverage in `tests/` using UniDec's unittest style. Run
   existing and extended native numerical tests, `tests.test_unichrom_workflow`,
   applicable CD-stack/workflow tests, and native smoke coverage after builds.
   Run GUI-launch tests only if GUI construction/bindings change. Record exact
   dependency/fixture limitations rather than claiming unrun checks passed.

For an isolated executable, set `UNIDEC_TEST_EXECUTABLE` to its path and run
`python -m unittest discover -s tests -p test_unichrom_native.py -v` from
`public/UniDec` using the configured interpreter. Discovery avoids a possible
collision with an installed package named `tests`. Extend the existing test
instead of writing a second copy of its 3-D oracle. Its current independent
oracle disables regularizers; the broader 55-case regularizer comparisons
were temporary old/new runs, so add durable coverage for regularizers actually
changed next. The 29 cases are not comprehensive coverage of malformed inputs,
all empty/sparse patterns, all thread backends, or real chromatograms.

The final retained source was rebuilt after reverting the `FFTW_MEASURE`
experiment. The current isolated executable passes the native equivalence
regression; focused routing, UCCD binary, native executable and shared-library
smoke tests also pass. A build invoked from the source directory reported zero
discovered tests because the test paths are package-relative; the subsequent
package-root invocation is the authoritative test run.

Report `Tcandidate_positive/Tzero_matched_outputs` and `Tcurrent_2D/Tcandidate`
for each workload, alongside quality and memory. If the target remains out of
reach, show the measured limiting phase and estimated remaining headroom.
Stop adding mechanisms when the simplest validated implementation meets the
target; retain broader-kernel fallback only where results justify it.

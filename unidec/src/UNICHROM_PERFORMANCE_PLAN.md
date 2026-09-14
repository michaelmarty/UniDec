# UniChrom performance investigation and implementation handoff

## Reassessment after charge-smoothing tiling (2026-09-11)

Recommendation: retain the current solver and move from open-ended optimization
to validation plus at most two bounded performance experiments. The current
implementation is near a sensible stopping point for this SEC workload, but
the evidence does not establish a hardware or algorithmic speed limit. The
historical priorities below are an experiment archive; this assessment takes
precedence when selecting new work.

The latest unprofiled native-command medians are 1.305 s positive width and
0.979 s zero width. Reaching 1.25 times that baseline would require about
0.081 s more saving, or 6.2% of current positive time. Charge smoothing still
takes about 0.612 s: reducing it another 13% could theoretically close that
gap if other costs stayed constant. This is headroom arithmetic, not a forecast.
The last exact change saved only about 0.047 s of process time. Large additional
speedups have no demonstrated low-risk route under the existing numerical rules.

### Work worth doing before declaring completion

1. Preserve a small reproducible benchmark runner and measure the actual GUI
   consumer workflow on current sources, including required grid generation and
   import. The latest 1.33 ratio came from `-decon` alone; it is not the earlier
   matched `-decon` plus `-grids` comparison. Record build identity, input and
   settings, iteration counts, first-run latency and warmed distributions. Do
   not compare different widths or historical machine-load conditions as if
   they were one experiment. Add a second representative real chromatogram.
2. Add durable coverage that actually enters the tiled smoother: multiple scans
   with more than 65,536 cells per scan, a partial final tile and neighbor reads
   across tiles, tested with positive charge and mass smoothing together and
   more than one worker count. The current 39 oracle configurations use small
   grids and exercise the per-scan fallback; the SEC comparison is currently
   the main evidence for the tiled branch. Keep singleton and negative-floor
   coverage as separate cases.
3. Treat sparse-mask correctness as partially established. Sparse/dense native
   equality rules out the static-list optimization as the cause of the observed
   discrepancy. It does not prove independent correctness where whole m/z rows
   have no allowed charge. Reconstructing observations from masked initialization
   loses those rows; the attempted broad-mask reconstruction also failed.
   Resolve that reference/numerical-conditioning case before claiming universal
   sparse equivalence. No production change is justified by that discrepancy yet.

### Remaining performance candidates, in order

| Candidate | Why it remains useful | Limit / decision rule |
| --- | --- | --- |
| Sweep existing worker settings (1, 2, 4, 8 where supported) | Lowest implementation cost; current conclusions use four workers, while scan count, tiles and small FFTs may scale differently | Measure total latency and memory on both real workloads; do not infer a universal default from one machine |
| Fuse pre-forward cleanup with charge projection | Removes a cube traversal and parallel-region boundary while retaining the exact cleanup point | The combined cleanup/projection/extension phase was about 0.119 s before tiling; only part is removable. Preserve stored zeros, finiteness checks, charge summation order and sparse/dense branches |
| Further smoothing tuning | Largest remaining phase; bounded tile-size or scheduling tests may help a different scan/worker shape | Optional only after the worker sweep; retain SIMD loops and the log/gather barrier. Do not keep a tuning framework or another arithmetic implementation for noise-level gains |

Do not prioritize an HDF5 handoff rewrite: preprocessing, merge and output
writes together were roughly 0.013 s, and only a fraction is avoidable. Output
reconvolution was about 0.100 s including about 0.064 s of transforms; even a
large percentage improvement saves little absolute time. Iteration FFTs were
about 0.128 s combined. Direct/hybrid convolution remains a conditional option
for a demonstrably FFT-dominated workload, with the current padded forward and
folded adjoint preserved. The earlier direct experiment failed scientific
acceptance, so it should not be revived merely to reach a timing ratio.

Compressed latent storage, GPU execution, reduced precision, relaxed log/exp
math, fewer iterations and altered kernel tails are deferred. These require
substantial complexity or a changed accuracy contract. Reconsider them only if
larger real datasets expose unacceptable latency or memory use. UCCD remains a
separate workload, not an automatic continuation of this HDF5 effort.

Stop after the validation work and bounded experiments if no repeatable,
practically useful gain remains. A reasonable proposed retention gate is at
least 3% end-to-end improvement for an added mechanism, above run-to-run noise,
with numerical equivalence and no material regression on the second workload.
The 1.25 ratio is a planning target, not a reason to accumulate complexity:
positive width performs a different coupled solve. If about 1.3 s and 129 MB
are acceptable in the actual workflow, the present architecture is a good
place to stop.

### Validation and bounded experiments completed (2026-09-14)

The completed benchmark used fresh temporary input copies and measured the
full GUI consumer path: configuration export, coupled native solve, per-scan
HDF5 import, grid refresh, peak picking/import and peak-object construction.
It recorded first and warmed runs, iteration counts, phase distributions and
peak native working set on Windows.

The worker sweep used unprofiled executable SHA-256
`0c43a9d7b5b18488f688514ab869834116a9b8a28b8ba87cac9afd49bf121a62`,
50 iterations, one first run and five warmed runs per point. Times below are
warmed medians in seconds; memory is the median peak native working set.

| Input | Workers | Full GUI consumer path | Coupled native solve | Peak MiB |
| --- | ---: | ---: | ---: | ---: |
| Bispecific, 10 scans, 16,230 m/z, 50 charges | 1 | 3.631 | 3.351 | 128.9 |
|  | 2 | 2.070 | 1.885 | 129.0 |
|  | 4 | 1.468 | 1.291 | 129.3 |
|  | 8 | **1.342** | **1.106** | 129.5 |
| Herceptin, 21 scans, 6,760 m/z, 50 charges | 1 | 3.252 | 3.006 | 154.9 |
|  | 2 | 1.938 | 1.714 | 155.0 |
|  | 4 | 1.361 | 1.124 | 155.2 |
|  | 8 | **1.147** | **0.897** | 155.4 |

Eight workers was fastest for both fixtures on this machine. This is a local
tuning result, so no application-wide default was changed. At eight workers,
the warmed `-grids` refresh took about 0.049--0.083 s and the Python work around
the solve and peak refresh added about 0.21--0.24 s. First full-path runs ranged
from 1.29 to 4.18 s across the sweep and were retained separately rather than
mixed into the warmed distributions.

The native oracle now includes a dedicated two-scan case whose merged grid has
65,540 m/z/charge cells per scan. It therefore enters the tiled positive-floor
smoother, leaves a partial second tile, and makes mass-neighbor reads across the
tile boundary. Positive charge and mass smoothing run together. Results match
the independent full-cube reference at one and four workers, and the two native
worker results are bit-for-bit identical. Existing singleton and negative-floor
cases remain separate.

The cleanup/projection fusion was implemented, built and tested, then rejected.
All native oracle tests passed, but seven alternating fresh-copy runs at eight
workers improved the full GUI path only 1.95% on Bispecific and 1.40% on
Herceptin. Native solve medians improved 0.44% and 2.94%; peak memory was
unchanged. This misses the 3% end-to-end gate on both workloads, so the fused
implementation was removed. These results support stopping performance work on
the current solver until a materially slower workload or changed accuracy
contract justifies another experiment.

## Implementation history

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
the native executable. The sparse list path now has a dedicated positive-mass-
window fixture in which every m/z row retains at least one valid charge; the
existing broad-window suite continues to exercise the dense fallback.

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

The earlier sparse-mask oracle failure was traced to the fixture rather than
the sparse implementation. Its wide m/z range contained rows with no valid
charge, while the oracle reconstructed the observed signal from a masked zero-
iteration output and therefore lost those rows. Forced sparse and dense native
paths were bit-identical across the expanded workload matrix. A revised narrow
positive-mass-window case keeps at least one valid charge in every m/z row and
now gives durable 3-D oracle coverage of the sparse path.

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

Point smoothing now writes from the active cube directly into the existing
scratch cube and swaps the two pointers. This removes its full-cube input copy
without changing the rolling-sum order or mask behavior. Matched IntelLLVM/MKL
Release builds produced bit-identical m/z and mass grids, axes and sums on the
SEC workload. In five alternating post-warmup runs with default worker settings,
median `-decon` process time fell from 1.001 s to 0.945 s (5.6%). Three profiled
runs reduced median point-smoothing CPU time from 0.202 s to 0.161 s (20%) and
solve time from 0.844 s to 0.798 s. The native 39-case numerical regression and
focused routing tests pass with the isolated candidate executable.

The positive-width command audit is now complete for the native commands.
`-proc`, `-extract`, and `-peaks` use their MetaUniDec operations without
running UniChrom. `-all` runs the coupled solve once and then refreshes merged
peaks through the existing `-grids` path. `-newgrids` likewise rebuilds the
coupled grids once before peak refresh. Commands that require unavailable
charge-resolved per-scan outputs (`-ultraextract`, `-charges`, and `-scanpeaks`)
return status 12 before doing a solve. The native regression covers processing,
the single-solve `-all` path, generated grids/peaks, and unsupported commands.

On the SEC workload with `OMP_NUM_THREADS=4`, five post-warmup command runs put
the corrected `-proc` median at 0.040 s versus 1.340 s for the previous
full-solve behavior. Corrected `-all` measured 1.583 s versus 1.605 s for
separate `-decon` and `-grids` processes. The single-process and two-process
outputs were bit-identical for both grids, both axes, both sums, peak data and
extracts. Command-specific timings use `OMP_NUM_THREADS` because the native
`-nthreads` option is accepted only as the primary command.

A matched current-build comparison now measures the complete `-decon` plus
`-grids` workflow for both widths. With four workers, five post-warmup runs gave
a 2.032 s positive-width median and 1.728 s zero-width median, a ratio of 1.18.
The ranges were 1.630-2.422 s and 1.633-1.850 s respectively, so machine load
was material. Peak working set was stable near 207 MB for positive width and
62 MB for zero width. The positive run reached its configured 50-iteration cap.
This workload meets the proposed 1.25-times speed target, while memory remains
the clearest gap.

VTune on the four-worker positive run attributes 0.591 CPU-s to vectorized
`logf`, 0.193 CPU-s to vectorized `expf`, and 0.722 CPU-s to OpenMP fork
barriers across the process. Charge smoothing remains the largest isolated
regularizer phase, but it already evaluates one log per cell and uses vector
math. Do not trade exact output for relaxed transcendental math. A future exact
candidate should target parallel-region overhead or memory only when a broader
workload misses the accepted target.

The next exact memory pass is complete. Charge, mass, harmonic-suppression and
point-smoothing tables are now allocated only when their corresponding mode is
active; the temporary coordinate grids are skipped when both charge and mass
smoothing are disabled. On the SEC configuration this removes the two inactive
mass-neighbor tables, about 6.5 MB during the solve and output phases.

Output reconvolution now reuses one FFT plan for batches of at most four charge
planes instead of allocating all 50 planes together. Relative to the matched
full-charge batch, peak working set on the SEC workload fell from about 207 MB
to 127 MB (39%). An eight-plane intermediate measured 134 MB. Five post-warmup
runs put the isolated output phase at 0.095 s for the full batch, 0.098-0.112 s
for the eight-plane batches, and 0.103-0.147 s for four planes; the overlapping
ranges and bimodal process timings do not establish a total-time change. The
four-plane outputs were bit-identical to the full and eight-plane results for
both grids, both axes and both sums, including a partial final batch. The native
expanded numerical and command regression passes. Retain four planes for the
measured memory reduction; the established per-charge allocation/plan-failure
fallback remains unchanged.

An unprofiled IntelLLVM/MKL Release build was then measured on seven alternating
post-warmup SEC runs with four workers. Positive `dtsig` had a 1.353 s median
(1.333-1.951 s) and 132.8 MB median peak working set; zero width had a 0.988 s
median (0.924-1.059 s) and 67.2 MB. The 1.37 timing ratio meets the 1.5 interim
target but not the proposed 1.25 target on this run. Repeated outputs were
bit-identical within each width, and the unprofiled executable emitted no phase
diagnostics.

The follow-up allocation and workload matrix also passes. None, charge, mass,
point, harmonic, and all-regularizer allocation modes were bit-identical to the
pre-memory-change executable across every output dataset. Narrow and broad
Gaussian kernels, a broad Lorentzian kernel, dense and sparse masks, an empty
scan gap, m/z-edge signals, and a scan-edge signal were likewise bit-identical.
The independent 3-D oracle checked the gap and both boundary-signal cases with
worst relative L2 error `7.41e-7`.

The optional profile now uses monotonic wall time when OpenMP is available and
separates setup, every regularizer, projection, forward/ratio/adjoint/update,
output planning/gather/transform/scatter, mass mapping and HDF5 writes. On the
SEC fixture's actual positive-width settings (`dtsig=1`, `zsig=1`, `psig=1`,
10 scans, 50 charges, `rawflag=0`), five post-warmup fixed 50-iteration runs
with four workers had a 1.296 s median profiled total and 1.327 s median process
latency. Peak working set was 129.2 MB. Normal positive `numit=50` also executed
all 50 iterations and measured 1.281 s, so early convergence provides no saving
for this fixture. A separately rebuilt unprofiled executable measured a 1.314 s
median (1.275-1.351 s); its outputs were bit-identical to the profiled build and
it emitted no diagnostics.

Charge smoothing is the next measured bottleneck: its 0.647 s median plus
0.128 s point smoothing accounts for 0.775 s, or about 69% of the 1.116 s
solve. Projection took 0.119 s, forward and adjoint transforms 0.065 s and
0.063 s, and output reconvolution 0.100 s, of which the transforms took
0.064 s. Setup was 0.039 s; preprocessing, merge, mass mapping and HDF5 output
were each at or below 0.015 s. Isolated 50-iteration medians were 0.109 s for
beta, 0.130 s for point smoothing, 0.638 s for charge smoothing, 0.575 s for
mass smoothing, and 0.121 s for harmonic suppression. The next exact candidate
should therefore reduce synchronization or cube traffic in the shared charge/
mass smoothing path while preserving its cached-log arithmetic. Direct
convolution, output changes and an in-memory HDF5 handoff are lower priorities
on this workload.

The first charge-smoothing parallel experiment collapsed every scan and cell
into two whole-cube OpenMP loops. Although numerically exact, it disrupted the
compiler's efficient per-scan vector-math loop and increased the charge phase
to about 1.01 s, so it was reverted. The retained implementation instead splits
each large scan into 65,536-cell tiles, distributes those tiles within one parallel
region, and keeps the existing SIMD log and neighbor/exp loops inside each
tile. Single-scan, small-scan and nonpositive-floor calls retain the established
per-scan helper.

Seven alternating profiled SEC comparisons reduced median charge smoothing
from 0.667 s to 0.612 s (9%) and profiled process latency from 1.351 s to
1.302 s (4%). Seven alternating unprofiled runs measured 1.352 s before and
1.305 s after, a 3.6% end-to-end improvement; the matched zero-width median was
0.979 s, giving a new positive/zero ratio of 1.33. All output datasets were
bit-identical across the SEC comparison and the 14-case allocation/kernel/mask/
gap/boundary matrix. Charge smoothing remains the largest measured phase, but
future work must retain the vectorized tile loops and justify more complexity
against this smaller remaining gap.

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
- `tests/test_unichrom_native.py`: 39 configurations checked against the
  original 3-D formulation, including singleton scan/charge dimensions,
  zero m/z width, narrow chromatography and masked mass outputs.
- Existing UniChrom routing and UCCD binary tests passed. Packaged binaries
  were not replaced; select the actual newly built executable for future work.
- The post-adjoint correction now checks the allowed mask, finiteness and
  nonnegative result while multiplying, removing the separate full-cube cleanup
  pass. The native regression, UniChrom routing test and UCCD binary test pass
  after this change.
- Point smoothing now allocates `scan_count * numz` scratch sums, runs scans in
  one OpenMP region, and uses a serial charge-block loop per scan. It writes
  directly into the existing scratch cube and swaps cube pointers, avoiding
  both nested parallel regions and the former full-cube input copy. Other
  callers retain the original parallel wrapper and in-place contract.
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
count was not separately recorded. This synthetic baseline did not include a
real-data or peak-memory measurement; the later SEC comparison above does.

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

The temporary comparison harness, profiler result and raw timing files were not
added to the repository; their results are summarized above. Make representative
benchmarks reproducible before building on these numbers. Extend the existing
regression source at `tests/test_unichrom_native.py` when changing this path.

## What the source establishes

| Location | Observation and implication |
| --- | --- |
| `unidec/src/UniDec.c`, `main` | Meta mode routes positive `dtsig` to `run_chromatogram`; zero goes to `run_metaunidec`. CDMS mode 2 routes separately to `run_unidec_UCCD`. |
| `unidec/src/UniChrom_Main.c`, `make_kernel_fft_UniChrom` | Builds the equivalent charge-independent 2-D kernel directly using `PeakDist` and the periodic scan peak formula. |
| Same file, `run_chromatogram` | Each iteration does four Q-cell FFT executions, plus dense charge projection and latent updates. There is no cube-sized ratio broadcast. |
| Same file, FFT initialization | Plans and kernel spectrum are already reused. Real transforms already exist. Recommending either as a new optimization would miss the actual opportunity. Plans use `FFTW_ESTIMATE`; the execution path now applies the available MKL/FFTW thread configuration around the reduced plans. |
| Same file, allocation/update loops | `blur`, `scratch`, `oldblur` still scale with the cube; FFT buffers scale with the measured grid. Pre-forward cleanup remains a full cube pass; post-adjoint correction now fuses cleanup into multiplication. Copies, smoothing and reductions still traverse the cube. |
| Same file, output reconvolution | For rawflag 0/2, charge planes are copied into contiguous temporary batches of at most four and transformed with one reused `plan_many` pair, then scattered back and masked. Allocation/plan failure retains the per-charge fallback. |
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
   iteration FFT executions and another `2*ceil(Z/min(Z,4))` batched output
   FFT executions for rawflag 0/2. The earlier per-charge path required 2Z;
   counts alone are not a timing
   breakdown. The retained path uses bounded contiguous packing; compare it
   against separable direct passes if this phase remains significant. Keep the
   per-charge fallback for allocation/plan failure. Do not allocate a full
   cube FFT again by default or parallelize the fallback charge loop using its
   shared mutable FFT workspace.
5. Three exact traffic/parallel candidates are complete: post-adjoint correction
   and cleanup are fused, point smoothing uses one scan-level parallel region
   with private sums, and that smoother now writes to the alternate cube instead
   of copying its input first. Reprofile before attempting another pass fusion.
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

First construct full-length axis arrays that reproduce the current UniChrom
sampled kernel, including origin/end image contributions. Preserve periodic
m/z convolution and the configured mirrored or zero scan extension, crop,
folded adjoint and sensitivity. Do not revert to the historical unpadded scan
operator or substitute a minimum-distance Gaussian. The split
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
- Positive-floor smoothing already caches one sanitized log per cell, then
  gathers immutable scratch values while writing results into data. UniChrom
  uses tiles for large multi-scan grids. Preserve the barrier between log
  generation and neighbor gathers; log caching is completed work.
- `softargmax` has its own OpenMP loop inside HDF5's scan-parallel call; avoid
  introducing another nested level there. Point smoothing's UniChrom path is
  scan-parallel with private sums and swaps its two cube buffers after writing
  the result; keep its other callers unchanged and never reintroduce shared
  scratch counters under an outer pragma.
- Inactive z/m regularizer tables, harmonic charge tables and point-smoothing
  sums are now allocated only when used. Mass-bin lower indexes and interpolation
  weights are precomputed after mass-axis
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
oracle includes charge and point smoothing; the broader 55-case regularizer
comparisons were temporary old/new runs, so add durable coverage for other
regularizers actually changed next. The 39 oracle cases are not comprehensive coverage of malformed inputs,
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

## Nonlinear processed-axis mode

Implemented after the FFT optimization work was closed. `UClineardecon=1`
remains the default and selects the common-grid FFT solver. Setting it to `0`
keeps each scan's regenerated `processed_data` axis and runs the coupled
Richardson-Lucy update with an explicit sparse direct forward/adjoint operator.
The GUI exposes this as **Linearize Before Coupled Deconvolution**, checked by
default. Missing HDF5 attributes also select the linear path; other values are
rejected.

The direct operator combines the existing peak-shape function in m/z with the
existing mirrored or zero-padded chromatographic response. It uses discrete
processed points, matching nonlinear UniDec's weighting convention. Charge,
mass, point and suppression regularizers reuse the established UCCD helpers or
their existing indexing rules. Empty scans keep their positions in the scan
kernel and produce zero output rows. Only completed results are interpolated to
common m/z and mass axes for the established HDF5 and GUI import contract.

`UCtype=0` keeps the default interpretation of `dtsig` in scans. `UCtype=1`
interprets it in retention-time units, forces the nonlinear direct solver and
builds the temporal response from the strictly increasing `retention_time`
attribute stored on each spectrum group. The GUI disables linearization while
time units are selected. Files without `UCtype` retain the scan-based default.

Full 50-iteration comparisons on the two repository examples met the accepted
0.98 normalized cosine threshold after interpolating results to the linear
reference axes:

| Example | m/z grid | mass grid | Dominant mass bin |
| --- | ---: | ---: | --- |
| SEC Native Herceptin | 0.99933 | 0.99975 | 148222 Da in both modes |
| SEC Native Bispecific Special | 0.99970 | 0.99972 | 193550/195900 Da in both modes |

These values use `dtsig=1` and Reconvolved/Profile output. Direct-to-linear
total intensity ratios were 0.99693/0.99690 for Herceptin and 0.99679/0.99080
(m/z/mass) for Bispecific. Against `dtsig=0`, the linear/direct m/z totals were
0.99941/0.99634 for Herceptin and 0.99948/0.99627 for Bispecific. The configured
Smart mass transform gave linear/direct mass totals of 0.99260/0.98952 and
0.98958/0.98048 relative to `dtsig=0`, respectively. Regression coverage
also checks unequal axes, output shape and finiteness, regularizers, a zero-width
m/z response, an empty middle scan, singleton charge, zero padding, HDF5 config
round trips, default and invalid dispatch, and the UniChrom GUI control.

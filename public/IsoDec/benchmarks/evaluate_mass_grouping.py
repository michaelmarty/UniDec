"""Compare mass grouping policies against fragment masses and observed centroids.

Run from the IsoDec package directory:
    python -m benchmarks.evaluate_mass_grouping spectrum.dat seq.fasta
The spectrum must already be centroided. No input files are modified.
"""

import argparse
import copy
import json
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import numpy as np

from isodec.brute_force_seq_match import brute_force_pep_match
from isodec.config import IsoDecConfig
from isodec.match import MatchedCollection, MatchedMass


def error_summary(errors, tolerance):
    errors = np.abs(np.asarray(errors, dtype=float))
    if not len(errors):
        return {"count": 0, "within_tolerance": 0, "median_abs_ppm": None,
                "p95_abs_ppm": None, "max_abs_ppm": None}
    return {
        "count": len(errors),
        "within_tolerance": int(np.sum(errors <= tolerance)),
        "median_abs_ppm": float(np.median(errors)),
        "p95_abs_ppm": float(np.percentile(errors, 95)),
        "max_abs_ppm": float(np.max(errors)),
    }


POLICIES = ("original", "mass_order", "matched_intensity", "apex_intensity", "exact_mass")


def group_peaks(peaks, config, policy):
    if policy == "mass_order":
        peaks.sort(key=lambda p: p.monoiso)
    elif policy == "matched_intensity":
        peaks.sort(key=lambda p: -p.matchedintensity)
    elif policy == "apex_intensity":
        peaks.sort(key=lambda p: -p.peakint)
    collection = MatchedCollection()
    if policy == "exact_mass":
        groups = {}
        for peak in peaks:
            if peak.monoiso in groups:
                groups[peak.monoiso].merge_in_pk(peak)
            else:
                groups[peak.monoiso] = MatchedMass(peak)
        collection.masses = list(groups.values())
    else:
        for peak in peaks:
            collection.add_pk_to_masses(peak, config)
    return collection


def benchmark(source, config, repeats):
    """Time warm grouping including sorting, excluding hit copying and scoring."""
    for policy in POLICIES:
        group_peaks(copy.deepcopy(source), config, policy)
    samples = {policy: [] for policy in POLICIES}
    rng = np.random.default_rng(0)
    for _ in range(repeats):
        for policy in rng.permutation(POLICIES):
            peaks = copy.deepcopy(source)
            start = perf_counter()
            result = group_peaks(peaks, config, policy)
            samples[policy].append((perf_counter() - start) * 1000)
            del result
    return {policy: {"repeats": repeats, "median_ms": float(np.median(times)),
                     "min_ms": min(times), "max_ms": max(times)}
            for policy, times in samples.items()}


def evaluate(source, theoretical, config, policy):
    # Keep replays independent, including when evaluating older merging code.
    collection = group_peaks(copy.deepcopy(source), config, policy)

    representative, candidates, fitted, observed = [], [], [], []
    groups_with_candidate = 0
    mixed_groups = 0
    for group in collection.masses:
        masses = np.unique([theoretical[p.sequence_match] for p in group.clusters])
        representative.append(float(np.max(np.abs(group.monoiso - masses) / masses * 1e6)))
        candidate_errors = [float(np.min(np.abs(m - masses) / masses * 1e6))
                            for m in group.monoisos]
        candidates.extend(candidate_errors)
        groups_with_candidate += any(e <= config.matchtol for e in candidate_errors)
        if len(masses) != 1:
            mixed_groups += 1
            continue
        mass = masses[0]
        # Recover the implicit monoisotopic origin of the fitted mass axis.
        seed = group.clusters[0]
        first_isotope = round((seed.massdist[0, 0] - mass) / config.mass_diff_c)
        fitted_origin = group.massdist[0, 0] - first_isotope * config.mass_diff_c
        fitted.append((fitted_origin - mass) / mass * 1e6)
        # Estimate mass from measured m/z, charge, and assigned isotope position.
        # This uses the existing assignments; it is not an independent ID test.
        residuals, weights = [], []
        for peak in group.clusters:
            positions = np.asarray(peak.isomatches, dtype=int)
            isotope_offsets = peak.massdist[positions, 0] - mass
            measured = (peak.matchedcentroids[:, 0] - config.adductmass) * peak.z
            residuals.extend(measured - isotope_offsets - mass)
            weights.extend(peak.matchedcentroids[:, 1])
        observed.append(np.average(residuals, weights=weights) / mass * 1e6)

    return {
        "groups": len(collection.masses),
        "mixed_theoretical_masses": mixed_groups,
        "representative": error_summary(representative, config.matchtol),
        "candidate_entries": error_summary(candidates, config.matchtol),
        "groups_with_valid_candidate": groups_with_candidate,
        "fitted_axis_origin": error_summary(fitted, config.matchtol),
        "observed_mass_estimate": error_summary(observed, config.matchtol),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("sequence", type=Path)
    parser.add_argument("--fragmentation", default="ETD")
    parser.add_argument("--timing-repeats", type=int, default=0,
                        help="Also measure warm grouping; excludes hit generation and copying")
    args = parser.parse_args()
    if args.timing_repeats < 0:
        parser.error("--timing-repeats must be nonnegative")
    sequence = "".join(line.strip() for line in args.sequence.read_text().splitlines()
                       if line.strip() and not line.startswith(">"))
    spectrum = np.loadtxt(args.spectrum)
    config = IsoDecConfig()
    config.activescan = 1
    config.activescanrt = 0
    config.activescanorder = 2
    # Obtain pristine accepted hits without the normal per-hit grouping.
    with patch("isodec.brute_force_seq_match._group_fragment_peaks", lambda pks, *a, **k: pks):
        hits = brute_force_pep_match(sequence, spectrum, centroided=True, config=config,
                                    fragmentation_type=args.fragmentation)
    results = {policy: evaluate(hits.peaks, hits.fragment_theoretical, config, policy)
               for policy in POLICIES}
    output = {"hits": len(hits.peaks), "tolerance_ppm": config.matchtol, "results": results}
    if args.timing_repeats:
        output["warm_grouping_timings"] = benchmark(hits.peaks, config, args.timing_repeats)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

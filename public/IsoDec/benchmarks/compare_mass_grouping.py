"""Compare grouping policies on accepted fragment hits and their isotope spectra.

Run from the IsoDec package directory:
    python -m benchmarks.compare_mass_grouping spectrum.dat seq.fasta
Use --profile for a profile-mode input that IsoDec must centroid first.
"""

import argparse
import copy
import json
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import numpy as np

from benchmarks.evaluate_mass_grouping import POLICIES, group_peaks
from isodec.brute_force_seq_match import brute_force_pep_match
from isodec.config import IsoDecConfig
from isodec.datatools import get_all_centroids
from isodec.match import MatchedCollection


def grouped(source, config, policy, rng=None):
    peaks = copy.deepcopy(source)
    index = {id(peak): i for i, peak in enumerate(peaks)}
    if rng is not None:
        rng.shuffle(peaks)
    collection = group_peaks(peaks, config, policy)
    members = [frozenset(index[id(peak)] for peak in mass.clusters)
               for mass in collection.masses]
    return collection, members


def projected(collection):
    """Align projected isotope intensities by fragment label and isotope number."""
    by_label = defaultdict(lambda: defaultdict(float))
    for group in collection.masses:
        label = group.clusters[0].sequence_match
        source = group.clusters[0].massdist[:, 0]
        first = round((source[0] - group.clusters[0].monoiso) / 1.0033)
        for i, row in enumerate(group.massdist):
            by_label[label][first + i] += float(row[1])
    return by_label


def similarity(reference, candidate):
    labels = sorted(set(reference) | set(candidate))
    total_a, total_b = [], []
    per_label = []
    differing = 0
    for label in labels:
        isotope_ids = sorted(set(reference[label]) | set(candidate[label]))
        a = np.asarray([reference[label].get(i, 0.0) for i in isotope_ids])
        b = np.asarray([candidate[label].get(i, 0.0) for i in isotope_ids])
        total_a.extend(a)
        total_b.extend(b)
        norm = np.linalg.norm(a) * np.linalg.norm(b)
        score = float(np.dot(a, b) / norm) if norm else float(not np.any(a) and not np.any(b))
        per_label.append((score, label))
        differing += not np.allclose(a, b, rtol=1e-6, atol=1e-10)
    a, b = np.asarray(total_a), np.asarray(total_b)
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return {
        "global_cosine": float(np.dot(a, b) / norm) if norm else None,
        "normalized_l1": float(np.sum(np.abs(a - b)) / np.sum(np.abs(a))) if np.sum(a) else None,
        "labels_with_changed_projection": differing,
        "minimum_label_cosine": min(per_label) if per_label else None,
    }


def comparison(reference_groups, other_groups, total_hits):
    reference = set(reference_groups)
    other = set(other_groups)
    reference_peers = {i: members for members in reference_groups for i in members}
    other_peers = {i: members for members in other_groups for i in members}
    return {
        "groups": len(other_groups),
        "identical_groups": len(reference & other),
        "hits_with_changed_group_membership": sum(
            reference_peers[i] != other_peers[i] for i in range(total_hits)),
    }


def group_field_differences(reference, reference_members, other, other_members):
    left = dict(zip(reference_members, reference.masses))
    right = dict(zip(other_members, other.masses))
    changed_charge_order = changed_charge_set = changed_total = 0
    axis_ppm = []
    for members in left.keys() & right.keys():
        a, b = left[members], right[members]
        changed_charge_order += not np.array_equal(a.zs, b.zs)
        changed_charge_set += set(a.zs) != set(b.zs)
        changed_total += not np.isclose(a.totalintensity, b.totalintensity,
                                        rtol=1e-10, atol=1e-10)
        if a.massdist.shape == b.massdist.shape:
            axis_ppm.append(float(np.max(np.abs(a.massdist[:, 0] - b.massdist[:, 0])
                                         / a.massdist[:, 0] * 1e6)))
    return {"shared_groups_with_changed_charge_order": changed_charge_order,
            "shared_groups_with_changed_charge_set": changed_charge_set,
            "shared_groups_with_changed_total_intensity": changed_total,
            "maximum_shared_group_axis_difference_ppm": max(axis_ppm, default=0.0)}


def run(spectrum, sequence, fragmentation, centroided, random_trials):
    config = IsoDecConfig()
    config.activescan = 1
    config.activescanrt = 0
    config.activescanorder = 2
    if not centroided:
        spectrum = get_all_centroids(spectrum, window=config.peakwindow,
                                     threshold=config.peakthresh * 0.1)
    spectrum = spectrum[np.argsort(spectrum[:, 0])]
    with patch("isodec.brute_force_seq_match._group_fragment_peaks", lambda pks, *a, **k: pks):
        hits = brute_force_pep_match(sequence, spectrum, centroided=True,
                                    config=config, fragmentation_type=fragmentation)
    source = hits.peaks
    baseline, baseline_members = grouped(source, config, "original")
    baseline_projected = projected(baseline)
    policies = {}
    for policy in POLICIES[1:]:
        groups, members = grouped(source, config, policy)
        policies[policy] = {**comparison(baseline_members, members, len(source)),
                            **group_field_differences(baseline, baseline_members, groups, members),
                            **similarity(baseline_projected, projected(groups))}
    random_results = []
    rng = np.random.default_rng(20260928)
    for _ in range(random_trials):
        groups, members = grouped(source, config, "original", rng)
        random_results.append({**comparison(baseline_members, members, len(source)),
                               **similarity(baseline_projected, projected(groups))})
    result = {"centroids": len(spectrum), "hits": len(source),
              "original_groups": len(baseline_members), "policies": policies,
              "random_replays": random_results}
    if len(source) <= 20:
        result["accepted_fragment_charges"] = [
            [peak.sequence_match, peak.z] for peak in source]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("sequence", type=Path)
    parser.add_argument("--fragmentation", default="ETD")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--random-trials", type=int, default=20)
    args = parser.parse_args()
    if args.random_trials < 0:
        parser.error("--random-trials must be nonnegative")
    sequence = "".join(line.strip() for line in args.sequence.read_text().splitlines()
                       if line.strip() and not line.startswith(">"))
    result = run(np.loadtxt(args.spectrum), sequence, args.fragmentation,
                 not args.profile, args.random_trials)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

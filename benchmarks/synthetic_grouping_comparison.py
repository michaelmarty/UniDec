"""Compare grouping policies on small ETD mixtures with known fragment sources.

Run from the IsoDec package directory:
    python -m benchmarks.synthetic_grouping_comparison
"""

import json

import isogen
import numpy as np

from benchmarks.compare_mass_grouping import grouped, run
from isodec.config import IsoDecConfig
from isodec.match import MatchedPeak


def spectrum_for_case(batch, labels, case):
    rng = np.random.default_rng(183)
    rows = []
    for label in labels:
        fragment = batch.labels.index(label)
        mass = batch.masses[fragment]
        isotope_intensity = batch.intensities[fragment]
        for charge in (2, 3):
            for isotope in range(5):
                if case == "missing_isotope" and isotope == 3:
                    continue
                if isotope_intensity[isotope] < 0.005:
                    continue
                mz = (mass + isotope * 1.0033) / charge + 1.007276467
                if case in ("noise", "nearby_shoulders"):
                    mz *= 1 + rng.normal(0, 0.35e-6)
                rows.append((mz, float(isotope_intensity[isotope] *
                                       (100 if label == "c4" else 30))))
    if case == "nearby_shoulders":
        rows.extend((mz + 0.015, intensity * 0.04)
                    for mz, intensity in rows[:4])
    if case in ("noise", "nearby_shoulders"):
        low, high = min(mz for mz, _ in rows), max(mz for mz, _ in rows)
        rows.extend(zip(rng.uniform(low, high, 80), rng.uniform(0.005, 0.03, 80)))
    return np.asarray(sorted(rows), dtype=float)


def close_mass_decoy(second_mass):
    """Two known sources with broad isotope envelopes."""
    config = IsoDecConfig()
    peaks = []
    for label, mass in (("A", 1000.0), ("B", second_mass)):
        for charge in (2, 3):
            dist = np.column_stack((mass + np.arange(5) * config.mass_diff_c,
                                    np.array([0.6, 0.8, 1.0, 0.8, 0.6]) * 100))
            isodist = dist.copy()
            isodist[:, 0] = isodist[:, 0] / charge + config.adductmass
            peak = MatchedPeak(charge, isodist[2, 0], centroids=isodist,
                               isodist=isodist, matchedindexes=list(range(5)),
                               isomatches=list(range(5)), config=config)
            peak.monoiso = mass
            peak.monoisos = [mass]
            peak.massdist = dist
            peak.peakint = float(isodist[:, 1].max())
            peak.avgmass = float(np.average(dist[:, 0], weights=dist[:, 1]))
            peak.sequence_match = label
            peaks.append(peak)
    return {policy: {"groups": len(groups.masses),
                     "mixed_identity_groups": sum(
                         len({peaks[i].sequence_match for i in members}) > 1
                         for members in member_sets)}
            for policy in ("original", "mass_order", "matched_intensity",
                           "apex_intensity", "exact_mass")
            for groups, member_sets in [grouped(peaks, config, policy)]}


def main():
    sequence = "PEPTIDEK"
    labels = ("c4", "z'4")
    batch = isogen.calc_pep_fragment_isodists(sequence, fragmentation_type="ETD")
    results = {}
    for case in ("clean", "missing_isotope", "noise", "nearby_shoulders"):
        result = run(spectrum_for_case(batch, labels, case), sequence, "ETD", True, 5)
        assert result["hits"] == 4 and result["original_groups"] == 2, case
        assert {tuple(item) for item in result["accepted_fragment_charges"]} == {
            (label, charge) for label in labels for charge in (2, 3)}, case
        assert all(item["hits_with_changed_group_membership"] == 0
                   for item in result["policies"].values()), case
        results[case] = {"centroids": result["centroids"],
                         "hits": result["hits"], "groups": result["original_groups"],
                         "minimum_cosine": min(item["global_cosine"]
                                               for item in result["policies"].values())}
    decoy = close_mass_decoy(1001.0033)
    isobaric = close_mass_decoy(1000.0)
    assert all(decoy[policy]["mixed_identity_groups"] == 1
               for policy in ("original", "mass_order", "matched_intensity",
                              "apex_intensity"))
    assert decoy["exact_mass"]["mixed_identity_groups"] == 0
    assert isobaric["exact_mass"]["mixed_identity_groups"] == 1
    print(json.dumps({"fragment_mixtures": results, "isotope_shift_decoy": decoy,
                      "isobaric_decoy": isobaric}, indent=2))


if __name__ == "__main__":
    main()

from pathlib import Path
from types import SimpleNamespace

import isogen
import numpy as np
import pytest

from isodec import match_fragments
from isodec.match import MatchedCollection, MatchedPeak

CA_SEQUENCE = (
    "S[Acetylation]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEAT"
    "SRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLV"
    "HWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNV"
    "LDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNR"
    "QVRGFPK"
)
CA_MASSES = tuple(
    float(value)
    for value in (
        Path(__file__).parents[1] / "extern/IsoGen/tests/data/ca_etd_masses.txt"
    ).read_text(encoding="utf-8-sig").split()
)


def collection(*masses):
    peaks = [SimpleNamespace(monoiso=mass, monoisos=[mass]) for mass in masses]
    return SimpleNamespace(peaks=peaks)


def test_carbonic_anhydrase_fragment_matching():
    peaks = []
    for mass in CA_MASSES:
        peak = MatchedPeak(1, mass)
        peak.monoiso = mass
        peak.monoisos = [mass]
        peaks.append(peak)
    pks = MatchedCollection().add_peaks(peaks)

    match_fragments(
        pks,
        CA_SEQUENCE,
        fragmentation_type="ETD",
        ppm_tolerance=20,
    )
    print("Sequence coverage:", pks.sequence_coverage)
    assert pks.sequence_coverage == pytest.approx(155 / 258)
    assert pks.fragment_match_percent == pytest.approx(166 / 605 * 100)
    assert any(
        peak.sequence_match.startswith("z'")
        for peak in pks.peaks
        if peak.sequence_match
    )


def test_matches_default_by_fragments_and_builds_summary_table():
    sequence = "PEPTIDE"
    fragments = isogen.calc_pep_fragments(sequence)
    pks = collection(fragments["b2"], fragments["y3"], fragments["b4"] + 1)

    result = match_fragments(pks, sequence)

    assert result is pks
    assert pks.peaks[0].sequence_match == "b2"
    assert pks.peaks[1].sequence_match == "y3"
    assert pks.peaks[2].sequence_match is None
    assert pks.fragment_match_percent == pytest.approx(200 / 3)
    assert pks.sequence_coverage == pytest.approx(2 / (len(sequence) - 1))
    assert pks.fragment_matches.loc[2, "b_match"] == pytest.approx(fragments["b2"])
    assert pks.fragment_matches.loc[4, "y_match"] == pytest.approx(fragments["y3"])
    assert pks.fragment_matches["a_mass"].isna().all()
    assert pks.fragment_matches.loc[2, "match_count"] == 1


def test_uses_all_peak_monoisotopic_candidates():
    sequence = "PEPTIDE"
    b3 = isogen.calc_pep_fragments(sequence)["b3"]
    peak = SimpleNamespace(monoiso=b3 + 10, monoisos=[b3 + 10, b3])
    pks = SimpleNamespace(peaks=[peak])

    match_fragments(pks, sequence)

    assert peak.sequence_match == "b3"
    assert pks.fragment_matches.loc[3, "b_match"] == pytest.approx(b3)


def test_same_fragment_can_annotate_peaks_from_multiple_scans():
    sequence = "PEPTIDE"
    b3 = isogen.calc_pep_fragments(sequence)["b3"]
    pks = collection(b3 + 0.0001, b3)

    match_fragments(pks, sequence)

    assert [peak.sequence_match for peak in pks.peaks] == ["b3", "b3"]
    assert pks.fragment_match_percent == 100
    assert pks.fragment_matches.loc[3, "b_match"] == pytest.approx(b3)


def test_duplicate_assignments_are_optional(monkeypatch):
    monkeypatch.setattr(
        isogen,
        "calc_pep_fragments",
        lambda *args, **kwargs: {"b1": 100.0, "y1": 100.0001},
    )

    pks = collection(100.00004)
    match_fragments(pks, "PEP", ppm_tolerance=5)
    assert pks.peaks[0].sequence_match == "b1"
    assert pks.fragment_matches.loc[1, "match_count"] == 1

    match_fragments(
        pks,
        "PEP",
        ppm_tolerance=5,
        allow_duplicate_assignments=True,
    )
    assert pks.peaks[0].sequence_match == ["b1", "y1"]
    assert pks.fragment_matches["match_count"].sum() == 2


def test_empty_collection_and_invalid_tolerance():
    pks = collection()
    match_fragments(pks, "PEPTIDE")
    assert pks.fragment_match_percent == 0
    assert pks.sequence_coverage == 0

    with pytest.raises(ValueError, match="non-negative"):
        match_fragments(pks, "PEPTIDE", ppm_tolerance=-1)

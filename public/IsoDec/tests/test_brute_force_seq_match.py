import isogen
import numpy as np
import pytest

from isodec import IsoDecRuntime, brute_force_pep_match


def _fragment_spectrum(sequence="PEPTIDE", label="b6", charge=2):
    fragments = isogen.calc_pep_fragment_isodists(sequence)
    index = fragments.labels.index(label)
    mass = fragments.masses[index]
    values = fragments.intensities[index]
    positions = np.flatnonzero(values > values.max() * 0.01)
    spectrum = np.column_stack((mass / charge + positions * 1.0033 / charge + 1.007276467,
                                values[positions] * 100))
    return mass, spectrum


def test_finds_fragment_and_charge_from_isotope_envelope():
    mass, spectrum = _fragment_spectrum()
    pks = brute_force_pep_match("PEPTIDE", spectrum, centroided=True)

    matches = [peak for peak in pks if peak.sequence_match == "b6" and peak.z == 2]
    assert len(matches) == 1
    peak = matches[0]
    assert peak.monoiso == pytest.approx(mass)
    assert peak.matchedintensity == pytest.approx(spectrum[:, 1].sum())
    assert pks.masses
    assert peak.peakint == pytest.approx(spectrum[:, 1].max())
    np.testing.assert_allclose(peak.isodist[:, 1], spectrum[:, 1], rtol=1e-6)
    np.testing.assert_allclose(peak.massdist[:, 1], spectrum[:, 1], rtol=1e-6)
    np.testing.assert_allclose(pks.masses[0].isodists[:, 1], spectrum[:, 1], rtol=1e-6)
    assert peak.match_score >= peak.config.css_thresh
    assert len(peak.matchedindexes) >= peak.config.minpeaks


def test_runtime_and_charge_limit():
    _, spectrum = _fragment_spectrum()
    runtime = IsoDecRuntime()
    pks = runtime.brute_force_pep_match("PEPTIDE", spectrum, centroided=True)
    assert pks is runtime.pks
    assert any(peak.sequence_match == "b6" and peak.z == 2 for peak in pks)

    limited = brute_force_pep_match("PEPTIDE", spectrum, centroided=True, max_charge=1)
    assert not any(peak.sequence_match == "b6" for peak in limited)


def test_rejects_mismatched_pattern_and_invalid_spectrum():
    _, spectrum = _fragment_spectrum()
    spectrum[:, 1] = [1, 1, 100, 1]
    pks = brute_force_pep_match("PEPTIDE", spectrum, centroided=True)
    assert not any(peak.sequence_match == "b6" and peak.z == 2 for peak in pks)

    with pytest.raises(ValueError, match="m/z and intensity"):
        brute_force_pep_match("PEPTIDE", np.ones((3, 3)))
    with pytest.raises(ValueError, match="positive integer"):
        brute_force_pep_match("PEPTIDE", spectrum, max_charge=0)
    with pytest.raises(ValueError, match="monoisotopic"):
        brute_force_pep_match("PEPTIDE", spectrum, monoisotopic=False)


def test_native_and_python_sequence_match_agree_for_modified_fragments():
    sequence = "S[Acetylation]HHS"
    _, spectrum = _fragment_spectrum(sequence, label="b3", charge=2)
    python = brute_force_pep_match(sequence, spectrum, centroided=True, native=False)
    native = brute_force_pep_match(sequence, spectrum, centroided=True, native=True)
    assert [(p.sequence_match, p.z) for p in native] == [
        (p.sequence_match, p.z) for p in python
    ]
    for first, second in zip(native, python):
        assert first.match_score == pytest.approx(second.match_score, abs=1e-9)
        np.testing.assert_allclose(first.isodist, second.isodist, rtol=1e-9)
        np.testing.assert_allclose(first.massdist, second.massdist, rtol=1e-9)
        assert first.matchedintensity == pytest.approx(second.matchedintensity)


def test_native_and_python_choose_same_centroids_with_nearby_peaks():
    _, spectrum = _fragment_spectrum()
    interferer = spectrum[1].copy()
    interferer[0] += interferer[0] * 2e-6
    interferer[1] *= 0.5
    spectrum = np.vstack((spectrum, interferer, [spectrum[-1, 0] + 0.2, 2.0]))
    native = brute_force_pep_match("PEPTIDE", spectrum, centroided=True)
    python = brute_force_pep_match("PEPTIDE", spectrum, centroided=True, native=False)
    assert [(p.sequence_match, p.z) for p in native] == [
        (p.sequence_match, p.z) for p in python
    ]
    for first, second in zip(native, python):
        assert first.matchedindexes == second.matchedindexes
        assert first.isomatches == second.isomatches

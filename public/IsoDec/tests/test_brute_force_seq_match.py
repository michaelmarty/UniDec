import isogen
import numpy as np
import pytest

from isodec import IsoDecRuntime, brute_force_pep_match
from isodec.isotope import calc_isotope_dist_dual


def _fragment_spectrum(sequence="PEPTIDE", label="b6", charge=2):
    mass = isogen.calc_pep_fragments(sequence)[label]
    _, distribution = calc_isotope_dist_dual(mass)
    spectrum = distribution.copy()
    spectrum[:, 0] = spectrum[:, 0] / charge + 1.007276467
    spectrum[:, 1] *= 100
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

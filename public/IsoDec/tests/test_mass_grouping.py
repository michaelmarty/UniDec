import numpy as np
import pytest

from isodec.config import IsoDecConfig
from isodec.match import MatchedMass, MatchedPeak, merge_massdist


def _peak(mass, intensity):
    config = IsoDecConfig()
    massdist = np.column_stack((mass + np.arange(3) * config.mass_diff_c,
                                intensity * np.array([0.6, 0.3, 0.1])))
    spectrum = massdist.copy()
    spectrum[:, 0] = spectrum[:, 0] / 2 + config.adductmass
    peak = MatchedPeak(2, spectrum[0, 0], centroids=spectrum, isodist=spectrum,
                       matchedindexes=[0, 1, 2], isomatches=[0, 1, 2], config=config)
    peak.monoiso = mass
    peak.monoisos = [mass]
    peak.massdist = massdist
    peak.peakint = spectrum[0, 1]
    peak.avgmass = np.average(massdist[:, 0], weights=massdist[:, 1])
    return peak


def test_mass_group_uses_matched_centroids_and_keeps_local_window():
    peak = _peak(1000.0, 2.0)
    original = peak.centroids.copy()
    distractor = np.array([[1000.004 / peak.z + peak.config.adductmass, 100.0]])
    peak.centroids = np.vstack((peak.centroids, distractor))
    group = MatchedMass(peak)
    np.testing.assert_allclose(group.decon_centroids[:, 0],
                               original[:, 0] * peak.z - peak.config.adductmass * peak.z)
    assert len(group.decon_centroids) == 3
    assert len(peak.centroids) == 4


def test_mass_group_recovers_matched_centroids_without_match_indexes():
    peak = _peak(1000.0, 2.0)
    peak.matchedcentroids = None
    peak.centroids = np.vstack((peak.centroids,
                                [[1000.04 / peak.z + peak.config.adductmass, 100.0]]))
    group = MatchedMass(peak)
    assert len(group.decon_centroids) == 3
    assert len(peak.centroids) == 4


def test_mass_group_uses_full_window_when_no_isotope_assignment_exists():
    peak = _peak(1000.0, 2.0)
    peak.matchedcentroids = None
    peak.isodist = None
    group = MatchedMass(peak)
    assert len(group.decon_centroids) == len(peak.centroids)


@pytest.mark.parametrize("offset", [0.0, 0.002])
def test_candidate_mass_is_weighted_average_without_mutating_peak(offset):
    first = _peak(1000.0, 2.0)
    second = _peak(1000.0 + offset, 6.0)
    group = MatchedMass(first)
    group.merge_in_pk(second)
    expected = 1000.0 + 0.75 * offset
    np.testing.assert_allclose(group.monoisos, [expected], rtol=0, atol=1e-10)
    assert group.monoiso == pytest.approx(expected)
    assert first.monoisos == [1000.0]
    assert second.monoisos == [1000.0 + offset]
    assert group.totalintensity == pytest.approx(8.0)


@pytest.mark.parametrize("shift", [-0.003, 0.0, 0.003])
def test_mass_axis_moves_toward_observed_centroid_without_repeated_drift(shift):
    theory = np.array([[1000.0, 10.0], [1001.0033, 5.0], [1002.0066, 2.0]])
    observed = theory.copy()
    observed[:, 0] += shift
    observed[:, 1] *= 2
    fitted = merge_massdist(theory.copy(), observed, 5.0)
    np.testing.assert_allclose(fitted, observed, rtol=0, atol=1e-10)
    for _ in range(3):
        fitted = merge_massdist(fitted, observed, 5.0)
    np.testing.assert_allclose(fitted, observed, rtol=0, atol=1e-10)


@pytest.mark.parametrize("shift", [-0.01, 0.01])
def test_mass_axis_does_not_follow_out_of_tolerance_centroid(shift):
    theory = np.array([[1000.0, 10.0]])
    observed = np.array([[1000.0 + shift, 20.0]])
    np.testing.assert_array_equal(merge_massdist(theory.copy(), observed, 5.0), theory)

"""Parity checks for the optional native mass-grouping ABI."""

import ctypes
import gzip
import json
import os
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import isogen
import numpy as np
import pytest

from isodec.brute_force_seq_match import brute_force_pep_match
from isodec.c_interface import IsoDecWrapper, MassGroupResultStruct
from isodec.config import IsoDecConfig
from isodec.match import MatchedCollection, MatchedPeak, find_matches


FROZEN_SNAPSHOT_RTOL = 1e-6
FROZEN_SNAPSHOT_ATOL = 1e-7


def _peak(mass, intensity):
    config = IsoDecConfig()
    massdist = np.column_stack((mass + np.arange(3) * config.mass_diff_c,
                                intensity * np.array([.6, .3, .1])))
    envelope = massdist.copy()
    envelope[:, 0] = envelope[:, 0] / 2 + config.adductmass
    peak = MatchedPeak(2, envelope[0, 0], centroids=envelope, isodist=envelope,
                       matchedindexes=[0, 1, 2], isomatches=[0, 1, 2], config=config)
    peak.monoiso, peak.monoisos = mass, [mass]
    peak.massdist, peak.peakint = massdist, envelope[0, 1]
    peak.avgmass = np.average(massdist[:, 0], weights=massdist[:, 1])
    return peak


def _native_wrapper():
    path = os.environ.get("ISODEC_TEST_NATIVE_LIBRARY")
    wrapper = IsoDecWrapper(path) if path else IsoDecWrapper()
    if wrapper._mass_group is None or wrapper._group_centroid_match is None:
        pytest.skip("native mass-grouping ABI is unavailable")
    return wrapper


def _reference(peaks, config, order):
    peaks = deepcopy(peaks)
    collection = MatchedCollection()
    collection.peaks = peaks
    insertion = (sorted(peaks, key=lambda p: -p.matchedintensity)
                 if order == "matched_intensity" else peaks)
    for peak in insertion:
        collection.add_pk_to_masses(peak, config)
    return collection


@pytest.mark.parametrize("order,expected_groups", [
    ("original", 197), ("matched_intensity", 196),
])
def test_native_grouping_matches_python_on_supplied_spectrum(spectrum, order, expected_groups):
    wrapper = _native_wrapper()
    config = wrapper.config
    config.mass_group_order = order
    native = wrapper.process_spectrum(spectrum, config=config)
    fallback = _native_wrapper()
    fallback._mass_group = None
    fallback._group_centroid_match = None
    reference = fallback.process_spectrum(spectrum, config=deepcopy(config))
    assert len(native.peaks) == 277
    assert len(native.masses) == len(reference.masses) == expected_groups
    np.testing.assert_array_equal(native.monoisos, reference.monoisos)
    for actual, expected in zip(native.masses, reference.masses):
        assert [native.peaks.index(p) for p in actual.clusters] == [
            reference.peaks.index(p) for p in expected.clusters]
        for field in ("monoiso", "monoisos", "massdist", "decon_centroids",
                      "totalintensity", "apexintensity", "mzs", "zs",
                      "isodists", "mzints", "scans", "avgmass", "minrt", "maxrt",
                      "apexrt", "minscan", "maxscan", "apexscan"):
            np.testing.assert_allclose(getattr(actual, field), getattr(expected, field),
                                       rtol=1e-12, atol=1e-9, err_msg=field)
        assert actual.scan_intensities == expected.scan_intensities
        assert actual.totalpeaks == expected.totalpeaks
    assert native.to_df().equals(reference.to_df())
    np.testing.assert_array_equal(native.to_mass_spectrum(), reference.to_mass_spectrum())


def test_native_rematch_uses_exported_window(spectrum):
    wrapper = _native_wrapper()
    peaks = wrapper.process_spectrum(spectrum)
    for peak in peaks:
        matched, _ = find_matches(peak.centroids, peak.isodist, peak.config.matchtol)
        expected = peak.centroids[np.unique(matched)] if matched else peak.centroids
        np.testing.assert_array_equal(peak.matchedcentroids, expected)
    selected = next(p for p in peaks if p.z == 1 and abs(p.mz - 649.2872) < 0.001)
    np.testing.assert_allclose(selected.matchedcentroids[:, 0],
                               [649.2871704, 650.2905884], rtol=0, atol=1e-4)


def test_wrapper_ignores_stale_native_isotope_tail(monkeypatch):
    wrapper = _native_wrapper()
    def native_result(mz, intensity, count, model, output, settings, kind):
        peak = output[0]
        peak.z, peak.mz, peak.monoiso, peak.avgmass = 1, 1001., 1000., 1000.5
        peak.peakint, peak.startindex, peak.endindex = 10., 0, 1
        peak.monoisos[0] = 1000.
        peak.realisolength = 2
        for i in range(3):
            peak.isomz[i] = 1001. + i
            peak.isomass[i] = 1000. + i
            peak.isodist[i] = 10. if i < 2 else 1000.
        return 1
    monkeypatch.setattr(wrapper.c_lib, "process_spectrum", native_result)
    result = wrapper.process_spectrum(np.array([[1001., 10.], [1002., 10.]]))
    assert len(result.peaks[0].isodist) == 2
    assert result.peaks[0].matchedintensity == 20.


def test_native_grouping_preserves_fragment_match_output():
    wrapper = _native_wrapper()
    batch = isogen.calc_pep_fragment_isodists("PEPTIDE")
    index = batch.labels.index("b6")
    mass, values = batch.masses[index], batch.intensities[index]
    positions = np.flatnonzero(values > values.max() * 0.01)
    spectrum = np.column_stack((mass / 2 + positions * 1.0033 / 2 + 1.007276467,
                                values[positions] * 100))
    native = brute_force_pep_match("PEPTIDE", spectrum, centroided=True,
                                   native_wrapper=wrapper)
    reference = brute_force_pep_match("PEPTIDE", spectrum, centroided=True,
                                      native=False)
    assert [(p.sequence_match, p.z) for p in native] == [
        (p.sequence_match, p.z) for p in reference]
    assert len(native.masses) == len(reference.masses)
    for actual, expected in zip(native.masses, reference.masses):
        np.testing.assert_allclose(actual.monoiso, expected.monoiso, rtol=1e-9)
        np.testing.assert_allclose(actual.massdist, expected.massdist, rtol=1e-9)
        np.testing.assert_allclose(actual.decon_centroids, expected.decon_centroids,
                                   rtol=1e-9)


def test_native_rematch_prefers_strongest_nearby_centroid_and_skips_zero():
    wrapper = _native_wrapper()
    spectrum = np.array([[100.0, 10.0], [100.0001, 20.0],
                         [101.0, 5.0], [101.0001, 0.0], [102.0, 2.0]])
    isotope = np.array([[100.00005, 1.0], [101.0, 1.0], [101.0005, 1.0]])
    peak = SimpleNamespace(startindex=0, endindex=4, isodist=isotope)
    expected, _ = find_matches(spectrum, isotope, 2.0)
    matched = wrapper.match_group_centroids_batch(spectrum, [peak], 2.0)
    assert matched[0].tolist() == expected == [1, 2]
    with pytest.raises(ValueError, match="centroid matching failed"):
        wrapper.match_group_centroids_batch(spectrum, [peak], -1.0)
    groups, lookup = wrapper.group_mass_peaks_batch([], wrapper.config)
    assert groups == [] and len(lookup) == 0


def test_native_grouping_abi_handles_empty_input_and_invalid_order():
    wrapper = _native_wrapper()
    groups = ctypes.POINTER(MassGroupResultStruct)()
    count = ctypes.c_int(-1)
    ids = ctypes.POINTER(ctypes.c_int)()
    arguments = (None, 0, 5.0, 3, 1.0033, 0.7, 100)
    assert wrapper._mass_group(*arguments, 0, ctypes.byref(groups),
                               ctypes.byref(count), ctypes.byref(ids)) == 0
    assert count.value == 0 and not groups and not ids
    assert wrapper._mass_group(*arguments, 2, ctypes.byref(groups),
                               ctypes.byref(count), ctypes.byref(ids)) == -1


def _assert_groups(native, reference):
    np.testing.assert_array_equal(native.monoisos, reference.monoisos)
    assert len(native.masses) == len(reference.masses)
    for actual, expected in zip(native.masses, reference.masses):
        assert vars(actual).keys() == vars(expected).keys()
        for field, value in vars(expected).items():
            observed = getattr(actual, field)
            if field == "clusters":
                assert [native.peaks.index(p) for p in observed] == [
                    reference.peaks.index(p) for p in value]
            elif field == "scan_intensities":
                assert observed == value
            elif value is None:
                assert observed is None
            elif field in ("zs", "scans", "minscan", "maxscan", "apexscan", "totalpeaks"):
                np.testing.assert_array_equal(observed, value, err_msg=field)
            elif field in ("monoiso", "monoisos", "mzs", "avgmass"):
                np.testing.assert_allclose(observed, value, rtol=0, atol=1e-9, err_msg=field)
            elif field in ("massdist", "decon_centroids", "isodists"):
                np.testing.assert_allclose(observed[:, 0], value[:, 0], rtol=0, atol=1e-9,
                                           err_msg=field + " masses")
                np.testing.assert_allclose(observed[:, 1], value[:, 1], rtol=1e-12, atol=1e-12,
                                           err_msg=field + " intensities")
            else:
                np.testing.assert_allclose(observed, value, rtol=1e-12, atol=1e-12, err_msg=field)


def _native_collection(wrapper, peaks, order="original"):
    collection = MatchedCollection()
    collection.peaks = deepcopy(peaks)
    collection.masses, collection.monoisos = wrapper.group_mass_peaks_batch(
        collection.peaks, wrapper.config, order)
    return collection


def _snapshot_arrays(collection):
    """Portable, pickle-free golden record of hits and every group field."""
    arrays = {"lookup": collection.monoisos}
    for index, peak in enumerate(collection.peaks):
        for field in ("sequence_match", "z", "mz", "monoiso", "monoisos", "matchedintensity",
                      "avgmass", "scan", "rt", "massdist", "isodist", "matchedcentroids"):
            arrays[f"peak_{index}_{field}"] = np.asarray(getattr(peak, field))
    for index, group in enumerate(collection.masses):
        for field, value in vars(group).items():
            if field == "clusters":
                value = [collection.peaks.index(p) for p in value]
            elif field == "scan_intensities":
                value = sorted(value.items())
            arrays[f"group_{index}_{field}"] = np.asarray(value)
    return arrays


@pytest.mark.parametrize("envelope", [True, False])
def test_native_recovers_centroids_and_handles_missing_envelope(envelope):
    wrapper = _native_wrapper()
    peak = _peak(1000., 2.)
    peak.matchedcentroids = None
    peak.centroids = np.vstack((peak.centroids, [[501.027276467, 100.]]))
    peak.centroids = peak.centroids[np.argsort(peak.centroids[:, 0])]
    if not envelope:
        peak.isodist = None
    # A stale cache must not override the current accepted centroids.
    peak.decon_centroids = np.array([[2000., 999.]])
    actual = _native_collection(wrapper, [peak])
    reference = _reference([peak], wrapper.config, "original")
    _assert_groups(actual, reference)
    assert len(actual.masses[0].decon_centroids) == (3 if envelope else 4)


def test_native_zero_intensity_centroids_remain_finite():
    wrapper = _native_wrapper()
    peak = _peak(1000., 2.)
    peak.matchedcentroids = np.vstack(([[500., 0.]], peak.matchedcentroids))
    peaks = [peak, deepcopy(peak)]
    actual = _native_collection(wrapper, peaks)
    _assert_groups(actual, _reference(peaks, wrapper.config, "original"))
    assert len(actual.masses) == 1
    assert np.isfinite(actual.masses[0].decon_centroids).all()


@pytest.mark.parametrize("field,value", [
    ("massdist", np.ones(4)), ("massdist", np.ones((3, 1))),
    ("monoisos", np.ones((2, 2))), ("centroids", np.ones((3, 3))),
    ("massdist", np.array([[1000., np.nan]])), ("massdist", np.empty((0, 2))),
    ("scan", 2 ** 32), ("z", 2 ** 32), ("scan", np.inf), ("z", 2.5),
])
def test_native_rejects_malformed_buffers_before_call(field, value, monkeypatch):
    wrapper = _native_wrapper()
    peak = _peak(1000., 2.)
    setattr(peak, field, value)
    def unexpected(*args):
        pytest.fail("invalid input reached C")
    monkeypatch.setattr(wrapper, "_mass_group", unexpected)
    with pytest.raises(ValueError):
        wrapper.group_mass_peaks_batch([peak], wrapper.config)


@pytest.mark.parametrize("value", [-1, 65, 2 ** 32, 1.5])
def test_native_rejects_invalid_shift_count(value):
    wrapper = _native_wrapper()
    wrapper.config.maxshift = value
    with pytest.raises(ValueError, match="maxshift"):
        wrapper.group_mass_peaks_batch([_peak(1000., 2.)], wrapper.config)


@pytest.mark.parametrize("order", ["original", "matched_intensity"])
def test_native_synthetic_merge_branches_and_scan_accounting(order):
    wrapper = _native_wrapper()
    peaks = []
    for i, (offset, intensity, scan, charge) in enumerate([
        (0., 2., 1, 2), (.002, 6., 1, 3), (1.0033, 4., 2, 2),
        (-1.0033, 4., 3, 4), (.001, 8., 104, 2), (0., 3., 105, 3),
        (10., 2., 1, 2), (10.003, 5., 101, 2),
    ]):
        peak = _peak(1000. + offset, intensity)
        peak.scan, peak.rt = scan, i * .25
        peak.z = charge
        peak.isodist = peak.massdist.copy()
        peak.isodist[:, 0] = peak.massdist[:, 0] / charge + peak.config.adductmass
        peak.centroids = peak.isodist.copy()
        peak.matchedcentroids = peak.centroids.copy()
        peak.mz = peak.centroids[0, 0]
        peaks.append(peak)
    for seed in range(5):
        shuffled = [peaks[i] for i in np.random.default_rng(seed).permutation(len(peaks))]
        _assert_groups(_native_collection(wrapper, shuffled, order),
                       _reference(shuffled, wrapper.config, order))


@pytest.mark.parametrize("offset,expected", [(0.004999, 1), (0.005001, 2)])
def test_native_mass_tolerance_boundary(offset, expected):
    wrapper = _native_wrapper()
    wrapper.config.matchtol = 5.
    peaks = [_peak(1000., 2.), _peak(1000. + offset, 6.)]
    actual = _native_collection(wrapper, peaks)
    _assert_groups(actual, _reference(peaks, wrapper.config, "original"))
    assert len(actual.masses) == expected


@pytest.mark.parametrize("order", ["original", "matched_intensity"])
def test_native_first_scan_and_python_accumulation_match_fallback(spectrum, order):
    wrapper = _native_wrapper()
    fallback = _native_wrapper()
    fallback._mass_group = fallback._group_centroid_match = None
    actual, expected = MatchedCollection(), MatchedCollection()
    for scan in (1, 2, 103):
        config = deepcopy(wrapper.config)
        config.activescan, config.activescanrt = scan, scan / 10
        config.mass_group_order = order
        actual = wrapper.process_spectrum(spectrum, pks=actual, config=config)
        expected = fallback.process_spectrum(spectrum, pks=expected, config=deepcopy(config))
        _assert_groups(actual, expected)


@pytest.mark.parametrize("order", ["original", "matched_intensity"])
def test_repository_etd_native_python_and_exports(order, tmp_path):
    wrapper = _native_wrapper()
    tests = Path(__file__).parent
    data = tests / "ca_etd_unidecfiles"
    spectrum = np.loadtxt(tests / "ca_etd.dat")
    sequence = "".join(line.strip() for line in (data / "seq.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    wrapper.config.mass_group_order = order
    wrapper.config.activescan = 1
    wrapper.config.activescanrt = 0
    wrapper.config.activescanorder = 2
    kwargs = dict(centroided=True, fragmentation_type="ETD")
    actual = brute_force_pep_match(sequence, spectrum, config=wrapper.config,
                                  native_wrapper=wrapper, **kwargs)
    expected = brute_force_pep_match(sequence, spectrum, config=deepcopy(wrapper.config),
                                    native=False, **kwargs)
    # IsoGen 1.1.3 produces seven additional ETD hits.
    assert len(actual) == len(expected) == 801
    assert len(actual.masses) == (300 if order == "original" else 299)
    assert [(p.sequence_match, p.z) for p in actual] == [(p.sequence_match, p.z) for p in expected]
    _assert_groups(actual, expected)
    if order == "original":
        observed = _snapshot_arrays(actual)
        with gzip.open(data / "reference_original.json.gz", "rt", encoding="utf-8") as stream:
            frozen = json.load(stream)
        assert set(observed) == set(frozen)
        for field, value in observed.items():
            if value.dtype.kind in "iUS":
                np.testing.assert_array_equal(value, frozen[field], err_msg=field)
            else:
                np.testing.assert_allclose(value, frozen[field], rtol=FROZEN_SNAPSHOT_RTOL,
                                           atol=FROZEN_SNAPSHOT_ATOL,
                                           err_msg=field)
    np.testing.assert_allclose(actual.to_mass_spectrum(), expected.to_mass_spectrum(),
                               rtol=1e-12, atol=1e-9)
    outputs = []
    reader = SimpleNamespace(get_isolation_mz_width=lambda scan: (1000., 2000.))
    for label, collection in (("native", actual), ("python", expected)):
        tsv = tmp_path / (label + ".tsv")
        collection.export_tsv(tsv)
        prefix = str(tmp_path / label)
        collection.export_msalign(wrapper.config, reader, filename=prefix, act_type="ETD")
        ms2 = Path(prefix + "_ms2.msalign").read_text()
        rows = [line for line in ms2.splitlines()
                if line and not line.startswith("#") and not line.startswith(("FILE NAME=", "FILE_NAME="))]
        assert "BEGIN IONS" in rows and "END IONS" in rows
        outputs.append((tsv.read_text(), rows))
    assert outputs[0] == outputs[1]

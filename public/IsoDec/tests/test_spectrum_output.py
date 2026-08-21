from pathlib import Path

import numpy as np
import pytest

from isodec import IsoDecRuntime
from isodec.io import TextSpectrumImporter


def test_supplied_spectrum_has_expected_shape_and_range(spectrum):
    assert spectrum.shape == (3537, 2)
    assert np.all(np.diff(spectrum[:, 0]) > 0)
    assert spectrum[0, 0] == pytest.approx(500.0240173)
    assert spectrum[-1, 0] == pytest.approx(1989.4266357)
    assert np.all(spectrum[:, 1] >= 0)


def test_deconvolution_output_regression(processed_spectrum):
    peaks = processed_spectrum.peaks
    assert len(peaks) == 277
    assert len(processed_spectrum.masses) == 195
    assert all(1 <= peak.z <= 50 for peak in peaks)
    assert all(np.isfinite(peak.monoiso) for peak in peaks)
    assert all(peak.matchedintensity > 0 for peak in peaks)

    fingerprint = [(peak.z, peak.mz, peak.monoiso) for peak in peaks[:8]]
    expected = [
        (1, 647.2921, 646.2848),
        (8, 734.2964, 5863.3037),
        (7, 801.0299, 5597.1484),
        (7, 839.0521, 5863.3047),
        (7, 877.2175, 6130.4619),
        (7, 915.0968, 6395.6172),
        (6, 934.3674, 5597.1514),
        (6, 978.7267, 5863.3076),
    ]
    assert [item[0] for item in fingerprint] == [item[0] for item in expected]
    np.testing.assert_allclose(
        [[item[1], item[2]] for item in fingerprint],
        [[item[1], item[2]] for item in expected],
        rtol=0,
        atol=0.03,
    )


def test_builtin_text_reader_and_file_pipeline(tmp_path):
    source = Path(__file__).with_name("test_spectrum.txt")
    reader = TextSpectrumImporter(source)
    assert reader.scans.tolist() == [1]
    assert reader.get_single_scan(1).shape == (3537, 2)

    runtime = IsoDecRuntime()
    returned_reader = runtime.process_file(
        str(source), assume_centroided=True, verbose=False
    )
    assert isinstance(returned_reader, TextSpectrumImporter)
    assert len(runtime.pks) == 277

    output = tmp_path / "assignments.tsv"
    runtime.pks.export_tsv(output)
    rows = output.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("Charge\tMost Abundant m/z\tMonoisotopic Mass")
    assert len(rows) > len(runtime.pks)

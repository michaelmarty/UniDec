import numpy as np
import pytest

from isodec.fwhm import _sigma


def test_sigma_integrates_gaussian_without_deprecated_numpy_api():
    x = np.linspace(-8., 8., 2001)
    data = np.column_stack((x, np.exp(-x * x / 2)))
    assert _sigma(data) == pytest.approx(1., abs=1e-10)
    assert _sigma(np.array([[1., 2.]])) == -1.

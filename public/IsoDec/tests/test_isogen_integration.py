import numpy as np

import isogen
from isodec.isotope import (
    calc_isotope_dist,
    calc_isotope_dist_dual,
    fft_gen_isodist,
    nn_gen_isodist,
)


def test_fft_adapter_is_exactly_isogen():
    expected = isogen.isodist(5000.0, isolen=64, method="FFT", dist_only=True)
    actual = fft_gen_isodist(5000.0, isolen=64)
    np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=1e-7)


def test_nn_adapter_is_exactly_isogen():
    expected = isogen.isodist(5000.0, isolen=64, method="NN", dist_only=True)
    actual = nn_gen_isodist(5000.0, isolen=64)
    np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=1e-7)


def test_averagine_compatibility_shapes_mass_and_mz_axes():
    mz_distribution, mass_distribution = calc_isotope_dist_dual(
        5000.0, charge=5, isotopethresh=0.001
    )
    assert mz_distribution.shape == mass_distribution.shape
    assert mz_distribution.shape[1] == 2
    np.testing.assert_allclose(mz_distribution[:, 1], mass_distribution[:, 1])
    np.testing.assert_allclose(
        mz_distribution[:, 0],
        (mass_distribution[:, 0] + 5 * 1.007276467) / 5,
    )

    direct = calc_isotope_dist(5000.0, charge=5, threshold=0.001)
    np.testing.assert_allclose(direct, mz_distribution)

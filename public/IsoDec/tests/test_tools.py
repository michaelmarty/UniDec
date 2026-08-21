import numpy as np

from isodec import tools


def test_data_helpers_do_not_mutate_inputs():
    data = np.array([[1.0, 1.0], [2.0, 3.0], [3.0, 2.0]])
    original = data.copy()
    chopped = tools.datachop(data, 1.5, 2.5)
    removed = tools.dataremove(data, 1.5, 2.5)
    np.testing.assert_array_equal(data, original)
    np.testing.assert_array_equal(chopped, [[2.0, 3.0]])
    np.testing.assert_array_equal(removed[:, 1], [1.0, 0.0, 2.0])


def test_safe_division_and_ppm_boundaries():
    np.testing.assert_array_equal(
        tools.safedivide(np.array([2.0, 2.0]), np.array([2.0, 0.0])),
        [1.0, 0.0],
    )
    assert tools.within_ppm(1000.0, 1000.004, 5)
    assert not tools.within_ppm(1000.0, 1000.006, 5)

import unittest

import numpy as np

from unidec import tools
from unidec.modules import MassSpecBuilder


class TestMassSpectrumBuilder(unittest.TestCase):
    @staticmethod
    def _charge_areas(resolution, fix_undersampling):
        mass = 10_000
        (spectrum, charges), params = MassSpecBuilder.simple_spectrum2(
            [mass], resolution=resolution, zwidth=1, fix_undersampling=fix_undersampling,
        )
        centers = (mass + 1.00727647 * charges) / charges
        areas = np.array([
            np.sum(spectrum[np.abs(spectrum[:, 0] - center) < 5, 1]) for center in centers
        ])
        return areas / np.sum(areas), tools.ndis_std(charges, params[0, 2], 1, norm_area=True)

    def test_undersampling_fix_preserves_high_resolution_charge_areas(self):
        for resolution in (5_000, 100_000):
            with self.subTest(resolution=resolution):
                observed, expected = self._charge_areas(resolution, True)
                np.testing.assert_allclose(observed, expected / np.sum(expected), atol=1e-12)

        legacy, expected = self._charge_areas(5_000, False)
        self.assertGreater(np.max(np.abs(legacy - expected / np.sum(expected))), 0.1)


if __name__ == "__main__":
    unittest.main()

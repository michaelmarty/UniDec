"""IsoGen-backed isotope-distribution helpers."""

import numpy as np
from isogen import isodist


mass_diff_c = 1.0033


def _intensities(mass, isolen, analyte_type, method, dtype=float):
    return np.asarray(
        isodist(
            mass,
            type=analyte_type.upper(),
            isolen=isolen,
            method=method,
            dist_only=True,
        ),
        dtype=dtype,
    )


def nn_gen_isodist(mass, type="PEPTIDE", isolen=64):
    """Return IsoGen neural-network isotope intensities."""
    if type is None:
        return None
    return _intensities(mass, isolen, type, "NN", np.float32)


def fft_gen_isodist(mass, type="PEPTIDE", isolen=128):
    """Return IsoGen FFT isotope intensities."""
    if type is None:
        return None
    return _intensities(mass, isolen, type, "FFT", np.float32)


def _mass_distribution(mass, isolen, analyte_type):
    intensities = _intensities(mass, isolen, analyte_type, "FFT")
    masses = mass + np.arange(len(intensities)) * mass_diff_c
    return np.column_stack((masses, intensities))


def calc_isotope_dist(
    mass, charge=1, adductmass=1.007276467, isolen=128, threshold=0.001,
):
    """Return a two-column m/z and intensity distribution using IsoGen."""
    distribution = _mass_distribution(mass, isolen, "PEPTIDE")
    if distribution.size:
        distribution = distribution[
            distribution[:, 1] > distribution[:, 1].max() * threshold
        ]
    if abs(charge) >= 1:
        distribution[:, 0] = (distribution[:, 0] + charge * adductmass) / abs(charge)
    return distribution


def calc_isotope_dist_dual(
    mass,
    charge=1,
    adductmass=1.007276467,
    isotopethresh=0.01,
    type="PEPTIDE",
):
    """Return paired m/z and neutral-mass distributions using IsoGen."""
    analyte_type = "PEPTIDE" if type is None else type
    distribution = _mass_distribution(mass, 128, analyte_type)
    if distribution.size:
        distribution = distribution[
            distribution[:, 1] > distribution[:, 1].max() * isotopethresh
        ]
    mass_distribution = distribution.copy()
    if abs(charge) >= 1:
        distribution[:, 0] = (distribution[:, 0] + charge * adductmass) / abs(charge)
    return distribution, mass_distribution


__all__ = [
    "calc_isotope_dist",
    "calc_isotope_dist_dual",
    "fft_gen_isodist",
    "mass_diff_c",
    "nn_gen_isodist",
]

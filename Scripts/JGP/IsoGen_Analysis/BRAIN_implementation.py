import numpy as np
import matplotlib as mpl
from unidec.IsoDec.IsoGen.isogenc import *
from unidec.modules.isotopetools import makemass
from unidec.IsoDec.IsoGen.isogen_tools import pepmass_to_isolen

PC = [0.9893, 0.0107]                # 12C, 13C
PH = [0.999885, 0.000115]            # 1H, 2H
PN = [0.99632, 0.00368]              # 14N, 15N
PO = [0.99757, 0.00038, 0.00205]     # 16O, 17O, 18O
PS = [0.9499, 0.0075, 0.0425, 0.0001]  # 32S, 33S, 34S, 36S
PP = [1,0]

C_root = -92.4579
H_root = -8694.6521
N_root = -270.7391
rho_O = 22.0595
theta_O = 1.5750

S_real = -424.8761
rho_S = 4.728
theta_S = 1.583


def element_roots(probabilities):
    """
    Given isotope probabilities ordered by increasing neutron number,
    return the roots of the element polynomial.

    Example:
        Carbon: [P(12C), P(13C)]
        Oxygen: [P(16O), P(17O), P(18O)]
    """
    # Polynomial coefficients must be highest degree first for numpy.roots
    coeffs = list(reversed(probabilities))
    return np.roots(coeffs)


def compute_psi(isolist, l):
    psi = 0
    psi += isolist[0] * (C_root**-l) #Add C portion
    psi += isolist[1] * (H_root**-l) #Add H portion
    psi += isolist[2] * (N_root**-l) #Add N portion

    psi += isolist[3] * 2 * (rho_O ** -l) * np.cos(l*theta_O) #Add the O portion

    psi += isolist[4] * (S_real**-l) #Add the real portion of S
    psi += isolist[4] * 2 * (rho_S ** -l) * np.cos(l * theta_S) #Add the complex portion of S

    return psi


def calc_brain_dist(mass, isolist):
    #n = isolist[0] + isolist[1] + isolist[2] + 2*isolist[3] + 4*isolist[4]
    n = pepmass_to_isolen(mass)
    # Initialize distribution
    q = np.zeros(n + 1, dtype=np.complex128)
    psi = np.zeros(n+1, dtype=np.complex128)

    q[0] = PC[0] * PH[0] * PN[0] * PO[0] * PS[0]

    for l in range(1, n+1):
        psi[l] = compute_psi(isolist, l)

    # Newton–Girard recursion
    for j in range(1, n + 1):
        s = 0.0 + 0.0j

        for l in range(1, j + 1):
            s += q[j - l] * psi[l]

        q[j] = -s / float(j)

    # Convert to real numbers if imaginary part negligible
    q_real = [qj.real for qj in q]

    return q_real





if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    protein_mass = 58960.38

    formula, minmassint, isolist = makemass(protein_mass)

    c_dist = fft_pep_mass_to_dist(protein_mass)
    brain_dist = calc_brain_dist(protein_mass, isolist)
    brain_dist = brain_dist / max(brain_dist)

    plt.plot(c_dist, label="FFT-C")
    plt.plot(brain_dist, label="BRAIN")
    plt.legend()
    plt.xlim(0)
    plt.ylim(0)
    plt.show()
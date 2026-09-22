"""Match deconvolved IsoDec peaks to theoretical protein fragments."""

import numpy as np
import pandas as pd

import isogen


ION_TYPES = tuple("abcxyz")


def _fragment_parts(label, sequence_length):
    """Return the ion series and cleavage site for an IsoGen fragment label."""
    base = label[0]
    remainder = label[1:].partition("#")[0]
    if remainder.startswith("'"):
        ion_type = base + "'"
        remainder = remainder[1:]
    elif remainder.endswith(("+1", "-1")):
        ion_type = base + remainder[-2:]
        remainder = remainder[:-2]
    else:
        ion_type = base
    length = int(remainder)
    row = length if base in "abc" else sequence_length - length
    return ion_type, row


def _peak_masses(peak):
    """Return the usable monoisotopic mass candidates for one peak."""
    masses = getattr(peak, "monoisos", None)
    if masses is None or len(masses) == 0:
        masses = [getattr(peak, "monoiso", np.nan)]
    return [float(mass) for mass in masses if np.isfinite(mass) and mass > 0]


def match_fragments(
    pks,
    sequence,
    ion_types=None,
    monoisotopic=True,
    ppm_tolerance=5,
    allow_duplicate_assignments=False,
    **isogen_kwargs,
):
    """Match a sequence's theoretical fragments to an IsoDec peak collection.

    The function annotates each peak's ``sequence_match`` and stores the summary
    results on ``pks`` as ``fragment_match_percent``, ``sequence_coverage``, and
    ``fragment_matches``. In duplicate mode, ``sequence_match`` is a list;
    otherwise it is the best label or ``None``.

    Args:
        pks: Collection whose peaks contain ``monoisos`` or ``monoiso`` masses.
        sequence: Amino acid sequence with optional ProForma annotations.
        ion_types: Ion series accepted by IsoGen, or ``None`` for its default.
        monoisotopic: Use monoisotopic rather than average fragment masses.
        ppm_tolerance: Maximum absolute mass error in ppm.
        allow_duplicate_assignments: Retain every valid label for each peak.
        **isogen_kwargs: Additional ``calc_pep_fragments`` options, including
            ``fragmentation_type`` and ``ambiguous_rule``.

    Returns:
        The annotated ``pks`` collection.
    """
    if ppm_tolerance < 0:
        raise ValueError("ppm_tolerance must be non-negative")

    plain_sequence = isogen.strip_proforma(sequence)
    theoretical = isogen.calc_pep_fragments(
        sequence,
        ion_types=ion_types,
        monoisotopic=monoisotopic,
        **isogen_kwargs,
    )

    fragment_parts = {
        label: _fragment_parts(label, len(plain_sequence))
        for label in theoretical
    }
    table_ion_types = dict.fromkeys(
        (*ION_TYPES, *(ion_type for ion_type, _ in fragment_parts.values()))
    )

    table = pd.DataFrame(
        index=pd.RangeIndex(1, len(plain_sequence), name="residue")
    )
    for ion_type in table_ion_types:
        table[f"{ion_type}_mass"] = np.nan
        table[f"{ion_type}_match"] = np.nan
    for label, mass in theoretical.items():
        ion_type, row = fragment_parts[label]
        table.loc[row, f"{ion_type}_mass"] = mass

    for peak in pks.peaks:
        peak.sequence_match = [] if allow_duplicate_assignments else None

    candidates = []
    for peak_index, peak in enumerate(pks.peaks):
        for observed_mass in _peak_masses(peak):
            for label, theoretical_mass in theoretical.items():
                ppm_error = abs(observed_mass - theoretical_mass) / theoretical_mass * 1e6
                if ppm_error <= ppm_tolerance:
                    candidates.append((ppm_error, peak_index, label, observed_mass))
    candidates.sort()

    matched_peaks = set()
    for _, peak_index, label, observed_mass in candidates:
        if not allow_duplicate_assignments and peak_index in matched_peaks:
            continue

        peak = pks.peaks[peak_index]
        if allow_duplicate_assignments:
            if label not in peak.sequence_match:
                peak.sequence_match.append(label)
        else:
            peak.sequence_match = label
            matched_peaks.add(peak_index)

        ion_type, row = fragment_parts[label]
        column = f"{ion_type}_match"
        if pd.isna(table.loc[row, column]):
            table.loc[row, column] = observed_mass

    match_columns = [f"{ion_type}_match" for ion_type in table_ion_types]
    table["match_count"] = table[match_columns].notna().sum(axis=1)

    assigned_peaks = sum(bool(peak.sequence_match) for peak in pks.peaks)
    pks.fragment_match_percent = (
        100.0 * assigned_peaks / len(pks.peaks) if pks.peaks else 0.0
    )
    pks.sequence_coverage = (
        float((table["match_count"] > 0).mean()) if len(table) else 0.0
    )
    pks.fragment_matches = table
    return pks

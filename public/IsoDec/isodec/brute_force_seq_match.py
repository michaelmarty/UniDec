"""Search a spectrum for the isotope envelopes of sequence fragments."""

import numpy as np
import isogen

from .config import IsoDecConfig
from .datatools import datacompsub, get_all_centroids
from .isotope import calc_isotope_dist_dual
from .match import MatchedCollection, MatchedPeak, calculate_cosinesimilarity, find_matches


def brute_force_pep_match(
    sequence,
    spectrum,
    fragmentation_type=None,
    ion_types=None,
    *,
    centroided=False,
    config=None,
    max_charge=None,
    monoisotopic=True,
    **isogen_kwargs,
):
    """Find theoretical sequence fragments directly in an m/z spectrum.

    ``spectrum`` has m/z and intensity columns. Profile data are centroided
    using IsoDec's peak detector unless ``centroided`` is true. IsoGen's
    ``calc_pep_fragments`` accepts the fragmentation options and any additional
    keyword arguments. One peak is returned per matched fragment and charge.
    Charge states are bounded by the observed m/z range; ``max_charge`` can
    impose a further bound.
    """
    config = IsoDecConfig() if config is None else config
    if not monoisotopic:
        raise ValueError("isotope-envelope matching requires monoisotopic fragment masses")
    if max_charge is not None and (not isinstance(max_charge, int) or max_charge < 1):
        raise ValueError("max_charge must be a positive integer")

    spectrum = np.asarray(spectrum, dtype=float)
    if spectrum.ndim != 2 or spectrum.shape[1] != 2:
        raise ValueError("spectrum must have m/z and intensity columns")
    if len(spectrum) and (not np.isfinite(spectrum).all() or
                          np.any(spectrum[:, 0] <= 0) or np.any(spectrum[:, 1] < 0)):
        raise ValueError("spectrum must contain finite, positive m/z and nonnegative intensity")

    pks = MatchedCollection()
    if len(spectrum) == 0:
        return pks
    spectrum = spectrum[np.argsort(spectrum[:, 0])]
    if not centroided:
        if config.background_subtraction > 0:
            spectrum = datacompsub(spectrum.copy(), config.background_subtraction)
        spectrum = get_all_centroids(spectrum, window=config.peakwindow,
                                     threshold=config.peakthresh * 0.1)
    if len(spectrum) == 0:
        return pks
    spectrum = spectrum[np.argsort(spectrum[:, 0])]
    spectrum = spectrum[spectrum[:, 1] > 0]
    if len(spectrum) == 0:
        return pks

    fragments = isogen.calc_pep_fragments(
        sequence, fragmentation_type=fragmentation_type, ion_types=ion_types,
        monoisotopic=monoisotopic, **isogen_kwargs,
    )
    min_mz, max_mz = spectrum[0, 0], spectrum[-1, 0]
    adduct = config.adductmass
    for label, mass in fragments.items():
        mass = float(mass)
        if not np.isfinite(mass) or mass <= 0:
            continue
        # Generate the neutral isotope axis once, then convert it for each z.
        _, massdist = calc_isotope_dist_dual(
            mass, isotopethresh=config.isotopethreshold,
        )
        if len(massdist) == 0 or not np.any(massdist[:, 1] > 0):
            continue
        first_mass, last_mass = massdist[0, 0], massdist[-1, 0]
        if max_mz <= adduct:
            continue
        first_z = max(1, int(np.ceil(first_mass / (max_mz - adduct))))
        if min_mz > adduct:
            last_z = int(np.floor(last_mass / (min_mz - adduct)))
        elif max_charge is not None:
            last_z = max_charge
        else:
            raise ValueError("max_charge is required when the spectrum reaches the adduct m/z")
        if max_charge is not None:
            last_z = min(last_z, max_charge)

        for z in range(first_z, last_z + 1):
            isodist = massdist.copy()
            isodist[:, 0] = massdist[:, 0] / z + adduct
            # Match within an isotope envelope plus IsoDec's m/z tolerance.
            margin = isodist[-1, 0] * config.matchtol * 1e-6
            left = np.searchsorted(spectrum[:, 0], isodist[0, 0] - margin)
            right = np.searchsorted(spectrum[:, 0], isodist[-1, 0] + margin, side="right")
            local = spectrum[left:right]
            if len(local) < config.minpeaks:
                continue
            matched, isotopes = find_matches(local, isodist, config.matchtol)
            if len(set(matched)) < config.minpeaks:
                continue

            observed = np.zeros(len(isodist))
            observed[np.asarray(isotopes, dtype=int)] = local[np.asarray(matched, dtype=int), 1]
            score = calculate_cosinesimilarity(observed, isodist[:, 1], 0, 0,
                                               minusoneaszero=config.minusoneaszero)
            if score < config.css_thresh:
                continue
            matched_iso = isodist[np.asarray(isotopes, dtype=int)]
            scale = np.max(observed) / np.max(isodist[:, 1])
            area_covered = np.sum(matched_iso[:, 1] * scale) / np.sum(local[:, 1])
            top_matched = np.sort(local[np.asarray(matched, dtype=int), 1])[-config.minpeaks:]
            top_local = np.sort(local[:, 1])[-config.minpeaks:]
            if area_covered <= config.minareacovered and not np.array_equal(top_matched, top_local):
                continue

            isodist[:, 1] *= scale
            matched_massdist = massdist.copy()
            matched_massdist[:, 1] *= scale
            peak = MatchedPeak(z, float(isodist[np.argmax(isodist[:, 1]), 0]),
                               centroids=local, isodist=isodist,
                               matchedindexes=matched, isomatches=isotopes, config=config)
            peak.monoiso = mass
            peak.monoisos = [mass]
            peak.peakint = float(np.max(local[np.asarray(matched, dtype=int), 1]))
            peak.massdist = matched_massdist
            peak.peakmass = float(np.average(massdist[:, 0], weights=massdist[:, 1]))
            peak.avgmass = peak.peakmass
            peak.sequence_match = label
            peak.match_score = float(score)
            peak.scan = config.activescan
            peak.rt = config.activescanrt
            peak.ms_order = config.activescanorder
            pks.add_peak(peak)
            pks.add_pk_to_masses(peak, config)
    return pks

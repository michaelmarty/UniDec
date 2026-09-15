#!/usr/bin/env python3
"""Annotate HEPES/sodium clusters in the averaged Apex MS2 spectrum.

The script averages the complete scan range with UniDecImporter, assigns two
singly charged HEPES/sodium cluster ladders, and writes an annotated figure
plus a peak table beside the RAW file.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.signal import find_peaks

from UniDecImporter import ImporterFactory


RAW_FILE = Path(
    r"Z:\Group Share\BHT\Data\Apex Data\Tarjani Lab Nanodisc\20260903"
    r"\20260903_BHT_ADCK4_MSP1D1_POPC_nanodisc_MS2_8000-3000_HCDAbs_60V.raw"
)
OUTPUT_BASENAME = RAW_FILE.parent / f"{RAW_FILE.stem}_HEPES_Na_annotated"

FULL_MZ_RANGE = (500.0, 4000.0)
ANNOTATION_MZ_RANGE = (500.0, 1100.0)
PEAK_MINIMUM_PERCENT = 0.05
MASS_TOLERANCE_PPM = 5.0
SHOW_PLOT = True

# Monoisotopic masses. The sodium and proton values are ionic masses because
# the assigned formulas are positive ions.
HEPES_MASS = 238.09872824168  # C8H18N2O4S
SODIUM_ION_MASS = 22.9897692820 - 0.000548579909065
PROTON_MASS = 1.007276466621
HEPES_NA_MINUS_H_MASS = HEPES_MASS + SODIUM_ION_MASS - PROTON_MASS

FAMILY_COLORS = {
    "balanced": "#0072B2",
    "sodium-rich": "#D55E00",
}

# (family, HEPES count, sodium count, removed-proton count)
ASSIGNMENTS = (
    ("balanced", 2, 2, 1),
    ("balanced", 3, 3, 2),
    ("sodium-rich", 2, 3, 2),
    ("sodium-rich", 3, 4, 3),
    ("sodium-rich", 4, 5, 4),
)


def theoretical_mz(n_hepes: int, n_sodium: int, n_protons: int) -> float:
    """Return the theoretical m/z of a singly charged HEPES/Na cluster."""
    return (
        n_hepes * HEPES_MASS
        + n_sodium * SODIUM_ION_MASS
        - n_protons * PROTON_MASS
    )


def plain_label(n_hepes: int, n_sodium: int, n_protons: int) -> str:
    """Return a plain-text ion label for the CSV and console."""
    proton_count = "" if n_protons == 1 else f"{n_protons} "
    return f"[{n_hepes} HEPES + {n_sodium} Na - {proton_count}H]+"


def plot_label(n_hepes: int, n_sodium: int, n_protons: int) -> str:
    """Return a math-text ion label for the figure."""
    proton_count = "" if n_protons == 1 else str(n_protons)
    return (
        rf"$[{n_hepes}\,\mathrm{{HEPES}}+{n_sodium}\,\mathrm{{Na}}"
        rf"-{proton_count}\,\mathrm{{H}}]^+$"
    )


def read_average_spectrum() -> tuple[np.ndarray, dict[str, object]]:
    """Read and average every scan in the Thermo RAW file."""
    if not RAW_FILE.is_file():
        raise FileNotFoundError(f"Thermo RAW file was not found: {RAW_FILE}")

    with ImporterFactory.create_importer(RAW_FILE, silent=True) as importer:
        first_scan = int(importer.scans[0])
        last_scan = int(importer.scans[-1])
        ms_orders = {int(importer.get_ms_order(scan)) for scan in importer.scans}
        scan_event = str(
            importer.msrun.source.GetScanEventStringForScanNumber(first_scan)
        )
        spectrum = np.asarray(
            importer.get_avg_scan(scan_range=(first_scan, last_scan)), dtype=float
        )
        metadata = {
            "first_scan": first_scan,
            "last_scan": last_scan,
            "scan_count": len(importer.scans),
            "ms_orders": ms_orders,
            "scan_event": scan_event,
        }

    if spectrum.ndim != 2 or spectrum.shape[1] != 2 or len(spectrum) == 0:
        raise ValueError(f"Expected a non-empty N x 2 spectrum; got {spectrum.shape}")
    spectrum = spectrum[np.all(np.isfinite(spectrum), axis=1)]
    spectrum = spectrum[np.argsort(spectrum[:, 0])]
    if not np.any(spectrum[:, 1] > 0):
        raise ValueError("The averaged spectrum contains no positive intensities")
    return spectrum, metadata


def centroid_profile(spectrum: np.ndarray) -> np.ndarray:
    """Pick local maxima above the configured relative-intensity threshold."""
    mz_mask = (
        (spectrum[:, 0] >= ANNOTATION_MZ_RANGE[0])
        & (spectrum[:, 0] <= ANNOTATION_MZ_RANGE[1])
    )
    region = spectrum[mz_mask]
    minimum = spectrum[:, 1].max() * PEAK_MINIMUM_PERCENT / 100.0
    indices, _ = find_peaks(region[:, 1], height=minimum)
    return region[indices]


def assign_clusters(spectrum: np.ndarray) -> list[dict[str, object]]:
    """Match the configured formulas to local maxima within the ppm tolerance."""
    peaks = centroid_profile(spectrum)
    base_peak = spectrum[:, 1].max()
    matches = []

    for family, n_hepes, n_sodium, n_protons in ASSIGNMENTS:
        calculated_mz = theoretical_mz(n_hepes, n_sodium, n_protons)
        tolerance = calculated_mz * MASS_TOLERANCE_PPM / 1e6
        candidates = peaks[np.abs(peaks[:, 0] - calculated_mz) <= tolerance]
        if len(candidates) == 0:
            raise RuntimeError(
                f"No peak matched {plain_label(n_hepes, n_sodium, n_protons)} "
                f"at {calculated_mz:.6f} m/z within {MASS_TOLERANCE_PPM:g} ppm"
            )

        observed_mz, intensity = candidates[np.argmax(candidates[:, 1])]
        matches.append(
            {
                "family": family,
                "ion": plain_label(n_hepes, n_sodium, n_protons),
                "plot_label": plot_label(n_hepes, n_sodium, n_protons),
                "n_hepes": n_hepes,
                "n_sodium": n_sodium,
                "n_protons_removed": n_protons,
                "theoretical_mz": calculated_mz,
                "observed_mz": observed_mz,
                "error_ppm": (observed_mz - calculated_mz) / calculated_mz * 1e6,
                "intensity_percent": intensity / base_peak * 100.0,
            }
        )

    return matches


def maximum_preserving_downsample(data: np.ndarray, maximum_points: int) -> np.ndarray:
    """Reduce plotting cost while retaining the highest point in each m/z bin."""
    if len(data) <= maximum_points:
        return data

    bin_size = int(np.ceil(len(data) / maximum_points))
    starts = np.arange(0, len(data), bin_size)
    maxima = [
        start + int(np.argmax(data[start:start + bin_size, 1]))
        for start in starts
    ]
    return data[maxima]


def annotate_mass_step(
    axis: plt.Axes,
    first: dict[str, object],
    second: dict[str, object],
    y_position: float,
) -> None:
    """Draw and label one observed HEPES + Na - H mass step."""
    first_mz = float(first["observed_mz"])
    second_mz = float(second["observed_mz"])
    observed_step = second_mz - first_mz
    step_error_ppm = (
        (observed_step - HEPES_NA_MINUS_H_MASS)
        / HEPES_NA_MINUS_H_MASS
        * 1e6
    )
    color = FAMILY_COLORS[str(first["family"])]

    axis.annotate(
        "",
        xy=(second_mz, y_position),
        xytext=(first_mz, y_position),
        arrowprops={"arrowstyle": "<->", "color": color, "linewidth": 1.6},
    )
    axis.text(
        (first_mz + second_mz) / 2.0,
        y_position + 1.8,
        f"{observed_step:.4f} Da ({step_error_ppm:+.1f} ppm)",
        ha="center",
        va="bottom",
        color=color,
        fontsize=9,
        weight="bold",
    )


def make_figure(
    spectrum: np.ndarray,
    matches: list[dict[str, object]],
    metadata: dict[str, object],
) -> plt.Figure:
    """Build a full-spectrum overview and an annotated HEPES/Na zoom."""
    normalized = spectrum.copy()
    normalized[:, 1] = normalized[:, 1] / normalized[:, 1].max() * 100.0

    full_mask = (
        (normalized[:, 0] >= FULL_MZ_RANGE[0])
        & (normalized[:, 0] <= FULL_MZ_RANGE[1])
    )
    zoom_mask = (
        (normalized[:, 0] >= ANNOTATION_MZ_RANGE[0])
        & (normalized[:, 0] <= ANNOTATION_MZ_RANGE[1])
    )
    full = maximum_preserving_downsample(normalized[full_mask], 25_000)
    zoom = maximum_preserving_downsample(normalized[zoom_mask], 35_000)

    figure, (overview_axis, detail_axis) = plt.subplots(
        2,
        1,
        figsize=(13, 8.5),
        gridspec_kw={"height_ratios": (1.0, 2.2)},
        constrained_layout=True,
    )
    figure.suptitle(
        "Apex MS2: HEPES + sodium cluster identification",
        fontsize=16,
        weight="bold",
    )

    overview_axis.plot(full[:, 0], full[:, 1], color="black", linewidth=0.65)
    overview_axis.axvspan(
        ANNOTATION_MZ_RANGE[0],
        ANNOTATION_MZ_RANGE[1],
        color="#E6E6E6",
        zorder=-1,
        label="annotated region",
    )
    overview_axis.set(
        xlim=FULL_MZ_RANGE,
        ylim=(0, 105),
        ylabel="Relative intensity (%)",
        title=(
            f"Average of scans {metadata['first_scan']}-{metadata['last_scan']} "
            f"({metadata['scan_count']} MS2 scans)"
        ),
    )
    overview_axis.legend(frameon=False, loc="upper right")

    detail_axis.plot(zoom[:, 0], zoom[:, 1], color="black", linewidth=0.8)
    detail_axis.set(
        xlim=ANNOTATION_MZ_RANGE,
        ylim=(0, 118),
        xlabel="m/z",
        ylabel="Relative intensity (%)",
    )
    detail_axis.text(
        0.01,
        0.97,
        (
            r"Repeated increment: HEPES + Na $-$ H = "
            f"{HEPES_NA_MINUS_H_MASS:.4f} Da"
        ),
        transform=detail_axis.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )

    label_offsets = {
        ("balanced", 2): (10, 40),
        ("balanced", 3): (-70, 30),
        ("sodium-rich", 2): (24, 8),
        ("sodium-rich", 3): (22, 25),
        ("sodium-rich", 4): (-80, -22),
    }
    for match in matches:
        family = str(match["family"])
        observed_mz = float(match["observed_mz"])
        intensity = float(match["intensity_percent"])
        color = FAMILY_COLORS[family]
        offset = label_offsets[(family, int(match["n_hepes"]))]

        detail_axis.vlines(
            observed_mz,
            0,
            intensity,
            color=color,
            linewidth=1.5,
            alpha=0.9,
        )
        detail_axis.scatter(
            observed_mz,
            intensity,
            s=24,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        detail_axis.annotate(
            f"{match['plot_label']}\n{observed_mz:.4f} m/z",
            xy=(observed_mz, intensity),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] >= 0 else "right",
            va="top" if offset[1] < 0 else "bottom",
            color=color,
            fontsize=9,
            arrowprops={
                "arrowstyle": "-",
                "color": color,
                "linewidth": 0.9,
                "shrinkA": 3,
                "shrinkB": 3,
            },
        )

    family_matches = {
        family: sorted(
            (match for match in matches if match["family"] == family),
            key=lambda match: int(match["n_hepes"]),
        )
        for family in FAMILY_COLORS
    }
    annotate_mass_step(
        detail_axis,
        family_matches["balanced"][0],
        family_matches["balanced"][1],
        24,
    )
    annotate_mass_step(
        detail_axis,
        family_matches["sodium-rich"][0],
        family_matches["sodium-rich"][1],
        89,
    )
    annotate_mass_step(
        detail_axis,
        family_matches["sodium-rich"][1],
        family_matches["sodium-rich"][2],
        66,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=FAMILY_COLORS["balanced"],
            linewidth=2,
            label=r"$[n\,\mathrm{HEPES}+n\,\mathrm{Na}-(n-1)\,\mathrm{H}]^+$",
        ),
        Line2D(
            [0],
            [0],
            color=FAMILY_COLORS["sodium-rich"],
            linewidth=2,
            label=(
                r"$[n\,\mathrm{HEPES}+(n+1)\,\mathrm{Na}"
                r"-n\,\mathrm{H}]^+$"
            ),
        ),
    ]
    detail_axis.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        fontsize=9,
    )

    for axis in (overview_axis, detail_axis):
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(direction="out")

    return figure


def write_match_table(matches: list[dict[str, object]]) -> Path:
    """Write the calculated and observed assignment values to CSV."""
    output = OUTPUT_BASENAME.with_suffix(".csv")
    columns = (
        "family",
        "ion",
        "n_hepes",
        "n_sodium",
        "n_protons_removed",
        "theoretical_mz",
        "observed_mz",
        "error_ppm",
        "intensity_percent",
    )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(matches)
    return output


def main() -> None:
    spectrum, metadata = read_average_spectrum()
    matches = assign_clusters(spectrum)

    print(f"Scan event: {metadata['scan_event']}")
    print(f"HEPES + Na - H increment: {HEPES_NA_MINUS_H_MASS:.6f} Da")
    print("\nAssigned HEPES/sodium clusters:")
    for match in matches:
        print(
            f"  {match['ion']:<31} "
            f"theory {match['theoretical_mz']:.6f}, "
            f"observed {match['observed_mz']:.6f}, "
            f"{match['error_ppm']:+.3f} ppm, "
            f"{match['intensity_percent']:.3f}%"
        )

    table_path = write_match_table(matches)
    figure = make_figure(spectrum, matches, metadata)
    png_path = OUTPUT_BASENAME.with_suffix(".png")
    svg_path = OUTPUT_BASENAME.with_suffix(".svg")
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(svg_path, bbox_inches="tight", facecolor="white")

    print(f"\nWrote {png_path}")
    print(f"Wrote {svg_path}")
    print(f"Wrote {table_path}")

    if SHOW_PLOT and "agg" not in plt.get_backend().lower():
        plt.show()
    else:
        plt.close(figure)


if __name__ == "__main__":
    main()

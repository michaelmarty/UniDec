import os
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb
from unidec import tools

matplotlib.use(os.environ.get("MPLBACKEND", "TkAgg"))
import matplotlib.pyplot as plt


OUTPUT_DIR = Path(os.environ.get(
    "ARTIFACT_SUPPRESSION_FIGURE_DIR",
    str(Path.home() / "The University of Texas at Austin" / "MartyLab - General"
        / "Papers" / "Artifact Suppression" / "Figures"),
))
MASSES = np.geomspace(10_000, 1_000_000, 10)
RESOLUTIONS = np.geomspace(20, 20000, 10)
MZ_BIN_SIZE = 1
INTENSITY_THRESHOLD = 0.001
BASELINE_RESOLUTION_THRESHOLD = 0.01
PLOT_PERCENT_THRESHOLD = 95
TOPX_VALUES = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.99]
USE_SECOND_PEAK_METRIC = True
QUALITY_LABEL = "Second-peak difference (%)" if USE_SECOND_PEAK_METRIC else "Correct peak area (%)"


def normalize(intensity):
    maximum = np.max(intensity)
    return intensity / maximum * 100 if maximum else np.zeros_like(intensity)


def invalid_resolution_reason(spectrum, mass, resolution, charge_states):
    peaks = tools.peakdetect(spectrum, window=10, threshold=INTENSITY_THRESHOLD)
    peak_count = len(peaks)
    if peak_count <= 2:
        return f"low resolution ({peak_count} peaks)"
    peak_indexes = np.searchsorted(spectrum[:, 0], peaks[:, 0])
    baseline_cutoff = BASELINE_RESOLUTION_THRESHOLD * np.max(spectrum[:, 1])
    if any(
        np.min(spectrum[left:right + 1, 1]) >= baseline_cutoff
        for left, right in zip(peak_indexes[:-1], peak_indexes[1:])
    ):
        return "not baseline resolved"
    if resolution > mass / (charge_states[-1] * MZ_BIN_SIZE):
        return "high resolution (m/z FWHM below one bin)"
    return None


def deconvolution_quality(engine, expected_mass, resolution):
    mass_data = engine.data.massdat
    if USE_SECOND_PEAK_METRIC:
        heights = sorted((peak.height for peak in engine.pks.peaks), reverse=True)
        if not heights or heights[0] <= 0:
            quality = 0
        else:
            second_height = heights[1] if len(heights) > 1 else 0
            quality = 100 * (1 - second_height / heights[0])
    else:
        full_area = np.trapezoid(mass_data[:, 1], mass_data[:, 0])
        half_width = 2 * expected_mass / resolution
        correct_peak = np.abs(mass_data[:, 0] - expected_mass) <= half_width
        correct_area = np.trapezoid(mass_data[correct_peak, 1], mass_data[correct_peak, 0])
        quality = 100 * correct_area / full_area if full_area > 0 else 0
    apex_mass = mass_data[np.argmax(mass_data[:, 1]), 0]
    return quality, apex_mass - expected_mass


def run_deconvolution(spectrum, mass, charge_states, topx=0):
    engine = ud.UniDec()
    engine.pass_data_in(spectrum)
    # engine.config.mzbins = 1
    # engine.config.linflag = 0
    # Normalize the data
    engine.config.datanorm = 1
    engine.config.intthresh = INTENSITY_THRESHOLD
    engine.config.numit = 10
    engine.config.startz = max(1, int(charge_states[0])-10)
    engine.config.endz = min(100, int(charge_states[-1])+10)
    engine.config.psig = 10
    engine.config.subbuff = 0
    engine.config.beta = 0
    engine.config.suppression_satellite = 0
    engine.config.suppression_harmonic = 0
    engine.config.suppression_topn = 0
    engine.config.suppression_topx = topx
    engine.config.masslb = mass * 0.25
    engine.config.massub = mass * 2.5
    engine.autorun(auto_peak_width=True, silent=True)
    return engine


def generate_results():
    results = []
    skipped = {"low": 0, "high": 0, "unresolved": 0}
    start = time.perf_counter()
    for mass in MASSES:
        for resolution in RESOLUTIONS:
            (spectrum, charge_states), params = msb.simple_spectrum2(
                [mass], resolution=resolution, psfun=0, baseline=0, noise=0.0,
                zwidth=1.0, mz_bin_size=MZ_BIN_SIZE,
            )
            # Normalize the spectrum intensity to a maximum of 1
            spectrum[:, 1] = normalize(spectrum[:, 1])/100.
            reason = invalid_resolution_reason(spectrum, mass, resolution, charge_states)
            if reason is not None:
                key = "low" if reason.startswith("low") else "high" if reason.startswith("high") else "unresolved"
                skipped[key] += 1
                print(f"Skipping mass={mass:g}, resolution={resolution:g}: {reason}")
                continue
            # plt.plot(spectrum[:, 0], spectrum[:, 1], label=f"Mass={mass:g}, R={resolution:g}")
            # plt.legend()
            # plt.show()
            engine = run_deconvolution(spectrum, mass, charge_states)
            assert engine.data.massdat.size, f"Deconvolution failed for mass={mass:g}, resolution={resolution:g}"
            quality, apex_error = deconvolution_quality(engine, mass, resolution)
            results.append((mass, resolution, params[0, 2], spectrum, engine.data.massdat,
                            quality, apex_error, charge_states))
    elapsed = time.perf_counter() - start
    assert len(results) + sum(skipped.values()) == len(MASSES) * len(RESOLUTIONS)
    print(
        f"Deconvolved {len(results)} spectra in {elapsed:.2f} s "
        f"({elapsed / len(results):.2f} s each); skipped {skipped['low']} low and "
        f"{skipped['high']} high resolution pairs and {skipped['unresolved']} "
        f"pairs that were not baseline resolved"
    )
    return results


def find_critical_topx(results):
    critical_results = []
    reruns = 0
    start = time.perf_counter()
    for mass, resolution, average_charge, spectrum, _, quality, _, charge_states in results:
        critical_topx = 0
        if quality <= PLOT_PERCENT_THRESHOLD:
            critical_topx = 1
            for topx in TOPX_VALUES:
                engine = run_deconvolution(spectrum, mass, charge_states, topx=topx)
                reruns += 1
                assert engine.data.massdat.size, (
                    f"Top X deconvolution failed for mass={mass:g}, resolution={resolution:g}, topx={topx}"
                )
                quality, _ = deconvolution_quality(engine, mass, resolution)
                if quality > PLOT_PERCENT_THRESHOLD:
                    critical_topx = topx
                    break
        critical_results.append((mass, resolution, average_charge, critical_topx))
    assert len(critical_results) == len(results)
    print(f"Ran {reruns} Top X deconvolutions in {time.perf_counter() - start:.2f} s")
    return critical_results


def plot_map(results):
    plotted_results = [result for result in results if result[5] < PLOT_PERCENT_THRESHOLD]
    assert plotted_results, f"No deconvolution results have correct peak area below {PLOT_PERCENT_THRESHOLD}%"
    fig, axis = plt.subplots(figsize=(6.5, 4.5))
    points = axis.scatter(
        [result[1] for result in plotted_results],
        [result[2] for result in plotted_results],
        c=[result[5] for result in plotted_results],
        s=22,
        cmap="viridis_r",
        norm=matplotlib.colors.PowerNorm(
            gamma=2, vmin=min(result[5] for result in plotted_results), vmax=100,
        ),
    )
    fig.colorbar(points, ax=axis, label=QUALITY_LABEL)
    axis.set_xscale("log")
    axis.set_yscale("function", functions=(np.sqrt, np.square))
    axis.set(xlabel="Resolution", ylabel="Average charge state")
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "resolution_map.png", dpi=300)
    return fig


def plot_topx_map(results):
    successful = [result for result in results if result[3] < 1]
    failed = [result for result in results if result[3] == 1]
    color_indexes = [
        0 if result[3] == 0 else TOPX_VALUES.index(result[3]) + 1
        for result in successful
    ]
    color_count = len(TOPX_VALUES) + 1
    color_map = matplotlib.colors.ListedColormap(
        plt.get_cmap("viridis")(np.linspace(0, 1, color_count))
    )
    fig, axis = plt.subplots(figsize=(6.5, 4.5))
    points = axis.scatter(
        [result[1] for result in successful],
        [result[2] for result in successful],
        c=color_indexes,
        s=22,
        cmap=color_map,
        norm=matplotlib.colors.BoundaryNorm(np.arange(-0.5, color_count + 0.5), color_count),
    )
    if failed:
        axis.scatter(
            [result[1] for result in failed],
            [result[2] for result in failed],
            marker="x",
            color="gray",
            s=28,
        )
    colorbar = fig.colorbar(points, ax=axis, label="Critical Top X", ticks=range(color_count))
    colorbar.ax.set_yticklabels(["0", *(f"{topx:g}" for topx in TOPX_VALUES)])
    axis.set_xscale("log")
    axis.set_yscale("function", functions=(np.sqrt, np.square))
    axis.set(xlabel="Resolution", ylabel="Average charge state")
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "resolution_map_topx.png", dpi=300)
    return fig


def plot_examples(results):
    selected = np.random.choice(len(results), 3, replace=False)
    fig, axes = plt.subplots(3, 2, figsize=(8, 8))
    for row, index in enumerate(selected):
        mass, resolution, average_charge, spectrum, mass_data, quality, apex_error, _ = results[index]
        axes[row, 0].plot(spectrum[:, 0], normalize(spectrum[:, 1]), color="darkblue")
        axes[row, 1].plot(mass_data[:, 0] / 1000, normalize(mass_data[:, 1]), color="darkblue")
        axes[row, 0].set_title(
            f"{mass / 1000:.0f} kDa, R={resolution:.0f}, z={average_charge:.1f}\n"
            f"{QUALITY_LABEL.removesuffix(' (%)')}={quality:.1f}%, apex error={apex_error:+.1f} Da",
            fontsize=10,
        )
        axes[row, 0].set_ylabel("Intensity (%)")
        axes[row, 1].set_ylabel("Intensity (%)")
        for axis in axes[row]:
            axis.set_ylim(0, 105)
            axis.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel(r"$m/z$")
    axes[-1, 1].set_xlabel("Mass (kDa)")
    axes[0, 0].text(0.98, 0.9, "Simulated spectrum", transform=axes[0, 0].transAxes,
                    ha="right", va="top")
    axes[0, 1].text(0.98, 0.9, "Deconvolved result", transform=axes[0, 1].transAxes,
                    ha="right", va="top")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "resolution_map_examples.png", dpi=300)
    return fig


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = generate_results()
    plot_map(results)
    plot_topx_map(find_critical_topx(results))
    plot_examples(results)
    plt.show()

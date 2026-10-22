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
RESOLUTIONS = np.geomspace(20, 1000, 10)
MZ_BIN_SIZE = 1
INTENSITY_THRESHOLD = 0.001


def normalize(intensity):
    maximum = np.max(intensity)
    return intensity / maximum * 100 if maximum else np.zeros_like(intensity)


def invalid_resolution_reason(spectrum, mass, resolution, charge_states):
    peak_count = len(tools.peakdetect(spectrum, window=10, threshold=INTENSITY_THRESHOLD))
    if peak_count <= 2:
        return f"low resolution ({peak_count} peaks)"
    if resolution > mass / (charge_states[-1] * MZ_BIN_SIZE):
        return "high resolution (m/z FWHM below one bin)"
    return None


def deconvolution_quality(mass_data, expected_mass, resolution):
    full_area = np.trapezoid(mass_data[:, 1], mass_data[:, 0])
    half_width = 2 * expected_mass / resolution
    correct_peak = np.abs(mass_data[:, 0] - expected_mass) <= half_width
    correct_area = np.trapezoid(mass_data[correct_peak, 1], mass_data[correct_peak, 0])
    area_percent = 100 * correct_area / full_area if full_area > 0 else 0
    apex_mass = mass_data[np.argmax(mass_data[:, 1]), 0]
    return area_percent, apex_mass - expected_mass


def run_deconvolution(spectrum, mass, charge_states):
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
    engine.config.masslb = mass * 0.25
    engine.config.massub = mass * 2.5
    engine.autorun(auto_peak_width=True, silent=True)
    return engine


def generate_results():
    results = []
    skipped = {"low": 0, "high": 0}
    start = time.perf_counter()
    for mass in MASSES:
        for resolution in RESOLUTIONS:
            (spectrum, charge_states), params = msb.simple_spectrum2(
                [mass], resolution=resolution, psfun=0, baseline=0, noise=0,
                zwidth=1.0, mz_bin_size=MZ_BIN_SIZE,
            )
            # Normalize the spectrum intensity to a maximum of 1
            spectrum[:, 1] = normalize(spectrum[:, 1])/100.
            reason = invalid_resolution_reason(spectrum, mass, resolution, charge_states)
            if reason is not None:
                skipped["low" if reason.startswith("low") else "high"] += 1
                print(f"Skipping mass={mass:g}, resolution={resolution:g}: {reason}")
                continue
            # plt.plot(spectrum[:, 0], spectrum[:, 1], label=f"Mass={mass:g}, R={resolution:g}")
            # plt.legend()
            # plt.show()
            engine = run_deconvolution(spectrum, mass, charge_states)
            assert engine.data.massdat.size, f"Deconvolution failed for mass={mass:g}, resolution={resolution:g}"
            quality, apex_error = deconvolution_quality(engine.data.massdat, mass, resolution)
            results.append((mass, resolution, params[0, 2], spectrum, engine.data.massdat,
                            quality, apex_error))
    elapsed = time.perf_counter() - start
    assert len(results) + sum(skipped.values()) == len(MASSES) * len(RESOLUTIONS)
    print(
        f"Deconvolved {len(results)} spectra in {elapsed:.2f} s "
        f"({elapsed / len(results):.2f} s each); skipped {skipped['low']} low and "
        f"{skipped['high']} high resolution pairs"
    )
    return results


def plot_map(results):
    fig, axis = plt.subplots(figsize=(6.5, 4.5))
    points = axis.scatter(
        [result[1] for result in results],
        [result[2] for result in results],
        c=[result[5] for result in results],
        s=22,
        cmap="viridis_r",
        norm=matplotlib.colors.PowerNorm(
            gamma=2, vmin=min(result[5] for result in results), vmax=100,
        ),
    )
    fig.colorbar(points, ax=axis, label="Correct peak area (%)")
    axis.set_xscale("log")
    axis.set(xlabel="Resolution", ylabel="Average charge state")
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "resolution_map.png", dpi=300)
    return fig


def plot_examples(results):
    selected = np.random.choice(len(results), 3, replace=False)
    fig, axes = plt.subplots(3, 2, figsize=(8, 8))
    for row, index in enumerate(selected):
        mass, resolution, average_charge, spectrum, mass_data, quality, apex_error = results[index]
        axes[row, 0].plot(spectrum[:, 0], normalize(spectrum[:, 1]), color="darkblue")
        axes[row, 1].plot(mass_data[:, 0] / 1000, normalize(mass_data[:, 1]), color="darkblue")
        axes[row, 0].set_title(
            f"{mass / 1000:.0f} kDa, R={resolution:.0f}, z={average_charge:.1f}\n"
            f"Area={quality:.1f}%, apex error={apex_error:+.1f} Da",
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
    # plot_examples(results)
    plt.show()

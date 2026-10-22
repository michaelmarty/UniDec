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

matplotlib.use(os.environ.get("MPLBACKEND", "TkAgg"))
import matplotlib.pyplot as plt


OUTPUT_DIR = Path(os.environ.get("ARTIFACT_SUPPRESSION_FIGURE_DIR", Path(__file__).parent))
MASSES = np.geomspace(10_000, 1_000_000, 10)
RESOLUTIONS = np.geomspace(20, 10000, 10)


def normalize(intensity):
    maximum = np.max(intensity)
    return intensity / maximum * 100 if maximum else np.zeros_like(intensity)


def run_deconvolution(spectrum, mass, charge_states):
    engine = ud.UniDec()
    engine.pass_data_in(spectrum)
    # engine.config.mzbins = 1
    # engine.config.linflag = 0
    # Normalize the data
    engine.config.datanorm = 1
    engine.config.intthresh = 0.001
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
    start = time.perf_counter()
    for mass in MASSES:
        for resolution in RESOLUTIONS:
            (spectrum, charge_states), params = msb.simple_spectrum2(
                [mass], resolution=resolution, psfun=0, baseline=0, noise=0,
                zwidth=1.0, mzbins=1,
            )
            # Normalize the spectrum intensity to a maximum of 1
            spectrum[:, 1] = normalize(spectrum[:, 1])/100.
            # plt.plot(spectrum[:, 0], spectrum[:, 1], label=f"Mass={mass:g}, R={resolution:g}")
            # plt.legend()
            # plt.show()
            engine = run_deconvolution(spectrum, mass, charge_states)
            assert engine.data.massdat.size, f"Deconvolution failed for mass={mass:g}, resolution={resolution:g}"
            results.append((mass, resolution, params[0, 2], spectrum, engine.data.massdat))
    elapsed = time.perf_counter() - start
    assert len(results) == len(MASSES) * len(RESOLUTIONS)
    print(f"Ran {len(results)} spectra in {elapsed:.2f} s ({elapsed / len(results):.2f} s each)")
    return results


def plot_map(results):
    fig, axis = plt.subplots(figsize=(6.5, 4.5))
    axis.scatter(
        [result[1] for result in results],
        [result[2] for result in results],
        s=22,
        color="darkblue",
    )
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
        mass, resolution, average_charge, spectrum, mass_data = results[index]
        axes[row, 0].plot(spectrum[:, 0], normalize(spectrum[:, 1]), color="darkblue")
        axes[row, 1].plot(mass_data[:, 0] / 1000, normalize(mass_data[:, 1]), color="darkblue")
        axes[row, 0].set_title(
            f"{mass / 1000:.0f} kDa, R={resolution:.0f}, z={average_charge:.1f}",
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
    plot_examples(results)
    plt.show()

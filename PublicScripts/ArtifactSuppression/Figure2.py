import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb
from PublicScripts.ArtifactSuppression.metrics import (
    description_with_uniscore, panel_metrics, write_metrics_workbook,
)

matplotlib.use(os.environ.get("MPLBACKEND", "WxAgg"))
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
})

FIGURE_DIR = Path(os.environ.get(
    "ARTIFACT_SUPPRESSION_FIGURE_DIR",
    str(Path.home() / "The University of Texas at Austin" / "MartyLab - General"
        / "Papers" / "Artifact Suppression" / "Figures"),
))


def normalize(data):
    """Normalize the intensity column to a maximum of 100."""
    maximum = np.max(data[:, 1])
    if maximum == 0:
        return np.zeros_like(data[:, 1])
    return data[:, 1] / maximum * 100


def clean_axes(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.set_ylim(0, 100)
    axis.set_yticks([0, 50, 100], labels=['0', '%', '100'])


def label_primary_charge_states(axis, spectrum, params, charge_states,
                                relative_threshold=0.1):
    normalized_intensity = normalize(spectrum)
    mass, _, charge_center, charge_width, _ = params[0]
    charge_weights = np.exp(
        -0.5 * ((charge_states - charge_center) / charge_width) ** 2
    )
    primary_charges = charge_states[
        charge_weights >= np.max(charge_weights) * relative_threshold
    ]

    for charge in primary_charges:
        mz = (mass + 1.00727647 * charge) / charge
        peak_index = np.argmin(np.abs(spectrum[:, 0] - mz))
        axis.annotate(
            f'{int(charge)}+',
            (spectrum[peak_index, 0], normalized_intensity[peak_index]),
            xytext=(0, 3),
            textcoords='offset points',
            ha='center',
            va='bottom',
        )


def run_deconvolution(spectrum, massrange=[5000, 50000], sat=0, harm=0, topn=0, topx=0,
                      baseline_subtraction=0, beta=0, psig=10):
    eng = ud.UniDec()
    eng.pass_data_in(spectrum)
    # eng.config.mzsig = 5
    eng.config.psig = psig
    eng.config.subbuff = baseline_subtraction
    eng.config.masslb = massrange[0]
    eng.config.massub = massrange[1]
    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = harm
    eng.config.suppression_topn = topn
    eng.config.suppression_topx = topx
    eng.config.beta = beta
    eng.autorun(auto_peak_width=True)

    return eng

# Set noise seed for reproducibility
np.random.seed(42)

# Create a sipmle spectrum with a mass of 20000 and a resolution of 500 plot it
spectrum, params = msb.simple_spectrum2([20000], resolution=50, psfun=0, baseline=0.03, noise=0.01, zwidth=0.75)

spectrum, ztab = spectrum

engines = [
    run_deconvolution(spectrum),
    run_deconvolution(spectrum, sat=1),
    run_deconvolution(spectrum, sat=2),
    run_deconvolution(spectrum, harm=1),
    run_deconvolution(spectrum, harm=1, sat=1),
    run_deconvolution(spectrum, topn=2),
    run_deconvolution(spectrum, topn=1),
    run_deconvolution(spectrum, topx=0.01),
    run_deconvolution(spectrum, topx=0.05),
]

fig = plt.figure(figsize=(6.5, 6.5))
fig.canvas.manager.resize(*map(round, fig.bbox.size))
outer_grid = fig.add_gridspec(
    2,
    1,
    height_ratios=[1, 4],
    hspace=0.2,
)
top_grid = outer_grid[0].subgridspec(1, 2, wspace=0.25)
mass_grid = outer_grid[1].subgridspec(4, 2, hspace=0.08, wspace=0.25)
axes = np.empty((5, 2), dtype=object)
for column in range(2):
    axes[0, column] = fig.add_subplot(top_grid[0, column])
    for row in range(1, 5):
        axes[row, column] = fig.add_subplot(mass_grid[row - 1, column])
flat_axes = axes.flat

axes[0, 0].plot(spectrum[:, 0], normalize(spectrum), color="darkblue")
axes[0, 0].set_xlabel(r'$\mathit{m/z}$')


mass_datasets = [
    engine.data.massdat for engine in engines
]
for axis, data in zip(list(flat_axes)[1:], mass_datasets):
    axis.plot(data[:, 0] / 1000, normalize(data), color="darkblue")
    axis.set_ylim(0, 1)

ypos = 15
for axis in (axes[0, 1], axes[2, 0], axes[3, 0], axes[4, 0]):
    axis.hlines(ypos, 16, 19, color='red', linewidth=3)
    axis.hlines(ypos, 21, 24, color='red', linewidth=3)

for axis in (axes[0, 1], axes[1, 0], axes[1, 1]):
    axis.hlines(ypos, 38.5, 41.5, color='magenta', linewidth=3)

plot_labels = [
    'Simulated Data',
    'Normal Settings',
    'Satellite 1',
    'Satellite 2',
    'Harmonic',
    'Harmonic + Satellite 1',
    'Top 2',
    'Top 1',
    'Top 0.01',
    'Top 0.05',
]
for panel_index, (axis, label) in enumerate(zip(axes.flat, plot_labels)):
    if panel_index > 0:
        label = description_with_uniscore(label, engines[panel_index - 1].pks.uniscore)
    axis.text(0.03, 0.95, chr(ord('A') + panel_index),
              transform=axis.transAxes, va='top', ha='left')
    axis.text(0.95, 0.95, label, transform=axis.transAxes, va='top', ha='right')

axes[-1, 0].set_xlabel('Mass (kDa)')
axes[-1, 1].set_xlabel('Mass (kDa)')

for axis in axes.flat:
    clean_axes(axis)

for axis in axes[:, 1:].flat:
    axis.tick_params(axis='y', labelleft=False)

# Panel A retains its m/z tick labels. Hide numeric x labels on the other
# intermediate mass panels, while preserving their axes and tick marks.
axes[0, 1].tick_params(axis='x', labelbottom=False)
for axis in axes[1:-1].flat:
    axis.tick_params(axis='x', labelbottom=False)
    axis.set_yticks([0, 50, 100], labels=['', '%', '100'])

axes[0, 0].set_ylim(0, 100)
label_primary_charge_states(axes[0, 0], spectrum, params, ztab)

fig.subplots_adjust(left=0.1, right=0.98, bottom=0.08, top=0.96)

os.chdir(FIGURE_DIR)
plt.savefig("Figure2.png", dpi=300)
plt.savefig("Figure2.pdf")
common_settings = {
    "reference_basis": "synthetic input mass",
    "simulation_resolution": 50,
    "simulation_baseline": 0.03,
    "simulation_noise": 0.01,
    "simulation_zwidth": 0.75,
    "random_seed": 42,
}
metric_rows = [panel_metrics('A', plot_labels[0], panel_type="input m/z spectrum",
                             reference_mass=20000, settings=common_settings)]
metric_rows.extend(
    panel_metrics(chr(ord('B') + index), label, engine, reference_mass=20000,
                  settings=common_settings)
    for index, (label, engine) in enumerate(zip(plot_labels[1:], engines))
)
write_metrics_workbook("Figure2.xlsx", metric_rows)


plt.show()

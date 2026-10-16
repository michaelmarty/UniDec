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
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_ylim(0, 100)
    axis.set_yticks([0, 50, 100], labels=['0', '%', '100'])


def run_deconvolution(spectrum, massrange, sat=0, topn=0, psig=10):
    eng = ud.UniDec()
    eng.pass_data_in(spectrum)
    eng.config.endz = 100
    eng.config.psig = psig
    eng.config.subbuff = 0
    eng.config.masslb = massrange[0]
    eng.config.massub = massrange[1]
    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = 0
    eng.config.suppression_topn = topn
    eng.config.suppression_topx = 0
    eng.config.beta = 0
    eng.autorun(auto_peak_width=True)

    return eng


np.random.seed(42)

masses = [150000, 600000]
resolutions = [150, 250]
satellite_1_engines = []
top_2_engines = []
top_1_engines = []

for mass, resolution in zip(masses, resolutions):
    spectrum, params = msb.simple_spectrum2(
        [mass],
        resolution=resolution,
        psfun=0,
        baseline=0.03,
        noise=0.01,
        zwidth=1,
    )
    spectrum, ztab = spectrum
    massrange = [mass * 0.25, mass * 2.5]

    satellite_1_engines.append(
        run_deconvolution(spectrum, massrange, sat=1)
    )
    top_2_engines.append(
        run_deconvolution(spectrum, massrange, topn=2)
    )
    top_1_engines.append(
        run_deconvolution(spectrum, massrange, topn=1)
    )

satellite_1_mass_datasets = [engine.data.massdat for engine in satellite_1_engines]
top_2_mass_datasets = [engine.data.massdat for engine in top_2_engines]
top_1_mass_datasets = [engine.data.massdat for engine in top_1_engines]

fig, axes = plt.subplots(3, 2, figsize=(6.5, 4.8))

for column, mass in enumerate(masses):
    axes[0, column].plot(
        top_1_mass_datasets[column][:, 0] / 1000,
        normalize(top_1_mass_datasets[column]),
        color="darkblue",
    )
    axes[1, column].plot(
        top_2_mass_datasets[column][:, 0] / 1000,
        normalize(top_2_mass_datasets[column]),
        color="darkblue",
    )
    axes[2, column].plot(
        satellite_1_mass_datasets[column][:, 0] / 1000,
        normalize(satellite_1_mass_datasets[column]),
        color="darkblue",
    )

    axes[0, column].set_title(f'{mass / 1000:g} kDa')
    axes[-1, column].set_xlabel('Mass (kDa)')

row_labels = ['Top 1', 'Top 2', 'Satellite 1']
method_engines = [top_1_engines, top_2_engines, satellite_1_engines]
for row, label in enumerate(row_labels):
    for column, axis in enumerate(axes[row]):
        axis.text(
            0.95,
            0.95,
            description_with_uniscore(label, method_engines[row][column].pks.uniscore),
            transform=axis.transAxes,
            va='top',
            ha='right',
        )

for panel_index, axis in enumerate(axes.flat):
    axis.text(
        0.03,
        0.95,
        chr(ord('A') + panel_index),
        transform=axis.transAxes,
        va='top',
        ha='left',
    )
    clean_axes(axis)

for axis in axes[:, 1:].flat:
    axis.tick_params(axis='y', labelleft=False)

ypos = 90
red_bar_panels = [
    (1, 0),  # C
    (1, 1),  # D
    (2, 1),  # F
]
for row, column in red_bar_panels:
    central_mass = masses[column] / 1000
    axes[row, column].hlines(
        ypos, central_mass * 0.80, central_mass * 0.98,
        color='red', linewidth=3
    )
    axes[row, column].hlines(
        ypos, central_mass * 1.02, central_mass * 1.20,
        color='red', linewidth=3
    )

plt.tight_layout()

os.chdir(FIGURE_DIR)
plt.savefig('FigureS2.png', dpi=300)
plt.savefig('FigureS2.pdf')
metric_rows = []
for row, (label, engines_for_method) in enumerate(zip(row_labels, method_engines)):
    for column, (mass, engine) in enumerate(zip(masses, engines_for_method)):
        settings = {
            "reference_basis": "synthetic input mass",
            "simulation_resolution": resolutions[column],
            "simulation_baseline": 0.03,
            "simulation_noise": 0.01,
            "simulation_zwidth": 1,
            "random_seed": 42,
        }
        metric_rows.append(panel_metrics(chr(ord('A') + row * 2 + column), label, engine,
                                         reference_mass=mass, settings=settings))
write_metrics_workbook('FigureS2.xlsx', metric_rows)

plt.show()

import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[3] / "public" / "UniDec"
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
})

POPC_DATA_FILE = (
    r"POPC_Nanodiscs.txt"
)
GROEL_DATA_FILE = (
    r"GroEL_SID1.RAW"
)
FIGURE_DIR = (
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)

os.chdir(FIGURE_DIR)

GROEL_CONFIG = {
    'endz': 60,
    'minmz': 990,
    'maxmz': 35500,
    'molig': 0,
    'msig': 0,
    'psig': 1,
    'mzsig': 25,
    'masslb': 1000,
    'massub': 1000000,
    'massbins': 10,
    'subbuff': 0,
}

groel_monomer_mass = 57300.0
groel_masses = groel_monomer_mass * np.arange(1, 15)

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


def label_most_intense_charge(axis, spectrum, mass):
    normalized_intensity = normalize(spectrum)
    peak_index = np.argmax(spectrum[:, 1])
    peak_mz = spectrum[peak_index, 0]
    charge = int(round(mass / (peak_mz - 1.00727647)))
    axis.annotate(
        f'{charge}+',
        (peak_mz, normalized_intensity[peak_index]),
        xytext=(-5, -6),
        textcoords='offset points',
        ha='right',
        va='top',
    )

def mark_groel_oligomers(axis, masses, label=True):
    labelcolor='goldenrod'

    # Put a label at the top of the axis for each mass
    for i, mass in enumerate(masses):
        axis.axvline(
            x=mass / 1000,
            color=labelcolor,
            linestyle='--',
            linewidth=1,
            alpha=0.5,
        )
        if label:
            axis.text(
                mass / 1000,
                95,
                f'{i + 1}',
                rotation=90,
                verticalalignment='top',
                horizontalalignment='right',
                fontsize=8,
                color=labelcolor,
            )

def run_deconvolution(file_path, sat=0, topx=0, topn=0):
    eng = ud.UniDec()
    eng.open_file(file_path, refresh=True)
    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = 0
    eng.config.suppression_topn = topn
    eng.config.suppression_topx = topx
    eng.config.beta = 0
    eng.config.minmz = 8000
    eng.config.maxmz = 13000
    eng.config.molig = 760
    eng.config.msig=1
    eng.config.psig=1
    eng.config.mzsig = 6
    eng.config.masslb = 90000
    eng.config.massub = 200000
    eng.autorun(auto_peak_width=True)
    return eng.data.data2.copy(), eng.data.massdat.copy()


def run_deconvolution_with_config(file_path, config_options, sat=0, topx=0, topn=0, beta=0):
    """Run deconvolution with file-specific UniDec configuration values."""
    eng = ud.UniDec()
    eng.open_file(file_path, refresh=True)
    for option, value in config_options.items():
        if not hasattr(eng.config, option):
            raise ValueError(f"Unknown UniDec config option: {option}")
        setattr(eng.config, option, value)

    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = 0
    eng.config.suppression_topn = topn
    eng.config.suppression_topx = topx
    eng.config.beta = beta
    eng.autorun(auto_peak_width=True)
    return eng.data.data2.copy(), eng.data.massdat.copy()


popc_spectrum, popc_default_mass_data = run_deconvolution(POPC_DATA_FILE)
_, popc_top_01_mass_data = run_deconvolution(POPC_DATA_FILE, topx=0.1)
_, popc_top_05_mass_data = run_deconvolution(POPC_DATA_FILE, topx=0.5)
_, popc_top_2_mass_data = run_deconvolution(POPC_DATA_FILE, topn=2)
_, popc_satellite_2_mass_data = run_deconvolution(POPC_DATA_FILE, sat=2)

groel_spectrum, groel_default_mass_data = run_deconvolution_with_config(
    GROEL_DATA_FILE, GROEL_CONFIG
)
_, groel_top_01_mass_data = run_deconvolution_with_config(
    GROEL_DATA_FILE, GROEL_CONFIG, topx=0.1
)
_, groel_satellite_1_mass_data = run_deconvolution_with_config(
    GROEL_DATA_FILE, GROEL_CONFIG, sat=1
)
_, groel_satellite_2_mass_data = run_deconvolution_with_config(
    GROEL_DATA_FILE, GROEL_CONFIG, sat=2
)
_, groel_beta_50_mass_data = run_deconvolution_with_config(
    GROEL_DATA_FILE, GROEL_CONFIG, beta=50
)

column_datasets = [
    (
        'Nanodiscs',
        popc_spectrum,
        [
            popc_default_mass_data,
            popc_top_01_mass_data,
            popc_top_05_mass_data,
            popc_top_2_mass_data,
            popc_satellite_2_mass_data,
        ],
    ),
    (
        'GroEL SID',
        groel_spectrum,
        [
            groel_default_mass_data,
            groel_top_01_mass_data,
            groel_satellite_1_mass_data,
            groel_satellite_2_mass_data,
            groel_beta_50_mass_data,
        ],
    ),
]

fig = plt.figure(figsize=(6.5, 7))
outer_grid = fig.add_gridspec(
    2,
    1,
    height_ratios=[1, 5],
    hspace=0.2,
)
spectrum_grid = outer_grid[0].subgridspec(1, 2, wspace=0.25)
mass_grid = outer_grid[1].subgridspec(5, 2, hspace=0.18, wspace=0.25)
axes = np.empty((6, 2), dtype=object)
for column in range(2):
    axes[0, column] = fig.add_subplot(spectrum_grid[0, column])
    for row in range(1, 6):
        axes[row, column] = fig.add_subplot(mass_grid[row - 1, column])

panel_labels = [
    ['Data', 'Data'],
    ['Normal', 'Normal'],
    ['Top 0.1', 'Top 0.1'],
    ['Top 0.5', 'Satellite 1'],
    ['Top 2', 'Satellite 2'],
    ['Satellite 2', 'SoftMax 50'],
]

for column, (title, spectrum, mass_datasets) in enumerate(column_datasets):
    axes[0, column].plot(spectrum[:, 0], normalize(spectrum), color="darkblue")
    axes[0, column].set_title(title)
    axes[0, column].set_xlabel(r'$\mathit{m/z}$')
    for row, mass_data in enumerate(mass_datasets, start=1):
        axes[row, column].plot(mass_data[:, 0] / 1000, normalize(mass_data), color="darkblue")
    axes[-1, column].set_xlabel('Mass (kDa)')

    # dominant_mass = mass_datasets[0][np.argmax(mass_datasets[0][:, 1]), 0]
    # label_most_intense_charge(axes[0, column], spectrum, dominant_mass)

for row, labels in enumerate(panel_labels):
    for axis, label in zip(axes[row], labels):
        axis.text(
            0.95,
            0.95,
            label,
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

# Keep the intermediate x axes and tick marks, but show numeric tick labels
# only on the bottom mass row.
for axis in axes[1:-1].flat:
    axis.tick_params(axis='x', labelbottom=False)
    axis.set_yticks([0, 50, 100], labels=['', '%', '100'])

ypos = 20
for axis in (axes[1, 0], axes[5, 0]):
    axis.hlines(
        ypos, 187, 200,
        color='orange', linewidth=3
    )

# GroEL mass panels D, F, H, J, and L.
for axis in axes[1:, 1]:
    labelflag = False
    # Turn on label for panel D (row 1, column 1) and turn off for the rest.
    if axis == axes[1, 1]:
        labelflag = True
    mark_groel_oligomers(axis, groel_masses, label=labelflag)

fig.subplots_adjust(left=0.06, right=0.98, bottom=0.07, top=0.96)

os.chdir(FIGURE_DIR)
plt.savefig('Figure4.png', dpi=300)
plt.savefig('Figure4.pdf')

plt.show()

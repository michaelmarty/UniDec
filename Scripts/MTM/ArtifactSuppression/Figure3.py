import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[3] / "public" / "UniDec"
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
})


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


def run_deconvolution(spectrum, massrange, sat=0, topx=0, psig=10):
    eng = ud.UniDec()
    print(massrange)
    eng.pass_data_in(spectrum)
    eng.config.endz = 100
    eng.config.psig = psig
    eng.config.subbuff = 0
    eng.config.masslb = massrange[0]
    eng.config.massub = massrange[1]
    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = 0
    eng.config.suppression_topn = 0
    eng.config.suppression_topx = topx
    eng.config.beta = 0
    eng.autorun(auto_peak_width=True)

    return eng.data.massdat


np.random.seed(42)

masses = [150000, 600000]
resolutions = [150, 250]
spectra = []
default_mass_datasets = []
top_05_mass_datasets = []
top_09_mass_datasets = []
satellite_2_mass_datasets = []

for i, mass in enumerate(masses):
    spectrum, params = msb.simple_spectrum2(
        [mass],
        resolution=resolutions[i],
        psfun=0,
        baseline=0.03,
        noise=0.01,
        zwidth=1,
    )
    spectrum, ztab = spectrum
    massrange = [mass * 0.25, mass * 2.5]

    spectra.append(spectrum)
    default_mass_datasets.append(run_deconvolution(spectrum, massrange))
    top_05_mass_datasets.append(
        run_deconvolution(spectrum, massrange, topx=0.5)
    )
    top_09_mass_datasets.append(
        run_deconvolution(spectrum, massrange, topx=0.9)
    )
    satellite_2_mass_datasets.append(
        run_deconvolution(spectrum, massrange, sat=2)
    )

fig = plt.figure(figsize=(6.5, 8))
outer_grid = fig.add_gridspec(
    2,
    1,
    height_ratios=[1, 4],
    hspace=0.15,
)
spectrum_grid = outer_grid[0].subgridspec(1, 2, wspace=0.25)
mass_grid = outer_grid[1].subgridspec(4, 2, hspace=0.08, wspace=0.25)
axes = np.empty((5, 2), dtype=object)
for column in range(2):
    axes[0, column] = fig.add_subplot(spectrum_grid[0, column])
    for row in range(1, 5):
        axes[row, column] = fig.add_subplot(mass_grid[row - 1, column])

for column, datasets in enumerate(
        zip(
            masses,
            spectra,
            default_mass_datasets,
            top_05_mass_datasets,
            top_09_mass_datasets,
            satellite_2_mass_datasets,
        )
):
    (
        mass,
        spectrum,
        default_data,
        top_05_data,
        top_09_data,
        satellite_2_data,
    ) = datasets
    axes[0, column].plot(spectrum[:, 0], normalize(spectrum), color="darkblue")
    axes[1, column].plot(
        default_data[:, 0] / 1000, normalize(default_data), color="darkblue"
    )
    axes[2, column].plot(
        top_05_data[:, 0] / 1000, normalize(top_05_data), color="darkblue"
    )
    axes[3, column].plot(
        top_09_data[:, 0] / 1000, normalize(top_09_data), color="darkblue"
    )
    axes[4, column].plot(
        satellite_2_data[:, 0] / 1000, normalize(satellite_2_data), color="darkblue"
    )

    axes[0, column].set_title(f'{mass / 1000:g} kDa')

    axes[0, column].set_xlabel(r'$\mathit{m/z}$')
    axes[-1, column].set_xlabel('Mass (kDa)')

row_labels = [
    'Data',
    'Normal',
    'Top 0.5',
    'Top 0.9',
    'Satellite 2',
]
for row, label in enumerate(row_labels):
    for axis in axes[row]:
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

for axis in axes[:, 1:].flat:
    axis.tick_params(axis='y', labelleft=False)

for axis in axes[1:-1].flat:
    axis.tick_params(axis='x', labelbottom=False)
    axis.set_yticks([0, 50, 100], labels=['', '%', '100'])

for column, (mass, spectrum) in enumerate(zip(masses, spectra)):
    label_most_intense_charge(axes[0, column], spectrum, mass)

ypos = 90
red_bar_panels = [
    (1, 0),
    (1, 1),
    (2, 1),
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

fig.subplots_adjust(left=0.1, right=0.98, bottom=0.08, top=0.96)

os.chdir(
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)
plt.savefig('Figure3.png', dpi=300)
plt.savefig('Figure3.pdf')

plt.show()

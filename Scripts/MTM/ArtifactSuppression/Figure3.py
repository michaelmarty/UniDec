import os

import matplotlib
import numpy as np
import unidec.engine as ud

matplotlib.use('WxAgg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
})

DATA_FILE = (
    r"POPC_Nanodiscs.txt"
)
FIGURE_DIR = (
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)

os.chdir(FIGURE_DIR)


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


def run_deconvolution(file_path, sat=0, topx=0, topn=0):
    eng = ud.UniDec()
    eng.open_file(file_path, refresh=False)
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


spectrum, default_mass_data = run_deconvolution(DATA_FILE)
_, top_01_mass_data = run_deconvolution(DATA_FILE, topx=0.1)
_, top_05_mass_data = run_deconvolution(DATA_FILE, topx=0.5)
_, top_2_mass_data = run_deconvolution(DATA_FILE, topn=2)
_, satellite_2_mass_data = run_deconvolution(DATA_FILE, sat=2)

dominant_mass = default_mass_data[np.argmax(default_mass_data[:, 1]), 0]

fig, axes = plt.subplots(6, 1, figsize=(3.5, 9.5))

mass_datasets = [
    default_mass_data,
    top_01_mass_data,
    top_05_mass_data,
    top_2_mass_data,
    satellite_2_mass_data,
]
row_labels = [
    'Data',
    'Default',
    'Top 0.1',
    'Top 0.5',
    'Top 2',
    'Satellite 2',
]

axes[0].plot(spectrum[:, 0], normalize(spectrum))
axes[0].set_title('Nanodiscs')
axes[0].set_xlabel(r'$\mathit{m/z}$')

for row, mass_data in enumerate(mass_datasets, start=1):
    axes[row].plot(mass_data[:, 0] / 1000, normalize(mass_data))

axes[-1].set_xlabel('Mass (kDa)')

for axis, label in zip(axes, row_labels):
    axis.text(
        0.95,
        0.95,
        label,
        transform=axis.transAxes,
        va='top',
        ha='right',
    )

for panel_index, axis in enumerate(axes):
    axis.text(
        0.03,
        0.95,
        chr(ord('A') + panel_index),
        transform=axis.transAxes,
        va='top',
        ha='left',
    )
    clean_axes(axis)

label_most_intense_charge(axes[0], spectrum, dominant_mass)

ypos = 20
bar_set = [1, 5]
for axis in axes[bar_set]:
    axis.hlines(
        ypos, 187, 200,
        color='orange', linewidth=3
    )

plt.tight_layout()

os.chdir(FIGURE_DIR)
plt.savefig('Figure3.png', dpi=300)
plt.savefig('Figure3.pdf')

plt.show()

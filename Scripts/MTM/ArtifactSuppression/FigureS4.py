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

FIGURE_DIR = Path(
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)

POPC_CONFIG = {
    'minmz': 8000,
    'maxmz': 13000,
    'molig': 760,
    'msig': 1,
    'psig': 1,
    'mzsig': 6,
    'masslb': 90000,
    'massub': 200000,
    'endz': 16,
    'startz': 7,
}

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


def run_default_deconvolution(file_path, config_options):
    engine = ud.UniDec()
    engine.open_file(file_path, refresh=True)
    for option, value in config_options.items():
        setattr(engine.config, option, value)

    engine.config.suppression_satellite = 0
    engine.config.suppression_harmonic = 0
    engine.config.suppression_topn = 0
    engine.config.suppression_topx = 0
    engine.config.beta = 0
    engine.autorun(auto_peak_width=True)
    return engine.data.mzgrid.copy(), engine.config.cmap


def plot_grid(axis, mzgrid, cmap, panel, title):
    mz = np.unique(mzgrid[:, 0])
    charge = np.unique(mzgrid[:, 1])
    intensity = mzgrid[:, 2].reshape(len(mz), len(charge))
    intensity = intensity / np.max(intensity) * 100

    contour = axis.contourf(
        mz,
        charge,
        intensity.T,
        levels=np.linspace(0, 100, 101),
        cmap=cmap,
    )
    axis.set_xlabel(r'$\mathit{m/z}$')
    axis.set_ylabel('Charge')
    axis.text(
        0.02, 0.95, panel,
        transform=axis.transAxes,
        va='top',
        ha='left',
        color='white',
    )
    axis.text(
        0.98, 0.95, title,
        transform=axis.transAxes,
        va='top',
        ha='right',
        color='white',
    )
    return contour



if __name__ == '__main__':
    datasets = [
        (
            'A',
            'Nanodiscs',
            FIGURE_DIR / 'POPC_Nanodiscs.txt',
            POPC_CONFIG,
        ),
        (
            'B',
            'GroEL SID',
            FIGURE_DIR / 'GroEL_SID1.RAW',
            GROEL_CONFIG,
        ),
    ]

    results = [
        (*dataset[:2], *run_default_deconvolution(dataset[2], dataset[3]))
        for dataset in datasets
    ]

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7))
    for axis, (panel, title, mzgrid, cmap) in zip(axes, results):
        contour = plot_grid(axis, mzgrid, cmap, panel, title)
        colorbar = fig.colorbar(contour, ax=axis, ticks=[0, 50, 100])
        colorbar.set_label('Relative intensity (%)')

    fig.subplots_adjust(left=0.11, right=0.93, bottom=0.08, top=0.98, hspace=0.28)
    fig.savefig(FIGURE_DIR / 'FigureS4.png', dpi=300)
    fig.savefig(FIGURE_DIR / 'FigureS4.pdf')

    plt.show()


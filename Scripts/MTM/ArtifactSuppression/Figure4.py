import os

import matplotlib
import numpy as np
from unidec.modules import CDEng

matplotlib.use('WxAgg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
})

FIGURE_DIR = (
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)
DATA_FILE = os.path.join(FIGURE_DIR, "GroEL_CDMS_1.RAW")

# Analysis settings from GroEL_CDMS_1_unidecfiles/conf.dat.
GROEL_CDMS_CONFIG = {
    'numit': 50,
    'startz': 50,
    'endz': 100,
    'zzsig': 1e-6,
    'psig': 5.0,
    'mzsig': 5.0,
    'csig': 5.0,
    'psfun': 2,
    'psfunz': 2,
    'rawflag': 0,
    'masslb': 600000.0,
    'massub': 1000000.0,
    'massbins': 10.0,
    'minmz': 9992.689903846154,
    'maxmz': 12573.891826923078,
    'subbuff': 0.0,
    'subtype': 0,
    'smooth': 0.0,
    'mzbins': 10.0,
    'intthresh': 0.0,
    'reductionpercent': 0.0,
    'adductmass': 1.007276467,
    'nativezlb': -100.0,
    'nativezub': 100.0,
    'poolflag': 1,
    'linflag': 2,
    'CDslope': 0.2074,
    'CDzbins': 1,
    'CDres': 0.0,
    'CDScanStart': -1,
    'CDScanEnd': -1,
    'CDprethresh': 0.0,
    'CDScanCompress': 0,
    'smoothdt': 0.0,
    'subbufdt': 0.0,
    'suppression_topn': 0,
    'suppression_percent': 0.0,
    'suppression_startit': 0,
    'suppression_harmonic': 0,
    'suppression_satellite': 0,
}


def normalize(data):
    maximum = np.max(data)
    if maximum == 0:
        return np.zeros_like(data)
    return data / maximum * 100


def apply_config(engine):
    for option, value in GROEL_CDMS_CONFIG.items():
        setattr(engine.config, option, value)


def run_deconvolution(engine, beta=0, topx=0):
    engine.config.beta = beta
    engine.config.suppression_topx = topx
    engine.run_deconvolution()
    return engine.harray.copy(), engine.data.massdat.copy()


engine = CDEng.UniDecCD()
engine.open_file(DATA_FILE)

# Keep generated files beside the selected RAW file even if a copied config
# contains paths from another machine or working directory.
engine.config.outfname = os.path.join(engine.config.udir, "")
engine.config.default_file_names(s="")
engine.config.confname = os.path.join(engine.config.udir, "conf.dat")
apply_config(engine)

engine.process_data()
raw_grid = engine.harray.copy()
raw_mass_data = engine.data.massdat.copy()
raw_mz = engine.mz.copy()
raw_charge = engine.ztab.copy()

default_grid, default_mass_data = run_deconvolution(engine, beta=0)
beta_grid, beta_mass_data = run_deconvolution(engine, beta=1)
top_grid, top_mass_data = run_deconvolution(engine, beta=0, topx=0.1)

rows = [
    ('Raw Data', raw_grid, raw_mass_data),
    ('Default', default_grid, default_mass_data),
    ('Beta 1', beta_grid, beta_mass_data),
    ('Top 0.1', top_grid, top_mass_data),
]

fig, axes = plt.subplots(
    4,
    2,
    figsize=(7, 8),
    gridspec_kw={'width_ratios': [1.6, 1], 'hspace': 0.3, 'wspace': 0.3},
)

for row, (label, grid, mass_data) in enumerate(rows):
    axis_2d = axes[row, 0]
    axis_mass = axes[row, 1]

    contour = axis_2d.contourf(
        raw_mz,
        raw_charge,
        normalize(grid),
        levels=np.linspace(0, 100, 101),
        cmap=engine.config.cmap,
    )
    colorbar = fig.colorbar(contour, ax=axis_2d, ticks=[0, 50, 100])
    colorbar.ax.set_yticklabels(['0', '%', '100'])
    axis_2d.set_xlabel(r'$\mathit{m/z}$')
    axis_2d.set_ylabel('Charge')

    axis_mass.plot(
        mass_data[:, 0] / 1000,
        normalize(mass_data[:, 1]),
        color="darkblue",
    )
    axis_mass.set_xlabel('Mass (kDa)')
    axis_mass.set_ylim(0, 100)
    axis_mass.set_yticks([0, 50, 100], labels=['0', '%', '100'])

    for column, axis in enumerate((axis_2d, axis_mass)):
        label_color = 'white' if column == 0 else 'black'
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.text(
            0.03,
            0.95,
            chr(ord('A') + row * 2 + column),
            transform=axis.transAxes,
            va='top',
            ha='left',
            color=label_color,
        )
        axis.text(
            0.97,
            0.95,
            label,
            transform=axis.transAxes,
            va='top',
            ha='right',
            color=label_color,
        )

for axis in axes[:-1].flat:
    axis.set_xlabel('')
    axis.tick_params(axis='x', labelbottom=False)

fig.subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.98)

plt.savefig(os.path.join(FIGURE_DIR, 'Figure4.png'), dpi=300)
plt.savefig(os.path.join(FIGURE_DIR, 'Figure4.pdf'))
plt.show()

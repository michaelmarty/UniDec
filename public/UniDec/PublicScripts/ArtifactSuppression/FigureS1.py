import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[3] / "public" / "UniDec"
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb

matplotlib.use('WxAgg')
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
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.set_ylim(0, 100)
    axis.set_yticks([0, 50, 100], labels=['0', '%', '100'])


def run_deconvolution(spectrum, massrange=[5000, 50000], sat=0, harm=0,
                      topn=0, topx=0, beta=0, psig=10,
                      baseline_subtraction=0):
    eng = ud.UniDec()
    eng.pass_data_in(spectrum)
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

    return eng.data.massdat


np.random.seed(42)

spectrum, params = msb.simple_spectrum2(
    [20000], resolution=50, psfun=0, baseline=0.03, noise=0.01, zwidth=0.75
)
spectrum, ztab = spectrum

massdat_topn2_sat1 = run_deconvolution(spectrum, topn=2, sat=1)
massdat_topx001_sat1 = run_deconvolution(
    spectrum, topx=0.01, sat=1, harm=1
)
massdat_topn2_psig1 = run_deconvolution(spectrum, topn=2, psig=1)
massdat_topn1_psig1 = run_deconvolution(spectrum, topn=1, psig=1)
massdat_topx005_psig1 = run_deconvolution(spectrum, topx=0.05, psig=1)
massdat_topx01_psig1 = run_deconvolution(spectrum, topx=0.1, psig=1)
massdat_beta_some = run_deconvolution(spectrum, beta=50)
massdat_beta_lots = run_deconvolution(spectrum, beta=500)

fig, axes = plt.subplots(4, 2, figsize=(6.5, 6.4))

mass_datasets = [
    massdat_topn2_sat1,
    massdat_topx001_sat1,
    massdat_topn2_psig1,
    massdat_topn1_psig1,
    massdat_topx005_psig1,
    massdat_topx01_psig1,
    massdat_beta_some,
    massdat_beta_lots,
]
plot_labels = [
    'Top 2 + Satellite 1',
    'Top 0.01 + Satellite 1\n+ Harmonics',
    'Top 2 + PS 1',
    'Top 1 + PS 1',
    'Top 0.05 + PS 1',
    'Top 0.1 + PS 1',
    'SoftMax 50',
    'SoftMax 500',
]

for panel_index, (axis, data, label) in enumerate(
        zip(axes.flat, mass_datasets, plot_labels)
):
    axis.plot(data[:, 0] / 1000, normalize(data), color="darkblue")
    axis.text(0.03, 0.95, chr(ord('A') + panel_index),
              transform=axis.transAxes, va='top', ha='left')
    axis.text(0.95, 0.95, label, transform=axis.transAxes,
              va='top', ha='right')
    clean_axes(axis)

for axis in axes[:, 1:].flat:
    axis.tick_params(axis='y', labelleft=False)

ypos = 15
for axis in (axes[1, 0], axes[2, 0]):
    axis.hlines(ypos, 16, 19, color='red', linewidth=3)
    axis.hlines(ypos, 21, 24, color='red', linewidth=3)

axes[2, 0].hlines(ypos, 38.5, 41.5, color='magenta', linewidth=3)
axes[0, 0].hlines(ypos, 38.5, 41.5, color='magenta', linewidth=3)

axes[-1, 0].set_xlabel('Mass (kDa)')
axes[-1, 1].set_xlabel('Mass (kDa)')

plt.tight_layout()

os.chdir(
    r"C:\Users\mm96978\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Artifact Suppression\Figures"
)
plt.savefig("FigureS1.png", dpi=300)
plt.savefig("FigureS1.pdf")

plt.show()

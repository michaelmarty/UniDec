import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb
import matplotlib
import numpy as np

matplotlib.use('WxAgg')
import matplotlib.pyplot as plt


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


def run_deconvolution(spectrum, massrange=[5000, 50000], sat=0, harm=0, topn=0, topx=0,
                      baseline_subtraction=0):
    eng = ud.UniDec()
    eng.pass_data_in(spectrum)
    # eng.config.mzsig = 5
    eng.config.psfun = 1
    eng.config.subbuff = baseline_subtraction
    eng.config.masslb = massrange[0]
    eng.config.massub = massrange[1]
    eng.config.suppression_satellite = sat
    eng.config.suppression_harmonic = harm
    eng.config.suppression_topn = topn
    eng.config.suppression_topx = topx
    eng.autorun(auto_peak_width=True)

    return eng.data.massdat

# Set noise seed for reproducibility
np.random.seed(42)

# Create a sipmle spectrum with a mass of 20000 and a resolution of 500 plot it
spectrum, params = msb.simple_spectrum2([20000], resolution=50, psfun=0, baseline=0.03, noise=0.01, zwidth=1)

spectrum, ztab = spectrum

massdat = run_deconvolution(spectrum)
massdat_sat = run_deconvolution(spectrum, sat=1)
massdat_harm = run_deconvolution(spectrum, harm=1)
massdat_topn = run_deconvolution(spectrum, topn=1)
massdat_topx = run_deconvolution(spectrum, topx=0.5)

fig, axes = plt.subplots(6, 1, figsize=(5, 12))

axes[0].plot(spectrum[:, 0], normalize(spectrum))
axes[0].set_xlabel('m/z')

mass_datasets = [massdat, massdat_sat, massdat_harm, massdat_topn, massdat_topx]
for axis, data in zip(axes[1:], mass_datasets):
    axis.plot(data[:, 0] / 1000, normalize(data))

ypos = 15
for axis in (axes[1], axes[3]):
    axis.hlines(ypos, 16, 19, color='red', linewidth=3)
    axis.hlines(ypos, 21, 24, color='red', linewidth=3)

for axis in (axes[1], axes[2]):
    axis.hlines(ypos, 38.5, 41.5, color='magenta', linewidth=3)

plot_labels = [
    'Raw Data',
    'Default Settings',
    'Satellite 1',
    'Harmonic',
    'Top 1',
    'Top 0.5',
]
for axis, label in zip(axes, plot_labels):
    axis.text(0.95, 0.95, label, transform=axis.transAxes, va='top', ha='right')

axes[-1].set_xlabel('Mass (kDa)')
fig.supylabel('Relative Intensity (%)')

for axis in axes:
    clean_axes(axis)

plt.tight_layout()
plt.show()

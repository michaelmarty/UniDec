import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb

matplotlib.use(os.environ.get("MPLBACKEND", "WxAgg"))
import matplotlib.pyplot as plt


np.random.seed(42)
(spectrum, _), _ = msb.simple_spectrum2(
    [20000], resolution=50, psfun=0, baseline=0.03, noise=0.01, zwidth=0.75
)

engine = ud.UniDec()
engine.pass_data_in(spectrum)
engine.config.psig = 10
engine.config.subbuff = 0
engine.config.masslb = 5000
engine.config.massub = 50000
engine.config.suppression_topn = 1
engine.config.beta = 0
engine.config.rescale_newblur = 1
# engine.config.mzsig = 30
# engine.config.endz=11
# engine.config.startz =5
engine.autorun(auto_peak_width=True)


peak = max(engine.pks.peaks, key=lambda item: item.height)
print(f"Most abundant peak: {peak.mass:.2f} Da")
print(f"UScore: {peak.uscore:.6f}")
print(f"MScore: {peak.mscore:.6f}")

mass = np.asarray(peak.zstack[0], dtype=float)
mass_by_charge = np.asarray(peak.zstack[1:], dtype=float)
assert len(peak.mzstack) == len(engine.data.ztab) == len(mass_by_charge)
summed_mass = np.sum(mass_by_charge, axis=0)
summed_mass /= np.max(summed_mass)

charge_weights = np.array([
    np.sum(np.asarray(stack[:, 2], dtype=float)) for stack in peak.mzstack
])
# shortcut: hide charge states below 1% of maximum fit area, lower this threshold to inspect trace states.
active = charge_weights > np.max(charge_weights) * 0.01
colors = plt.cm.viridis(np.linspace(0, 1, np.count_nonzero(active)))
mz_values = np.unique(engine.data.mzgrid[:, 0])
charge_values = np.unique(engine.data.mzgrid[:, 1])
grid_intensity = engine.data.mzgrid[:, 2].reshape(len(mz_values), len(charge_values))

fig, ((mz_axis, mass_axis), (fit_axis, grid_axis)) = plt.subplots(2, 2, figsize=(12, 10))
for color, index in zip(colors, np.flatnonzero(active)):
    charge = engine.data.ztab[index]

    mz_stack = np.asarray(peak.mzstack[index], dtype=float)
    mz = mz_stack[:, 0]
    data = mz_stack[:, 1]
    fit = mz_stack[:, 2]
    if engine.config.orbimode == 1:
        fit *= charge
    uscore = 1 - np.sum(np.abs(data - fit)) / np.sum(data)

    mass_data = mass_by_charge[index]
    mass_fit = summed_mass * np.sum(mass_data) / np.sum(summed_mass)
    mscore = 1 - np.sum(np.abs(mass_data - mass_fit)) / np.sum(mass_fit)

    print(f"Charge {charge:g}+: UScore={uscore:.6f}, MScore={mscore:.6f}")
    mz_axis.plot(mz, data, color=color, alpha=0.55)
    mz_axis.plot(mz, fit, color=color, linestyle="--", label=f"{charge:g}+")
    mass_axis.plot(mass, mass_data, color=color, alpha=0.55)
    mass_axis.plot(mass, mass_fit, color=color, linestyle="--", label=f"{charge:g}+")
    fit_axis.plot(mz_values, grid_intensity[:, index], color=color, label=f"{charge:g}+")

mz_axis.set(
    xlabel=r"$m/z$",
    ylabel="Intensity",
    title=f"UScore = {peak.uscore:.4f}\nSolid: data; dashed: charge-state fit",
)
mass_axis.set(
    xlabel="Mass (Da)",
    ylabel="Intensity",
    title=f"MScore = {peak.mscore:.4f}\nSolid: charge-state data; dashed: scaled summed fit",
)
fit_axis.plot(engine.data.data2[:, 0], engine.data.data2[:, 1], label="data2", color="black")
fit_axis.plot(engine.data.data2[:, 0], engine.data.fitdat, label="fitdat", color="red")
fit_axis.plot(mz_values, np.sum(grid_intensity, axis=1), "--", color="blue", label="charge sum")
fit_axis.set(xlabel=r"$m/z$", ylabel="Intensity", title="Full data, fit, and charge-state rows")

grid_plot = grid_axis.pcolormesh(
    mz_values, charge_values, grid_intensity.T, shading="auto", cmap=engine.config.cmap
)
grid_axis.set(xlabel=r"$m/z$", ylabel="Charge", title=r"$m/z$ vs. charge grid")
active_charges = engine.data.ztab[active]
grid_axis.set_ylim(active_charges[0] - 0.5, active_charges[-1] + 0.5)
fig.colorbar(grid_plot, ax=grid_axis, label="Intensity")

for axis in (mz_axis, mass_axis, fit_axis, grid_axis):
    axis.spines[["top", "right"]].set_visible(False)
for axis in (mz_axis, mass_axis):
    axis.legend(title="Charge", fontsize=8)
fit_axis.legend()

fig.suptitle(f"Top 1 suppression: most abundant peak at {peak.mass:.2f} Da")
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".png"), dpi=300)
plt.show()

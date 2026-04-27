import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# ----------------------------
# USER SETTINGS
# ----------------------------

# 4 text files
files = [
r"Z:\Group Share\Estevan\Q Exactive Data\20250902\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_unidecfiles\msp1d1_minus_1_to_80_DMPC_CDMS_2025-09-02-01-20-17_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251106\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_CDMS_2025-11-06-11-40-25_unidecfiles\20251106_EC_D1sfCheryshort_1to80_DMPC_rep2_numba2_CDMS_2025-11-06-11-40-25_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20251105\20251105_EC_D1sfGFP1to80_DMPC_CDMSlowerspray_2025-11-05-11-20-47_unidecfiles\20251105_EC_D1sfGFP1to80_DMPC_CDMSlowerspray_2025-11-05-11-20-47_mass.txt",
r"Z:\Group Share\Estevan\Q Exactive Data\20251111\20251111_EC_D1mblueberry_1to80_DMPC_CDMS_2025-11-11-11-03-38_unidecfiles\20251111_EC_D1mblueberry_1to80_DMPC_CDMS_2025-11-11-11-03-38_mass.txt"
]

# Base masses (kDa) for each subplot
base_masses = [44.088, 95.492, 97.722, 97.486]   # <-- EDIT THESE

# Colors for each replicate
colors = ["C0", "C1", "C2"]

# ----------------------------
# PROCESSING + PLOTTING
# ----------------------------

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for i, ax in enumerate(axes):

    file = files[i]
    base_mass = base_masses[i]

    # Load file (A4/B4 = skip first 3 rows)
    df = pd.read_csv(file, sep=None, engine="python", skiprows=3, header=None)
    df.columns = ["x", "y"]

    # Assume triplicates: split into 3 equal chunks
    n = len(df) // 3
    reps = [df.iloc[j*n:(j+1)*n] for j in range(3)]

    # ----------------------------
    # PLOT EACH REPLICATE
    # ----------------------------
    peak_positions = []
    fwhm_values = []

    for r, rep in enumerate(reps):

        x = rep["x"].values
        y = rep["y"].values

        # Find peak
        peak_idx = np.argmax(y)
        peak_x = x[peak_idx]
        peak_positions.append(peak_x)

        # Compute FWHM
        half_max = y[peak_idx] / 2
        left = np.where(y[:peak_idx] < half_max)[0]
        right = np.where(y[peak_idx:] < half_max)[0]

        if len(left) > 0 and len(right) > 0:
            left_idx = left[-1]
            right_idx = right[0] + peak_idx
            fwhm = x[right_idx] - x[left_idx]
        else:
            fwhm = np.nan

        fwhm_values.append(fwhm)

        ax.plot(x, y, color=colors[r], label=f"Rep {r+1}")

    # ----------------------------
    # USE AVERAGE PEAK + FWHM
    # ----------------------------
    peak_avg = np.nanmean(peak_positions)
    fwhm_avg = np.nanmean(fwhm_values)

    # ----------------------------
    # COMPUTE REMAINDER
    # remainder = (peak - base_mass) / 678
    # ----------------------------
    remainder_kDa = (peak_avg - base_mass) / 678

    # Build remainder string
    remainder_str = (
        f"Remainder: {remainder_kDa:.2f} ± {int(round(fwhm_avg))} kDa"
    )

    # Add text annotation
    ax.text(
        0.98, 0.85, remainder_str,
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=9,
        color="black"
    )

    # Subplot title
    ax.set_title(f"File: {file}\nBase mass: {base_mass} kDa")
    ax.set_xlabel("Mass (kDa)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend()

plt.tight_layout()
plt.show()

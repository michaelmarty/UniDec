# plot_stacked_excitation.py
# Requirements: pandas, numpy, matplotlib
# Paste your excitation Excel file paths into `sample_files` below and run.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# -----------------------
# rcParams (match your reference)
# -----------------------
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12
rcParams['font.weight'] = 200
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 700
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

# -----------------------
# USER: Paste your excitation file paths here
# Each dict must have keys: "label", "exc"
# Example:
# sample_files = [
#   {"label":"Rep1", "exc": r"C:\path\Rep1_exc.xlsx"},
#   {"label":"Rep2", "exc": r"C:\path\Rep2_exc.xlsx"},
# ]
sample_files = [
    {"label":"MSP1D1-mBlueberry2", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mBlueberry2\MSP1D1_mBlueberry2_1uM_excitation.xlsx"},
    {"label":"MSP1D1-mCherry", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mCherry\MSP1D1_mCherry_1uM_excitation.xlsx"},
    # add more...
]

# Output
outdir = os.path.join(os.getcwd(), "excitation_plots")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "stacked_excitation_spectra.pdf")

# Plot settings
color = "#2C7FB8"
peak_threshold = 0.10  # fraction of max to consider peaks
ymax_plot = 1.08

# -----------------------
# Helper functions
# -----------------------
def read_two_col_excel(path):
    """
    Reads Excel expected to contain wavelength and intensity columns.
    Tries to detect columns named like 'Wavelength'/'nm' and 'Intensity'/'a.u.'.
    Falls back to first two numeric columns.
    Returns (wavelengths, intensities) numpy arrays.
    """
    df = pd.read_excel(path)
    cols = list(df.columns)
    wl_col = None
    inten_col = None
    for c in cols:
        lc = str(c).lower()
        if 'wavelength' in lc or 'nm' in lc:
            wl_col = c
        if 'intensity' in lc or 'mau' in lc or 'a.u.' in lc or 'au' in lc:
            inten_col = c
    if wl_col is None or inten_col is None:
        numeric_cols = [c for c in cols if np.issubdtype(df[c].dtype, np.number)]
        if len(numeric_cols) >= 2:
            wl_col, inten_col = numeric_cols[0], numeric_cols[1]
        else:
            raise ValueError(f"Could not detect wavelength/intensity columns in {path}. Columns: {cols}")
    df[wl_col] = pd.to_numeric(df[wl_col], errors='coerce')
    df[inten_col] = pd.to_numeric(df[inten_col], errors='coerce')
    df = df.dropna(subset=[wl_col, inten_col])
    return df[wl_col].values.astype(float), df[inten_col].values.astype(float)

def normalize_by_abs_max(y):
    """Normalize by absolute maximum: $$I_{\\text{norm}} = I / \\max(|I|)$$"""
    y = np.array(y, dtype=float)
    m = np.nanmax(np.abs(y))
    if m == 0 or np.isnan(m):
        return y
    return y / m

def find_peaks_simple(x, y, min_rel_height=0.10):
    """
    Simple local-maxima detector: point greater than both neighbors + relative height threshold.
    Returns indices of peaks.
    """
    y = np.array(y, dtype=float)
    if y.size < 3:
        return np.array([], dtype=int)
    left = y > np.roll(y, 1)
    right = y > np.roll(y, -1)
    local_max = left & right
    local_max[0] = False
    local_max[-1] = False
    y_max = np.nanmax(np.abs(y))
    if y_max == 0 or np.isnan(y_max):
        return np.array([], dtype=int)
    thresh = min_rel_height * y_max
    peaks_idx = np.where(local_max & (y >= thresh))[0]
    return peaks_idx

# -----------------------
# Validate files exist
# -----------------------
missing = []
for s in sample_files:
    p = s.get('exc')
    if not p or not os.path.exists(p):
        missing.append((s.get('label'), p))
if missing:
    msg = "Missing excitation files:\n" + "\n".join([f"{m[0]}: {m[1]}" for m in missing])
    raise FileNotFoundError(msg)

# -----------------------
# Plot stacked excitation spectra (one row per sample)
# -----------------------
n = len(sample_files)
fig, axes = plt.subplots(n, 1, figsize=(8, 2.7*n), dpi=300, sharex=False)
if n == 1:
    axes = [axes]

for i, s in enumerate(sample_files):
    label = s.get('label', f"Sample {i+1}")
    path = s['exc']
    wl, inten = read_two_col_excel(path)
    if wl.size == 0:
        print(f"Warning: no data in {path}. Skipping.")
        continue
    sort_idx = np.argsort(wl)
    wl = wl[sort_idx]; inten = inten[sort_idx]
    y = normalize_by_abs_max(inten)
    ax = axes[i]
    ax.plot(wl, y, color=color, linewidth=1.8, label=label)
    ax.fill_between(wl, y, 0, color=color, alpha=0.9, linewidth=0)
    # Peak detect and annotate tallest peak above threshold
    peaks = find_peaks_simple(wl, y, min_rel_height=peak_threshold)
    if len(peaks) > 0:
        tallest = peaks[np.argmax(y[peaks])]
        px = wl[tallest]; py = y[tallest]
        ax.plot(px, py, marker='o', color='k', markersize=4)
        ax.axvline(px, color=color, linestyle='--', linewidth=0.8, alpha=0.6)
        ax.text(px, py + 0.04, f"{px:.1f} nm", ha='center', va='bottom', fontsize=9,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    # Formatting (similar to your reference)
    ax.set_ylim(0, ymax_plot)
    ax.margins(y=0)
    ax.set_xlim(np.min(wl), np.max(wl))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([0, 1])
    ax.tick_params(axis="y", length=3, width=0.8)
    ax.tick_params(axis="x", length=3, width=0.8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.5)
    if i == 0:
        ax.set_title("Fluorescence Excitation Spectra")

# Shared labels
fig.text(0.04, 0.5, "Normalized intensity", va="center", rotation="vertical")
axes[-1].set_xlabel("Wavelength (nm)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
plt.savefig(outfile)
print(f"Saved stacked excitation plot to: {outfile}")
plt.show()
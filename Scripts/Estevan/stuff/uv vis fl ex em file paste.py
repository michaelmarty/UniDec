# paste_paths_plot_spectra.py
# Requirements: pandas, numpy, matplotlib
# Edit the sample_files list below to include your Excel file paths,
# then run this script.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# -----------------------
# Styling (match your reference rcParams)
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
# USER: Paste your file paths here (Option A)
# Each sample dict must include keys: "label", "uv", "exc", "em"
# Example:
# sample_files = [
#     {"label":"Rep1",
#      "uv": r"C:\path\Rep1_uv.xlsx",
#      "exc": r"C:\path\Rep1_exc.xlsx",
#      "em" : r"C:\path\Rep1_em.xlsx"},
#     {"label":"Rep2",
#      "uv": r"C:\path\Rep2_uv.xlsx",
#      "exc": r"C:\path\Rep2_exc.xlsx",
#      "em" : r"C:\path\Rep2_em.xlsx"},
# ]
# Replace the example entries below with your actual file paths.
sample_files = [
    {"label":"Rep1",
     "uv": r"C:\path\to\Rep1_uv.xlsx",
     "exc": r"C:\path\to\Rep1_exc.xlsx",
     "em" : r"C:\path\to\Rep1_em.xlsx"},
    {"label":"Rep2",
     "uv": r"C:\path\to\Rep2_uv.xlsx",
     "exc": r"C:\path\to\Rep2_exc.xlsx",
     "em" : r"C:\path\to\Rep2_em.xlsx"},
    # Add more samples as needed...
]

# -----------------------
# Output settings
outdir = os.path.join(os.getcwd(), "spectra_plots")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "samples_uv_exc_em_layout.pdf")

# Plot settings
panel_titles = ['UV-Vis', 'Fluorescence Excitation', 'Fluorescence Emission']
colors = ["#333333", "#C76BA0", "#2C7FB8"]  # column colors
peak_threshold = 0.10  # fraction of max for peak detection
ymax_plot = 1.08

# -----------------------
# Helper functions
def read_two_col_excel(path):
    """
    Reads an Excel file expected to have headers like:
      'Wavelength (nm)' and 'Intensity (a.u.)'
    If header names differ, attempts to infer the first two numeric columns.
    Returns (wavelengths, intensity) as numpy arrays.
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
    """
    Normalize by absolute maximum:
    $$I_{\\text{norm}}=\\frac{I}{\\max(|I|)}$$
    """
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
missing = []
for s in sample_files:
    for key in ('uv','exc','em'):
        path = s.get(key)
        if not path or not os.path.exists(path):
            missing.append((s.get('label'), key, path))
if missing:
    msg = "Missing files detected:\n" + "\n".join([f"Sample '{m[0]}', key '{m[1]}': {m[2]}" for m in missing])
    raise FileNotFoundError(msg)

# -----------------------
# Plotting layout
n_samples = len(sample_files)
fig, axes = plt.subplots(n_samples, 3, figsize=(10, 3.0*n_samples), dpi=300, sharex='col')
if n_samples == 1:
    axes = axes[np.newaxis, :]

for i, s in enumerate(sample_files):
    label = s.get("label", f"Sample {i+1}")
    for j, key in enumerate(['uv', 'exc', 'em']):
        ax = axes[i, j]
        path = s[key]
        wl, inten = read_two_col_excel(path)
        if wl.size == 0:
            continue
        # Sort by wavelength
        sort_idx = np.argsort(wl)
        wl = wl[sort_idx]
        inten = inten[sort_idx]
        y = normalize_by_abs_max(inten)
        ax.plot(wl, y, color=colors[j % len(colors)], linewidth=1.8, label=panel_titles[j])
        ax.fill_between(wl, y, 0, color=colors[j % len(colors)], alpha=0.9, linewidth=0)
        # Peak detection: annotate tallest peak above threshold
        peaks = find_peaks_simple(wl, y, min_rel_height=peak_threshold)
        if len(peaks) > 0:
            tallest = peaks[np.argmax(y[peaks])]
            px = wl[tallest]; py = y[tallest]
            ax.plot(px, py, marker='o', color='k', markersize=4)
            ax.axvline(px, color=colors[j % len(colors)], linestyle='--', linewidth=0.8, alpha=0.6)
            ax.text(px, py + 0.04, f"{px:.1f} nm", ha='center', va='bottom', fontsize=9,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        # Formatting
        ax.set_ylim(0, ymax_plot)
        ax.margins(y=0)
        ax.set_xlim(np.min(wl), np.max(wl))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_yticks([0, 1])
        ax.tick_params(axis="y", length=3, width=0.8)
        ax.tick_params(axis="x", length=3, width=0.8)
        if j == 0:
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.5)
        if i == 0:
            ax.set_title(panel_titles[j])
    # Row label at left
    axes[i,0].text(-0.12, 0.5, label, va='center', ha='right', transform=axes[i,0].transAxes, fontsize=12, fontweight=600)

# Shared labels
fig.text(0.04, 0.5, "Normalized intensity", va="center", rotation="vertical")
axes[-1,1].set_xlabel("Wavelength (nm)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
plt.savefig(outfile)
print(f"Saved plot to: {outfile}")
plt.show()
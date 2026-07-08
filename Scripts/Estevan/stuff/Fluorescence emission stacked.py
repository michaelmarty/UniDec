# stacked_emission_colors_list_only_strict.py
# Requires: pandas, numpy, matplotlib
# Edit USER SETTINGS (sample_files, colors, GLOBAL_WL_RANGE, per_sample_wl, FILL) and run.

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
# USER SETTINGS (edit)
# -----------------------
# Paste your emission files here. Each dict must include: "label", "em".
sample_files = [
    {"label": "MSP1D1-mBlueberry2 Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mBlueberry2\MSP1D1_mBlueberry2_nanodisc_emission.xlsx"},
    {"label": "MSP1D1-sfGFP Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfGFP\MSP1D1_sfGFP_nanodisc_emission.xlsx"},
    {"label": "MSP1E3D1-sfGFP Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfGFP\MSP1E3D1_sfGFP_nanodisc_emission.xlsx"},
    #{"label": "AnnexinV-sfGFP", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\AnnexinV\Annexin_V_sfGFP_emission.xlsx"},
    {"label": "MSP1D1-mCherry Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mCherry\MSP1D1_mCherry_nanodisc_emission.xlsx"},
    {"label": "MSP1D1-sfCherryS Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryS\MSP1D1_sfCherryS_nanodisc_emission.xlsx"},
    {"label": "MSP1D1-sfCherryL Nanodiscs", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryL\MSP1D1_sfCherryL_nanodisc_emission.xlsx"},
    #{"label": "MSP1E3D1-sfCherry3C", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfCherry3C\MSP1E3D1 sfCherry3C emission.xlsx"},
    # add more...
]

# === CUSTOM HEX COLORS ===
# Edit this list exactly like your example: comment/uncomment hex codes as needed.
colors = [
    "#3895D3",
    "#4C9A2A",
    "#4C9A2A",
    #"#4C9A2A",
    "#9043B7",
    "#9043B7",
    "#9043B7",
    #"#9043B7",
]

# Global wavelength range (set to tuple (min,max) to force same x-limits for all, or None to auto-scale)
GLOBAL_WL_RANGE = None  # e.g. (300, 700) or None

# Per-sample wavelength ranges (overrides GLOBAL_WL_RANGE). Keys must match sample "label".
per_sample_wl = {
    # "Sample1": (450, 650),
}

# Fill area under curve? Set False to disable.
FILL = False

# Save per-sample PNG files as well?
SAVE_INDIVIDUAL_PNG = False

# Output directory and filename
outdir = os.path.join(os.getcwd(), "emission_plots")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "stacked_emission_spectra.pdf")

# Plot settings
PEAK_THRESHOLD = 0.10  # fraction of max for peak detection
YMAX_PLOT = 1.08
ANNOTATE_TALLEST_ONLY = True  # set False to annotate all detected peaks

# -----------------------
# Helper functions
# -----------------------
def read_two_col_excel(path):
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
    $$I_{\\text{norm}} = \\dfrac{I}{\\max(|I|)}$$
    """
    y = np.array(y, dtype=float)
    m = np.nanmax(np.abs(y))
    if m == 0 or np.isnan(m):
        return y
    return y / m

def find_peaks_simple(x, y, min_rel_height=0.10):
    """
    Simple local-maxima detector: strictly greater than neighbors + relative height threshold.
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
    p = s.get('em')
    if not p or not os.path.exists(p):
        missing.append((s.get('label'), p))
if missing:
    msg = "Missing emission files:\n" + "\n".join([f"{m[0]}: {m[1]}" for m in missing])
    raise FileNotFoundError(msg)

# -----------------------
# Validate colors list is non-empty after user edits
if not any(isinstance(c, str) and c.strip() for c in colors):
    raise ValueError("The 'colors' list is empty or contains no active hex codes. Please uncomment/add at least one hex color.")

# Use colors exactly as provided (uncommented strings)
color_list = [c for c in colors if isinstance(c, str) and c.strip()]

# -----------------------
# Plot stacked emission spectra (colors only from the list)
n = len(sample_files)
fig, axes = plt.subplots(n, 1, figsize=(8, 2.7*n), dpi=300, sharex=True)
if n == 1:
    axes = [axes]

# Compute global min/max if GLOBAL_WL_RANGE not set
global_min = None
global_max = None
if GLOBAL_WL_RANGE is None:
    wl_mins = []
    wl_maxs = []
    for s in sample_files:
        wl, _ = read_two_col_excel(s['em'])
        if wl.size > 0:
            wl_mins.append(np.nanmin(wl))
            wl_maxs.append(np.nanmax(wl))
    if wl_mins and wl_maxs:
        global_min = min(wl_mins)
        global_max = max(wl_maxs)
else:
    global_min, global_max = GLOBAL_WL_RANGE

for i, s in enumerate(sample_files):
    label = s.get('label', f"Sample {i+1}")
    # Choose color by cycling through the provided colors list only
    sample_color = color_list[i % len(color_list)]

    path = s['em']
    wl, inten = read_two_col_excel(path)
    if wl.size == 0:
        print(f"Warning: no data in {path}. Skipping.")
        continue
    idx = np.argsort(wl)
    wl = wl[idx]; inten = inten[idx]
    y = normalize_by_abs_max(inten)
    ax = axes[i]
    ax.plot(wl, y, color=sample_color, linewidth=1.8, label=label)
    if FILL:
        ax.fill_between(wl, y, 0, color=sample_color, alpha=0.9, linewidth=0)
    # Peak detect and annotate
    peaks = find_peaks_simple(wl, y, min_rel_height=PEAK_THRESHOLD)
    if peaks.size > 0:
        if ANNOTATE_TALLEST_ONLY:
            tallest_idx = peaks[np.argmax(y[peaks])]
            px = wl[tallest_idx]; py = y[tallest_idx]
            ax.plot(px, py, marker='o', color='k', markersize=4)
            ax.axvline(px, color=sample_color, linestyle='--', linewidth=0.8, alpha=0.6)
            ax.text(px, py + 0.04, f"{px:.1f} nm", ha='center', va='bottom', fontsize=9,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        else:
            for p_idx in peaks:
                px = wl[p_idx]; py = y[p_idx]
                ax.plot(px, py, marker='o', color='k', markersize=4)
                ax.axvline(px, color=sample_color, linestyle='--', linewidth=0.6, alpha=0.6)
                ax.text(px, py + 0.03, f"{px:.1f} nm", ha='center', va='bottom', fontsize=8,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    # Determine x-limits: per-sample override > GLOBAL_WL_RANGE > combined data-range > sample data
    if label in per_sample_wl:
        xmin, xmax = per_sample_wl[label]
    elif GLOBAL_WL_RANGE is not None:
        xmin, xmax = GLOBAL_WL_RANGE
    elif global_min is not None and global_max is not None:
        xmin, xmax = global_min, global_max
    else:
        xmin, xmax = np.min(wl), np.max(wl)
    ax.set_xlim(xmin, xmax)
    # Formatting
    ax.set_ylim(0, YMAX_PLOT)
    ax.margins(y=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([0, 1])
    ax.tick_params(axis="y", length=3, width=0.8)
    ax.tick_params(axis="x", length=3, width=0.8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.5)
    if i == 0:
        ax.set_title("Fluorescence Emission Spectra")

# Shared labels
fig.text(0.04, 0.5, "Normalized intensity", va="center", rotation="vertical")
axes[-1].set_xlabel("Wavelength (nm)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
plt.savefig(outfile)
print(f"Saved stacked emission plot to: {outfile}")

# Optionally save individual PNGs
if SAVE_INDIVIDUAL_PNG:
    for i, s in enumerate(sample_files):
        label = s.get('label', f"Sample{i+1}")
        png_path = os.path.join(outdir, f"{label}_emission.png")
        fig_i, ax_i = plt.subplots(figsize=(6,3), dpi=300)
        wl, inten = read_two_col_excel(s['em'])
        if wl.size == 0:
            continue
        idx = np.argsort(wl)
        wl = wl[idx]; inten = inten[idx]
        y = normalize_by_abs_max(inten)
        color_i = color_list[i % len(color_list)]
        ax_i.plot(wl, y, color=color_i, linewidth=1.8, label=label)
        if FILL:
            ax_i.fill_between(wl, y, 0, color=color_i, alpha=0.9, linewidth=0)
        # set x-limits same as stacked plot's for consistency
        ax_i.set_xlim(xmin, xmax)
        ax_i.set_ylim(0, YMAX_PLOT)
        ax_i.set_xlabel("Wavelength (nm)")
        ax_i.set_ylabel("Normalized intensity")
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)
        fig_i.tight_layout()
        fig_i.savefig(png_path)
        plt.close(fig_i)
    print("Saved individual PNGs to:", outdir)

plt.show()
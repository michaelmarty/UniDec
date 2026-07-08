# stacked_uvvis_colors_list_only_no_normalize_with0tick.py
# Requires: pandas, numpy, matplotlib
# Edit USER SETTINGS (sample_files, colors, GLOBAL_WL_RANGE, per_sample_wl, FILL) and run.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator

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
rcParams['axes.linewidth'] = 1
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

# -----------------------
# USER SETTINGS (edit)
# -----------------------
sample_files = [
    {"label": "1", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1D1 mBlueberry2\MSP1D1 mBlueberry2 UVvis Avg.xlsx"},
    {"label": "2", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1D1 sfGFP\MSP1D1 sfGFP UVvis Avg.xlsx"},
    {"label": "3", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1E3D1 sfGFP\MSP1E3D1 sfGFP UVvis Avg.xlsx"},
    {"label": "4", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\Annexin V sfGFP\AnnexinV sfGFP UV Vis Avg.xlsx"},
    {"label": "5", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1D1 mCherry\MSP1D1 mCherry UVvis Avg.xlsx"},
    {"label": "6", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1D1 sfCherryS\MSP1D1 sfCherryS UVvis Avg.xlsx"},
    {"label": "7", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1D1 sfCherryL\MSP1D1 sfCherryL UVvis Avg.xlsx"},
    {"label": "8", "uv": r"D:\UT Summer 2026\FLuoro ND paper revision\UV Vis\UV VIS\MSP1E3D1 sfCherry3C\MSP1E3D1 sfCherry 3C UV Vis Avg.xlsx"},
]

colors = [
    "#3895D3",
    "#4C9A2A",
    "#4C9A2A",
    "#4C9A2A",
    "#9043B7",
    "#9043B7",
    "#9043B7",
    "#9043B7",
]

GLOBAL_WL_RANGE = (250, 650)
per_sample_wl = {}
FILL = False
SAVE_INDIVIDUAL_PNG = False

outdir = os.path.join(os.getcwd(), "uvvis_plots")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "stacked_uvvis_spectra.pdf")

PEAK_THRESHOLD = 0.10
ANNOTATE_TALLEST_ONLY = True

# -----------------------
# Helpers
def read_two_col_excel(path):
    df = pd.read_excel(path)
    cols = list(df.columns)
    wl_col = None
    inten_col = None
    for c in cols:
        lc = str(c).lower()
        if 'wavelength' in lc or 'nm' in lc:
            wl_col = c
        if 'intensity' in lc or 'absorb' in lc or 'a.u.' in lc or 'au' in lc:
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
    y = np.array(y, dtype=float)
    m = np.nanmax(np.abs(y))
    if m == 0 or np.isnan(m):
        return y
    return y / m

def find_peaks_simple(x, y, min_rel_height=0.10):
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

# Ensure files exist and colors valid
missing = []
for s in sample_files:
    p = s.get('uv')
    if not p or not os.path.exists(p):
        missing.append((s.get('label'), p))
if missing:
    msg = "Missing UV-Vis files:\n" + "\n".join([f"{m[0]}: {m[1]}" for m in missing])
    raise FileNotFoundError(msg)

if not any(isinstance(c, str) and c.strip() for c in colors):
    raise ValueError("The 'colors' list is empty or contains no active hex codes.")

color_list = [c for c in colors if isinstance(c, str) and c.strip()]

# -----------------------
# Plot stacked UV-Vis spectra
n = len(sample_files)
fig, axes = plt.subplots(n, 1, figsize=(8, 2.7*n), dpi=300, sharex=True)
if n == 1:
    axes = [axes]

# Global wl range if needed
global_min = None
global_max = None
if GLOBAL_WL_RANGE is None:
    wl_mins, wl_maxs = [], []
    for s in sample_files:
        wl, _ = read_two_col_excel(s['uv'])
        if wl.size > 0:
            wl_mins.append(np.nanmin(wl)); wl_maxs.append(np.nanmax(wl))
    if wl_mins and wl_maxs:
        global_min = min(wl_mins); global_max = max(wl_maxs)
else:
    global_min, global_max = GLOBAL_WL_RANGE

for i, s in enumerate(sample_files):
    label = s.get('label', f"Sample {i+1}")
    sample_color = color_list[i % len(color_list)]

    path = s['uv']
    wl, inten = read_two_col_excel(path)
    if wl.size == 0:
        print(f"Warning: no data in {path}. Skipping.")
        continue
    idx = np.argsort(wl)
    wl = wl[idx]; inten = inten[idx]

    # Apply wavelength cutoff
    if label in per_sample_wl:
        xmin_cut, xmax_cut = per_sample_wl[label]
    elif GLOBAL_WL_RANGE is not None:
        xmin_cut, xmax_cut = GLOBAL_WL_RANGE
    else:
        xmin_cut, xmax_cut = 250, np.max(wl)

    mask = (wl >= xmin_cut) & (wl <= xmax_cut)
    wl_cut = wl[mask]; inten_cut = inten[mask]

    if wl_cut.size == 0:
        print(f"Warning: no data for sample '{label}' in wavelength range {xmin_cut}–{xmax_cut}. Skipping this sample.")
        continue

    # Raw intensity
    y = inten_cut

    ax = axes[i]
    ax.plot(wl_cut, y, color=sample_color, linewidth=1.8, label=label)
    if FILL:
        ax.fill_between(wl_cut, y, 0, color=sample_color, alpha=0.9, linewidth=0)

    # Peak detection & annotate
    peaks = find_peaks_simple(wl_cut, y, min_rel_height=PEAK_THRESHOLD)
    if peaks.size > 0:
        if ANNOTATE_TALLEST_ONLY:
            tallest_idx = peaks[np.argmax(y[peaks])]
            px = wl_cut[tallest_idx]; py = y[tallest_idx]
            ax.plot(px, py, marker='o', color='k', markersize=4)
            ax.axvline(px, color=sample_color, linestyle='--', linewidth=0.8, alpha=0.6)
            ax.text(px, py + 0.04*(np.nanmax(y)-np.nanmin(y)), f"{px:.1f} nm", ha='center', va='bottom', fontsize=9,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        else:
            for p_idx in peaks:
                px = wl_cut[p_idx]; py = y[p_idx]
                ax.plot(px, py, marker='o', color='k', markersize=4)
                ax.axvline(px, color=sample_color, linestyle='--', linewidth=0.6, alpha=0.6)
                ax.text(px, py + 0.03*(np.nanmax(y)-np.nanmin(y)), f"{px:.1f} nm", ha='center', va='bottom', fontsize=8,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))

    # x-limits
    if label in per_sample_wl:
        xmin, xmax = per_sample_wl[label]
    elif GLOBAL_WL_RANGE is not None:
        xmin, xmax = GLOBAL_WL_RANGE
    elif global_min is not None and global_max is not None:
        xmin, xmax = global_min, global_max
    else:
        xmin, xmax = np.min(wl_cut), np.max(wl_cut)
    ax.set_xlim(xmin, xmax)

    # y-limits: start at 0, top = 1.08 * max(y) (per-sample)
    y_maxdata = np.nanmax(y)
    if np.isnan(y_maxdata) or y_maxdata == 0:
        y_top = 0.01
    else:
        y_top = y_maxdata * 1.08
    ax.set_ylim(0.0, y_top)

    # Ensure 0 is included among ticks while keeping nice automatic ticks
    locator = MaxNLocator(nbins=4, prune='both')
    ticks = locator.tick_values(0.0, y_top)
    # include 0 explicitly if not present
    if 0.0 not in ticks:
        ticks = np.concatenate(([0.0], ticks))
    # keep only ticks within [0, y_top]
    ticks = np.array([t for t in ticks if 0.0 <= t <= y_top])
    ticks = np.unique(np.sort(ticks))
    ax.set_yticks(ticks)

    # Formatting
    ax.margins(y=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=3, width=0.8)
    ax.tick_params(axis="x", length=3, width=0.8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.5)
    if i == 0:
        ax.set_title("UV-Vis Absorbance Spectra (raw)")

# Shared labels
fig.text(0.04, 0.5, "Absorbance", va="center", rotation="vertical")
axes[-1].set_xlabel("Wavelength (nm)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
plt.savefig(outfile)
print(f"Saved stacked UV-Vis plot to: {outfile}")

# Individual PNGs
if SAVE_INDIVIDUAL_PNG:
    for i, s in enumerate(sample_files):
        label = s.get('label', f"Sample{i+1}")
        png_path = os.path.join(outdir, f"{label}_uvvis.png")
        fig_i, ax_i = plt.subplots(figsize=(6,3), dpi=300)
        wl, inten = read_two_col_excel(s['uv'])
        if wl.size == 0:
            continue
        idx = np.argsort(wl)
        wl = wl[idx]; inten = inten[idx]

        if label in per_sample_wl:
            xmin_cut, xmax_cut = per_sample_wl[label]
        elif GLOBAL_WL_RANGE is not None:
            xmin_cut, xmax_cut = GLOBAL_WL_RANGE
        else:
            xmin_cut, xmax_cut = 250, np.max(wl)

        mask = (wl >= xmin_cut) & (wl <= xmax_cut)
        wl_cut = wl[mask]; inten_cut = inten[mask]
        if wl_cut.size == 0:
            continue

        y = inten_cut
        color_i = color_list[i % len(color_list)]
        ax_i.plot(wl_cut, y, color=color_i, linewidth=1.8, label=label)
        if FILL:
            ax_i.fill_between(wl_cut, y, 0, color=color_i, alpha=0.9, linewidth=0)

        ax_i.set_xlim(xmin, xmax)

        # y-limits start at 0, top = 1.08 * max(y)
        y_maxdata = np.nanmax(y)
        if np.isnan(y_maxdata) or y_maxdata == 0:
            y_top = 0.01
        else:
            y_top = y_maxdata * 1.08
        ax_i.set_ylim(0.0, y_top)

        # Ensure 0 is a tick
        locator = MaxNLocator(nbins=4, prune='both')
        ticks = locator.tick_values(0.0, y_top)
        if 0.0 not in ticks:
            ticks = np.concatenate(([0.0], ticks))
        ticks = np.array([t for t in ticks if 0.0 <= t <= y_top])
        ticks = np.unique(np.sort(ticks))
        ax_i.set_yticks(ticks)

        ax_i.set_xlabel("Wavelength (nm)")
        ax_i.set_ylabel("Absorbance")
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)
        fig_i.tight_layout()
        fig_i.savefig(png_path)
        plt.close(fig_i)
    print("Saved individual PNGs to:", outdir)

plt.show()
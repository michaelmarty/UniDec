# merged_excitation_emission_stacked.py
# Requires: pandas, numpy, matplotlib
# This script merges your excitation and emission lists (pairing by label)
# and plots excitation (solid) + emission (dashed) on the same subplot per sample.
# Edit USER SETTINGS (exc_samples, em_samples, colors, GLOBAL_WL_RANGE, per_sample_wl, FILL) and run.

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
rcParams['axes.linewidth'] = 1
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

# -----------------------
# USER SETTINGS (edit)
# Paste your excitation entries here (label, exc)
exc_samples = [
    {"label": "i", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mBlueberry2\MSP1D1_mBlueberry2_1uM_excitation.xlsx"},
    {"label": "ii", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfGFP\MSP1D1_sfGFP_1uM_excitation.xlsx"},
{"label": "iii", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfGFP\MSP1E3D1_sfGFP_1uM_excitation.xlsx"},
    {"label": "iv", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\AnnexinV\Annexin_V_sfGFP excitation.xlsx"},
    {"label": "v", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mCherry\MSP1D1_mCherry_1uM_excitation.xlsx"},
    {"label": "vi", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryS\MSP1D1_sfCherryS_1uM_excitation.xlsx"},
    {"label": "vii", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryL\MSP1D1_sfCherryL_1uM_excitation.xlsx"},
{"label": "viii", "exc": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfCherry3C\MSP1E3D1 sfCherry3C excitation.xlsx"},
]

# Paste your emission entries here (label, em)
em_samples = [
    {"label": "i", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mBlueberry2\MSP1D1_mBlueberry2_1uM_emission.xlsx"},
    {"label": "ii", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfGFP\MSP1D1_sfGFP_1uM_emission.xlsx"},
{"label": "iii", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfGFP\MSP1E3D1_sfGFP_1uM_emission.xlsx"},
    {"label": "iv", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\AnnexinV\Annexin_V_sfGFP_emission.xlsx"},
    {"label": "v", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_mCherry\MSP1D1_mCherry_1uM_emission.xlsx"},
    {"label": "vi", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryS\MSP1D1_sfCherryS_1uM_emission.xlsx"},
    {"label": "vii", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1D1_sfCherryL\MSP1D1_sfCherryL_1uM_emission.xlsx"},
{"label": "viii", "em": r"D:\UT Summer 2026\FLuoro ND paper revision\Fluorescence\MSP1E3D1 sfCherry3C\MSP1E3D1 sfCherry3C emission.xlsx"},
]

# === CUSTOM HEX COLORS ===
# Edit this list exactly like your example: comment/uncomment hex codes as needed.
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

# Global wavelength range (set to tuple (min,max) to force same x-limits for all, or None to auto-scale)
GLOBAL_WL_RANGE = None  # e.g. (250, 800) or None

# Per-sample wavelength ranges (overrides GLOBAL_WL_RANGE). Keys must match sample "label".
per_sample_wl = {
    # "MSP1D1-mBlueberry2 Nanodiscs": (300, 650),
}

# Fill area under curve? Set False to disable.
FILL = False

# Save per-sample PNG files as well?
SAVE_INDIVIDUAL_PNG = False

# Output directory and filename
outdir = os.path.join(os.getcwd(), "fluorescence_plots")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "stacked_excitation_and_emission_spectra.pdf")

# Plot settings
PEAK_THRESHOLD = 0.10  # fraction of max for peak detection
YMAX_PLOT = 1.08
ANNOTATE_TALLEST_ONLY = True  # set False to annotate all detected peaks
# Styles
EXC_STYLE = dict(linestyle='-', linewidth=1.8)   # excitation line style
EM_STYLE  = dict(linestyle='--', linewidth=1.6)  # emission line style

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

# -----------------------
# Build dicts and pair samples by label
exc_map = {s['label']: s['exc'] for s in exc_samples}
em_map  = {s['label']: s['em']  for s in em_samples}

# build paired list of samples present in both
paired = []
missing_labels = []
for lbl in sorted(set(list(exc_map.keys()) + list(em_map.keys()))):
    if lbl in exc_map and lbl in em_map:
        paired.append({"label": lbl, "exc": exc_map[lbl], "em": em_map[lbl]})
    else:
        missing_labels.append(lbl)

if missing_labels:
    print("Warning: some labels are not present in both excitation and emission lists. They will be skipped unless you add matching entries:")
    for m in missing_labels:
        print("  -", m)

if not paired:
    raise RuntimeError("No paired samples found. Ensure labels match between exc_samples and em_samples.")

# -----------------------
# Validate colors list non-empty
if not any(isinstance(c, str) and c.strip() for c in colors):
    raise ValueError("The 'colors' list is empty or contains no active hex codes. Please uncomment/add at least one hex color.")
color_list = [c for c in colors if isinstance(c, str) and c.strip()]

# -----------------------
# Plot stacked paired spectra
n = len(paired)
fig, axes = plt.subplots(n, 1, figsize=(9, 2.6*n), dpi=300, sharex=True)
if n == 1:
    axes = [axes]

# Determine combined global xlimits if desired
global_min = None
global_max = None
if GLOBAL_WL_RANGE is None:
    wl_mins = []
    wl_maxs = []
    for s in paired:
        # check both exc and em files
        wl_e, _ = read_two_col_excel(s['exc'])
        wl_m, _ = read_two_col_excel(s['em'])
        if wl_e.size>0:
            wl_mins.append(np.nanmin(wl_e)); wl_maxs.append(np.nanmax(wl_e))
        if wl_m.size>0:
            wl_mins.append(np.nanmin(wl_m)); wl_maxs.append(np.nanmax(wl_m))
    if wl_mins and wl_maxs:
        global_min = min(wl_mins); global_max = max(wl_maxs)
else:
    global_min, global_max = GLOBAL_WL_RANGE

for i, s in enumerate(paired):
    label = s['label']
    color = color_list[i % len(color_list)]
    # read excitation
    wl_exc, inten_exc = read_two_col_excel(s['exc'])
    # read emission
    wl_em, inten_em = read_two_col_excel(s['em'])

    # sort
    if wl_exc.size>0:
        idx_exc = np.argsort(wl_exc); wl_exc = wl_exc[idx_exc]; inten_exc = inten_exc[idx_exc]
    if wl_em.size>0:
        idx_em = np.argsort(wl_em); wl_em = wl_em[idx_em]; inten_em = inten_em[idx_em]

    # apply wavelength cutoff logic (can be different per curve if you want)
    # Here we use per_sample_wl or GLOBAL_WL_RANGE; if neither set, keep full ranges
    if label in per_sample_wl:
        xmin_cut, xmax_cut = per_sample_wl[label]
    elif GLOBAL_WL_RANGE is not None:
        xmin_cut, xmax_cut = GLOBAL_WL_RANGE
    elif global_min is not None and global_max is not None:
        xmin_cut, xmax_cut = global_min, global_max
    else:
        xmin_cut, xmax_cut = None, None

    def apply_cut(wl, inten, xmin, xmax):
        if xmin is None and xmax is None:
            return wl, inten
        mask = np.ones_like(wl, dtype=bool)
        if xmin is not None:
            mask &= (wl >= xmin)
        if xmax is not None:
            mask &= (wl <= xmax)
        return wl[mask], inten[mask]

    wl_exc_cut, inten_exc_cut = apply_cut(wl_exc, inten_exc, xmin_cut, xmax_cut)
    wl_em_cut, inten_em_cut   = apply_cut(wl_em,   inten_em,   xmin_cut, xmax_cut)

    if wl_exc_cut.size == 0 and wl_em_cut.size == 0:
        print(f"Warning: no data in selected wavelength range for sample '{label}'. Skipping.")
        continue

    # normalize each curve independently
    y_exc = normalize_by_abs_max(inten_exc_cut) if wl_exc_cut.size>0 else np.array([])
    y_em  = normalize_by_abs_max(inten_em_cut)  if wl_em_cut.size>0 else np.array([])

    ax = axes[i]
    # plot excitation (solid)
    if wl_exc_cut.size>0:
        ax.plot(wl_exc_cut, y_exc, color=color, label=f"{label} (exc)", **EXC_STYLE)
        if FILL:
            ax.fill_between(wl_exc_cut, y_exc, 0, color=color, alpha=0.12, linewidth=0)
    # plot emission (dashed)
    if wl_em_cut.size>0:
        ax.plot(wl_em_cut, y_em, color=color, label=f"{label} (em)", **EM_STYLE)

    # Peak detection and annotate tallest per curve
    if wl_exc_cut.size>0:
        peaks_exc = find_peaks_simple(wl_exc_cut, y_exc, min_rel_height=PEAK_THRESHOLD)
    else:
        peaks_exc = np.array([], dtype=int)
    if wl_em_cut.size>0:
        peaks_em = find_peaks_simple(wl_em_cut, y_em, min_rel_height=PEAK_THRESHOLD)
    else:
        peaks_em = np.array([], dtype=int)

    if ANNOTATE_TALLEST_ONLY:
        if peaks_exc.size>0:
            t_exc = peaks_exc[np.argmax(y_exc[peaks_exc])]
            px, py = wl_exc_cut[t_exc], y_exc[t_exc]
            ax.plot(px, py, 'o', color='k', markersize=4)
            ax.axvline(px, color=color, linestyle='--', linewidth=0.6, alpha=0.6)
            ax.text(px, py+0.04, f"{px:.1f} nm", ha='center', va='bottom', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        if peaks_em.size>0:
            t_em = peaks_em[np.argmax(y_em[peaks_em])]
            px, py = wl_em_cut[t_em], y_em[t_em]
            ax.plot(px, py, 'o', color='k', markersize=4)
            ax.axvline(px, color=color, linestyle=':', linewidth=0.6, alpha=0.6)
            ax.text(px, py+0.04, f"{px:.1f} nm", ha='center', va='bottom', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    else:
        for p in peaks_exc:
            px, py = wl_exc_cut[p], y_exc[p]
            ax.plot(px, py, 'o', color='k', markersize=4)
            ax.axvline(px, color=color, linestyle='--', linewidth=0.6, alpha=0.6)
            ax.text(px, py+0.03, f"{px:.1f} nm", ha='center', va='bottom', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
        for p in peaks_em:
            px, py = wl_em_cut[p], y_em[p]
            ax.plot(px, py, 'o', color='k', markersize=4)
            ax.axvline(px, color=color, linestyle=':', linewidth=0.6, alpha=0.6)
            ax.text(px, py+0.03, f"{px:.1f} nm", ha='center', va='bottom', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))

    # Determine x-limits for this subplot:
    if label in per_sample_wl:
        xmin, xmax = per_sample_wl[label]
    elif GLOBAL_WL_RANGE is not None:
        xmin, xmax = GLOBAL_WL_RANGE
    elif global_min is not None and global_max is not None:
        xmin, xmax = global_min, global_max
    else:
        # combine ranges available
        ranges = []
        if wl_exc_cut.size>0:
            ranges.append((np.min(wl_exc_cut), np.max(wl_exc_cut)))
        if wl_em_cut.size>0:
            ranges.append((np.min(wl_em_cut), np.max(wl_em_cut)))
        xs = [r[0] for r in ranges] + [r[1] for r in ranges]
        xmin, xmax = (min(xs), max(xs)) if xs else (0,1)
    ax.set_xlim(xmin, xmax)

    # Formatting
    ax.set_ylim(0, YMAX_PLOT)
    ax.margins(y=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([0, 1])
    ax.tick_params(axis="y", length=3, width=0.8)
    ax.tick_params(axis="x", length=3, width=0.8)
    # legend: show which is exc/em
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.5)
    if i == 0:
        ax.set_title("Fluorescence: Excitation (solid) and Emission (dashed)")

# Shared labels
fig.text(0.04, 0.5, "Normalized intensity", va="center", rotation="vertical")
axes[-1].set_xlabel("Wavelength (nm)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
plt.savefig(outfile)
print(f"Saved stacked excitation+emission plot to: {outfile}")

# Optionally save individual PNGs for each paired sample
if SAVE_INDIVIDUAL_PNG:
    for i, s in enumerate(paired):
        label = s['label']
        png_path = os.path.join(outdir, f"{label}_exc_em.png")
        fig_i, ax_i = plt.subplots(figsize=(6,3), dpi=300)
        color = color_list[i % len(color_list)]
        wl_exc, inten_exc = read_two_col_excel(s['exc'])
        wl_em, inten_em = read_two_col_excel(s['em'])
        idx_exc = np.argsort(wl_exc); wl_exc = wl_exc[idx_exc]; inten_exc = inten_exc[idx_exc]
        idx_em = np.argsort(wl_em); wl_em = wl_em[idx_em]; inten_em = inten_em[idx_em]
        wl_exc_cut, inten_exc_cut = apply_cut(wl_exc, inten_exc, xmin_cut, xmax_cut)
        wl_em_cut, inten_em_cut   = apply_cut(wl_em, intent_em,   xmin_cut, xmax_cut) if 'apply_cut' in globals() else (wl_em, inten_em)
        if wl_exc_cut.size>0:
            y_exc = normalize_by_abs_max(inten_exc_cut)
            ax_i.plot(wl_exc_cut, y_exc, color=color, **EXC_STYLE)
        if wl_em_cut.size>0:
            y_em = normalize_by_abs_max(inten_em_cut)
            ax_i.plot(wl_em_cut, y_em, color=color, **EM_STYLE)
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
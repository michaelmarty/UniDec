# dls_malvern_number_plot_fixed_v3.py
# Fixed to properly interpret Distribution Fit Data as log10(d.nm)
# and handle length mismatches by constructing standard Malvern grid
# Requires: pandas, numpy, matplotlib

import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import LogLocator, FuncFormatter
from io import StringIO

rcParams['font.size'] = 12
rcParams['font.sans-serif'] = 'Arial'
rcParams['axes.linewidth'] = 1

# ---------- USER SETTINGS ----------
FILES = [
   r"D:\UT Summer 2026\FLuoro ND paper revision\DLS\mcherry.txt",
    r"D:\UT Summer 2026\FLuoro ND paper revision\DLS\MSP try 2 stuff.txt",
    r"D:\UT Summer 2026\FLuoro ND paper revision\DLS\sfGFP.txt",
]
FILE_GLOB = r"./*.txt"

OUTDIR = "./dls_plots_debug"
os.makedirs(OUTDIR, exist_ok=True)
OUT_OVERLAY_PDF = os.path.join(OUTDIR, "dls_number_distribution_overlay.pdf")
OUT_OVERLAY_PNG = os.path.join(OUTDIR, "dls_number_distribution_overlay.png")

COLORS = ["#d62728", "#C76BA0", "#00A9A5", "#333333", "#FF7F0E"]

PLOT_OVERLAY = True
SAVE_INDIVIDUAL = True

# Match Malvern: log X axis, 0.1–10000 nm display range
LOG_X = True
X_MIN_DISPLAY = 0.1
X_MAX_DISPLAY = 10000.0

# If the file provides Intensities[...] (Percent) but not Numbers[...],
# set COMPUTE_NUMBER_FROM_INTENSITY=True to estimate Number% via I/d^6.
COMPUTE_NUMBER_FROM_INTENSITY = False
RAYLEIGH_POWER = 6  # exponent (d^6)
MIN_SIZE_FOR_DIV = 1e-6  # nm guard


# ---------- Helpers ----------
def try_read_csv(path):
    encodings = ['utf-8', 'utf-8-sig', 'latin1', 'cp1252']
    for enc in encodings:
        try:
            df = pd.read_csv(path, header=0, encoding=enc, engine='python')
            print(f"Read {os.path.basename(path)} with encoding {enc}")
            return df
        except Exception:
            continue
    # fallback: read bytes and decode replacing invalid bytes
    with open(path, 'rb') as f:
        raw = f.read()
    txt = raw.decode('utf-8', errors='replace')
    df = pd.read_csv(StringIO(txt), header=0, engine='python')
    print(f"Read {os.path.basename(path)} with fallback replace-decode")
    return df


def find_cols(cols, patterns):
    """Return list of columns whose name contains any substring from patterns (case-insensitive)"""
    out = []
    for c in cols:
        lc = c.lower()
        for p in patterns:
            if p in lc:
                out.append(c)
                break
    return out


def construct_standard_malvern_grid(n_bins):
    """
    Construct the standard Malvern log-spaced size grid.
    Malvern instruments typically use 0.1 to 10000 nm.
    """
    return np.logspace(np.log10(X_MIN_DISPLAY), np.log10(X_MAX_DISPLAY), n_bins)


def read_malvern_file_with_diag(path):
    """Return list of measurements (dicts) AND print diagnostics."""
    df = try_read_csv(path)
    # strip column names
    df.columns = [str(c).strip() for c in df.columns]
    cols = list(df.columns)

    # Detect columns
    dist_cols = find_cols(cols, ['distribution fit data', 'distribution fit'])  # size bin centers
    numbers_cols = find_cols(cols, ['numbers[', ' numbers', 'numbers (percent)'])
    intensity_cols = find_cols(cols, ['intensity[', 'intensities', 'intensity (percent)'])

    # Fallbacks
    if not numbers_cols:
        # also try columns starting with 'numbers'
        numbers_cols = [c for c in cols if c.lower().startswith('numbers')]
    if not dist_cols:
        dist_cols = [c for c in cols if 'distribution fit data' in c.lower()]

    # Build sizes array from dist_cols if available (take first numeric value per column)
    sizes_from_dist = None
    if dist_cols:
        vals = []
        for c in dist_cols:
            v = np.nan
            for r in range(len(df)):
                try:
                    vv = df.iloc[r][c]
                    if isinstance(vv, str):
                        vv = vv.strip()
                        if vv == '':
                            continue
                    v = float(vv)
                    break
                except Exception:
                    continue
            vals.append(v)
        sizes_all = np.array(vals, dtype=float)
        # drop NaNs
        valid = ~np.isnan(sizes_all)
        sizes_log = sizes_all[valid]
        dist_cols = [c for j, c in enumerate(dist_cols) if valid[j]]

        # KEY FIX: Transform from log10(d.nm) to d.nm
        sizes_from_dist = 10 ** sizes_log

        print(f"\n  Distribution Fit Data interpretation:")
        print(f"    Raw values (first 8): {sizes_log[:8]}")
        print(f"    Converted to nm (10^x): {sizes_from_dist[:8]}")

    # Diagnostics: print detected columns counts
    print(f"\nFile: {os.path.basename(path)}")
    print(f"  columns total: {len(cols)}")
    print(f"  found Distribution columns: {len(dist_cols)}")
    print(f"  found Numbers columns: {len(numbers_cols)}")
    print(f"  found Intensities columns: {len(intensity_cols)}")

    measurements = []
    for row in range(len(df)):
        # label
        label = None
        for cand in ['Sample Name', 'SampleName', 'Sample', 'Filename', 'Type']:
            if cand in df.columns:
                label = str(df.iloc[row][cand])
                break
        if label is None:
            label = f"{os.path.splitext(os.path.basename(path))[0]}_row{row}"

        # extract numbers
        nums = []
        for c in numbers_cols:
            try:
                v = df.iloc[row][c]
                if isinstance(v, str):
                    v = v.strip().replace('%', '')
                    if v == '':
                        v = 0.0
                nums.append(float(v))
            except Exception:
                nums.append(np.nan)
        nums = np.array(nums, dtype=float) if nums else np.array([], dtype=float)

        # extract intensities
        ints = []
        for c in intensity_cols:
            try:
                v = df.iloc[row][c]
                if isinstance(v, str):
                    v = v.strip().replace('%', '')
                    if v == '':
                        v = 0.0
                ints.append(float(v))
            except Exception:
                ints.append(np.nan)
        ints = np.array(ints, dtype=float) if ints else np.array([], dtype=float)

        # KEY FIX: Determine proper size grid
        sizes_use = None
        n_numbers = len(nums) if nums.size > 0 else len(ints)

        if sizes_from_dist is not None:
            # Check if lengths match
            if len(sizes_from_dist) == n_numbers:
                # Perfect match - use the converted Distribution Fit Data
                sizes_use = sizes_from_dist.copy()
                print(f"  Row {row}: Using Distribution Fit Data (converted from log), {len(sizes_use)} bins")
            else:
                # Length mismatch - construct standard Malvern grid instead
                print(f"  Row {row}: Length mismatch (dist={len(sizes_from_dist)}, numbers={n_numbers})")
                print(
                    f"           Constructing standard {n_numbers}-bin grid from {X_MIN_DISPLAY} to {X_MAX_DISPLAY} nm")
                sizes_use = construct_standard_malvern_grid(n_numbers)
        else:
            # No Distribution Fit Data - construct standard grid
            if n_numbers == 0:
                print(f"  Warning: no numeric distribution columns for row {row}. Skipping.")
                continue
            print(f"  Row {row}: No Distribution Fit Data found, constructing standard {n_numbers}-bin grid")
            sizes_use = construct_standard_malvern_grid(n_numbers)

        # Now align arrays (should already match but be safe)
        L = min(len(sizes_use), len(nums) if nums.size > 0 else len(ints))
        if nums.size > 0:
            nums = nums[:L]
        if ints.size > 0:
            ints = ints[:L]
        sizes_use = sizes_use[:L]

        # Replace NaNs
        nums = np.nan_to_num(nums, nan=0.0)
        ints = np.nan_to_num(ints, nan=0.0)

        # If Numbers not present but Intensities present and user wants compute option:
        numbers_pct = nums
        if numbers_pct.size == 0 and ints.size > 0 and COMPUTE_NUMBER_FROM_INTENSITY:
            # Estimate number from intensity via d^power
            sizes_safe = np.maximum(sizes_use, MIN_SIZE_FOR_DIV)
            Nraw = ints / (sizes_safe ** RAYLEIGH_POWER)
            Nraw = np.nan_to_num(Nraw, nan=0.0, posinf=0.0, neginf=0.0)
            s = np.sum(Nraw)
            numbers_pct = (Nraw / s) * 100.0 if s > 0 else np.zeros_like(Nraw)
            note = "(computed from intensities)"
        elif numbers_pct.size > 0:
            s = np.sum(numbers_pct)
            note = "(numbers provided; sum={:.3f})".format(s)
        else:
            numbers_pct = np.zeros_like(sizes_use)
            note = "(no numbers or intensities) "

        # Diagnostics per measurement
        print(f"  Measurement: {label} row {row}: bins={len(sizes_use)}")
        print(f"    Size range: {sizes_use[0]:.4f} to {sizes_use[-1]:.4f} nm")
        print(f"    First 8 sizes (nm): {np.round(sizes_use[:8], 3)}")
        print(f"    numbers columns read: {len(nums)}, sum(numbers)={np.nansum(nums):.4f}")
        if ints.size > 0:
            print(f"    intensities cols read: {len(ints)}, sum(intensities)={np.nansum(ints):.4f}")
        print(f"    using note: {note}")
        print(f"    first numbers (pct): {np.round(numbers_pct[:8], 4)}")
        print(f"    Peak position (approx): {sizes_use[np.argmax(numbers_pct)]:.2f} nm\n")

        measurements.append(
            {
                'label': label,
                'sizes_nm': sizes_use.copy(),
                'numbers_pct': numbers_pct.copy(),
                'intensities_pct': ints.copy(),
                'source_file': path,
                'row_index': row,
            }
        )
    return measurements


def apply_malvern_axes(ax, sizes, numbers):
    """Configure axes so they visually match the Malvern 'Size Distribution by Number' plot."""
    ax.set_xlabel("Size (d.nm)")
    ax.set_ylabel("Number (Percent)")

    if LOG_X:
        ax.set_xscale('log')
        ax.set_xlim(X_MIN_DISPLAY, X_MAX_DISPLAY)

        # Use decade ticks but label them as plain nm values (0.1, 1, 10, 100, 1000, 10000)
        ax.xaxis.set_major_locator(LogLocator(base=10.0))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))

        # Keep minor ticks for visual reference but no labels
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda x, pos: ""))
    else:
        # if linear is requested, still ensure non-negative lower limit
        ax.set_xlim(left=0.0)

    # Y range: from 0 up to a rounded max based on the data
    if numbers.size:
        ymax = numbers.max()
        if not np.isfinite(ymax) or ymax <= 0:
            ymax = 25.0
        else:
            # round up to nearest 5 to mimic Malvern ticks (0, 5, 10, ...)
            ymax = 5.0 * np.ceil(ymax / 5.0)
        ax.set_ylim(0.0, 30.0)
    else:
        ax.set_ylim(0.0, 30.0)


# ---------- Main ----------
files = FILES[:] if FILES else sorted(glob.glob(FILE_GLOB))
if not files:
    raise RuntimeError("No files found.")

all_measurements = []
for f in files:
    try:
        m = read_malvern_file_with_diag(f)
        all_measurements.extend(m)
    except Exception as e:
        print(f"Error reading {f}: {e}")

if not all_measurements:
    raise RuntimeError("No measurements parsed.")

# Plot overlay
if PLOT_OVERLAY:
    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)

    for i, m in enumerate(all_measurements):
        sizes = m['sizes_nm']
        numbers = m['numbers_pct']
        color = COLORS[i % len(COLORS)]

        # Filled area + outline to match Malvern style
        ax.fill_between(sizes, numbers, 0.0, color=color, alpha=0.35, linewidth=0)
        ax.plot(sizes, numbers, color=color, linewidth=1.5, label=m['label'])

    apply_malvern_axes(ax, all_measurements[0]['sizes_nm'], all_measurements[0]['numbers_pct'])
    ax.set_title("Size Distribution by Number")
    ax.set_xlabel("Size (r.nm)")
    ax.legend(frameon=False, fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(OUT_OVERLAY_PDF)
    plt.savefig(OUT_OVERLAY_PNG, dpi=300)
    plt.close(fig)
    print("Saved overlay:", OUT_OVERLAY_PDF, OUT_OVERLAY_PNG)

# Individual plots
if SAVE_INDIVIDUAL:
    for i, m in enumerate(all_measurements):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        sizes = m['sizes_nm']
        numbers = m['numbers_pct']
        color = COLORS[i % len(COLORS)]

        ax.fill_between(sizes, numbers, 0.0, color=color, alpha=0.6, linewidth=0)
        ax.plot(sizes, numbers, color=color, linewidth=1.6)

        apply_malvern_axes(ax, sizes, numbers)
        ax.set_title("Size Distribution by Number")

        plt.tight_layout()
        fname = os.path.join(
            OUTDIR,
            f"{os.path.splitext(os.path.basename(m['source_file']))[0]}_row{m['row_index']}_number.png",
        )
        plt.savefig(fname, dpi=300)
        plt.close(fig)
    print("Saved individual PNGs to:", OUTDIR)

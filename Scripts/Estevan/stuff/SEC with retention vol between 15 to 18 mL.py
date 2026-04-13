import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import numpy as np

# -----------------------
# Settings
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

# === FILES TO PLOT ===
files = [
    r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep1.csv",
    r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep2.csv",
    r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep3.csv",
]

# === LEGEND LABELS ===
sample_labels = [
    "MSP1E3D1-sfGFP",
    "MSP1E3D1-sfGFP",
    "MSP1E3D1-sfGFP",
]

# === CUSTOM HEX COLORS ===
colors = [
    "#333333",
    "#333333",
    "#333333",
]

# -----------------------
# SEC STANDARD VOLUMES (mL)
# -----------------------
standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_color = "#8C8C8C"

# -----------------------
# Peak detection / grouping parameters
# -----------------------
rv_window = (15.0, 18.0)     # retention volume window to consider (mL)
group_rv_tol = 0.15          # group candidate peaks within this mL distance as one peak
min_rel_height = 0.0         # optional: minimum relative height (0..1) to keep candidates; 0 disables

# -----------------------
# Helper: detect likely delimiter from file start
# -----------------------
def detect_delimiter(file_path, encodings=("utf-8", "utf-8-sig", "latin1", "utf-16")):
    for enc in encodings:
        try:
            with open(file_path, "r", encoding=enc, errors="strict") as f:
                sample = f.read(2048)
            if "\t" in sample:
                return "\t"
            if "," in sample:
                return ","
            return None
        except Exception:
            continue
    return None

# -----------------------
# Read Excel/CSV starting at A4/B4 (robust to encodings/delimiters)
# -----------------------
def read_sec_file(file_path, csv_sep=None):
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    if ext in (".xls", ".xlsx"):
        df = pd.read_excel(file_path, skiprows=3, header=None, usecols="A:B")
    elif ext == ".csv":
        base_kwargs = {"skiprows": 3, "header": None, "usecols": [0, 1]}
        encodings_to_try = ["utf-8", "utf-8-sig", "latin1", "utf-16"]
        df = None
        last_exc = None

        if csv_sep is None:
            csv_sep = detect_delimiter(file_path, encodings=encodings_to_try)

        if csv_sep is not None:
            for enc in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, sep=csv_sep, encoding=enc, **base_kwargs)
                    break
                except Exception as e:
                    last_exc = e

        if df is None:
            for enc in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, encoding=enc, **base_kwargs)
                    break
                except Exception as e:
                    last_exc = e
                    continue

        if df is None:
            try:
                df = pd.read_csv(file_path, engine="python", sep=None, encoding="latin1", **base_kwargs)
            except Exception as e:
                raise RuntimeError(f"Failed to read CSV file '{file_path}'. Last error: {e}") from e
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Supported:.xls.xlsx.csv")

    df.columns = ["ml", "mAU"]
    df = df.dropna()
    return df

# -----------------------
# Local maxima detector (neighborhood comparison)
# -----------------------
def find_local_maxima(y):
    if len(y) < 3:
        return np.array([], dtype=int)
    left = y > np.roll(y, 1)
    right = y > np.roll(y, -1)
    local_max = left & right
    local_max[0] = False
    local_max[-1] = False
    return np.where(local_max)[0]

# -----------------------
# Collapse nearby candidate peaks into single distinguished peaks
# -----------------------
def pick_distinguished_peaks(x, y, candidate_idx, group_tol):
    """
    Group candidate indices that are within group_tol (in x units) and select
    the index with the highest y in each group.
    Returns an array of picked indices (sorted).
    """
    if len(candidate_idx) == 0:
        return np.array([], dtype=int)

    # Sort candidates by x value
    sorted_idx = candidate_idx[np.argsort(x[candidate_idx])]
    groups = []
    current_group = [sorted_idx[0]]

    for idx in sorted_idx[1:]:
        if x[idx] - x[current_group[-1]] <= group_tol:
            current_group.append(idx)
        else:
            groups.append(current_group)
            current_group = [idx]
    groups.append(current_group)

    # From each group, pick index with max y
    picked = []
    for g in groups:
        g = np.array(g, dtype=int)
        max_pos = g[np.argmax(y[g])]
        picked.append(int(max_pos))

    return np.array(sorted(picked), dtype=int)

# -----------------------
# Plot setup
# -----------------------
fig, axes = plt.subplots(len(files), 1, figsize=(8, 4), sharex=True, dpi=300)
if len(files) == 1:
    axes = [axes]

for i, (file, label) in enumerate(zip(files, sample_labels)):
    print(f"\n--- Processing file: {os.path.basename(file)} ---")
    df = read_sec_file(file)

    # Ensure numeric and drop non-numeric rows
    df["ml"] = pd.to_numeric(df["ml"], errors="coerce")
    df["mAU"] = pd.to_numeric(df["mAU"], errors="coerce")
    df = df.dropna().reset_index(drop=True)

    # Check and sort ml if needed
    if not df["ml"].is_monotonic_increasing:
        print("Warning: ml column not strictly increasing. Sorting by ml.")
        df = df.sort_values("ml").reset_index(drop=True)

    # Exclude anything after 20 mL to avoid air peaks
    df = df[df["ml"] <= 20.0].reset_index(drop=True)
    if df.empty:
        print(f"Warning: no data <= 20 mL for file: {file}. Skipping plot for this file.")
        continue

    # Normalize by absolute max of the filtered data
    max_abs = df["mAU"].abs().max()
    if max_abs == 0 or np.isnan(max_abs):
        norm_mAU = df["mAU"]
    else:
        norm_mAU = df["mAU"] / max_abs
    df["mAU"] = norm_mAU

    x = df["ml"].values
    y = df["mAU"].values

    # Plot line and fill
    axes[i].plot(x, y, color=colors[i % len(colors)], linewidth=1.8, label=label)
    axes[i].fill_between(x, y, 0, color=colors[i % len(colors)], alpha=0.9, linewidth=0)

    # Candidate local maxima (no grouping yet)
    candidate_idx = find_local_maxima(y)

    # Optionally filter candidates by minimum relative height to suppress tiny noise
    if min_rel_height > 0.0:
        height_mask = y[candidate_idx] >= (min_rel_height * np.nanmax(y))
        candidate_idx = candidate_idx[height_mask]

    # Diagnostics: print all candidates (index, ml, mAU)
    if len(candidate_idx) > 0:
        print("Local maxima candidates (idx, ml, mAU):")
        for idx in candidate_idx:
            print(f"  {idx:4d}  {x[idx]:6.3f} mL  {y[idx]:6.4f}")
    else:
        print("No local maxima candidates found.")

    # Keep only candidates inside the requested rv window
    in_win_mask = (x[candidate_idx] >= rv_window[0]) & (x[candidate_idx] <= rv_window[1])
    candidates_in_window = candidate_idx[in_win_mask]

    print(f"Candidates in window {rv_window[0]}-{rv_window[1]} mL: {[round(float(x[idx]),3) for idx in candidates_in_window]}")

    # Collapse nearby candidates and pick the max of each distinguished peak
    picked_peak_idx = pick_distinguished_peaks(x, y, candidates_in_window, group_rv_tol)

    # Diagnostics: show groups/picks
    if len(picked_peak_idx) > 0:
        print("Picked distinguished peaks (ml, mAU):")
        for p_idx in picked_peak_idx:
            print(f"  {x[p_idx]:6.3f} mL  {y[p_idx]:6.4f}")
    else:
        print("No distinguished peaks picked in window.")

    # Annotate the picked peaks
    if len(picked_peak_idx) > 0:
        picked_rvs = []
        for p_idx in picked_peak_idx:
            peak_x = x[p_idx]
            peak_y = y[p_idx]
            picked_rvs.append(peak_x)

            # Marker
            axes[i].plot(peak_x, peak_y, marker='o', color='k', markersize=5)
            # Vertical dashed line
            axes[i].axvline(peak_x, color=colors[i % len(colors)], linestyle='--', linewidth=0.9, alpha=0.75)
            # Text annotation slightly above the peak
            axes[i].text(
                peak_x,
                peak_y + 0.04,
                f"{peak_x:.2f} mL",
                ha="center",
                va="bottom",
                fontsize=9,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1)
            )

        print(f"Final picked peaks in {rv_window[0]}-{rv_window[1]} mL: {', '.join(f'{v:.2f}' for v in picked_rvs)}")
    else:
        print(f"No final picks in {rv_window[0]}-{rv_window[1]} mL.")

    # Also mark all local maxima candidates (for debugging) in red small dots
    if len(candidate_idx) > 0:
        axes[i].plot(x[candidate_idx], y[candidate_idx], 'ro', markersize=3, alpha=0.6)

    # Axis formatting
    axes[i].set_ylim(0, 1.08)
    axes[i].margins(y=0)
    axes[i].set_xlim(10, 20)

    axes[i].spines["top"].set_visible(False)
    axes[i].spines["right"].set_visible(False)

    axes[i].set_yticks([0, 1])
    axes[i].tick_params(axis="y", length=3, width=0.8)
    axes[i].tick_params(axis="x", length=3, width=0.8)

    axes[i].legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        handlelength=1.5
    )

    if i == 0:
        top_ax = axes[i].secondary_xaxis('top')
        top_ax.set_xticks(standard_volumes)
        top_ax.set_xticklabels([])
        top_ax.tick_params(axis='x', length=5, width=1.0, colors=standard_color)
        top_ax.spines['top'].set_visible(False)
        top_ax.set_xlabel("SEC standards", color=standard_color, labelpad=6)

# -----------------------
# Shared axis labels
# -----------------------
fig.text(0.04, 0.5, "Normalized Absorbance", va="center", rotation="vertical")
axes[-1].set_xlabel("Retention Volume (mL)")

# -----------------------
# Layout and save
# -----------------------
plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
outdir = r"D:\UT Spring 2026\general lab data\FluorND"
os.makedirs(outdir, exist_ok=True)
plt.savefig(os.path.join(outdir, "sec_sec_distinguished_peaks_15_18mL.pdf"))

plt.show()
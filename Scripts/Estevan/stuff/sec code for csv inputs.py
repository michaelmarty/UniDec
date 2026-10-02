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
rcParams['axes.linewidth'] = 1
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

# === FILES TO PLOT ===
files = [
     r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep1.csv",
     r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep1.csv",
    r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep3.csv",
    #r"D:\UT Summer 2026\shawn chromatograms\E3T0_1to160_DMPC_SC 001.csv",
    #r"D:\UT Summer 2026\shawn chromatograms\E3T0_1to170_DMPC 001.csv",
    #r"E:\FLND AKTA files\XLSX chromatograms\D1mCheery_1to80DMPC.xlsx",
    #r"E:\FLND AKTA files\XLSX chromatograms\d1bLUEBERRY_1TO80_dmpc.xlsx",
    #r"E:\FLND AKTA files\XLSX chromatograms\D1HaloTag_1to90_DMPC.xlsx",
#r"E:\FLND AKTA files\XLSX chromatograms\E3sfGFP_1to150_DMPC.xlsx",
]

# === LEGEND LABELS ===
sample_labels = [
     "MSP1D1",
    "MSP1D1-sfGFP",
    "MSP1D1-mCherry",
    #"MSP1E3D1 1:160  DMPC",
    #"MSP1E3D1 1:170  DMPC",
    #"MSP1D1-mCherry",
    #"MSP1D1-mBlueberry2",
    #"MSP1D1-HaloTag",
#"MSP1E3D1-sfGFP",
]

# === CUSTOM HEX COLORS ===
colors = [
    "#333333",
    "#00A9A5",
    "#C76BA0",
    #"#C76BA0",
    #"#C76BA0",
    #"#C76BA0",
    #"#2C7FB8",
    #"#7A7A7A",
#"#00A9A5",
]

# -----------------------
# SEC STANDARD VOLUMES (mL)
# -----------------------
standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_color = "#8C8C8C"

# -----------------------
# Peak detection threshold (fraction of max)
# -----------------------
# Default: $$\text{peak_threshold} = 0.10$$
peak_threshold = 0.10

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
# Simple peak detector (local maxima + min height)
# -----------------------
def find_peaks_simple(x, y, min_rel_height=0.10):
    if len(y) < 3:
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
# Plot setup
# -----------------------
fig, axes = plt.subplots(len(files), 1, figsize=(8, 4), sharex=True, dpi=300)
if len(files) == 1:
    axes = [axes]

for i, (file, label) in enumerate(zip(files, sample_labels)):
    df = read_sec_file(file)

    # Ensure numeric and drop non-numeric rows
    df["ml"] = pd.to_numeric(df["ml"], errors="coerce")
    df["mAU"] = pd.to_numeric(df["mAU"], errors="coerce")
    df = df.dropna()

    # Exclude anything after 20 mL to avoid air peaks
    df = df[df["ml"] <= 20.0]
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

    # Peak detection
    peak_idx = find_peaks_simple(x, y, min_rel_height=peak_threshold)

    # --- NEW: annotate only the tallest peak ---
    if len(peak_idx) > 0:
        # pick the peak index with max y value
        tallest_rel_idx = peak_idx[np.argmax(y[peak_idx])]
        peak_x = x[tallest_rel_idx]
        peak_y = y[tallest_rel_idx]

        # Marker
        axes[i].plot(peak_x, peak_y, marker='o', color='k', markersize=4)
        # Vertical dashed line
        axes[i].axvline(peak_x, color=colors[i % len(colors)], linestyle='--', linewidth=0.8, alpha=0.6)
        # Text annotation slightly above the peak
        axes[i].text(
            peak_x,
            peak_y + 0.04,
            f"{peak_x:.2f} mL",
            ha="center",
            va="bottom",
            fontsize=9,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1)
        )
    # -------------------------------------------------

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
plt.savefig(os.path.join(outdir, "sec_sec_std_with_tallest_peak.pdf"))

plt.show()
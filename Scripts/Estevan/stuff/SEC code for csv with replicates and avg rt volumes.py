import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import numpy as np

# -----------------------
# Settings (unchanged)
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
# FILES TO PLOT (groups of replicates)
# Replace these lists with your file paths. Each inner list is one subplot (replicates).
# Example: two samples with 3 and 2 replicates respectively.
files = [
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep1.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep2.csv", r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep3.csv"]
    # add more groups as needed...
]

# Labels for each group (must match len(files))
sample_labels = [
    "MSP1D1",
    "MSP1D1-mBlueberry",
    "MSP1D1-mCherry",
    "MSP1D1-sfCherryS",
    "MSP1D1-sfCherryL",
    "MSP1D1-sfGFP",
    "MSP1E3D1-sfGFP",
    "MSP1D1-HaloTag"
    # add more labels corresponding to groups above...
]

# Colors per group
colors = [
    "#333333",
    "#1261A0",
    "#B947B1",
    "#B947B1",
    "#B947B1",
    "#68BB59",
"#68BB59",
    "#ACADA8"
    #"#2C7FB8",
    #"#7A7A7A",
    #"#C76BA0",
]

# -----------------------
# SEC STANDARD VOLUMES (mL)
# -----------------------
standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_color = "#8C8C8C"

# Peak detection threshold (fraction of max)
peak_threshold = 0.10  # adjust if necessary

# -----------------------
# Helper: detect delimiter
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
# Read function (robust CSV/XLSX)
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

        if csv_sep is None:
            csv_sep = detect_delimiter(file_path, encodings=encodings_to_try)

        if csv_sep is not None:
            for enc in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, sep=csv_sep, encoding=enc, **base_kwargs)
                    break
                except Exception:
                    continue

        if df is None:
            for enc in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, encoding=enc, **base_kwargs)
                    break
                except Exception:
                    continue

        if df is None:
            df = pd.read_csv(file_path, engine="python", sep=None, encoding="latin1", **base_kwargs)
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Supported:.xls.xlsx.csv")

    df.columns = ["ml", "mAU"]
    df = df.dropna()
    return df

# -----------------------
# Simple peak detector (local maxima + min height)
# -----------------------
def find_peaks_simple(y, min_rel_height=0.10):
    """
    Return indices of peaks where y is a local maximum and y >= min_rel_height * y_max.
    Input: y = 1D numpy array
    """
    y = np.asarray(y)
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
n_groups = len(files)
fig, axes = plt.subplots(n_groups, 1, figsize=(8, 4), sharex=True, dpi=300)
if n_groups == 1:
    axes = [axes]

for i, (group_files, label) in enumerate(zip(files, sample_labels)):
    peak_times = []   # collect tallest-peak retention times for this group
    color = colors[i % len(colors)]
    plotted_any = False

    for rep_file in group_files:
        try:
            df = read_sec_file(rep_file)
        except Exception as e:
            print(f"Failed reading {rep_file}: {e}")
            continue

        # Ensure numeric, drop non-numeric
        df["ml"] = pd.to_numeric(df["ml"], errors="coerce")
        df["mAU"] = pd.to_numeric(df["mAU"], errors="coerce")
        df = df.dropna()

        # exclude anything after 20 mL
        df = df[df["ml"] <= 20.0]
        if df.empty:
            print(f"Warning: no data <=20 mL for file {rep_file}. Skipping replicate.")
            continue

        # normalize each replicate independently by its absolute max
        max_abs = df["mAU"].abs().max()
        if max_abs == 0 or np.isnan(max_abs):
            norm_mAU = df["mAU"].values
        else:
            norm_mAU = (df["mAU"] / max_abs).values

        x = df["ml"].values
        y = norm_mAU

        # Plot replicate trace with reduced alpha so overlaps show density
        axes[i].plot(x, y, color=color, linewidth=1.2, alpha=0.6)
        axes[i].fill_between(x, y, 0, color=color, alpha=0.25, linewidth=0)

        # find tallest peak of this replicate
        peaks_idx = find_peaks_simple(y, min_rel_height=peak_threshold)
        if len(peaks_idx) > 0:
            tallest_idx = peaks_idx[np.argmax(y[peaks_idx])]
            peak_time = float(x[tallest_idx])
            peak_times.append(peak_time)

        plotted_any = True

    if not plotted_any:
        axes[i].text(0.5, 0.5, "No valid replicates", transform=axes[i].transAxes, ha="center")
        axes[i].set_xlim(10, 20)
        axes[i].set_ylim(0, 1.08)
        continue

    # compute mean and std of the collected peak times
    if len(peak_times) == 0:
        mean_peak = None
        std_peak = None
    else:
        mean_peak = float(np.mean(peak_times))
        std_peak = float(np.std(peak_times, ddof=1)) if len(peak_times) > 1 else 0.0

    # Annotate mean +/- std on the subplot at the mean peak position (if we have a mean)
    if mean_peak is not None:
        # draw vertical dashed line at mean
        axes[i].axvline(mean_peak, color=color, linestyle='--', linewidth=1.0, alpha=0.9)

        # annotation string: mean ± std (2 decimal places)
        if std_peak is None:
            ann_text = f"{mean_peak:.2f} mL"
        else:
            ann_text = f"{mean_peak:.2f} ± {std_peak:.2f} mL"

        # place annotation slightly above the highest y at that x (approx)
        # We'll find y at the nearest x in last plotted replicate if available:
        # find y-value nearest to mean_peak among the last processed replicate's x (fallback)
        # If no good y, place annotation at y=0.9
        try:
            # get the last processed x,y from the loop above (they exist if plotted_any)
            idx_near = (np.abs(x - mean_peak)).argmin()
            y_at_mean = y[idx_near]
            text_y = min(y_at_mean + 0.06, 1.02)
        except Exception:
            text_y = 0.9

        # put the text
        axes[i].text(mean_peak, text_y, ann_text, ha="center", va="bottom",
                     fontsize=9, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

    # Axis formatting (same as before)
    axes[i].set_ylim(0, 1.08)
    axes[i].margins(y=0)
    axes[i].set_xlim(10, 20)

    axes[i].spines["top"].set_visible(False)
    axes[i].spines["right"].set_visible(False)

    axes[i].set_yticks([0, 1])
    axes[i].tick_params(axis="y", length=3, width=0.8)
    axes[i].tick_params(axis="x", length=3, width=0.8)

    axes[i].legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                   frameon=False, handlelength=1.5, labels=[label])

    # TOP ticks for standards on the first subplot
    if i == 0:
        top_ax = axes[i].secondary_xaxis('top')
        top_ax.set_xticks(standard_volumes)
        top_ax.set_xticklabels([])
        top_ax.tick_params(axis='x', length=5, width=1.0, colors=standard_color)
        top_ax.spines['top'].set_visible(False)
        top_ax.set_xlabel("SEC standards", color=standard_color, labelpad=6)

# Shared axis labels & save
fig.text(0.04, 0.5, "Normalized Absorbance", va="center", rotation="vertical")
axes[-1].set_xlabel("Retention Volume (mL)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)
outdir = r"D:\UT Spring 2026\general lab data\FluorND"
os.makedirs(outdir, exist_ok=True)
plt.savefig(os.path.join(outdir, "sec_sec_std_replicates_meanstd.pdf"))

plt.show()
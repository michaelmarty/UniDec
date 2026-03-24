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
files = [
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1 DMPC Nanodiscs\MSP1D1 DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mBlueberry DMPC Nanodiscs\MSP1D1 mBlueberry DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCherryS DMPC Nanodiscs\MSP1D1 sfCherryS DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfCheryL DMPC Nanodiscs\MSP1D1 sfCherryL DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1E3D1sfGFP DMPC Nanodiscs\MSP1E3D1 sfGFP DMPC rep3.csv"],
    [r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep1.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep2.csv",
     r"D:\UT Spring 2026\general lab data\FluorND\Chromatogram CSV files\MSP1D1HaloTag DMPC Nanodiscs\MSP1D1 HaloTag DMPC rep3.csv"]
    # add more groups...
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
    # add more labels...
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

    
]

# -----------------------
# SEC STANDARD VOLUMES (mL) AND CORRESPONDING STOKES DIAMETERS
# Provide these two lists (same length, corresponding pairs).
# Example placeholders (replace with your actual standards):
standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_stokes =    [17, 10.52,  5.62,  3.5,  1.7]  # units: e.g., nm

# Color for standard ticks
standard_color = "#8C8C8C"

# Peak detection threshold (fraction of max)
peak_threshold = 0.10  # adjust if necessary

# -----------------------
# Helper: delimiter detection
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
def find_peaks_simple(y, min_rel_height=0.10):
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
# Calibration: fit linear mapping from volume -> Stokes diameter
# We assume D = a*V + b (linear). If you want a different model, change here.
if len(standard_volumes) != len(standard_stokes):
    raise ValueError("standard_volumes and standard_stokes must have the same length.")
V = np.array(standard_volumes, dtype=float)
D = np.array(standard_stokes, dtype=float)
# Fit linear (degree=1)
coeffs = np.polyfit(V, D, 1)   # coeffs[0]=a, coeffs[1]=b
a_fit, b_fit = coeffs
# Optional: compute R^2 for fit quality
D_pred = a_fit * V + b_fit
ss_res = np.sum((D - D_pred) ** 2)
ss_tot = np.sum((D - D.mean()) ** 2)
r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

# Cal function
def ml_to_stokes(ml):
    """Convert retention volume (ml) to Stokes diameter using the linear calibration."""
    return a_fit * ml + b_fit

# -----------------------
# Plot setup
n_groups = len(files)
fig, axes = plt.subplots(n_groups, 1, figsize=(8, 4), sharex=True, dpi=300)
if n_groups == 1:
    axes = [axes]

# Prepare summary list for CSV output
summary_rows = []

for i, (group_files, label) in enumerate(zip(files, sample_labels)):
    color = colors[i % len(colors)]
    peak_times = []    # retention times (mL) of tallest peak per replicate
    peak_diams = []    # converted Stokes diameters per replicate
    plotted_any = False

    for rep_file in group_files:
        try:
            df = read_sec_file(rep_file)
        except Exception as e:
            print(f"Failed reading {rep_file}: {e}")
            continue

        # ensure numeric and drop non-numeric rows
        df["ml"] = pd.to_numeric(df["ml"], errors="coerce")
        df["mAU"] = pd.to_numeric(df["mAU"], errors="coerce")
        df = df.dropna()

        # exclude anything after 20 mL
        df = df[df["ml"] <= 20.0]
        if df.empty:
            print(f"Warning: no data <= 20 mL for file {rep_file}. Skipping replicate.")
            continue

        # normalize replicate independently
        max_abs = df["mAU"].abs().max()
        if max_abs == 0 or np.isnan(max_abs):
            y = df["mAU"].values
        else:
            y = (df["mAU"] / max_abs).values
        x = df["ml"].values

        # plot replicate
        axes[i].plot(x, y, color=color, linewidth=1.2, alpha=0.6)
        axes[i].fill_between(x, y, 0, color=color, alpha=0.25, linewidth=0)

        # find tallest peak in this replicate
        peaks_idx = find_peaks_simple(y, min_rel_height=peak_threshold)
        if len(peaks_idx) > 0:
            tallest_idx = peaks_idx[np.argmax(y[peaks_idx])]
            peak_time = float(x[tallest_idx])
            peak_times.append(peak_time)
            # convert retention time to Stokes diameter
            peak_diams.append(float(ml_to_stokes(peak_time)))

        plotted_any = True

    if not plotted_any:
        axes[i].text(0.5, 0.5, "No valid replicates", transform=axes[i].transAxes, ha="center")
        axes[i].set_xlim(10, 20)
        axes[i].set_ylim(0, 1.08)
        summary_rows.append({
            "label": label,
            "n_replicates": 0,
            "mean_mL": np.nan,
            "std_mL": np.nan,
            "mean_Stokes": np.nan,
            "std_Stokes": np.nan
        })
        continue

    # compute mean/std in retention time (mL)
    if len(peak_times) == 0:
        mean_ml = None
        std_ml = None
    else:
        mean_ml = float(np.mean(peak_times))
        std_ml = float(np.std(peak_times, ddof=1)) if len(peak_times) > 1 else 0.0

    # compute mean/std in Stokes diameter
    if len(peak_diams) == 0:
        mean_d = None
        std_d = None
    else:
        mean_d = float(np.mean(peak_diams))
        std_d = float(np.std(peak_diams, ddof=1)) if len(peak_diams) > 1 else 0.0

    # Annotate mean retention and mean Stokes on the subplot
    if mean_ml is not None:
        axes[i].axvline(mean_ml, color=color, linestyle='--', linewidth=1.0, alpha=0.9)
        if std_ml is None:
            ann_ml = f"{mean_ml:.2f} mL"
        else:
            ann_ml = f"{mean_ml:.2f} ± {std_ml:.2f} mL"
    else:
        ann_ml = "No peak"

    if mean_d is not None:
        if std_d is None:
            ann_d = f"{mean_d:.2f}"
        else:
            ann_d = f"{mean_d:.2f} ± {std_d:.2f}"
        ann_text = f"{ann_ml}\nStokes: {ann_d}"
    else:
        ann_text = f"{ann_ml}\nStokes: N/A"

    # place annotation: try to position near the top of the curve at mean_ml
    try:
        idx_near = (np.abs(x - (mean_ml if mean_ml is not None else 0))).argmin()
        y_at_mean = y[idx_near] if 'y' in locals() else 0.9
        text_y = min(y_at_mean + 0.06, 1.02)
    except Exception:
        text_y = 0.9

    # put the text box
    axes[i].text((mean_ml if mean_ml is not None else 11.0), text_y, ann_text,
                 ha="center", va="bottom", fontsize=9,
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=2))

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

    if i == 0:
        top_ax = axes[i].secondary_xaxis('top')
        top_ax.set_xticks(standard_volumes)
        top_ax.set_xticklabels([])
        top_ax.tick_params(axis='x', length=5, width=1.0, colors=standard_color)
        top_ax.spines['top'].set_visible(False)
        top_ax.set_xlabel("SEC standards", color=standard_color, labelpad=6)

    # Append summary row
    summary_rows.append({
        "label": label,
        "n_replicates": len(peak_times),
        "mean_mL": mean_ml,
        "std_mL": std_ml,
        "mean_Stokes": mean_d,
        "std_Stokes": std_d
    })

# -----------------------
# Shared axis labels & save
fig.text(0.04, 0.5, "Normalized Absorbance", va="center", rotation="vertical")
axes[-1].set_xlabel("Retention Volume (mL)")

plt.tight_layout(rect=[0.08, 0.05, 0.95, 0.97], h_pad=0.6)
plt.subplots_adjust(hspace=0.55)

outdir = r"D:\UT Spring 2026\general lab data\FluorND"
os.makedirs(outdir, exist_ok=True)
figpath = os.path.join(outdir, "sec_sec_std_replicates_stokes.pdf")
plt.savefig(figpath)

# Also save CSV summary of mean/stdevs
summary_df = pd.DataFrame(summary_rows)
summary_csv = os.path.join(outdir, "sec_peak_stokes_summary.csv")
summary_df.to_csv(summary_csv, index=False)

# Print calibration info
print(f"Calibration linear fit: D = a*V + b with a={a_fit:.6g}, b={b_fit:.6g}, R^2={r2:.4f}")
print(f"Saved figure to: {figpath}")
print(f"Saved summary CSV to: {summary_csv}")

plt.show()
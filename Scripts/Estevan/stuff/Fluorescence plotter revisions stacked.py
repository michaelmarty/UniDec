import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
figsize = (7, 6)   # width, height (height large enough for stacked plots)
dpi = 300

# Define samples: each is a list of file paths (one or more XLSX per sample)
sample_A_files = [
    r"path\to\sampleA_rep1.xlsx",
    # r"path\to\sampleA_rep2.xlsx",
]
sample_B_files = [
    r"path\to\sampleB_rep1.xlsx",
    # r"path\to\sampleB_rep2.xlsx",
]
# If you want more than 2 samples, add more lists and include them in `all_samples` below.

# Labels and colors for each sample (must match number of samples in all_samples)
sample_labels = ["Sample A", "Sample B"]
sample_colors = ["#1f77b4", "#ff7f0e"]

# Processing options
apply_smoothing = True
savgol_params = {"window_length": 11, "polyorder": 3}   # window_length must be odd
normalize_each = True     # normalize each replicate individually
offset_increment_within = 0.0   # progressive offsets applied within each sample (per replicate)

# Excel parsing
header_row = 0   # set None if sheet has no header row
wavelength_names = ["wavelength (nm)", "wavelength", "nm"]
intensity_names = ["intensity (a.u.)", "intensity", "a.u."]

# Axis/figure settings
xlabel = "Wavelength (nm)"
ylabel = "Intensity (a.u.)"
x_range = None   # e.g. (300, 700)
y_range_per_sample = None  # or list of tuples per sample e.g. [(0,1.2),(0,0.8)]
share_x = True

# -------------------- end of user settings --------------------

plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

def smooth_savgol(y, window_length=11, polyorder=3):
    try:
        from scipy.signal import savgol_filter
    except Exception:
        return y
    n = len(y)
    wl = int(window_length)
    if wl >= n:
        wl = max(3, n - (1 - n % 2))
    if wl % 2 == 0:
        wl += 1
    return savgol_filter(y, wl, polyorder)

def _ensure_df_from_read_excel(obj):
    if isinstance(obj, dict):
        for df in obj.values():
            if isinstance(df, pd.DataFrame) and df.shape[0] > 0:
                return df
        return next(iter(obj.values()))
    if isinstance(obj, pd.DataFrame):
        return obj
    raise TypeError("Unexpected read_excel return type")

def detect_columns(df):
    cols_l = [str(c).lower() for c in df.columns]
    wcol = None
    icol = None
    for name in wavelength_names:
        if name in cols_l:
            wcol = df.columns[cols_l.index(name)]
            break
    for name in intensity_names:
        if name in cols_l:
            icol = df.columns[cols_l.index(name)]
            break
    if wcol is None:
        wcol = df.columns[0]
    if icol is None:
        icol = df.columns[1] if df.shape[1] > 1 else df.columns[0]
    return wcol, icol

def read_spectrum(path, sheet=None, header=header_row):
    if sheet is None:
        raw = pd.read_excel(path, sheet_name=0, header=header)
    else:
        raw = pd.read_excel(path, sheet_name=sheet, header=header)
    df = _ensure_df_from_read_excel(raw)
    if df.shape[0] == 0:
        raise ValueError(f"No data in file: {path}")
    wcol, icol = detect_columns(df)
    x = pd.to_numeric(df[wcol], errors='coerce').values
    y = pd.to_numeric(df[icol], errors='coerce').values
    mask = ~np.isnan(x) & ~np.isnan(y)
    return x[mask], y[mask]

def plot_stacked_samples(samples_list, labels=None, colors=None,
                         normalize=normalize_each, smoothing=apply_smoothing,
                         savgol_kwargs=None, offset_inc_within=offset_increment_within,
                         figsize_local=figsize, dpi_local=dpi, share_x_local=share_x,
                         xlim=x_range, ylim_per_sample=y_range_per_sample):
    """
    samples_list: list of lists. Each element is a list of file paths for one sample (replicates).
    labels: list of sample labels
    colors: list of colors per sample
    """
    savgol_kwargs = savgol_kwargs or savol_kwargs = savgol_params  # fallback binding
    n_samples = len(samples_list)
    labels = labels or [f"Sample {i+1}" for i in range(n_samples)]
    colors = colors or plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, axes = plt.subplots(n_samples, 1, figsize=figsize_local, dpi=dpi_local,
                             sharex=share_x_local, squeeze=False)
    axes = axes.flatten()

    for si, sample_files in enumerate(samples_list):
        ax = axes[si]
        color = colors[si % len(colors)]
        offsets = [i * offset_inc_within for i in range(len(sample_files))]
        lines = []
        for rep_i, (fpath, off) in enumerate(zip(sample_files, offsets)):
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
            x, y = read_spectrum(fpath)
            if smoothing and len(y) > 3:
                y = smooth_savgol(y, **savgol_kwargs)
            if normalize:
                mv = np.nanmax(np.abs(y))
                if mv != 0:
                    y = y / mv
            y = y + off
            ln, = ax.plot(x, y, color=color, alpha=0.8, linewidth=1.4)
            lines.append(ln)

        # styling for this subplot
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.set_ylabel(labels[si], fontname=font_name, fontsize=font_size)  # label at y-axis

        if xlim is not None:
            ax.set_xlim(xlim)
        if isinstance(ylim_per_sample, (list, tuple)) and len(ylim_per_sample) == n_samples:
            if ylim_per_sample[si] is not None:
                ax.set_ylim(ylim_per_sample[si])
        elif isinstance(ylim_per_sample, (list, tuple)) and len(ylim_per_sample) == 2:
            ax.set_ylim(ylim_per_sample)

    # common x label on bottom plot
    axes[-1].set_xlabel(xlabel, fontname=font_name, fontsize=font_size)

    # optional common y label (commented out to keep individual labels). Uncomment to use:
    # fig.text(0.04, 0.5, ylabel + (" (normalized)" if normalize else ""), va='center', rotation='vertical', fontsize=font_size)

    plt.tight_layout()
    return fig, axes

# -------------------------
# Example usage:
if __name__ == "__main__":
    all_samples = [sample_A_files, sample_B_files]  # add more lists for more stacked plots
    fig, axes = plot_stacked_samples(all_samples, labels=sample_labels, colors=sample_colors,
                                     normalize=normalize_each, smoothing=apply_smoothing,
                                     savgol_kwargs=savgol_params, offset_inc_within=offset_increment_within,
                                     xlim=x_range, ylim_per_sample=y_range_per_sample)
    # Save if wanted:
    # fig.savefig("stacked_samples.png", dpi=300)
    plt.show()
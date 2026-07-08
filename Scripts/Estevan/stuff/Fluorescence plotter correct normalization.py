import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
figsize = (7, 4)
dpi = 300

# Axis limits (None to autoscale)
x_range = (500, 700)    # e.g. (300, 700)
y_range = None

# Smoothing / normalization / offsets
apply_smoothing = True
savgol_params = {"window_length": 11, "polyorder": 3}  # window_length must be odd

# Normalization mode:
# "each"  -> normalize each spectrum to its own max (original behavior)
# "global"-> normalize all spectra to the maximum intensity across all files
# "first" -> normalize all spectra to the maximum intensity of the first file in `files`
normalize_mode = "global"   # choose: "each", "global", or "first"

offset_increment = 0.0   # progressive offsets: 0.0 (no offset), or e.g. 0.2

# === FILES TO PLOT (one spectrum per XLSX) ===
files = [
    r"E:\UT Summer 2026\Fluorescence data\FRET\AnnexinV_1uM.xlsx",
    r"E:\UT Summer 2026\Fluorescence data\FRET\MSP1D1_DMPG_1uM_w_1uM_ANNEXIN.xlsx",
    r"E:\UT Summer 2026\Fluorescence data\FRET\MSP1D1_sfCherryS_DMPG_1uM.xlsx",
    r"E:\UT Summer 2026\Fluorescence data\FRET\MSP1D1_sfCherryS_DMPG_1uM_plus1uM_Annexin.xlsx",
r"E:\UT Summer 2026\Fluorescence data\FRET\MSP1D1_sfCherryS_DMPC_1uM_w_1uM_ANNEXIN.xlsx",
    # add more paths...
]

# Optional sheet names or indices (use None for first sheet)
sheets = [
    None,
    None,
    None,
    None,
]

# Legend labels (one per file)
labels = [
    "AnnexinV",
    "MSP1D1-DMPG + Annexin",
    "MSP1D1-sfCherryS-DMPG",
    "MSP1D1-sfCherryS-DMPG + Annexin",
"MSP1D1-sfCherryS-DMPC + Annexin",
]

# Custom colors (will cycle if fewer than files)
colors = [
    "#D4D4D4",
    "#a97510",
    "#B407DB",
    "#E366B1",
    "#91666e",
]

# Column detection: preferred header names (case-insensitive)
wavelength_names = ["Wavelength (nm)", "Wavelength", "wavelength (nm)", "wavelength", "nm"]
intensity_names = ["Intensity (a.u.)", "Intensity", "intensity (a.u.)", "intensity", "a.u."]

# If header_row is None, pandas will try to infer headers; set to 0 if first row is header
header_row = 0

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

def find_columns(df):
    """Try to find wavelength and intensity columns by name, else use first two columns."""
    cols_lower = {c.lower(): c for c in df.columns.astype(str)}
    # try wavelength
    for name in wavelength_names:
        if name.lower() in cols_lower:
            wcol = cols_lower[name.lower()]
            break
    else:
        # fallback to first column
        wcol = df.columns[0]
    # try intensity
    for name in intensity_names:
        if name.lower() in cols_lower:
            icol = cols_lower[name.lower()]
            break
    else:
        # fallback to second column if exists, else first
        icol = df.columns[1] if df.shape[1] > 1 else df.columns[0]
    return wcol, icol

# Validate lists
n_files = len(files)
if len(sheets) != n_files:
    sheets = [None] * n_files
if len(labels) != n_files:
    labels = [os.path.basename(f) for f in files]

# offsets: progressive offsets using offset_increment
offsets = [i * offset_increment for i in range(n_files)]

# ------------------------
# PRELOAD all spectra so we can compute reference max if needed
spectra = []  # list of dicts: {"file":..., "x":..., "y":..., "label":..., "sheet":...}
for file, sheet, label in zip(files, sheets, labels):
    if not os.path.exists(file):
        raise FileNotFoundError(f"File not found: {file}")
    try:
        if sheet is None:
            df = pd.read_excel(file, header=header_row)
        else:
            df = pd.read_excel(file, sheet_name=sheet, header=header_row)
    except Exception as e:
        raise RuntimeError(f"Failed reading {file}: {e}")

    if df.shape[0] == 0:
        raise ValueError(f"No data in file: {file}")

    wcol, icol = find_columns(df)
    x = pd.to_numeric(df[wcol], errors='coerce').values
    y = pd.to_numeric(df[icol], errors='coerce').values
    mask = ~np.isnan(x) & ~np.isnan(y)
    x = x[mask]
    y = y[mask]

    if len(x) == 0:
        raise ValueError(f"No numeric wavelength/intensity pairs found in {file}")

    # smoothing (pre-plot)
    if apply_smoothing and len(y) > 3:
        y = smooth_savgol(y, **savgol_params)

    spectra.append({"file": file, "x": x, "y": y, "label": label, "sheet": sheet})

# Determine normalization reference
if normalize_mode == "each":
    # per-file normalization on plotting step (no global reference needed)
    reference_max = None
elif normalize_mode == "global":
    # maximum across all loaded spectra
    all_maxes = [np.nanmax(np.abs(s["y"])) if len(s["y"])>0 else 0.0 for s in spectra]
    reference_max = np.nanmax(all_maxes) if all_maxes else 1.0
    if reference_max == 0:
        reference_max = 1.0
elif normalize_mode == "first":
    # normalize to the first file's max
    if len(spectra) == 0:
        raise ValueError("No spectra loaded to determine 'first' normalization reference")
    first_max = np.nanmax(np.abs(spectra[0]["y"])) if len(spectra[0]["y"])>0 else 0.0
    reference_max = first_max if first_max != 0 else 1.0
else:
    raise ValueError("normalize_mode must be one of 'each', 'global', or 'first'")

# ------------------------
# Plot
fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
handles = []
legend_labels = []

for i, s in enumerate(spectra):
    x = s["x"]
    y = s["y"].copy()

    # normalization
    if normalize_mode == "each":
        maxv = np.nanmax(np.abs(y)) if len(y)>0 else 0.0
        if maxv != 0:
            y = y / maxv
    else:
        # use shared reference_max computed above
        y = y / reference_max

    # apply offset
    y = y + offsets[i]

    color = colors[i % len(colors)]
    line, = ax.plot(x, y, color=color, linewidth=1.5, label=s["label"])
    handles.append(line)
    legend_labels.append(s["label"])

# styling: remove gridlines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
if x_range is not None:
    ax.set_xlim(x_range)
if y_range is not None:
    ax.set_ylim(y_range)

ax.set_xlabel("Wavelength (nm)", fontname=font_name, fontsize=font_size)
norm_label = f" (normalized to {('each spectrum' if normalize_mode=='each' else ('first file' if normalize_mode=='first' else 'global max'))})"
ax.set_ylabel("Intensity (a.u.)" + (norm_label if normalize_mode else ""), fontname=font_name, fontsize=font_size)

# Turn off gridlines explicitly
ax.grid(False)

ax.legend(handles, legend_labels, loc="upper right", frameon=False, fontsize=font_size)
plt.tight_layout()
plt.show()
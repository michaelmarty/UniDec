import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ==========================
# USER SETTINGS
# ==========================

# Font settings (to match existing lipid-count plotting style)
font_name = "Arial"
font_size = 10

# Mass window around the most intense peak (Da)
window_Da = 50.0  # e.g., ±50 Da

# Normalization option
# If True: normalize intensity within the window (each trace goes 0–100% in its own window)
# If False: normalize by the global max intensity in the full spectrum
normalize_within_window = True

# Optional: path to save the figure (set to None to skip saving)
save_path = "zero_aligned_native_ms.png"  # or None if you only want to show the plot

# FOLDER PATH: Paste the full path to the folder containing all your mass text files
# Example: r"Z:\\Group Share\\Estevan\\native_ms\\mass_files"
folder_path = r"C:\Users\ec45294\Downloads\Mass Texts"  # FILL THIS WITH YOUR FOLDER PATH

# Optional: Specify a file extension filter (e.g., ".txt" or ".dat")
# Set to None to process all files in the folder
file_extension = ".txt"

# Optional: Manual ordering/labeling
# If you want to control the order and labels, fill this list with filenames (without path).
# The script will process files in this order and use these labels.
# Example:
# file_order = [
#     "protein1_mass.txt",
#     "protein2_mass.txt",
#     "protein3_mass.txt",
# ]
# If empty, the script will process all files in alphabetical order and use filenames as labels.
file_order = []

# ==========================
# PLOTTING SETUP
# ==========================

plt.rcParams.update({
    "font.family": font_name,
    "font.size": font_size,
})

# Adjusted figure size to accommodate legend on the right
fig, ax = plt.subplots(figsize=(8, 3), dpi=300)

# Predefined color cycle (extend or modify as needed)
colors = [
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
]

# ==========================
# BUILD FILE LIST
# ==========================

if not folder_path:
    print("[ERROR] No folder_path specified. Please set folder_path in the USER SETTINGS section.")
    exit(1)

if not os.path.isdir(folder_path):
    print(f"[ERROR] The specified folder does not exist: {folder_path}")
    exit(1)

# Get all files in the folder
all_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

# Filter by extension if specified
if file_extension:
    all_files = [f for f in all_files if f.endswith(file_extension)]

# If file_order is specified, use that order; otherwise sort alphabetically
if file_order:
    # Use the specified order, but verify files exist
    files_to_process = []
    for fname in file_order:
        if fname in all_files:
            files_to_process.append(fname)
        else:
            print(f"[WARNING] File '{fname}' specified in file_order but not found in folder. Skipping.")
else:
    # Process all files in alphabetical order
    files_to_process = sorted(all_files)

if not files_to_process:
    print("[ERROR] No files to process. Check folder_path, file_extension, and file_order settings.")
    exit(1)

print(f"Found {len(files_to_process)} file(s) to process:")
for fname in files_to_process:
    print(f"  - {fname}")

# ==========================
# MAIN LOOP OVER FILES
# ==========================

for i, fname in enumerate(files_to_process):
    file_path = os.path.join(folder_path, fname)

    # Use filename (without extension) as label by default
    label = os.path.splitext(fname)[0]

    try:
        # Read the UniDec eMass-style file: two columns, whitespace-separated, optional '#' comments
        df = pd.read_csv(
            file_path,
            comment="#",
            delim_whitespace=True,
            header=None,
            names=["Mass_Da", "Intensity"],
        )
    except Exception as e:
        print(f"[ERROR] Failed to read file '{fname}': {e}")
        continue

    # Drop any rows with missing or non-numeric values
    df = df.dropna(subset=["Mass_Da", "Intensity"]).copy()

    if df.empty:
        print(f"[WARNING] File '{fname}' is empty after cleaning. Skipping.")
        continue

    # Ensure numeric types
    df["Mass_Da"] = pd.to_numeric(df["Mass_Da"], errors="coerce")
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    df = df.dropna(subset=["Mass_Da", "Intensity"]).copy()

    if df.empty:
        print(f"[WARNING] File '{fname}' has no valid numeric data. Skipping.")
        continue

    # Identify the most intense peak (global maximum intensity)
    max_idx = df["Intensity"].idxmax()
    center_mass = df.loc[max_idx, "Mass_Da"]

    # Compute relative mass axis
    df["Rel_Mass_Da"] = df["Mass_Da"] - center_mass

    # Apply window around zero
    mask = (df["Rel_Mass_Da"] >= -window_Da) & (df["Rel_Mass_Da"] <= window_Da)
    df_win = df.loc[mask].copy()

    if df_win.empty:
        print(
            f"[WARNING] File '{fname}' has no data points within ±{window_Da} Da of the main peak. Skipping."
        )
        continue

    # Normalization
    if normalize_within_window:
        max_intensity = df_win["Intensity"].max()
    else:
        max_intensity = df["Intensity"].max()

    if max_intensity <= 0:
        print(f"[WARNING] File '{fname}' has non-positive max intensity. Skipping.")
        continue

    df_win["Percent"] = (df_win["Intensity"] / max_intensity) * 100.0

    # Sort by relative mass for cleaner plotting (optional but nice)
    df_win = df_win.sort_values("Rel_Mass_Da")

    # Choose color
    color = colors[i % len(colors)]

    # Plot trace
    ax.plot(
        df_win["Rel_Mass_Da"],
        df_win["Percent"],
        label=label,
        color=color,
        linewidth=1.0,
    )

# ==========================
# STYLING AND OUTPUT
# ==========================

# Vertical dashed line at zero (alignment point)
ax.axvline(0, color="k", linestyle="--", linewidth=0.8)

# Clean up spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Axis limits and labels
ax.set_xlim(-window_Da, window_Da)
ax.set_ylim(0, 100)
ax.set_yticks([0, 50, 100])

ax.set_xlabel("Relative mass (Da)")
ax.set_ylabel("Relative Intensity (%)")

# Legend positioned outside on the right
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)

plt.tight_layout()

# Save figure if requested
if save_path is not None:
    try:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved to: {save_path}")
    except Exception as e:
        print(f"[ERROR] Failed to save figure to '{save_path}': {e}")

# Show plot
plt.show()

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === USER SETTINGS ===
font_name = "Arial"
font_size = 10
figsize = (8, 5)
dpi = 300

# Plot limits (None to autoscale)
x_range = None
y_range = None

# Processing options
apply_smoothing = True
savgol_params = {"window_length": 11, "polyorder": 3}  # must be odd
normalize_each = True          # normalize each spectrum to its own max
offset_increment = 0.0         # progressive offset per spectrum (within group)

# Visual options
group_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]  # colors used per group (cycled)
within_group_alpha = 0.6      # alpha for individual spectra
within_group_linestyles = None  # e.g. ['-', '--', ':'] or None to use solid lines

# Column detection (case-insensitive)
wavelength_names = ["wavelength (nm)", "wavelength", "nm"]
intensity_names = ["intensity (a.u.)", "intensity", "a.u."]

# Excel reading defaults
header_row = 0   # set to None if your file has no header row

# -------------------- end of user settings --------------------

plt.rcParams.update({'font.size': font_size, 'font.family': font_name})

def smooth_savgol(y, window_length=11, polyorder=3):
    try:
        from |
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

# Set directory where files live
data_dir = "Z:\\Group Share\\Annika\\Stellar"

chrom_file = os.path.join(data_dir, "CHROMATOGRAM747.txt")
spectrum_file = os.path.join(data_dir, "MS747.txt")

# === Read data ===
chrom = pd.read_csv(chrom_file, delim_whitespace=True, header=None, names=["time", "intensity"])
spec = pd.read_csv(spectrum_file, delim_whitespace=True, header=None, names=["mz", "intensity"])

# === Smooth MS intensity to reduce tiny noise spikes ===
mz = spec["mz"].values
intensity = spec["intensity"].values
window = 3  # simple moving average
intensity_smooth = np.convolve(intensity, np.ones(window)/window, mode='same')

# === Pick peaks using prominence instead of just height ===
peaks, _ = find_peaks(intensity_smooth, prominence=intensity_smooth.max()*0.01)
mz_centroid = mz[peaks]
intensity_centroid = intensity[peaks]

# === Plot ===
fig, axes = plt.subplots(2, 1, figsize=(9, 10), sharex=False)  # taller, less wide

# Chromatogram (top)
axes[0].plot(chrom["time"], chrom["intensity"], color="blue", linewidth=2)


plt.xlim(150, 750)
axes[0].set_ylim(bottom=0)
axes[0].set_xlim(left=0)

# Spectrum (bottom) - centroided peaks
axes[1].stem(mz_centroid, intensity_centroid, linefmt='b-', markerfmt=' ', basefmt='k-')



axes[1].set_ylim(bottom=0)
axes[1].set_xlim(left=0)

plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl

mpl.use("WxAgg")

def plot_spec_match(row, exp_col="MS/MS spectrum", ref_col="Ref Spec", name_col="Metabolite name"):
    # Stick plot of the experimental and reference spectra with ref in negative butterfly style
    exp_spec = row[exp_col]
    ref_spec = row[ref_col]

    # Parse the spectra into m/z and intensity
    exp_mz, exp_intensity = zip(*[map(float, peak.split(':')) for peak in exp_spec.split(' ')])
    ref_mz, ref_intensity = zip(*[map(float, peak.split(':')) for peak in ref_spec.split(' ')])

    # Normalize the intensities for better visualization
    exp_intensity = np.array(exp_intensity) / max(exp_intensity)
    ref_intensity = np.array(ref_intensity) / max(ref_intensity)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.bar(exp_mz, exp_intensity, width=1, color='blue', label='Experimental', alpha=1)
    plt.bar(ref_mz, [-i for i in ref_intensity], width=1, color='red', label='Reference', alpha=1)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(row[name_col])
    plt.legend()
    plt.show()

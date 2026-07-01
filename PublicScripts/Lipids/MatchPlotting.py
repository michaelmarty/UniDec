import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl

mpl.use("WxAgg")

def butterfly_plot(exp_data, ref_data, title="Butterfly Plot"):
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.bar(exp_data[:, 0], exp_data[:, 1], width=1, color='blue', label='Experimental', alpha=1)
    plt.bar(ref_data[:, 0], [-i for i in ref_data[:, 1]], width=1, color='red', label='Reference', alpha=1)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(title)
    plt.legend()
    plt.show()

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

    butterfly_plot(np.column_stack([exp_mz, exp_intensity]), np.column_stack([ref_mz, ref_intensity]), title=row[name_col])

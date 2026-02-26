import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---------- Publication Font Settings ----------
rcParams["font.family"] = "Arial"
rcParams["font.size"] = 12

def plot_stokes_with_error(
    cal_RV,
    cal_D,
    sample_RV_means,
    sample_RV_stds,
    sample_labels=None,
    sample_colors=None,
    fig_title="Stokes Diameter Calibration",
):
    """
    Plot calibration (Stokes diameter vs retention volume)
    and sample points with error bars (black).
    """

    cal_RV = np.array(cal_RV)
    cal_D = np.array(cal_D)
    sample_RV_means = np.array(sample_RV_means)
    sample_RV_stds = np.array(sample_RV_stds)

    # ----- Fit the calibration curve -----
    # Fit RV = m*D + b
    m, b = np.polyfit(cal_D, cal_RV, 1)

    # ----- Predict retention volume for calibration line -----
    D_fit = np.linspace(min(cal_D), max(cal_D), 200)
    RV_fit = m * D_fit + b

    # ----- Convert sample RV → Stokes diameter -----
    # D = (RV - b) / m
    sample_D = (sample_RV_means - b) / m
    sample_D_err = sample_RV_stds / abs(m)

    # ----- Plot -----
    fig, ax = plt.subplots(figsize=(6, 4.5))

    # Calibration line
    ax.plot(D_fit, RV_fit, "--", color="black", lw=1.5, label="Calibration fit")

    # Scatter samples + error bars
    if sample_labels is None:
        sample_labels = [f"Sample {i+1}" for i in range(len(sample_D))]

    if sample_colors is None:
        sample_colors = ["tab:blue"] * len(sample_D)

    for Dv, Rv, d_err, label, color in zip(
        sample_D, sample_RV_means, sample_D_err, sample_labels, sample_colors
    ):
        ax.errorbar(
            Dv,
            Rv,
            xerr=d_err,
            linestyle="none",
            marker="o",
            markersize=6,
            color=color,
            ecolor="black",    # black error bars
            elinewidth=1.2,
            capsize=4,
            label=label,
        )

    # ----- Axis labels -----
    ax.set_xlabel("Stokes Diameter (nm)")
    ax.set_ylabel("Retention Volume (mL)")
    ax.set_title(fig_title)

    # ----- Clean formatting -----
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.legend(frameon=False)

    plt.tight_layout()
    return fig, ax

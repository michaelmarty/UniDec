import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


DATASETS = [
    (
        r"Z:\Group Share\BHT\Data\Q Exactive HF Data\Walti HSP60\20251211\20251211_BHT_0.1mgmL_cHSP60_CDMS_long_mix_1_unidecfiles\Swoop_mz__4360_z__15_chrom.txt",
        "m/z 4360, z 15",
    ),
    (
        r"Z:\Group Share\BHT\Data\Q Exactive HF Data\Walti HSP60\20251211\20251211_BHT_0.1mgmL_cHSP60_CDMS_long_mix_1_unidecfiles\mz__8353-10104_z__41-58_chrom.txt",
        "m/z 8353-10104, z 41-58",
    ),
    (
        r"Z:\Group Share\BHT\Data\Q Exactive HF Data\Walti HSP60\20251211\20251211_BHT_0.1mgmL_cHSP60_CDMS_long_mix_1_unidecfiles\mz__11492-13926_z__63-84_chrom.txt",
        "m/z 11492-13926, z 63-84",
    ),
]

OUTPUT_NAME = "mz_extract_chromatogram_overlay_sum_normalized"


def load_trace(file_path, label):
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"{label} must contain at least two columns: time and intensity.")
    return data[:, 0], data[:, 1]


def normalize_traces_to_sum(raw_traces):
    reference_time = raw_traces[0][1]
    for label, time, _ in raw_traces[1:]:
        if time.shape != reference_time.shape or not np.allclose(time, reference_time):
            raise ValueError(f"Time axis for {label} does not match the first trace.")

    intensity_sum = np.sum([intensity for _, _, intensity in raw_traces], axis=0)
    normalized_traces = []
    for label, time, intensity in raw_traces:
        normalized_intensity = np.divide(
            intensity,
            intensity_sum,
            out=np.zeros_like(intensity, dtype=float),
            where=intensity_sum > 0,
        )
        normalized_traces.append((label, time, normalized_intensity))

    return normalized_traces


def main():
    raw_traces = [(label, *load_trace(file_path, label)) for file_path, label in DATASETS]
    traces = normalize_traces_to_sum(raw_traces)

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 28
    plt.rcParams["xtick.major.size"] = 12
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 12
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.major.pad"] = 8
    plt.rcParams["ytick.major.pad"] = 8
    plt.rcParams["axes.labelpad"] = 18
    plt.rcParams["axes.titlepad"] = 18
    plt.rcParams["axes.linewidth"] = plt.rcParams["xtick.major.width"]
    plt.rcParams["pdf.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(16, 9))

    for label, time, intensity in traces:
        ax.plot(time, intensity, linewidth=2.5, label=label)

    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("Fraction of Summed Intensity")
    ax.set_ylim(0, 1.05)

    all_times = np.concatenate([time for _, time, _ in traces])
    ax.set_xlim(np.min(all_times), np.max(all_times))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    formatter = ScalarFormatter(useOffset=False, useMathText=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon=False)
    fig.tight_layout()

    output_dir = os.path.dirname(DATASETS[0][0])
    output_path_png = os.path.join(output_dir, f"{OUTPUT_NAME}.png")
    output_path_pdf = os.path.join(output_dir, f"{OUTPUT_NAME}.pdf")

    fig.savefig(output_path_png, dpi=300, transparent=True)
    fig.savefig(output_path_pdf, dpi=300)

    print(f"Plot saved as PNG at: {output_path_png}")
    print(f"Plot saved as PDF at: {output_path_pdf}")


if __name__ == "__main__":
    main()

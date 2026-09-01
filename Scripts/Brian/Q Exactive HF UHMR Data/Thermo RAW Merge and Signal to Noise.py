"""Calculate measured and theoretical S/N for 1 through N merged Thermo RAW TICs."""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import UniDecImporter.ImporterFactory as importer_factory

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "public" / "UniDec"))

from unidec.modules.HTEng import UniChromCDEng


RAW_FOLDERS = [
    Path(r"Z:\Group Share\BHT\Temp\Paper\DMT\Waters mAb 21 Injections Set 1"),
    Path(r"Z:\Group Share\BHT\Temp\Paper\DMT\Waters mAb 21 Injections Set 2"),
    Path(r"Z:\Group Share\BHT\Temp\Paper\DMT\Waters mAb 21 Injections Set 3"),
]
OUTPUT_FOLDER_NAME = "SN Calculation"
NOISE_WINDOW_MIN = (4.0, 5.0)
SIGNAL_WINDOW_MIN = (1.65, 2.2)


def natural_sort_key(path: Path) -> list[object]:
    """Sort numbered RAW files as 1, 2, ..., 10 instead of 1, 10, ..., 2."""
    return [int(part) if part.isdigit() else part.casefold() for part in re.split(r"(\d+)", path.name)]


def find_raw_files(folder: Path) -> list[Path]:
    raw_files = sorted(
        (path for path in folder.iterdir() if path.is_file() and path.suffix.casefold() == ".raw"),
        key=natural_sort_key,
    )
    if not raw_files:
        raise FileNotFoundError(f"No Thermo .RAW files found in {folder}")
    return raw_files


def read_tic(raw_path: Path) -> np.ndarray:
    """Read and validate a two-column [time, intensity] TIC."""
    importer = importer_factory.create_importer(str(raw_path))
    tic = np.asarray(importer.get_tic(), dtype=float)
    if tic.ndim != 2 or tic.shape[1] < 2 or tic.shape[0] < 2:
        raise ValueError(f"Invalid TIC returned for {raw_path.name}: shape {tic.shape}")

    tic = tic[:, :2]
    tic = tic[np.all(np.isfinite(tic), axis=1)]
    tic = tic[np.argsort(tic[:, 0])]
    unique_times, unique_indices = np.unique(tic[:, 0], return_index=True)
    tic = np.column_stack((unique_times, tic[unique_indices, 1]))
    if len(tic) < 2 or np.any(np.diff(tic[:, 0]) <= 0):
        raise ValueError(f"TIC times are not strictly increasing for {raw_path.name}")
    return tic


def make_common_time_axis(tics: list[np.ndarray]) -> np.ndarray:
    """Create a regular time axis over the range shared by every acquisition."""
    start = max(tic[0, 0] for tic in tics)
    stop = min(tic[-1, 0] for tic in tics)
    time_steps = np.concatenate([np.diff(tic[:, 0]) for tic in tics])
    time_steps = time_steps[np.isfinite(time_steps) & (time_steps > 0)]
    step = float(np.median(time_steps))
    point_count = int(np.floor((stop - start) / step)) + 1
    # UniDec's MRS irfft path returns one fewer point for odd-length cycles.
    # Drop only the final shared endpoint so every cycle has an even length.
    if point_count % 2:
        point_count -= 1
    if point_count < 2:
        raise ValueError(f"TICs do not have a usable shared time range: {start:g} to {stop:g} min")
    return start + np.arange(point_count, dtype=float) * step


def demultiplex_mrs(time_axis: np.ndarray, aligned_tics: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate aligned runs and fold them with UniDec's MRS demultiplexer."""
    multiplexed_intensity = aligned_tics.reshape(-1)
    time_step = float(np.median(np.diff(time_axis)))
    multiplexed_time = time_axis[0] + np.arange(multiplexed_intensity.size) * time_step

    engine = UniChromCDEng()
    engine.config.demultiplexmode = "MRS"
    engine.config.HTtimepad = 0
    engine.fulltime = multiplexed_time
    engine.fullscans = np.arange(multiplexed_intensity.size)
    engine.fulltic = multiplexed_intensity
    engine.setup_demultiplex(cycleindex=len(time_axis))
    demultiplexed, _ = engine.run_demultiplex(multiplexed_intensity, normalize=False)
    demultiplexed = np.asarray(demultiplexed, dtype=float)

    # For a regular MRS kernel, the FFT demultiplexer must equal a direct sum
    # across acquisitions. Keep this assertion to catch ordering or padding errors.
    direct_sum = np.sum(aligned_tics, axis=0)
    if not np.allclose(demultiplexed, direct_sum, rtol=1e-10, atol=1e-4):
        max_error = float(np.max(np.abs(demultiplexed - direct_sum)))
        raise RuntimeError(f"MRS demultiplexing validation failed (maximum error {max_error:g})")
    return np.column_stack((multiplexed_time, multiplexed_intensity)), demultiplexed


def window_mask(time_axis: np.ndarray, bounds: tuple[float, float], label: str) -> np.ndarray:
    mask = (time_axis >= bounds[0]) & (time_axis <= bounds[1])
    if not np.any(mask):
        raise ValueError(f"No TIC points found in the {label} window {bounds[0]}-{bounds[1]} min")
    return mask


def calculate_signal_to_noise(time_axis: np.ndarray, intensity: np.ndarray) -> dict[str, float | int | list[float] | str]:
    noise_mask = window_mask(time_axis, NOISE_WINDOW_MIN, "noise")
    signal_mask = window_mask(time_axis, SIGNAL_WINDOW_MIN, "signal")
    noise_values = intensity[noise_mask]
    signal_times = time_axis[signal_mask]
    signal_values = intensity[signal_mask]

    noise_mean = float(np.mean(noise_values))
    noise_std = float(np.std(noise_values, ddof=0))
    if noise_std == 0:
        raise ZeroDivisionError("Noise standard deviation is zero; S/N cannot be calculated")

    peak_index = int(np.argmax(signal_values))
    peak_signal = float(signal_values[peak_index])
    baseline_corrected_peak = peak_signal - noise_mean
    return {
        "snr_definition": "(peak_signal - noise_mean) / noise_standard_deviation",
        "noise_window_min": list(NOISE_WINDOW_MIN),
        "signal_window_min": list(SIGNAL_WINDOW_MIN),
        "noise_point_count": int(noise_values.size),
        "signal_point_count": int(signal_values.size),
        "noise_mean": noise_mean,
        "noise_standard_deviation": noise_std,
        "noise_3sigma": 3.0 * noise_std,
        "peak_signal": peak_signal,
        "peak_time_min": float(signal_times[peak_index]),
        "baseline_corrected_peak_signal": baseline_corrected_peak,
        "signal_to_noise": baseline_corrected_peak / noise_std,
        "peak_over_noise_standard_deviation": peak_signal / noise_std,
        "baseline_corrected_peak_over_3sigma": baseline_corrected_peak / (3.0 * noise_std),
    }


def dataset_label(raw_folder: Path) -> str:
    match = re.search(r"Set\s+\d+$", raw_folder.name, flags=re.IGNORECASE)
    return match.group(0).title() if match else raw_folder.name


def analyze_dataset(raw_folder: Path) -> dict[str, object]:
    label = dataset_label(raw_folder)
    raw_files = find_raw_files(raw_folder)
    print(f"\n{label}: found {len(raw_files)} Thermo RAW files in {raw_folder}")

    tics: list[np.ndarray] = []
    for index, raw_path in enumerate(raw_files, start=1):
        print(f"[{label} {index:02d}/{len(raw_files):02d}] Reading {raw_path.name}")
        tics.append(read_tic(raw_path))

    time_axis = make_common_time_axis(tics)
    aligned_tics = np.vstack([np.interp(time_axis, tic[:, 0], tic[:, 1]) for tic in tics])
    subset_tics: list[np.ndarray] = []
    subset_results: list[dict[str, object]] = []
    for injection_count in range(1, len(raw_files) + 1):
        _, demultiplexed_intensity = demultiplex_mrs(time_axis, aligned_tics[:injection_count])
        result: dict[str, object] = calculate_signal_to_noise(time_axis, demultiplexed_intensity)
        result.update(
            {
                "injection_count": injection_count,
                "last_raw_file_added": raw_files[injection_count - 1].name,
                "raw_files_in_subset": [path.name for path in raw_files[:injection_count]],
            }
        )
        subset_tics.append(demultiplexed_intensity)
        subset_results.append(result)
        print(f"  {label}, {injection_count:02d} injection(s): S/N {float(result['signal_to_noise']):.6f}")

    return {
        "label": label,
        "raw_folder": raw_folder,
        "raw_files": raw_files,
        "time_axis": time_axis,
        "subset_tics": np.vstack(subset_tics),
        "subset_results": subset_results,
    }


def save_replicate_outputs(output_folder: Path, datasets: list[dict[str, object]]) -> dict[str, np.ndarray]:
    output_folder.mkdir(exist_ok=True)
    result_counts = {len(dataset["subset_results"]) for dataset in datasets}
    if len(result_counts) != 1:
        raise ValueError(f"Replicate sets contain different injection counts: {sorted(result_counts)}")

    injection_count = result_counts.pop()
    counts = np.arange(1, injection_count + 1, dtype=float)
    replicate_snr = np.array(
        [
            [float(result["signal_to_noise"]) for result in dataset["subset_results"]]
            for dataset in datasets
        ]
    )
    mean_snr = np.mean(replicate_snr, axis=0)
    sd_snr = np.std(replicate_snr, axis=0, ddof=1)
    sem_snr = sd_snr / np.sqrt(len(datasets))
    theoretical_snr = mean_snr[0] * np.sqrt(counts)

    for dataset in datasets:
        label = str(dataset["label"])
        safe_label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_").lower()
        time_axis = np.asarray(dataset["time_axis"])
        subset_tics = np.asarray(dataset["subset_tics"])
        tic_path = output_folder / f"demultiplexed_tics_{safe_label}.csv"
        with tic_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["Time_min", *[f"TIC_{count:02d}_merged_injections" for count in range(1, injection_count + 1)]]
            )
            writer.writerows(np.column_stack((time_axis, subset_tics.T)))

        merged_tic = subset_tics[-1]
        merged_result = dataset["subset_results"][-1]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(time_axis, merged_tic, color="#6A4C93", linewidth=1.8, label="Demultiplexed TIC")
        ax.axvspan(*NOISE_WINDOW_MIN, color="#4C956C", alpha=0.18, label="Noise window")
        ax.axvspan(*SIGNAL_WINDOW_MIN, color="#F2C14E", alpha=0.22, label="Signal window")
        ax.axvline(
            float(merged_result["peak_time_min"]),
            color="#D1495B",
            linestyle="--",
            linewidth=1.2,
            label="Peak",
        )
        ax.set(
            xlabel="Time (min)",
            ylabel="Summed TIC intensity",
            title=f"{label} Merged TIC: S/N = {float(merged_result['signal_to_noise']):.3f}",
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_folder / f"demultiplexed_merged_tic_{safe_label}.png", dpi=300)
        fig.savefig(output_folder / f"demultiplexed_merged_tic_{safe_label}.pdf")
        plt.close(fig)

    long_fields = [
        "replicate",
        "raw_folder",
        "injection_count",
        "experimental_signal_to_noise",
        "noise_mean",
        "noise_standard_deviation",
        "noise_3sigma",
        "peak_signal",
        "baseline_corrected_peak_signal",
        "peak_time_min",
        "last_raw_file_added",
    ]
    with (output_folder / "signal_to_noise_replicates_long.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=long_fields)
        writer.writeheader()
        for dataset in datasets:
            for result in dataset["subset_results"]:
                writer.writerow(
                    {
                        "replicate": dataset["label"],
                        "raw_folder": dataset["raw_folder"],
                        "injection_count": result["injection_count"],
                        "experimental_signal_to_noise": result["signal_to_noise"],
                        "noise_mean": result["noise_mean"],
                        "noise_standard_deviation": result["noise_standard_deviation"],
                        "noise_3sigma": result["noise_3sigma"],
                        "peak_signal": result["peak_signal"],
                        "baseline_corrected_peak_signal": result["baseline_corrected_peak_signal"],
                        "peak_time_min": result["peak_time_min"],
                        "last_raw_file_added": result["last_raw_file_added"],
                    }
                )

    summary_fields = [
        "injection_count",
        "replicate_count",
        "mean_experimental_signal_to_noise",
        "sample_standard_deviation",
        "standard_error_of_mean",
        "minimum_experimental_signal_to_noise",
        "maximum_experimental_signal_to_noise",
        "theoretical_signal_to_noise",
        "mean_over_theoretical",
        "mean_percent_of_theoretical",
    ]
    with (output_folder / "signal_to_noise_replicate_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        for index, count in enumerate(counts):
            writer.writerow(
                {
                    "injection_count": int(count),
                    "replicate_count": len(datasets),
                    "mean_experimental_signal_to_noise": mean_snr[index],
                    "sample_standard_deviation": sd_snr[index],
                    "standard_error_of_mean": sem_snr[index],
                    "minimum_experimental_signal_to_noise": np.min(replicate_snr[:, index]),
                    "maximum_experimental_signal_to_noise": np.max(replicate_snr[:, index]),
                    "theoretical_signal_to_noise": theoretical_snr[index],
                    "mean_over_theoretical": mean_snr[index] / theoretical_snr[index],
                    "mean_percent_of_theoretical": 100.0 * mean_snr[index] / theoretical_snr[index],
                }
            )

    report = {
        "output_folder": str(output_folder),
        "replicate_count": len(datasets),
        "subset_definition": "Cumulative natural file order within each replicate set",
        "snr_definition": "(peak_signal - noise_mean) / noise_standard_deviation",
        "noise_window_min": list(NOISE_WINDOW_MIN),
        "signal_window_min": list(SIGNAL_WINDOW_MIN),
        "theoretical_curve_definition": "mean one-injection S/N across replicate sets multiplied by sqrt(injection_count)",
        "theoretical_curve_equation": "S/N(N) = mean replicate S/N(1) * sqrt(N)",
        "variability_definition": "Sample standard deviation across replicate sets (ddof=1)",
        "merge_method": "interpolate each TIC to a replicate-specific shared time grid, concatenate cumulative cycles, and demultiplex with UniChromCDEng MRS mode",
        "datasets": [
            {
                "label": dataset["label"],
                "raw_folder": str(dataset["raw_folder"]),
                "raw_files_in_cumulative_order": [path.name for path in dataset["raw_files"]],
                "common_time_start_min": float(np.asarray(dataset["time_axis"])[0]),
                "common_time_stop_min": float(np.asarray(dataset["time_axis"])[-1]),
                "common_time_step_min": float(np.median(np.diff(np.asarray(dataset["time_axis"])))),
                "common_time_point_count": int(np.asarray(dataset["time_axis"]).size),
                "results": dataset["subset_results"],
            }
            for dataset in datasets
        ],
        "summary": [
            {
                "injection_count": int(counts[index]),
                "replicate_signal_to_noise": [float(value) for value in replicate_snr[:, index]],
                "mean_experimental_signal_to_noise": float(mean_snr[index]),
                "sample_standard_deviation": float(sd_snr[index]),
                "standard_error_of_mean": float(sem_snr[index]),
                "theoretical_signal_to_noise": float(theoretical_snr[index]),
                "mean_percent_of_theoretical": float(100.0 * mean_snr[index] / theoretical_snr[index]),
            }
            for index in range(injection_count)
        ],
    }
    with (output_folder / "signal_to_noise_replicates.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    colors = ["#6A4C93", "#2A9D8F", "#E9A23B"]
    fig, ax = plt.subplots(figsize=(10, 6))
    for index, dataset in enumerate(datasets):
        ax.plot(
            counts,
            replicate_snr[index],
            color=colors[index % len(colors)],
            marker="o",
            markersize=3.5,
            linewidth=1.3,
            alpha=0.75,
            label=str(dataset["label"]),
        )
    ax.fill_between(counts, mean_snr - sd_snr, mean_snr + sd_snr, color="#34495E", alpha=0.15, label="Mean ± SD")
    ax.plot(counts, mean_snr, color="#263238", marker="o", markersize=5, linewidth=2.5, label="Replicate mean")
    ax.plot(
        counts,
        theoretical_snr,
        color="#D1495B",
        linestyle="--",
        linewidth=2.2,
        label=r"Theoretical S/N = mean S/N$_1\sqrt{N}$",
    )
    ax.set(
        xlabel="Number of merged injections",
        ylabel="Signal-to-noise ratio",
        title="Signal-to-Noise vs. Number of Merged Injections",
        xticks=counts.astype(int),
        xlim=(0.5, injection_count + 0.5),
    )
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output_folder / "signal_to_noise_replicates_by_merged_injections.png", dpi=300)
    fig.savefig(output_folder / "signal_to_noise_replicates_by_merged_injections.pdf")
    plt.close(fig)

    return {
        "counts": counts,
        "replicate_snr": replicate_snr,
        "mean_snr": mean_snr,
        "sd_snr": sd_snr,
        "theoretical_snr": theoretical_snr,
    }


def main() -> None:
    datasets = [analyze_dataset(raw_folder) for raw_folder in RAW_FOLDERS]
    output_folder = RAW_FOLDERS[0] / OUTPUT_FOLDER_NAME
    summary = save_replicate_outputs(output_folder, datasets)

    mean_snr = summary["mean_snr"]
    sd_snr = summary["sd_snr"]
    theoretical_snr = summary["theoretical_snr"]
    print(f"\nMean one-injection S/N: {mean_snr[0]:.6f} ± {sd_snr[0]:.6f} SD")
    print(f"Mean twenty-one-injection S/N: {mean_snr[-1]:.6f} ± {sd_snr[-1]:.6f} SD")
    print(f"Theoretical twenty-one-injection S/N: {theoretical_snr[-1]:.6f}")
    print(f"Outputs saved in {output_folder}")


if __name__ == "__main__":
    main()

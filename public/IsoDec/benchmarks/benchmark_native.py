"""Measure end-to-end native IsoDec deconvolution performance."""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np

from isodec import IsoDecWrapper, __version__


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True, help="Compiler/build label")
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=7)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.warmups < 0 or args.iterations < 1:
        raise ValueError("warmups must be nonnegative and iterations must be positive")

    spectrum = np.loadtxt(args.spectrum)
    wrapper = IsoDecWrapper()

    expected_peak_count = None
    for _ in range(args.warmups):
        expected_peak_count = len(wrapper.process_spectrum(spectrum))

    samples = []
    for _ in range(args.iterations):
        started = time.perf_counter_ns()
        peak_count = len(wrapper.process_spectrum(spectrum))
        samples.append((time.perf_counter_ns() - started) / 1_000_000)
        if expected_peak_count is None:
            expected_peak_count = peak_count
        elif peak_count != expected_peak_count:
            raise RuntimeError(
                f"Non-deterministic peak count: {peak_count} != {expected_peak_count}"
            )

    median_ms = statistics.median(samples)
    result = {
        "label": args.label,
        "isodec_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "spectrum_points": len(spectrum),
        "peak_count": expected_peak_count,
        "warmups": args.warmups,
        "iterations": args.iterations,
        "median_ms": round(median_ms, 6),
        "mean_ms": round(statistics.mean(samples), 6),
        "min_ms": round(min(samples), 6),
        "stdev_ms": round(statistics.stdev(samples), 6) if len(samples) > 1 else 0.0,
        "spectra_per_second": round(1000 / median_ms, 6),
        "samples_ms": [round(sample, 6) for sample in samples],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    json.dump(result, sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    main()

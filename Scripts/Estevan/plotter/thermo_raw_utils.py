# thermo_raw_utils.py
from __future__ import annotations
import os
import shutil
import subprocess
import sys
from typing import List, Tuple, Optional, Dict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Optional imports
try:
    from pymsfilereader import MSFileReader  # Windows + Thermo MSFileReader
    PYMS_AVAILABLE = True
except Exception:
    PYMS_AVAILABLE = False

try:
    import pymzml  # for reading mzML after conversion
    PYMZML_AVAILABLE = True
except Exception:
    PYMZML_AVAILABLE = False

# Consistent sizing
DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN     = 1.35   # per subplot
MIN_HEIGHT_IN     = 3.8    # minimum total height


def _apply_style(font_name: str, font_size: int, tick_length: float = 5.0, tick_width: float = 1.0):
    mpl.rcParams.update({
        "font.family": font_name,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "axes.linewidth": 1.2,
        "xtick.major.size": tick_length,
        "ytick.major.size": tick_length,
        "xtick.major.width": tick_width,
        "ytick.major.width": tick_width,
        "legend.frameon": False,
        "legend.borderpad": 0.2,
        "legend.handlelength": 2.0,
        "legend.handletextpad": 0.6,
    })


def _find_msconvert(env_path: Optional[str] = None) -> Optional[str]:
    """
    Locate msconvert in this priority:
      1) Provided msconvert path (env_path)
      2) Environment variable PROTEOWIZARD_MSCONVERT
      3) PATH via shutil.which
      4) Common Windows install locations
    """
    candidates = []
    if env_path:
        candidates.append(env_path)

    env_var = os.environ.get("PROTEOWIZARD_MSCONVERT")
    if env_var:
        candidates.append(env_var)

    exe = shutil.which("msconvert")
    if exe:
        candidates.append(exe)

    candidates += [
        r"C:\Program Files\ProteoWizard\ProteoWizard 3.0\msconvert.exe",
        r"C:\Program Files\ProteoWizard\msconvert.exe",
        r"C:\Program Files (x86)\ProteoWizard\msconvert.exe",
    ]

    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def _convert_raw_to_mzml(raw_path: str, out_dir: Optional[str] = None, msconvert_path: Optional[str] = None) -> str:
    msconvert = _find_msconvert(msconvert_path)
    if not msconvert:
        raise RuntimeError(
            "msconvert not found. Install ProteoWizard and either:\n"
            "- add msconvert.exe to PATH,\n"
            "- set PROTEOWIZARD_MSCONVERT env var,\n"
            "- or specify msconvert path in the GUI."
        )
    if out_dir is None:
        out_dir = os.path.dirname(raw_path)
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(raw_path))[0]
    out_file = f"{base}.mzML"
    cmd = [msconvert, raw_path, "--mzML", "--outfile", out_file, "--outdir", out_dir]
    subprocess.run(cmd, check=True)
    mzml_path = os.path.join(out_dir, out_file)
    if not os.path.isfile(mzml_path):
        raise RuntimeError("msconvert finished but mzML file not found.")
    return mzml_path


def _read_mzml_scans(mzml_path: str, ms_level: int = 1, max_scans: Optional[int] = None) -> List[Dict[str, np.ndarray]]:
    if not PYMZML_AVAILABLE:
        raise RuntimeError("pymzml not installed. pip install pymzml")
    run = pymzml.run.Reader(mzml_path)
    scans: List[Dict[str, np.ndarray]] = []
    for spec in run:
        if getattr(spec, "ms_level", None) != ms_level:
            continue
        rt = getattr(spec, "scan_time_in_minutes", None)
        if rt is None:
            try:
                st = spec.get("scan start time")
                rt = float(st) if st is not None else np.nan
            except Exception:
                rt = np.nan
        mzs = np.asarray(spec.mz, dtype=float)
        ints = np.asarray(spec.i, dtype=float)
        if mzs.size and ints.size:
            scans.append({"rt": np.array(rt, dtype=float), "mz": mzs, "int": ints})
        if max_scans and len(scans) >= max_scans:
            break
    return scans


def _read_raw_windows_scans(raw_path: str, ms_level: int = 1, max_scans: Optional[int] = 100) -> List[Dict[str, np.ndarray]]:
    if not PYMS_AVAILABLE:
        raise RuntimeError("pymsfilereader not available. pip install pymsfilereader and install Thermo MSFileReader.")
    reader = MSFileReader(raw_path)
    n = int(reader.GetNumSpectra())
    scans: List[Dict[str, np.ndarray]] = []
    count = 0
    for sn in range(1, n + 1):
        try:
            hdr = reader.GetScanHeaderInfoForScanNum(sn)
        except Exception:
            continue
        ms_order = int(hdr[0]) if hdr and isinstance(hdr[0], (int, float)) else 1
        if ms_order != ms_level:
            continue
        rt = np.nan
        for method in ("RTFromScanNum", "GetRetentionTimeFromScanNum"):
            try:
                val = getattr(reader, method)(sn)
                rt = float(val)
                break
            except Exception:
                continue
        try:
            mz, intens = reader.GetMassListFromScanNum(sn)
        except Exception:
            continue
        mz = np.asarray(mz, dtype=float)
        intens = np.asarray(intens, dtype=float)
        if mz.size and intens.size:
            scans.append({"rt": np.array(rt, dtype=float), "mz": mz, "int": intens})
            count += 1
            if max_scans and count >= max_scans:
                break
    return scans


def load_raw_scans(
    raw_path: str,
    ms_level: int = 1,
    max_scans: Optional[int] = 200,
    prefer_native: bool = True,
    msconvert_path: Optional[str] = None,
) -> List[Dict[str, np.ndarray]]:
    """
    Load scans from Thermo.raw or.mzML file.
    prefer_native=True: on Windows, try native MSFileReader first; otherwise fallback to msconvert + pymzml.
    msconvert_path: optional explicit path to msconvert.exe (overrides PATH/env).
    Returns a list of dicts: {"rt": float array (minutes), "mz": array, "int": array}
    """
    ext = os.path.splitext(raw_path)[1].lower()
    if ext == ".mzml":
        return _read_mzml_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
    if ext == ".raw":
        if sys.platform.startswith("win") and prefer_native:
            if PYMS_AVAILABLE:
                try:
                    return _read_raw_windows_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
                except Exception:
                    # If native fails, fall back to msconvert
                    pass
        # Fallback: convert to mzML, then read
        mzml_path = _convert_raw_to_mzml(raw_path, msconvert_path=msconvert_path)
        return _read_mzml_scans(mzml_path, ms_level=ms_level, max_scans=max_scans)
    raise ValueError(f"Unsupported file extension: {ext}. Use.raw (Thermo) or.mzML.")


def compute_tic(scans: List[Dict[str, np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]:
    if not scans:
        return np.array([], dtype=float), np.array([], dtype=float)
    rt = np.array([float(s["rt"]) if np.size(s["rt"]) else np.nan for s in scans], dtype=float)
    tic = np.array([float(np.nan_to_num(s["int"]).sum()) for s in scans], dtype=float)
    mask = np.isfinite(rt)
    return rt[mask], tic[mask]


def pick_spectrum(scans: List[Dict[str, np.ndarray]], mode: str = "apex", index: Optional[int] = None) -> Dict[str, np.ndarray]:
    if not scans:
        return {"rt": np.array(np.nan), "mz": np.array([]), "int": np.array([])}
    if mode == "index" and index is not None:
        i = max(0, min(index, len(scans) - 1))
        return scans[i]
    rt, tic = compute_tic(scans)
    if rt.size == 0:
        return scans[0]
    apex_rt = rt[int(np.argmax(tic))]
    diffs = np.abs(np.array([float(s["rt"]) for s in scans]) - apex_rt)
    i = int(np.argmin(diffs))
    return scans[i]


def plot_raw_chrom_and_spectrum(
    raw_path: str,
    ms_level: int = 1,
    max_scans: Optional[int] = 200,
    spectrum_mode: str = "apex",
    spectrum_index: Optional[int] = None,
    normalize_spectrum: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    dpi: int = 300,
    font_name: str = "Arial",
    font_size: int = 10,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    save: Optional[str] = None,
    show: bool = False,
    prefer_native: bool = True,
    msconvert_path: Optional[str] = None,
):
    """
    Load Thermo RAW or mzML, compute TIC, pick a representative spectrum, and plot both.

    Styling and consistency:
      - Unified fonts and tick styling.
      - Spectrum normalization to percent (optional) with y-ticks $$[0, \%, 100]$$.
      - Fixed width; comfortable height for 2-row layout.
      - Clean spines.

    Returns (fig, (ax_chrom, ax_spec)).
    """
    _apply_style(font_name, font_size, tick_length=tick_length, tick_width=tick_width)

    scans = load_raw_scans(
        raw_path,
        ms_level=ms_level,
        max_scans=max_scans,
        prefer_native=prefer_native,
        msconvert_path=msconvert_path,
    )
    rt, tic = compute_tic(scans)
    spec = pick_spectrum(scans, mode=spectrum_mode, index=spectrum_index)

    # Respect provided width; ensure comfortable height for 2 rows
    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = max(MIN_HEIGHT_IN, (ROW_HEIGHT_IN * 2) + 0.8)

    fig = plt.figure(figsize=(width_in, target_h), dpi=dpi)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.45)

    # Chromatogram (TIC)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(rt, tic, color="#2C7FB8", lw=1.8)
    ax1.set_xlabel("Time (min)")
    ax1.set_ylabel("Total Ion Current (a.u.)")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.margins(y=0)

    # m/z spectrum
    ax2 = fig.add_subplot(gs[1, 0])
    mz = spec.get("mz", np.array([]))
    intens = spec.get("int", np.array([]))

    if normalize_spectrum and intens.size:
        m = float(np.nanmax(intens))
        y = (intens / m) * 100.0 if m > 0 else intens
        ax2.set_ylim(0, 100)
        ax2.set_yticks([0, 50, 100])
        ax2.set_yticklabels(["0", "%", "100"])
        ax2.set_ylabel("Relative Intensity (%)")
    else:
        y = intens
        if y.size:
            ymin, ymax = float(np.nanmin(y)), float(np.nanmax(y))
            ax2.set_ylim(max(0.0, ymin), max(1.0, ymax))
        ax2.set_ylabel("Intensity (a.u.)")

    ax2.plot(mz, y, color="#CC79A7", lw=1.8)
    ax2.set_xlabel("m/z")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.margins(y=0)

    # Optional, minimal title (kept; remove if you prefer no title)
    base = os.path.basename(raw_path)
    where = "apex" if spectrum_mode == "apex" else f"index {spectrum_index}"
    rt_val = float(spec["rt"]) if np.size(spec.get("rt", np.array(np.nan))) else np.nan
    fig.suptitle(f"{base} — TIC and m/z (scan: {where}, rt≈{rt_val:.2f} min)", fontsize=font_size + 2)

    # Layout
    plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96))

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, (ax1, ax2)
# thermo_raw_utils.py
from __future__ import annotations
import os
import shutil
import subprocess
import sys
from typing import List, Tuple, Optional, Dict

import numpy as np
import matplotlib.pyplot as plt

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


def _find_msconvert() -> Optional[str]:
    exe = shutil.which("msconvert")
    if exe:
        return exe
    candidates = [
        r"C:\Program Files\ProteoWizard\ProteoWizard 3.0.23117\msconvert.exe",
        r"C:\Program Files\ProteoWizard\msconvert.exe",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None

def _convert_raw_to_mzml(raw_path: str, out_dir: Optional[str] = None) -> str:
    msconvert = _find_msconvert()
    if not msconvert:
        raise RuntimeError("msconvert not found. Install ProteoWizard and add it to PATH.")
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

def load_raw_scans(raw_path: str, ms_level: int = 1, max_scans: Optional[int] = 200) -> List[Dict[str, np.ndarray]]:
    ext = os.path.splitext(raw_path)[1].lower()
    if ext == ".mzml":
        return _read_mzml_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
    if ext == ".raw":
        if sys.platform.startswith("win") and PYMS_AVAILABLE:
            return _read_raw_windows_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
        mzml_path = _convert_raw_to_mzml(raw_path)
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
    spectrum_mode: str = "apex",  # "apex" or "index"
    spectrum_index: Optional[int] = None,
    normalize_spectrum: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    dpi: int = 120,
    font_name: str = "Arial",
    font_size: int = 10,
    save: Optional[str] = None,
    show: bool = False,
):
    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    scans = load_raw_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
    rt, tic = compute_tic(scans)
    spec = pick_spectrum(scans, mode=spectrum_mode, index=spectrum_index)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.4)

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(rt, tic, color="#2C7FB8", lw=1.8)
    ax1.set_xlabel("Retention Time (min)")
    ax1.set_ylabel("TIC (sum intensity)")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2 = fig.add_subplot(gs[1, 0])
    mz = spec["mz"]; intens = spec["int"]
    if normalize_spectrum and intens.size:
        m = np.nanmax(intens)
        y = (intens / m) * 100.0 if m > 0 else intens
        ax2.set_ylabel("Relative Intensity (%)")
    else:
        y = intens
        ax2.set_ylabel("Intensity")
    ax2.plot(mz, y, color="#CC79A7", lw=1.5)
    ax2.set_xlabel("m/z")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    base = os.path.basename(raw_path)
    where = "apex" if spectrum_mode == "apex" else f"index {spectrum_index}"
    rt_val = float(spec["rt"]) if np.size(spec["rt"]) else np.nan
    fig.suptitle(f"{base} — TIC and m/z (scan: {where}, rt≈{rt_val:.2f} min)", fontsize=font_size + 2)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, (ax1, ax2)
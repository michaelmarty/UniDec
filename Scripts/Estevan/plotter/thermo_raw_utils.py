# thermo_raw_utils.py (updated - streaming reader + conversion + helpers)
from __future__ import annotations
import os
import shutil
import subprocess
import sys
from typing import List, Tuple, Optional, Dict, Generator

import numpy as np
import matplotlib.pyplot as plt

# Optional native reader (pymsfilereader)
try:
    from pymsfilereader import MSFileReader
    PYMS_AVAILABLE = True
except Exception:
    PYMS_AVAILABLE = False

# pymzml for mzML reading
try:
    import pymzml
    PYMZML_AVAILABLE = True
except Exception:
    PYMZML_AVAILABLE = False

from plot_helpers import style_context, compute_target_height, rc_for_style

# Sizing constants
DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN = 1.35
MIN_HEIGHT_IN = 3.8


def _find_msconvert(env_path: Optional[str] = None) -> Optional[str]:
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


def _convert_raw_to_mzml(
    raw_path: str,
    out_dir: Optional[str] = None,
    msconvert_path: Optional[str] = None,
    ms_level_filter: Optional[str] = None,
    peakpick: bool = True,
    compress: bool = True,
) -> str:
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
    if compress:
        cmd.append("--zlib")
    if peakpick:
        mslevel_arg = ms_level_filter if ms_level_filter is not None else "1-"
        cmd.extend(["--filter", f"peakPicking true {mslevel_arg}"])
    if ms_level_filter:
        cmd.extend(["--filter", f"msLevel {ms_level_filter}"])
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else ""
        stdout = exc.stdout.strip() if exc.stdout else ""
        msg = f"msconvert failed (exit {exc.returncode}).\nstdout:\n{stdout}\nstderr:\n{stderr}"
        raise RuntimeError(msg) from exc
    except Exception as exc:
        raise RuntimeError(f"Failed to run msconvert: {exc}") from exc
    mzml_path = os.path.join(out_dir, out_file)
    if not os.path.isfile(mzml_path):
        raise RuntimeError("msconvert finished but mzML file not found.")
    return mzml_path


def _read_mzml_scans_list(mzml_path: str, ms_level: int = 1, max_scans: Optional[int] = None) -> List[Dict[str, np.ndarray]]:
    """
    Non-streaming reader kept for compatibility (returns list).
    """
    if not PYMZML_AVAILABLE:
        raise RuntimeError("pymzml not installed. pip install pymzml")
    reader = pymzml.run.Reader(mzml_path)
    out = []
    count = 0
    for spec in reader:
        spec_ms = getattr(spec, "ms_level", None)
        if spec_ms is None:
            try:
                spec_ms = int(spec['ms level']) if 'ms level' in spec else 1
            except Exception:
                spec_ms = 1
        if int(spec_ms) != int(ms_level):
            continue
        rt = getattr(spec, "scan_time_in_minutes", None)
        if rt is None:
            try:
                rt_val = spec.get("scan start time") if hasattr(spec, "get") else None
                rt = float(rt_val) if rt_val is not None else np.nan
            except Exception:
                rt = np.nan
        try:
            mz_array = np.asarray(spec.mz, dtype=float)
            int_array = np.asarray(spec.i, dtype=float)
        except Exception:
            continue
        if mz_array.size and int_array.size:
            out.append({"rt": np.array(rt, dtype=float), "mz": mz_array, "int": int_array})
            count += 1
            if max_scans and count >= max_scans:
                break
    return out


def read_mzml_scans_stream(mzml_path: str, ms_level: int = 1) -> Generator[Dict[str, np.ndarray], None, None]:
    """
    Streaming generator that yields MS-level scans from an mzML using pymzml.
    Yields dicts: {"rt": float, "mz": np.ndarray, "int": np.ndarray}
    This avoids storing the whole file in memory.
    """
    if not PYMZML_AVAILABLE:
        raise RuntimeError("pymzml not installed. pip install pymzml")
    reader = pymzml.run.Reader(mzml_path)
    for spec in reader:
        spec_ms = getattr(spec, "ms_level", None)
        if spec_ms is None:
            try:
                spec_ms = int(spec['ms level']) if 'ms level' in spec else 1
            except Exception:
                spec_ms = 1
        if int(spec_ms) != int(ms_level):
            continue
        rt = getattr(spec, "scan_time_in_minutes", None)
        if rt is None:
            try:
                rt_val = spec.get("scan start time") if hasattr(spec, "get") else None
                rt = float(rt_val) if rt_val is not None else np.nan
            except Exception:
                rt = np.nan
        try:
            mz_array = np.asarray(spec.mz, dtype=float)
            int_array = np.asarray(spec.i, dtype=float)
        except Exception:
            continue
        if mz_array.size and int_array.size:
            yield {"rt": np.array(rt, dtype=float), "mz": mz_array, "int": int_array}


def _read_raw_windows_scans(raw_path: str, ms_level: int = 1, max_scans: Optional[int] = 100) -> List[Dict[str, np.ndarray]]:
    """
    Native Thermo MSFileReader path (unchanged).
    """
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
    peakpick_on_convert: bool = True,
    compress_mzml: bool = True,
) -> List[Dict[str, np.ndarray]]:
    """
    Load scans from.raw or.mzML, returning a list (compatibility wrapper).
    Uses native reader if available, otherwise converts with msconvert and reads mzML.
    """
    ext = os.path.splitext(raw_path)[1].lower()
    if ext == ".mzml":
        return _read_mzml_scans_list(raw_path, ms_level=ms_level, max_scans=max_scans)
    if ext == ".raw":
        if sys.platform.startswith("win") and prefer_native:
            if PYMS_AVAILABLE:
                try:
                    return _read_raw_windows_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
                except Exception:
                    pass
        mzml_path = _convert_raw_to_mzml(
            raw_path,
            msconvert_path=msconvert_path,
            ms_level_filter=str(ms_level) if ms_level is not None else None,
            peakpick=bool(peakpick_on_convert),
            compress=bool(compress_mzml),
        )
        return _read_mzml_scans_list(mzml_path, ms_level=ms_level, max_scans=max_scans)
    raise ValueError(f"Unsupported file extension: {ext}. Use.raw or.mzML.")


def compute_tic(scans: List[Dict[str, np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]:
    if not scans:
        return np.array([], dtype=float), np.array([], dtype=float)
    rt = np.array([float(s["rt"]) if np.size(s["rt"]) else np.nan for s in scans], dtype=float)
    tic = np.array([float(np.nan_to_num(s["int"]).sum()) for s in scans], dtype=float)
    mask = np.isfinite(rt)
    return rt[mask], tic[mask]


def compute_tic_stream(mzml_path: str, ms_level: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute TIC from an mzML by streaming (low memory).
    Returns (rt_array, tic_array).
    """
    rts = []
    tics = []
    for s in read_mzml_scans_stream(mzml_path, ms_level=ms_level):
        rts.append(float(s.get("rt", np.nan)))
        tics.append(float(np.nansum(s["int"])))
    if not rts:
        return np.array([], dtype=float), np.array([], dtype=float)
    rts_arr = np.array(rts, dtype=float)
    tics_arr = np.array(tics, dtype=float)
    mask = np.isfinite(rts_arr)
    return rts_arr[mask], tics_arr[mask]


def pick_spectrum_from_stream(mzml_path: str, ms_level: int = 1, method: str = "apex", index: Optional[int] = None) -> Dict[str, np.ndarray]:
    """
    Pick a representative spectrum from an mzML by streaming:
    - method 'apex' finds the scan with maximum TIC
    - method 'index' uses the given index
    Returns a single scan dict.
    """
    if method == "index" and index is not None:
        # iterate until index-th matching ms_level
        i = 0
        for s in read_mzml_scans_stream(mzml_path, ms_level=ms_level):
            if i == index:
                return s
            i += 1
        return {"rt": np.array(np.nan), "mz": np.array([]), "int": np.array([])}

    # Apex: compute TIC on the fly and keep the best scan (single-pass)
    best_tic = -np.inf
    best_scan = None
    for s in read_mzml_scans_stream(mzml_path, ms_level=ms_level):
        cur_tic = float(np.nansum(s["int"]))
        if cur_tic > best_tic:
            best_tic = cur_tic
            best_scan = s
    if best_scan is None:
        return {"rt": np.array(np.nan), "mz": np.array([]), "int": np.array([])}
    return best_scan


def plot_raw_chrom_and_spectrum(
    raw_path: str,
    ms_level: int = 1,
    max_scans: Optional[int] = 200,
    spectrum_mode: str = "apex",
    spectrum_index: Optional[int] = None,
    normalize_spectrum: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    # NEW: optional y max for the normalized spectrum axis
    spec_y_max: Optional[float] = None,
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
    peakpick_on_convert: bool = True,
    compress_mzml: bool = True,
):
    """
    High-level plot function that converts (if needed), reads, computes TIC and picks spectrum.
    Uses streaming helpers internally to avoid large memory usage.
    """
    rc = rc_for_style(font_name, font_size, tick_length=tick_length, tick_width=tick_width)

    ext = os.path.splitext(raw_path)[1].lower()
    if ext == ".raw" and (not (sys.platform.startswith("win") and prefer_native and PYMS_AVAILABLE)):
        # convert then stream
        mzml_path = _convert_raw_to_mzml(
            raw_path,
            msconvert_path=msconvert_path,
            ms_level_filter=str(ms_level) if ms_level is not None else None,
            peakpick=bool(peakpick_on_convert),
            compress=bool(compress_mzml),
        )
    elif ext == ".raw" and (sys.platform.startswith("win") and prefer_native and PYMS_AVAILABLE):
        # Use native first if available; otherwise convert
        try:
            scans = _read_raw_windows_scans(raw_path, ms_level=ms_level, max_scans=max_scans)
            rt, tic = compute_tic(scans)
            spec = pick_spectrum(scans, mode=spectrum_mode, index=spectrum_index)
            # Plot using same logic below (we have scans list)
            # Reuse below plotting block by setting mzml_path=None and using scans/spec directly
            use_native = True
        except Exception:
            mzml_path = _convert_raw_to_mzml(
                raw_path,
                msconvert_path=msconvert_path,
                ms_level_filter=str(ms_level) if ms_level is not None else None,
                peakpick=bool(peakpick_on_convert),
                compress=bool(compress_mzml),
            )
            use_native = False
    else:
        mzml_path = raw_path
        use_native = False

    # If we used native already and have scans, plot directly
    if 'use_native' in locals() and use_native:
        # Plot using rc context
        with plt.rc_context(rc):
            width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
            target_h = max(MIN_HEIGHT_IN, (ROW_HEIGHT_IN * 2) + 0.8)
            fig = plt.figure(figsize=(width_in, target_h), dpi=dpi)
            gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.45)

            ax1 = fig.add_subplot(gs[0, 0])
            if rt.size and tic.size:
                ax1.plot(rt, tic, color="#2C7FB8", lw=1.8)
            ax1.set_xlabel("Time (min)")
            ax1.set_ylabel("Total Ion Current (a.u.)")
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_visible(False)
            ax1.margins(y=0)

            ax2 = fig.add_subplot(gs[1, 0])
            mz = spec.get("mz", np.array([]))
            intens = spec.get("int", np.array([]))

            if normalize_spectrum and intens.size:
                m = float(np.nanmax(intens)) if intens.size else 0.0
                if m > 0:
                    if norm_mode == "percent":
                        y = (intens / m) * 100.0
                    else:
                        y = intens / m
                else:
                    y = intens if intens.size else np.array([])

                # Default normalized ymax per mode
                if norm_mode == "percent":
                    default_ymax = 100.0
                else:
                    default_ymax = 1.0

                if spec_y_max is None:
                    ymax = default_ymax
                else:
                    if norm_mode == "percent":
                        ymax = max(0.0, min(float(spec_y_max), 100.0))
                    else:
                        ymax = max(0.0, min(float(spec_y_max), 1.0))

                ax2.set_ylim(0.0, ymax)
                mid = ymax / 2.0
                if norm_mode == "percent":
                    ax2.set_yticks([0.0, mid, ymax])
                    ax2.set_yticklabels(["0", "%", f"{ymax:.0f}"])
                else:
                    ax2.set_yticks([0.0, mid, ymax])
                    ax2.set_yticklabels(["0", f"{mid:.2g}", f"{ymax:.2g}"])
                ax2.set_ylabel("Normalized Intensity")
            else:
                y = intens
                if y.size:
                    ymin, ymax = float(np.nanmin(y)), float(np.nanmax(y))
                    ax2.set_ylim(max(0.0, ymin), max(1.0, ymax))
                ax2.set_ylabel("Intensity (a.u.)")

            if mz.size and getattr(y, "size", None):
                ax2.plot(mz, y, color="#CC79A7", lw=1.8)
            ax2.set_xlabel("m/z")
            ax2.spines["top"].set_visible(False)
            ax2.spines["right"].set_visible(False)
            ax2.margins(y=0)

            base = os.path.basename(raw_path)
            where = "apex" if spectrum_mode == "apex" else f"index {spectrum_index}"
            rt_val = float(spec.get("rt", np.nan)) if np.size(spec.get("rt", np.array([]))) else np.nan
            fig.text(0.5, 0.98, f"{base} — TIC and m/z (scan: {where}, rt≈{rt_val:.2f} min)",
                     ha="center", va="top", fontsize=font_size + 2)

            plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96))

            if save:
                fig.savefig(save, dpi=dpi, bbox_inches="tight")
            if show:
                plt.show()

        return fig, (ax1, ax2)

    # Otherwise use mzml_path and streaming helpers for low-memory processing
    rt, tic = compute_tic_stream(mzml_path, ms_level=ms_level)
    spec = pick_spectrum_from_stream(mzml_path, ms_level=ms_level, method=spectrum_mode, index=spectrum_index)

    with plt.rc_context(rc):
        width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
        target_h = max(MIN_HEIGHT_IN, (ROW_HEIGHT_IN * 2) + 0.8)
        fig = plt.figure(figsize=(width_in, target_h), dpi=dpi)
        gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.45)

        ax1 = fig.add_subplot(gs[0, 0])
        if rt.size and tic.size:
            ax1.plot(rt, tic, color="#2C7FB8", lw=1.8)
        ax1.set_xlabel("Time (min)")
        ax1.set_ylabel("Total Ion Current (a.u.)")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.margins(y=0)

        ax2 = fig.add_subplot(gs[1, 0])
        mz = spec.get("mz", np.array([]))
        intens = spec.get("int", np.array([]))

        if normalize_spectrum and intens.size:
            m = float(np.nanmax(intens)) if intens.size else 0.0
            if m > 0:
                if norm_mode == "percent":
                    y = (intens / m) * 100.0
                else:
                    y = intens / m
            else:
                y = intens if intens.size else np.array([])

            if norm_mode == "percent":
                default_ymax = 100.0
            else:
                default_ymax = 1.0

            if spec_y_max is None:
                ymax = default_ymax
            else:
                if norm_mode == "percent":
                    ymax = max(0.0, min(float(spec_y_max), 100.0))
                else:
                    ymax = max(0.0, min(float(spec_y_max), 1.0))

            ax2.set_ylim(0.0, ymax)
            mid = ymax / 2.0
            if norm_mode == "percent":
                ax2.set_yticks([0.0, mid, ymax])
                ax2.set_yticklabels(["0", "%", f"{ymax:.0f}"])
            else:
                ax2.set_yticks([0.0, mid, ymax])
                ax2.set_yticklabels(["0", f"{mid:.2g}", f"{ymax:.2g}"])
            ax2.set_ylabel("Normalized Intensity")
        else:
            y = intens
            if y.size:
                ymin, ymax = float(np.nanmin(y)), float(np.nanmax(y))
                ax2.set_ylim(max(0.0, ymin), max(1.0, ymax))
            ax2.set_ylabel("Intensity (a.u.)")

        if mz.size and getattr(y, "size", None):
            ax2.plot(mz, y, color="#CC79A7", lw=1.8)
        ax2.set_xlabel("m/z")
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.margins(y=0)

        base = os.path.basename(raw_path)
        where = "apex" if spectrum_mode == "apex" else f"index {spectrum_index}"
        rt_val = float(spec.get("rt", np.nan)) if np.size(spec.get("rt", np.array([]))) else np.nan
        fig.text(0.5, 0.98, f"{base} — TIC and m/z (scan: {where}, rt≈{rt_val:.2f} min)",
                 ha="center", va="top", fontsize=font_size + 2)

        plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96))

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, (ax1, ax2)
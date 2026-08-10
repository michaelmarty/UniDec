"""
Interactive TIC peak integrator for UniChromCD TIC_decon.txt files.

The integration is intentionally simple: each peak area is the sum of the
y-values whose x-values fall inside the peak's left/right boundaries.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD

    TkBase = TkinterDnD.Tk
except Exception:
    DND_FILES = None
    TkBase = tk.Tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter


DEFAULT_TIC_PATH = ""


THEME_COLORS = {
    "Light": {
        "root_bg": "#f4f6f8",
        "panel_bg": "#f4f6f8",
        "frame_bg": "#f4f6f8",
        "text": "#17212b",
        "muted_text": "#4d5b68",
        "entry_bg": "#ffffff",
        "button_bg": "#e6ebf0",
        "button_active": "#d8e0e8",
        "select_bg": "#2f6f9f",
        "select_fg": "#ffffff",
        "tree_bg": "#ffffff",
        "tree_heading_bg": "#e6ebf0",
        "plot_bg": "#ffffff",
        "figure_bg": "#ffffff",
        "spine": "#17212b",
        "grid": "#d7dde4",
        "tic": "#17324d",
        "fill_selected": "#f08c2e",
        "line_selected": "#c4422f",
        "fill_unselected": "#708a9f",
        "line_unselected": "#708a9f",
    },
    "Dark": {
        "root_bg": "#1f2227",
        "panel_bg": "#252a31",
        "frame_bg": "#252a31",
        "text": "#edf0f3",
        "muted_text": "#b9c0c8",
        "entry_bg": "#15181d",
        "button_bg": "#343a43",
        "button_active": "#424a55",
        "select_bg": "#4d8ac8",
        "select_fg": "#ffffff",
        "tree_bg": "#181c22",
        "tree_heading_bg": "#343a43",
        "plot_bg": "#15181d",
        "figure_bg": "#15181d",
        "spine": "#d8dde3",
        "grid": "#343a43",
        "tic": "#8cc8ff",
        "fill_selected": "#f0a650",
        "line_selected": "#ff9b73",
        "fill_unselected": "#7f93a8",
        "line_unselected": "#7f93a8",
    },
}


def detect_system_theme() -> str:
    """Return the OS app theme on Windows, falling back to Light elsewhere."""

    try:
        import winreg

        key_path = r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize"
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, key_path) as key:
            apps_use_light_theme, _ = winreg.QueryValueEx(key, "AppsUseLightTheme")
        return "Light" if int(apps_use_light_theme) else "Dark"
    except Exception:
        return "Light"


def split_drop_paths(widget: tk.Widget, raw_data: str) -> list[str]:
    """Parse a TkDND file-list payload into regular filesystem paths."""

    try:
        paths = widget.tk.splitlist(raw_data)
    except tk.TclError:
        paths = raw_data.split()
    return [str(path).strip().strip("{}") for path in paths if str(path).strip()]


def normalize_tic_paths(paths: Iterable[str | os.PathLike[str]]) -> list[str]:
    """Return existing TIC file paths, accepting direct files or folders containing TIC_decon.txt."""

    normalized: list[str] = []
    seen: set[str] = set()
    for raw_path in paths:
        path = Path(os.fspath(raw_path).strip().strip('"').strip("'"))
        if path.is_dir():
            direct_tic = path / "TIC_decon.txt"
            if direct_tic.exists():
                path = direct_tic
            else:
                continue
        if not path.is_file():
            continue
        resolved = str(path)
        key = os.path.normcase(os.path.abspath(resolved))
        if key not in seen:
            normalized.append(resolved)
            seen.add(key)
    return normalized


@dataclass
class IntegrationWindow:
    """A fixed x-range used for single-file or batch TIC integration."""

    label: str
    left_x: float
    right_x: float

    def normalized_bounds(self) -> tuple[float, float]:
        left = float(min(self.left_x, self.right_x))
        right = float(max(self.left_x, self.right_x))
        return left, right

    def point_mask(self, x: np.ndarray) -> np.ndarray:
        left, right = self.normalized_bounds()
        return (x >= left) & (x <= right)

    def integrate(self, x: np.ndarray, y: np.ndarray) -> dict[str, float | int | None]:
        mask = self.point_mask(x)
        if not np.any(mask):
            return {"area": 0.0, "points": 0, "apex_x": None, "height": None}

        indices = np.flatnonzero(mask)
        local_y = y[mask]
        apex_index = int(indices[np.argmax(local_y)])
        return {
            "area": float(np.sum(local_y)),
            "points": int(len(local_y)),
            "apex_x": float(x[apex_index]),
            "height": float(y[apex_index]),
        }


@dataclass
class PeakIntegration:
    """Editable integration region for one TIC peak."""

    apex_index: int
    left_x: float
    right_x: float
    label: str
    enabled: bool = True

    def normalized_bounds(self) -> tuple[float, float]:
        left = float(min(self.left_x, self.right_x))
        right = float(max(self.left_x, self.right_x))
        return left, right

    def point_mask(self, x: np.ndarray) -> np.ndarray:
        left, right = self.normalized_bounds()
        return (x >= left) & (x <= right)

    def area(self, x: np.ndarray, y: np.ndarray) -> float:
        mask = self.point_mask(x)
        if not np.any(mask):
            return 0.0
        return float(np.sum(y[mask]))

    def height(self, y: np.ndarray) -> float:
        if 0 <= self.apex_index < len(y):
            return float(y[self.apex_index])
        return 0.0

    def apex_x(self, x: np.ndarray) -> float:
        if 0 <= self.apex_index < len(x):
            return float(x[self.apex_index])
        return 0.0


def load_tic_file(path: str | os.PathLike[str]) -> tuple[np.ndarray, np.ndarray]:
    """Load a UniChromCD TIC text file and return clean x/y arrays."""

    path = os.fspath(path)
    errors: list[str] = []
    for delimiter in (None, "\t", ","):
        try:
            data = np.genfromtxt(
                path,
                comments="#",
                delimiter=delimiter,
                invalid_raise=False,
                dtype=float,
            )
        except Exception as exc:
            errors.append(str(exc))
            continue

        if data.size == 0:
            continue

        data = np.atleast_2d(data)
        if data.shape[1] < 2:
            continue

        data = data[:, :2]
        data = data[np.isfinite(data).all(axis=1)]
        if len(data) < 2:
            continue

        order = np.argsort(data[:, 0])
        data = data[order]
        return data[:, 0].astype(float), data[:, 1].astype(float)

    detail = "; ".join(errors[-2:]) if errors else "no numeric two-column data found"
    raise ValueError(f"Could not read TIC data from {path!r}: {detail}")


def moving_average(y: np.ndarray, window: int) -> np.ndarray:
    """Return a centered moving average without changing the array length."""

    window = int(max(1, window))
    if window <= 1 or len(y) < 3:
        return y.copy()

    if window > len(y):
        window = len(y)
    if window % 2 == 0:
        window -= 1
    if window <= 1:
        return y.copy()

    kernel = np.ones(window, dtype=float) / float(window)
    padded = np.pad(y, window // 2, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def _fallback_find_peaks(
    x: np.ndarray,
    y: np.ndarray,
    min_prominence: float,
    min_distance_x: float,
) -> np.ndarray:
    """Small local-maxima fallback used if scipy is not available."""

    if len(y) < 3:
        return np.array([], dtype=int)

    baseline = float(np.nanpercentile(y, 5))
    candidates = np.flatnonzero((y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:])) + 1
    candidates = candidates[y[candidates] >= baseline + min_prominence]
    if len(candidates) == 0:
        return np.array([], dtype=int)

    ordered = sorted(candidates, key=lambda idx: y[idx], reverse=True)
    kept: list[int] = []
    for idx in ordered:
        if all(abs(float(x[idx] - x[other])) >= min_distance_x for other in kept):
            kept.append(int(idx))
    return np.array(sorted(kept), dtype=int)


def find_peak_indices(
    x: np.ndarray,
    y: np.ndarray,
    min_prominence: float,
    min_distance_x: float,
) -> np.ndarray:
    """Find peak apex indices using scipy when available."""

    if len(y) < 3:
        return np.array([], dtype=int)

    min_prominence = float(max(0.0, min_prominence))
    min_distance_x = float(max(0.0, min_distance_x))
    dx = float(np.nanmedian(np.diff(x))) if len(x) > 2 else 1.0
    if not np.isfinite(dx) or dx <= 0:
        dx = 1.0
    distance_points = max(1, int(round(min_distance_x / dx)))

    try:
        from scipy.signal import find_peaks

        peaks, _ = find_peaks(y, prominence=min_prominence, distance=distance_points)
        return peaks.astype(int)
    except Exception:
        return _fallback_find_peaks(x, y, min_prominence, min_distance_x)


def _initial_bounds_for_peak(
    y_for_bounds: np.ndarray,
    apex_index: int,
    boundary_fraction: float,
) -> tuple[int, int]:
    """Walk away from the apex until the trace drops below a simple threshold."""

    baseline = float(np.nanpercentile(y_for_bounds, 5))
    peak_height = float(y_for_bounds[apex_index] - baseline)
    if peak_height <= 0:
        return max(0, apex_index - 1), min(len(y_for_bounds) - 1, apex_index + 1)

    threshold = baseline + float(boundary_fraction) * peak_height
    left = int(apex_index)
    right = int(apex_index)

    while left > 0 and y_for_bounds[left] > threshold:
        left -= 1
    while right < len(y_for_bounds) - 1 and y_for_bounds[right] > threshold:
        right += 1

    if left == right:
        left = max(0, apex_index - 1)
        right = min(len(y_for_bounds) - 1, apex_index + 1)
    return left, right


def detect_peak_integrations(
    x: np.ndarray,
    y: np.ndarray,
    min_prominence: float,
    min_distance_x: float,
    boundary_fraction: float,
    smooth_points: int,
) -> list[PeakIntegration]:
    """Detect TIC peaks and create editable integration regions."""

    if len(x) != len(y) or len(x) < 3:
        return []

    y_smooth = moving_average(y, smooth_points)
    peak_indices = find_peak_indices(x, y_smooth, min_prominence, min_distance_x)
    if len(peak_indices) == 0:
        return []

    boundary_fraction = float(np.clip(boundary_fraction, 0.0, 1.0))
    left_indices: list[int] = []
    right_indices: list[int] = []

    for apex_index in peak_indices:
        left, right = _initial_bounds_for_peak(y_smooth, int(apex_index), boundary_fraction)
        left_indices.append(left)
        right_indices.append(right)

    for i in range(len(peak_indices) - 1):
        if right_indices[i] < left_indices[i + 1]:
            continue

        left_apex = int(peak_indices[i])
        right_apex = int(peak_indices[i + 1])
        if right_apex <= left_apex:
            continue

        valley_offset = int(np.argmin(y_smooth[left_apex : right_apex + 1]))
        valley_index = left_apex + valley_offset
        right_indices[i] = min(right_indices[i], valley_index)
        left_indices[i + 1] = max(left_indices[i + 1], valley_index)

    peaks: list[PeakIntegration] = []
    for number, (apex_index, left_index, right_index) in enumerate(
        zip(peak_indices, left_indices, right_indices),
        start=1,
    ):
        peaks.append(
            PeakIntegration(
                apex_index=int(apex_index),
                left_x=float(x[left_index]),
                right_x=float(x[right_index]),
                label=f"Peak {number}",
            )
        )
    return peaks


class TICTIntegratorApp(TkBase):
    """Tkinter/matplotlib GUI for TIC peak detection and manual integration."""

    def __init__(self) -> None:
        super().__init__()
        self.title("CDMS TIC Integrator")
        self.geometry("1450x900")
        self.minsize(1100, 700)

        self.x: np.ndarray | None = None
        self.y: np.ndarray | None = None
        self.file_path = tk.StringVar(value=DEFAULT_TIC_PATH)
        self.peaks: list[PeakIntegration] = []
        self.selected_peak_index: int | None = None
        self.dragging_boundary: str | None = None
        self.last_cursor_x: float | None = None
        self.batch_window: BatchProcessingWindow | None = None
        self.style = ttk.Style(self)

        self._build_variables()
        self.apply_theme(self.theme_var.get(), redraw=False)
        self._build_layout()
        self._bind_plot_events()
        self._enable_drag_and_drop()
        self.apply_theme(self.theme_var.get(), redraw=False)

        if DEFAULT_TIC_PATH and Path(DEFAULT_TIC_PATH).exists():
            self.load_file(DEFAULT_TIC_PATH)
        else:
            self._set_status("Open a TIC_decon.txt file to begin.")

    def _build_variables(self) -> None:
        self.theme_var = tk.StringVar(value=detect_system_theme())
        self.prominence_var = tk.DoubleVar(value=5.0)
        self.distance_var = tk.DoubleVar(value=0.20)
        self.boundary_fraction_var = tk.DoubleVar(value=0.10)
        self.smooth_points_var = tk.IntVar(value=3)
        self.left_boundary_var = tk.StringVar()
        self.right_boundary_var = tk.StringVar()
        self.selected_area_var = tk.StringVar(value="")
        self.status_var = tk.StringVar(value="")

    def _build_layout(self) -> None:
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        controls = ttk.Frame(self, padding=8)
        controls.grid(row=0, column=0, sticky="ns")
        controls.columnconfigure(0, weight=1)

        plot_area = ttk.Frame(self, padding=(0, 8, 8, 8))
        plot_area.grid(row=0, column=1, sticky="nsew")
        plot_area.columnconfigure(0, weight=1)
        plot_area.rowconfigure(0, weight=1)

        self._build_file_controls(controls)
        self._build_detection_controls(controls)
        self._build_peak_table(controls)
        self._build_boundary_controls(controls)
        self._build_plot(plot_area)

        status = ttk.Label(self, textvariable=self.status_var, anchor="w", padding=(8, 4))
        status.grid(row=1, column=0, columnspan=2, sticky="ew")

    def _build_file_controls(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="TIC File", padding=8)
        frame.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        frame.columnconfigure(0, weight=1)

        self.file_entry = ttk.Entry(frame, textvariable=self.file_path, width=48)
        self.file_entry.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 6))

        ttk.Button(frame, text="Browse", command=self.browse_file).grid(row=1, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(frame, text="Load", command=self.load_current_path).grid(row=1, column=1, sticky="ew")

        theme_frame = ttk.Frame(frame)
        theme_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        theme_frame.columnconfigure(1, weight=1)
        ttk.Label(theme_frame, text="Theme").grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.theme_combo = ttk.Combobox(
            theme_frame,
            textvariable=self.theme_var,
            values=("Light", "Dark"),
            state="readonly",
            width=10,
        )
        self.theme_combo.grid(row=0, column=1, sticky="ew")
        self.theme_combo.bind("<<ComboboxSelected>>", self.on_theme_change)

    def _build_detection_controls(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Peak Detection", padding=8)
        frame.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        frame.columnconfigure(1, weight=1)

        rows = [
            ("Prominence", self.prominence_var),
            ("Distance (min)", self.distance_var),
            ("Boundary Fraction", self.boundary_fraction_var),
            ("Smoothing Points", self.smooth_points_var),
        ]
        for row, (label, variable) in enumerate(rows):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=2)
            ttk.Entry(frame, textvariable=variable, width=12).grid(row=row, column=1, sticky="ew", pady=2)

        button_row = len(rows)
        ttk.Button(frame, text="Auto Scale", command=self.auto_scale_detection).grid(
            row=button_row,
            column=0,
            sticky="ew",
            padx=(0, 4),
            pady=(8, 0),
        )
        ttk.Button(frame, text="Detect Peaks", command=self.detect_peaks).grid(
            row=button_row,
            column=1,
            sticky="ew",
            pady=(8, 0),
        )

    def _build_peak_table(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Integrations", padding=8)
        frame.grid(row=2, column=0, sticky="nsew", pady=(0, 8))
        parent.rowconfigure(2, weight=1)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = ("apex", "left", "right", "area", "height")
        self.peak_table = ttk.Treeview(
            frame,
            columns=columns,
            show="tree headings",
            selectmode="browse",
            height=16,
        )
        self.peak_table.heading("#0", text="Peak")
        self.peak_table.heading("apex", text="Apex")
        self.peak_table.heading("left", text="Left")
        self.peak_table.heading("right", text="Right")
        self.peak_table.heading("area", text="Area")
        self.peak_table.heading("height", text="Height")
        self.peak_table.column("#0", width=70, stretch=False)
        self.peak_table.column("apex", width=70, anchor="e", stretch=False)
        self.peak_table.column("left", width=70, anchor="e", stretch=False)
        self.peak_table.column("right", width=70, anchor="e", stretch=False)
        self.peak_table.column("area", width=95, anchor="e", stretch=False)
        self.peak_table.column("height", width=80, anchor="e", stretch=False)
        self.peak_table.grid(row=0, column=0, sticky="nsew")
        self.peak_table.bind("<<TreeviewSelect>>", self.on_table_select)

        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=self.peak_table.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.peak_table.configure(yscrollcommand=scrollbar.set)

        buttons = ttk.Frame(frame)
        buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        buttons.columnconfigure((0, 1, 2), weight=1)
        ttk.Button(buttons, text="Add View", command=self.add_peak_from_view).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Delete", command=self.delete_selected_peak).grid(
            row=0,
            column=1,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Clear", command=self.clear_peaks).grid(row=0, column=2, sticky="ew")

        ttk.Button(frame, text="Export CSV", command=self.export_csv).grid(
            row=2,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(8, 0),
        )
        ttk.Button(frame, text="Batch Processing", command=self.open_batch_window).grid(
            row=3,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(6, 0),
        )

    def _build_boundary_controls(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Selected Peak", padding=8)
        frame.grid(row=3, column=0, sticky="ew")
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Left").grid(row=0, column=0, sticky="w", pady=2)
        ttk.Entry(frame, textvariable=self.left_boundary_var, width=12).grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Label(frame, text="Right").grid(row=1, column=0, sticky="w", pady=2)
        ttk.Entry(frame, textvariable=self.right_boundary_var, width=12).grid(row=1, column=1, sticky="ew", pady=2)
        ttk.Label(frame, text="Area").grid(row=2, column=0, sticky="w", pady=2)
        ttk.Label(frame, textvariable=self.selected_area_var).grid(row=2, column=1, sticky="e", pady=2)

        ttk.Button(frame, text="Apply Boundaries", command=self.apply_boundary_entries).grid(
            row=3,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(8, 0),
        )
        ttk.Button(frame, text="Use View Boundaries", command=self.use_view_boundaries).grid(
            row=4,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(6, 0),
        )

    def _build_plot(self, parent: ttk.Frame) -> None:
        colors = self.current_colors()
        self.figure = Figure(figsize=(9, 6), dpi=100, facecolor=colors["figure_bg"])
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel("Time (min)")
        self.ax.set_ylabel("Intensity (# of ions)")
        self._style_axes()
        self.figure.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")

        self.toolbar = NavigationToolbar2Tk(self.canvas, parent, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.grid(row=1, column=0, sticky="ew")

    def _bind_plot_events(self) -> None:
        self.canvas.mpl_connect("button_press_event", self.on_plot_press)
        self.canvas.mpl_connect("button_release_event", self.on_plot_release)
        self.canvas.mpl_connect("motion_notify_event", self.on_plot_motion)

    def _enable_drag_and_drop(self) -> None:
        enabled = False
        enabled |= self.register_file_drop_target(self, self.on_main_file_drop)
        enabled |= self.register_file_drop_target(self.file_entry, self.on_main_file_drop)
        enabled |= self.register_file_drop_target(self.canvas_widget, self.on_main_file_drop)
        if enabled:
            self._set_status("Open, browse, or drop a TIC_decon.txt file to begin.")

    def register_file_drop_target(self, widget: tk.Widget, callback: Callable[[list[str]], None]) -> bool:
        if DND_FILES is None:
            return False
        if not hasattr(widget, "drop_target_register") or not hasattr(widget, "dnd_bind"):
            return False

        try:
            widget.drop_target_register(DND_FILES)
            widget.dnd_bind("<<Drop>>", lambda event: callback(split_drop_paths(widget, event.data)))
        except tk.TclError:
            return False
        return True

    def on_main_file_drop(self, paths: list[str]) -> None:
        tic_paths = normalize_tic_paths(paths)
        if not tic_paths:
            self._set_status("Dropped item did not contain a readable TIC file.")
            return
        self.load_file(tic_paths[0])
        if len(tic_paths) > 1:
            self._set_status(f"Loaded {Path(tic_paths[0]).name}; ignored {len(tic_paths) - 1} additional files.")

    def get_integration_windows(self) -> list[IntegrationWindow]:
        windows: list[IntegrationWindow] = []
        for peak in self.peaks:
            left, right = peak.normalized_bounds()
            windows.append(IntegrationWindow(label=peak.label, left_x=left, right_x=right))
        return windows

    def open_batch_window(self) -> None:
        if self.batch_window is not None and self.batch_window.winfo_exists():
            self.batch_window.lift()
            self.batch_window.focus_force()
            self.batch_window.refresh_peak_windows()
            return

        self.batch_window = BatchProcessingWindow(self)

    def current_colors(self) -> dict[str, str]:
        theme = self.theme_var.get() if hasattr(self, "theme_var") else "Light"
        return THEME_COLORS.get(theme, THEME_COLORS["Light"])

    def apply_theme(self, theme_name: str | None = None, redraw: bool = True) -> None:
        if theme_name not in THEME_COLORS:
            theme_name = "Light"
        if hasattr(self, "theme_var"):
            self.theme_var.set(theme_name)

        colors = THEME_COLORS[theme_name]
        try:
            self.style.theme_use("clam")
        except tk.TclError:
            pass

        self.configure(background=colors["root_bg"])
        self.option_add("*TCombobox*Listbox.background", colors["tree_bg"])
        self.option_add("*TCombobox*Listbox.foreground", colors["text"])
        self.option_add("*TCombobox*Listbox.selectBackground", colors["select_bg"])
        self.option_add("*TCombobox*Listbox.selectForeground", colors["select_fg"])

        self.style.configure(".", background=colors["panel_bg"], foreground=colors["text"])
        self.style.configure("TFrame", background=colors["panel_bg"])
        self.style.configure("TLabelframe", background=colors["panel_bg"], foreground=colors["text"])
        self.style.configure("TLabelframe.Label", background=colors["panel_bg"], foreground=colors["text"])
        self.style.configure("TLabel", background=colors["panel_bg"], foreground=colors["text"])
        self.style.configure("TButton", background=colors["button_bg"], foreground=colors["text"])
        self.style.map(
            "TButton",
            background=[("active", colors["button_active"]), ("pressed", colors["button_active"])],
            foreground=[("disabled", colors["muted_text"])],
        )
        self.style.configure(
            "TEntry",
            fieldbackground=colors["entry_bg"],
            foreground=colors["text"],
            insertcolor=colors["text"],
        )
        self.style.configure(
            "TCombobox",
            fieldbackground=colors["entry_bg"],
            background=colors["button_bg"],
            foreground=colors["text"],
            arrowcolor=colors["text"],
        )
        self.style.map(
            "TCombobox",
            fieldbackground=[("readonly", colors["entry_bg"])],
            selectbackground=[("readonly", colors["entry_bg"])],
            selectforeground=[("readonly", colors["text"])],
        )
        self.style.configure(
            "Treeview",
            background=colors["tree_bg"],
            fieldbackground=colors["tree_bg"],
            foreground=colors["text"],
        )
        self.style.configure(
            "Treeview.Heading",
            background=colors["tree_heading_bg"],
            foreground=colors["text"],
        )
        self.style.map(
            "Treeview",
            background=[("selected", colors["select_bg"])],
            foreground=[("selected", colors["select_fg"])],
        )

        self._style_toolbar()
        if hasattr(self, "figure"):
            self.figure.patch.set_facecolor(colors["figure_bg"])
            self._style_axes()
        if hasattr(self, "canvas_widget"):
            self.canvas_widget.configure(background=colors["plot_bg"], highlightbackground=colors["plot_bg"])
        if (
            hasattr(self, "batch_window")
            and self.batch_window is not None
            and self.batch_window.winfo_exists()
        ):
            self.batch_window.apply_theme()

        if redraw and hasattr(self, "canvas"):
            self.refresh_plot(reset_limits=False)
        elif hasattr(self, "canvas"):
            self.canvas.draw_idle()

    def _style_toolbar(self) -> None:
        if not hasattr(self, "toolbar"):
            return

        colors = self.current_colors()
        self._configure_tk_widget(self.toolbar, background=colors["panel_bg"])
        for child in self.toolbar.winfo_children():
            self._configure_tk_widget(
                child,
                background=colors["panel_bg"],
                foreground=colors["text"],
                activebackground=colors["button_active"],
                activeforeground=colors["text"],
                highlightbackground=colors["panel_bg"],
                selectcolor=colors["button_active"],
            )

    @staticmethod
    def _configure_tk_widget(widget: tk.Widget, **options: str) -> None:
        for option, value in options.items():
            try:
                widget.configure(**{option: value})
            except tk.TclError:
                pass

    def on_theme_change(self, event: object | None = None) -> None:
        self.apply_theme(self.theme_var.get(), redraw=True)

    @staticmethod
    def _plain_number_formatter() -> ScalarFormatter:
        formatter = ScalarFormatter(useOffset=False, useMathText=False)
        formatter.set_scientific(False)
        return formatter

    def _style_axes(self) -> None:
        if not hasattr(self, "ax"):
            return

        colors = self.current_colors()
        self.ax.set_facecolor(colors["plot_bg"])
        self.ax.figure.patch.set_facecolor(colors["figure_bg"])
        self.ax.xaxis.set_major_formatter(self._plain_number_formatter())
        self.ax.yaxis.set_major_formatter(self._plain_number_formatter())
        self.ax.xaxis.get_offset_text().set_visible(False)
        self.ax.yaxis.get_offset_text().set_visible(False)
        self.ax.xaxis.label.set_color(colors["text"])
        self.ax.yaxis.label.set_color(colors["text"])
        self.ax.title.set_color(colors["text"])
        self.ax.tick_params(axis="both", colors=colors["text"])
        for spine in self.ax.spines.values():
            spine.set_color(colors["spine"])

    def browse_file(self) -> None:
        initial_path = Path(self.file_path.get())
        initial_dir = initial_path.parent if initial_path.parent.exists() else Path.cwd()
        path = filedialog.askopenfilename(
            title="Open TIC Text File",
            initialdir=str(initial_dir),
            filetypes=[
                ("TIC text files", "TIC*.txt"),
                ("Text files", "*.txt"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.file_path.set(path)
            self.load_file(path)

    def load_current_path(self) -> None:
        self.load_file(self.file_path.get())

    def load_file(self, path: str) -> None:
        path = path.strip()
        if not path:
            messagebox.showinfo("Load TIC File", "Choose a TIC_decon.txt file first.")
            self._set_status("No TIC file selected.")
            return

        try:
            x, y = load_tic_file(path)
        except Exception as exc:
            messagebox.showerror("Load TIC File", str(exc))
            self._set_status("File load failed.")
            return

        self.x = x
        self.y = y
        self.file_path.set(path)
        self.peaks.clear()
        self.selected_peak_index = None
        self.auto_scale_detection()
        self.refresh_table()
        self.refresh_plot(reset_limits=True)
        self._set_status(f"Loaded {len(x):,} points from {Path(path).name}.")

    def auto_scale_detection(self) -> None:
        if self.y is None or len(self.y) == 0:
            return

        y_range = float(np.nanmax(self.y) - np.nanmin(self.y))
        if np.isfinite(y_range) and y_range > 0:
            self.prominence_var.set(max(y_range * 0.05, 1.0))
        else:
            self.prominence_var.set(1.0)

        if self.x is not None and len(self.x) > 1:
            x_range = float(np.nanmax(self.x) - np.nanmin(self.x))
            self.distance_var.set(max(x_range * 0.01, 0.05))

    def detect_peaks(self) -> None:
        if self.x is None or self.y is None:
            messagebox.showinfo("Detect Peaks", "Load a TIC file first.")
            return

        try:
            min_prominence = float(self.prominence_var.get())
            min_distance_x = float(self.distance_var.get())
            boundary_fraction = float(self.boundary_fraction_var.get())
            smooth_points = int(self.smooth_points_var.get())
        except Exception:
            messagebox.showerror("Detect Peaks", "Detection settings must be numeric.")
            return

        self.peaks = detect_peak_integrations(
            self.x,
            self.y,
            min_prominence=min_prominence,
            min_distance_x=min_distance_x,
            boundary_fraction=boundary_fraction,
            smooth_points=smooth_points,
        )
        self.selected_peak_index = 0 if self.peaks else None
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)
        self._set_status(f"Detected {len(self.peaks)} peak integrations.")

    def add_peak_from_view(self) -> None:
        if self.x is None or self.y is None:
            messagebox.showinfo("Add View", "Load a TIC file first.")
            return

        left, right = self.ax.get_xlim()
        left, right = sorted((float(left), float(right)))
        data_min = float(np.nanmin(self.x))
        data_max = float(np.nanmax(self.x))
        left = max(left, data_min)
        right = min(right, data_max)
        if right <= left:
            messagebox.showerror("Add View", "The visible x-range does not overlap the TIC.")
            return

        mask = (self.x >= left) & (self.x <= right)
        if not np.any(mask):
            messagebox.showerror("Add View", "No TIC points are inside the visible x-range.")
            return

        visible_indices = np.flatnonzero(mask)
        apex_index = int(visible_indices[np.argmax(self.y[mask])])
        peak = PeakIntegration(
            apex_index=apex_index,
            left_x=left,
            right_x=right,
            label=f"Peak {len(self.peaks) + 1}",
        )
        self.peaks.append(peak)
        self.selected_peak_index = len(self.peaks) - 1
        self._renumber_peaks()
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)
        self._set_status("Added one integration from the visible x-range.")

    def delete_selected_peak(self) -> None:
        if self.selected_peak_index is None:
            return
        if not (0 <= self.selected_peak_index < len(self.peaks)):
            return

        del self.peaks[self.selected_peak_index]
        self._renumber_peaks()
        if self.peaks:
            self.selected_peak_index = min(self.selected_peak_index, len(self.peaks) - 1)
        else:
            self.selected_peak_index = None
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)

    def clear_peaks(self) -> None:
        self.peaks.clear()
        self.selected_peak_index = None
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._set_status("Cleared all integrations.")

    def export_csv(self) -> None:
        if self.x is None or self.y is None or not self.peaks:
            messagebox.showinfo("Export CSV", "There are no peak integrations to export.")
            return

        source_path = Path(self.file_path.get())
        default_name = f"{source_path.stem}_integrations.csv"
        output_path = filedialog.asksaveasfilename(
            title="Export Integrations",
            initialdir=str(source_path.parent) if source_path.parent.exists() else str(Path.cwd()),
            initialfile=default_name,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not output_path:
            return

        with open(output_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "label",
                    "apex_time",
                    "apex_intensity",
                    "left_time",
                    "right_time",
                    "points",
                    "sum_y_area",
                    "source_file",
                ]
            )
            for peak in self.peaks:
                left, right = peak.normalized_bounds()
                mask = peak.point_mask(self.x)
                writer.writerow(
                    [
                        peak.label,
                        peak.apex_x(self.x),
                        peak.height(self.y),
                        left,
                        right,
                        int(np.count_nonzero(mask)),
                        peak.area(self.x, self.y),
                        str(source_path),
                    ]
                )

        self._set_status(f"Exported integrations to {output_path}.")

    def apply_boundary_entries(self) -> None:
        peak = self.get_selected_peak()
        if peak is None:
            return

        try:
            left = float(self.left_boundary_var.get())
            right = float(self.right_boundary_var.get())
        except Exception:
            messagebox.showerror("Selected Peak", "Boundary values must be numeric.")
            return

        self.set_peak_bounds(peak, left, right)
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)

    def use_view_boundaries(self) -> None:
        peak = self.get_selected_peak()
        if peak is None:
            return

        left, right = self.ax.get_xlim()
        self.set_peak_bounds(peak, float(left), float(right))
        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)

    def set_peak_bounds(self, peak: PeakIntegration, left: float, right: float) -> None:
        if self.x is None or self.y is None:
            return

        data_min = float(np.nanmin(self.x))
        data_max = float(np.nanmax(self.x))
        left = float(np.clip(left, data_min, data_max))
        right = float(np.clip(right, data_min, data_max))
        if left == right:
            right = min(data_max, left + np.nanmedian(np.diff(self.x)))
        peak.left_x = min(left, right)
        peak.right_x = max(left, right)

        mask = peak.point_mask(self.x)
        if np.any(mask):
            indices = np.flatnonzero(mask)
            peak.apex_index = int(indices[np.argmax(self.y[mask])])

        self.update_selected_fields()

    def get_selected_peak(self) -> PeakIntegration | None:
        if self.selected_peak_index is None:
            return None
        if 0 <= self.selected_peak_index < len(self.peaks):
            return self.peaks[self.selected_peak_index]
        return None

    def on_table_select(self, event: object | None = None) -> None:
        selection = self.peak_table.selection()
        if not selection:
            return
        item = selection[0]
        try:
            self.selected_peak_index = int(item.replace("peak_", ""))
        except ValueError:
            return
        self.update_selected_fields()
        self.refresh_plot(reset_limits=False)

    def on_plot_press(self, event: object) -> None:
        if event.inaxes != self.ax or event.xdata is None:
            return
        if self._toolbar_is_active():
            return

        self.last_cursor_x = float(event.xdata)
        peak = self.get_selected_peak()
        if peak is not None and getattr(event, "button", None) == 1:
            boundary = self._nearest_boundary(float(event.xdata), peak)
            if boundary is not None:
                self.dragging_boundary = boundary
                return

        if getattr(event, "button", None) == 1:
            self.select_peak_at_x(float(event.xdata))

    def on_plot_motion(self, event: object) -> None:
        if event.inaxes == self.ax and event.xdata is not None:
            self.last_cursor_x = float(event.xdata)

        if self.dragging_boundary is None:
            return
        if event.inaxes != self.ax or event.xdata is None:
            return

        peak = self.get_selected_peak()
        if peak is None:
            return

        x_value = float(event.xdata)
        if self.dragging_boundary == "left":
            self.set_peak_bounds(peak, x_value, peak.right_x)
        elif self.dragging_boundary == "right":
            self.set_peak_bounds(peak, peak.left_x, x_value)

        self.refresh_table()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)

    def on_plot_release(self, event: object) -> None:
        if self.dragging_boundary is not None:
            self.dragging_boundary = None
            self.refresh_table()
            self.refresh_plot(reset_limits=False)
            self._select_table_row(self.selected_peak_index)

    def select_peak_at_x(self, x_value: float) -> None:
        if not self.peaks or self.x is None:
            return

        containing = [
            index
            for index, peak in enumerate(self.peaks)
            if peak.normalized_bounds()[0] <= x_value <= peak.normalized_bounds()[1]
        ]
        if containing:
            self.selected_peak_index = containing[0]
        else:
            self.selected_peak_index = min(
                range(len(self.peaks)),
                key=lambda index: abs(self.peaks[index].apex_x(self.x) - x_value),
            )
        self.update_selected_fields()
        self.refresh_plot(reset_limits=False)
        self._select_table_row(self.selected_peak_index)

    def _nearest_boundary(self, x_value: float, peak: PeakIntegration) -> str | None:
        xlim = self.ax.get_xlim()
        tolerance = abs(float(xlim[1] - xlim[0])) * 0.015
        if tolerance <= 0:
            return None

        left_distance = abs(x_value - peak.left_x)
        right_distance = abs(x_value - peak.right_x)
        if min(left_distance, right_distance) > tolerance:
            return None
        return "left" if left_distance <= right_distance else "right"

    def _toolbar_is_active(self) -> bool:
        mode = getattr(self.toolbar, "mode", "")
        return bool(mode)

    def refresh_plot(self, reset_limits: bool) -> None:
        old_xlim = self.ax.get_xlim()
        old_ylim = self.ax.get_ylim()

        self.ax.clear()
        self.ax.set_xlabel("Time (min)")
        self.ax.set_ylabel("Intensity (# of ions)")
        self._style_axes()
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["right"].set_visible(False)

        if self.x is not None and self.y is not None:
            colors = self.current_colors()
            self.ax.plot(self.x, self.y, color=colors["tic"], linewidth=1.25, label="TIC")
            self._draw_peak_overlays()
            legend = self.ax.legend(loc="upper right", frameon=False)
            for text in legend.get_texts():
                text.set_color(colors["text"])

        if self.x is not None and self.y is not None:
            if reset_limits:
                self.ax.set_xlim(float(np.nanmin(self.x)), float(np.nanmax(self.x)))
                ymin = float(np.nanmin(self.y))
                ymax = float(np.nanmax(self.y))
                yrange = ymax - ymin
                if yrange <= 0:
                    yrange = max(abs(ymax), 1.0)
                self.ax.set_ylim(ymin - yrange * 0.05, ymax + yrange * 0.10)
            else:
                self.ax.set_xlim(old_xlim)
                self.ax.set_ylim(old_ylim)

        self.figure.tight_layout()
        self.canvas.draw_idle()
        self.update_selected_fields()

    def _draw_peak_overlays(self) -> None:
        if self.x is None or self.y is None:
            return

        for index, peak in enumerate(self.peaks):
            left, right = peak.normalized_bounds()
            selected = index == self.selected_peak_index
            colors = self.current_colors()
            fill_color = colors["fill_selected"] if selected else colors["fill_unselected"]
            line_color = colors["line_selected"] if selected else colors["line_unselected"]
            alpha = 0.22 if selected else 0.10

            self.ax.axvspan(left, right, color=fill_color, alpha=alpha, linewidth=0)
            self.ax.axvline(left, color=line_color, linewidth=1.5 if selected else 0.8, linestyle="--")
            self.ax.axvline(right, color=line_color, linewidth=1.5 if selected else 0.8, linestyle="--")
            self.ax.plot(
                peak.apex_x(self.x),
                peak.height(self.y),
                marker="o",
                markersize=5 if selected else 4,
                color=line_color,
            )

    def refresh_table(self) -> None:
        if self.x is None or self.y is None:
            self.peak_table.delete(*self.peak_table.get_children())
            return

        previous_selection = self.selected_peak_index
        self.peak_table.delete(*self.peak_table.get_children())

        for index, peak in enumerate(self.peaks):
            left, right = peak.normalized_bounds()
            self.peak_table.insert(
                "",
                "end",
                iid=f"peak_{index}",
                text=peak.label,
                values=(
                    f"{peak.apex_x(self.x):.4f}",
                    f"{left:.4f}",
                    f"{right:.4f}",
                    f"{peak.area(self.x, self.y):.4g}",
                    f"{peak.height(self.y):.4g}",
                ),
            )

        self._select_table_row(previous_selection)
        self.update_selected_fields()

    def update_selected_fields(self) -> None:
        peak = self.get_selected_peak()
        if peak is None or self.x is None or self.y is None:
            self.left_boundary_var.set("")
            self.right_boundary_var.set("")
            self.selected_area_var.set("")
            return

        left, right = peak.normalized_bounds()
        self.left_boundary_var.set(f"{left:.6g}")
        self.right_boundary_var.set(f"{right:.6g}")
        self.selected_area_var.set(f"{peak.area(self.x, self.y):.6g}")

    def _select_table_row(self, index: int | None) -> None:
        if index is None:
            return
        item = f"peak_{index}"
        if self.peak_table.exists(item):
            self.peak_table.selection_set(item)
            self.peak_table.focus(item)
            self.peak_table.see(item)

    def _renumber_peaks(self) -> None:
        for index, peak in enumerate(self.peaks, start=1):
            peak.label = f"Peak {index}"

    def _set_status(self, message: str) -> None:
        self.status_var.set(message)


class BatchProcessingWindow(tk.Toplevel):
    """Batch integration window using peak windows from the main TIC viewer."""

    def __init__(self, app: TICTIntegratorApp) -> None:
        super().__init__(app)
        self.app = app
        self.title("Batch TIC Integration")
        self.geometry("1250x760")
        self.minsize(950, 620)

        self.dataset_paths: list[str] = []
        self.windows: list[IntegrationWindow] = []
        self.result_rows: list[dict[str, object]] = []
        self.status_var = tk.StringVar(value="")

        self.protocol("WM_DELETE_WINDOW", self.on_close)
        self._build_layout()
        self.refresh_peak_windows()
        self._enable_drag_and_drop()
        self.apply_theme()

    def _build_layout(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        top = ttk.Frame(self, padding=8)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(0, weight=2)
        top.columnconfigure(1, weight=1)

        self._build_dataset_panel(top)
        self._build_windows_panel(top)
        self._build_results_panel()

        status = ttk.Label(self, textvariable=self.status_var, anchor="w", padding=(8, 4))
        status.grid(row=2, column=0, sticky="ew")

    def _build_dataset_panel(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Datasets", padding=8)
        frame.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = ("file", "folder")
        self.dataset_table = ttk.Treeview(frame, columns=columns, show="headings", selectmode="extended", height=8)
        self.dataset_table.heading("file", text="File")
        self.dataset_table.heading("folder", text="Folder")
        self.dataset_table.column("file", width=320, stretch=False)
        self.dataset_table.column("folder", width=520, stretch=True)
        self.dataset_table.grid(row=0, column=0, sticky="nsew")

        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=self.dataset_table.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        self.dataset_table.configure(yscrollcommand=y_scroll.set)

        buttons = ttk.Frame(frame)
        buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        buttons.columnconfigure((0, 1, 2, 3), weight=1)
        ttk.Button(buttons, text="Add Files", command=self.add_files).grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(buttons, text="Add Current", command=self.add_current_dataset).grid(
            row=0,
            column=1,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Remove", command=self.remove_selected_datasets).grid(
            row=0,
            column=2,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Clear", command=self.clear_datasets).grid(row=0, column=3, sticky="ew")

    def _build_windows_panel(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Peak Windows", padding=8)
        frame.grid(row=0, column=1, sticky="nsew")
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = ("label", "left", "right")
        self.window_table = ttk.Treeview(frame, columns=columns, show="headings", selectmode="browse", height=8)
        self.window_table.heading("label", text="Peak")
        self.window_table.heading("left", text="Left")
        self.window_table.heading("right", text="Right")
        self.window_table.column("label", width=110, stretch=True)
        self.window_table.column("left", width=90, anchor="e", stretch=False)
        self.window_table.column("right", width=90, anchor="e", stretch=False)
        self.window_table.grid(row=0, column=0, sticky="nsew")

        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=self.window_table.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        self.window_table.configure(yscrollcommand=y_scroll.set)

        buttons = ttk.Frame(frame)
        buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        buttons.columnconfigure((0, 1), weight=1)
        ttk.Button(buttons, text="Refresh Windows", command=self.refresh_peak_windows).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Process Batch", command=self.process_batch).grid(row=0, column=1, sticky="ew")

    def _build_results_panel(self) -> None:
        frame = ttk.LabelFrame(self, text="Results", padding=8)
        frame.grid(row=1, column=0, sticky="nsew", padx=8)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = ("dataset", "peak", "left", "right", "points", "area", "apex", "height", "status")
        self.results_table = ttk.Treeview(frame, columns=columns, show="headings", selectmode="extended")
        headings = {
            "dataset": "Dataset",
            "peak": "Peak",
            "left": "Left",
            "right": "Right",
            "points": "Points",
            "area": "Area",
            "apex": "Apex",
            "height": "Height",
            "status": "Status",
        }
        widths = {
            "dataset": 300,
            "peak": 100,
            "left": 90,
            "right": 90,
            "points": 70,
            "area": 110,
            "apex": 90,
            "height": 110,
            "status": 360,
        }
        for column in columns:
            self.results_table.heading(column, text=headings[column])
            anchor = "e" if column in {"left", "right", "points", "area", "apex", "height"} else "w"
            self.results_table.column(column, width=widths[column], anchor=anchor, stretch=column in {"dataset", "status"})

        self.results_table.grid(row=0, column=0, sticky="nsew")
        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=self.results_table.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll = ttk.Scrollbar(frame, orient="horizontal", command=self.results_table.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.results_table.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        buttons = ttk.Frame(frame)
        buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        buttons.columnconfigure((0, 1), weight=1)
        ttk.Button(buttons, text="Export Results CSV", command=self.export_results_csv).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(buttons, text="Clear Results", command=self.clear_results).grid(row=0, column=1, sticky="ew")

    def _enable_drag_and_drop(self) -> None:
        enabled = False
        enabled |= self.app.register_file_drop_target(self, self.add_dataset_paths)
        enabled |= self.app.register_file_drop_target(self.dataset_table, self.add_dataset_paths)
        if enabled:
            self.set_status("Drop TIC files onto the dataset list or use Add Files.")

    def apply_theme(self) -> None:
        colors = self.app.current_colors()
        self.configure(background=colors["root_bg"])

    def add_files(self) -> None:
        initial_dir = self._initial_directory()
        paths = filedialog.askopenfilenames(
            title="Add TIC Files",
            initialdir=str(initial_dir),
            filetypes=[
                ("TIC text files", "TIC*.txt"),
                ("Text files", "*.txt"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        self.add_dataset_paths(list(paths))

    def add_current_dataset(self) -> None:
        current_path = self.app.file_path.get().strip()
        if not current_path:
            messagebox.showinfo("Add Current", "No TIC file is loaded in the main window.")
            return
        self.add_dataset_paths([current_path])

    def add_dataset_paths(self, paths: list[str]) -> None:
        tic_paths = normalize_tic_paths(paths)
        existing = {os.path.normcase(os.path.abspath(path)) for path in self.dataset_paths}
        added = 0
        for path in tic_paths:
            key = os.path.normcase(os.path.abspath(path))
            if key in existing:
                continue
            self.dataset_paths.append(path)
            existing.add(key)
            added += 1
        self.refresh_dataset_table()
        if added:
            self.set_status(f"Added {added} dataset(s).")
        elif paths:
            self.set_status("No new readable TIC files were added.")

    def remove_selected_datasets(self) -> None:
        selected = self.dataset_table.selection()
        if not selected:
            return
        remove_indices = sorted((int(item.replace("dataset_", "")) for item in selected), reverse=True)
        for index in remove_indices:
            if 0 <= index < len(self.dataset_paths):
                del self.dataset_paths[index]
        self.refresh_dataset_table()
        self.set_status(f"Removed {len(remove_indices)} dataset(s).")

    def clear_datasets(self) -> None:
        self.dataset_paths.clear()
        self.refresh_dataset_table()
        self.set_status("Cleared datasets.")

    def refresh_dataset_table(self) -> None:
        self.dataset_table.delete(*self.dataset_table.get_children())
        for index, path in enumerate(self.dataset_paths):
            file_path = Path(path)
            self.dataset_table.insert(
                "",
                "end",
                iid=f"dataset_{index}",
                values=(file_path.name, str(file_path.parent)),
            )

    def refresh_peak_windows(self) -> None:
        self.windows = self.app.get_integration_windows()
        self.window_table.delete(*self.window_table.get_children())
        for index, window in enumerate(self.windows):
            left, right = window.normalized_bounds()
            self.window_table.insert(
                "",
                "end",
                iid=f"window_{index}",
                values=(window.label, f"{left:.6g}", f"{right:.6g}"),
            )
        self.set_status(f"Using {len(self.windows)} peak window(s) from the main TIC viewer.")

    def process_batch(self) -> None:
        self.refresh_peak_windows()
        if not self.windows:
            messagebox.showinfo("Process Batch", "Define at least one integration peak in the main window first.")
            return
        if not self.dataset_paths:
            messagebox.showinfo("Process Batch", "Add at least one dataset first.")
            return

        rows: list[dict[str, object]] = []
        for path in self.dataset_paths:
            dataset_name = Path(path).name
            try:
                x, y = load_tic_file(path)
            except Exception as exc:
                rows.append(
                    {
                        "dataset": dataset_name,
                        "source_file": path,
                        "peak": "",
                        "left": "",
                        "right": "",
                        "points": "",
                        "area": "",
                        "apex": "",
                        "height": "",
                        "status": str(exc),
                    }
                )
                continue

            for window in self.windows:
                left, right = window.normalized_bounds()
                stats = window.integrate(x, y)
                rows.append(
                    {
                        "dataset": dataset_name,
                        "source_file": path,
                        "peak": window.label,
                        "left": left,
                        "right": right,
                        "points": stats["points"],
                        "area": stats["area"],
                        "apex": stats["apex_x"],
                        "height": stats["height"],
                        "status": "OK",
                    }
                )

        self.result_rows = rows
        self.refresh_results_table()
        self.set_status(f"Processed {len(self.dataset_paths)} dataset(s) with {len(self.windows)} peak window(s).")

    def refresh_results_table(self) -> None:
        self.results_table.delete(*self.results_table.get_children())
        for index, row in enumerate(self.result_rows):
            self.results_table.insert(
                "",
                "end",
                iid=f"result_{index}",
                values=(
                    row["dataset"],
                    row["peak"],
                    self._format_number(row["left"]),
                    self._format_number(row["right"]),
                    row["points"],
                    self._format_number(row["area"]),
                    self._format_number(row["apex"]),
                    self._format_number(row["height"]),
                    row["status"],
                ),
            )

    def clear_results(self) -> None:
        self.result_rows.clear()
        self.refresh_results_table()
        self.set_status("Cleared results.")

    def export_results_csv(self) -> None:
        if not self.result_rows:
            messagebox.showinfo("Export Results CSV", "There are no batch results to export.")
            return

        output_path = filedialog.asksaveasfilename(
            title="Export Batch Results",
            initialdir=str(self._initial_directory()),
            initialfile="tic_batch_integrations.csv",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not output_path:
            return

        columns = [
            "dataset",
            "source_file",
            "peak",
            "left",
            "right",
            "points",
            "area",
            "apex",
            "height",
            "status",
        ]
        with open(output_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in self.result_rows:
                writer.writerow({column: row.get(column, "") for column in columns})

        self.set_status(f"Exported batch results to {output_path}.")

    def _initial_directory(self) -> Path:
        if self.dataset_paths:
            first = Path(self.dataset_paths[0])
            if first.parent.exists():
                return first.parent
        current_path = self.app.file_path.get().strip()
        if current_path:
            parent = Path(current_path).parent
            if parent.exists():
                return parent
        return Path.cwd()

    @staticmethod
    def _format_number(value: object) -> str:
        if value is None or value == "":
            return ""
        try:
            return f"{float(value):.6g}"
        except (TypeError, ValueError):
            return str(value)

    def set_status(self, message: str) -> None:
        self.status_var.set(message)

    def on_close(self) -> None:
        self.app.batch_window = None
        self.destroy()


def main(argv: Iterable[str] | None = None) -> None:
    argv = list(argv or [])
    app = TICTIntegratorApp()
    if argv:
        app.file_path.set(argv[0])
        app.load_file(argv[0])
    app.mainloop()


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])

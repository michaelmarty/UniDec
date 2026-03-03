from __future__ import annotations
import os
import threading
from dataclasses import dataclass
from typing import List, Sequence, Tuple, Optional, Callable

import numpy as np
import tkinter as tk
from tkinter import filedialog
from tkinter.colorchooser import askcolor

import ttkbootstrap as tb
from ttkbootstrap.constants import *
from ttkbootstrap.dialogs import Messagebox
try:
    # Newer ttkbootstrap
    from ttkbootstrap.widgets.scrolled import ScrolledFrame
except Exception:  # fallback for older versions
    from ttkbootstrap.scrolled import ScrolledFrame

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Drag-and-drop root (optional)
try:
    from tkinterdnd2 import TkinterDnD, DND_FILES
    DND_AVAILABLE = True
except Exception:
    DND_AVAILABLE = False

# Plotting utilities
from chromatogram_plot_utils import plot_chrom_stack, read_chrom_file
from cdms_mass_plot_utils import (
    plot_cdms_mass_stack,
    plot_single_cdms_mass,
    _load_mass,
)
from mwd_chromatogram_utils import plot_mwd_chromatograms
from mz_spectrum_plot_utils import (
    plot_single_mz_spectrum,
    plot_mz_spectrum_stack,
    read_mz_spectrum_file,
)
from thermo_raw_utils import (
    plot_raw_chrom_and_spectrum, _convert_raw_to_mzml
)

# --------------------------- Defaults & Helpers ---------------------------

DEFAULT_CHROM_COLORS = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CDMS_COLORS  = ["#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CHROM_STANDARDS = "13.81, 16.82, 18.09, 19.05, 20.79"
DEFAULT_MWD_COLORS   = ["#0F0C08", "#3CAEA3"]
DEFAULT_MZ_COLORS    = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442"]
THEMES = ["flatly", "cosmo", "lumen", "minty", "sandstone", "superhero", "cyborg", "darkly", "morph"]

DISPLAY_DPI = 110  # preview DPI


def parse_float_list(text: str) -> List[float]:
    text = text.strip()
    if not text:
        return []
    out: List[float] = []
    for part in text.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(float(part))
        except ValueError:
            pass
    return out


def parse_str_list(text: str) -> List[str]:
    text = text.strip()
    if not text:
        return []
    return [s.strip() for s in text.replace(";", ",").split(",") if s.strip()]


def filenames_to_labels(paths: Sequence[str]) -> List[str]:
    return [os.path.splitext(os.path.basename(p))[0] for p in paths]


def clamp_range(val_min: float, val_max: float, lo: float, hi: float) -> Tuple[float, float]:
    vmin = max(lo, min(hi, val_min))
    vmax = max(lo, min(hi, val_max))
    if vmax < vmin:
        vmin, vmax = vmax, vmin
    return vmin, vmax


def normalize_y_values(
    y: np.ndarray,
    normalize: bool,
    mode: str
) -> np.ndarray:
    if not normalize:
        return y
    m = float(np.nanmax(y))
    if m <= 0:
        return y
    if mode == "percent":
        return (y / m) * 100.0
    return y / m


@dataclass
class CommonConfig:
    save: Optional[str]
    font_name: str
    font_size: int
    dpi: int
    figsize: Tuple[float, float]
    transparent: bool
    show: bool


# ------------------------------ Main App ------------------------------

class PlotApp:
    def __init__(self, root: tk.Tk, themename: str = "flatly"):
        self.root = root
        try:
            self.root.tk.call("package", "require", "msgcat")
        except tk.TclError:
            pass

        try:
            self.style = tb.Style(master=self.root, theme=themename)
        except TypeError:
            self.style = tb.Style(themename)
        try:
            self.style.theme_use(themename)
        except Exception:
            pass

        self.win_w = 1200
        self.win_h = 850
        self.root.geometry(f"{self.win_w}x{self.win_h}")
        try:
            self.root.configure(bg=self.style.colors.bg)
        except Exception:
            pass

        self.header = tb.Frame(self.root, bootstyle="secondary")
        self.header.pack(fill=X)
        self.centerer = tb.Frame(self.root, bootstyle="secondary")
        self.centerer.pack(fill=BOTH, expand=YES)
        self.footer = tb.Frame(self.root, bootstyle="secondary")
        self.footer.pack(fill=X, side=BOTTOM)

        self.centerer.grid_columnconfigure(0, weight=0)
        self.centerer.grid_columnconfigure(1, weight=1)
        self.centerer.grid_columnconfigure(2, weight=0)
        self.centerer.grid_rowconfigure(0, weight=1)

        self.body = tb.Frame(self.centerer, bootstyle="secondary")
        self.body.grid(row=0, column=1, sticky="nsew")

        self._build_header(self.header)
        self._build_body(self.body)
        self._build_footer(self.footer)

        # Bind automatic clamping for global Y-limits (after tabs are built
        # so that chrom_norm_mode exists)
        self._bind_global_y_clamping()

        self._center_on_screen()

        self.files: List[str] = []
        self._set_status("Ready. Drag files into the list or use the pickers, then click Plot.")

        for seq in ("<F5>", "<Control-Return>", "<Control-KP_Enter>"):
            self.root.bind(seq, lambda e: self._plot())

    # --------------------------- Header ---------------------------

    def _build_header(self, parent: tk.Widget):
        header = tb.Frame(parent, padding=10, bootstyle="secondary")
        header.pack(fill=X)
        tb.Label(
            header,
            text="Plotting GUI: Chromatograms, CD-MS Mass, MWD, m/z Spectra, and Thermo RAW",
            font=("Arial", 16, "bold"),
        ).pack(side=LEFT)

        theme_frame = tb.Frame(header, bootstyle="secondary")
        theme_frame.pack(side=RIGHT)
        tb.Label(theme_frame, text="Theme:", font=("Arial", 10, "bold")).pack(side=LEFT, padx=(0, 6))
        self.theme_var = tk.StringVar(value=self.style.theme_use())
        combo = tb.Combobox(
            theme_frame,
            values=THEMES,
            textvariable=self.theme_var,
            state="readonly",
            width=10,
        )
        combo.pack(side=LEFT)
        combo.bind("<<ComboboxSelected>>", self._on_theme_change)

        bar = tb.Frame(parent, bootstyle="secondary")
        bar.pack(fill=X)
        self.status_label = tb.Label(bar, text="Status: Ready", anchor=W)
        self.status_label.pack(side=LEFT, padx=8, pady=4)

        self.plot_btn = tb.Button(bar, text="Plot", bootstyle="success", command=self._plot)
        self._apply_plot_button_sizing()
        self.plot_btn.pack(side=RIGHT, padx=8, pady=4)

    def _on_theme_change(self, _=None):
        self.style.theme_use(self.theme_var.get())
        try:
            self.root.configure(bg=self.style.colors.bg)
        except Exception:
            pass
        self._apply_plot_button_sizing()

    def _apply_plot_button_sizing(self):
        if not hasattr(self, "plot_btn") or self.plot_btn is None:
            return
        try:
            scale = float(self.root.tk.call("tk", "scaling"))
        except Exception:
            scale = 1.0

        self.plot_btn.configure(width=12)
        base_pad_x, base_pad_y = 14, 10
        pad_x, pad_y = int(base_pad_x * scale), int(base_pad_y * scale)
        try:
            self.plot_btn.configure(padding=(pad_x, pad_y))
        except Exception:
            pass

    # --------------------------- Scrollable Left ---------------------------

    def _build_scrollable_left(self, paned: tb.Panedwindow):
        container = tb.Frame(paned, bootstyle="secondary")

        sf = ScrolledFrame(
            master=container,
            autohide=False,
            bootstyle="secondary"
        )
        sf.pack(fill=BOTH, expand=YES)

        if hasattr(sf, "frame"):
            inner = sf.frame
        elif hasattr(sf, "scrollable_frame"):
            inner = sf.scrollable_frame
        elif hasattr(sf, "innerframe"):
            inner = sf.innerframe
        else:
            inner = sf

        canvas = getattr(sf, "canvas", None)
        vbar = getattr(sf, "vscrollbar", None)

        def update_scrollbar_visibility(_=None):
            if not canvas or not vbar:
                return
            canvas.update_idletasks()
            bbox = canvas.bbox("all")
            if not bbox:
                return
            content_height = bbox[3] - bbox[1]
            view_height = canvas.winfo_height()
            if view_height > 1 and content_height <= view_height:
                canvas.yview_moveto(0.0)
                if vbar.winfo_manager():
                    vbar.pack_forget()
            else:
                if not vbar.winfo_manager():
                    vbar.pack(side=RIGHT, fill=Y)

        inner.bind("<Configure>", update_scrollbar_visibility)
        sf.bind("<Configure>", update_scrollbar_visibility)

        return container, inner

    # --------------------------- Body ---------------------------

    def _build_body(self, parent: tk.Widget):
        paned = tb.Panedwindow(parent, orient=HORIZONTAL, bootstyle="secondary")
        paned.pack(fill=BOTH, expand=YES, padx=10, pady=0)

        left_container, left_inner = self._build_scrollable_left(paned)
        paned.add(left_container, weight=1)

        left = tb.Frame(left_inner, padding=(0, 0, 16, 0), bootstyle="secondary")
        left.pack(fill=BOTH, expand=YES)

        right = tb.Frame(paned, padding=10, bootstyle="secondary")
        paned.add(right, weight=3)
        self.plot_area = right

        # Notebook (tabs)
        self.nb = tb.Notebook(left, bootstyle="primary")
        self.nb.pack(fill=BOTH, expand=YES)

        self._build_common_options(left)
        self._build_chrom_tab(self.nb)
        self._build_cdms_tab(self.nb)
        self._build_mwd_tab(self.nb)
        self._build_mz_tab(self.nb)
        self._build_raw_tab(self.nb)

        # Rebuild file panel after tabs
        self._build_file_panel(left)

        # Bind tab change to update per-tab default ranges
        self.nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)
        return  # ensure no extra calls below

    # --------------------------- Footer ---------------------------

    def _build_footer(self, parent: tk.Widget):
        frame = tb.Frame(parent, bootstyle="secondary")
        frame.pack(fill=X, pady=6)
        self.progress = tb.Progressbar(frame, mode="indeterminate", bootstyle="info")
        self.progress.pack(side=LEFT, fill=X, expand=YES, padx=6)

    def _center_on_screen(self):
        self.root.update_idletasks()
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        x = (sw - self.win_w) // 2
        y = (sh - self.win_h) // 2
        self.root.geometry(f"{self.win_w}x{self.win_h}+{x}+{y}")
        try:
            self.root.state("zoomed")
        except Exception:
            pass

    # --------------------------- Common Options ---------------------------

    def _build_common_options(self, parent: tk.Widget):
        lf = tb.Labelframe(parent, text="Common Options", padding=10, bootstyle="secondary")
        lf.pack(fill=X, pady=6)

        # Save path
        sp = tb.Frame(lf, bootstyle="secondary")
        sp.pack(fill=X, pady=4)
        tb.Label(sp, text="Save as:", width=12).pack(side=LEFT)
        self.save_var = tk.StringVar()
        tb.Entry(sp, textvariable=self.save_var).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(sp, text="Browse", bootstyle="secondary", command=self._browse_save).pack(side=LEFT)

        grid = tb.Frame(lf, bootstyle="secondary")
        grid.pack(fill=X, pady=4)

        self.font_name = tk.StringVar(value="Arial")
        self.font_size = tk.IntVar(value=10)
        self.dpi = tk.IntVar(value=300)
        self.fig_w = tk.DoubleVar(value=10.0)
        self.fig_h = tk.DoubleVar(value=6.0)

        def label_entry(row, col, text, var, width):
            tb.Label(grid, text=text).grid(row=row, column=col, sticky=W)
            tb.Entry(grid, textvariable=var, width=width).grid(row=row, column=col + 1, sticky=W, padx=4)

        label_entry(0, 0, "Font:", self.font_name, 14)
        label_entry(0, 2, "Size:", self.font_size, 6)
        label_entry(0, 4, "DPI:", self.dpi, 8)
        tb.Label(grid, text="Figsize (W,H):").grid(row=0, column=6, sticky=W)
        tb.Entry(grid, textvariable=self.fig_w, width=6).grid(row=0, column=7, sticky=W, padx=2)
        tb.Entry(grid, textvariable=self.fig_h, width=6).grid(row=0, column=8, sticky=W, padx=2)

        # Global axis / style controls
        axes_frame = tb.Frame(lf, bootstyle="secondary")
        axes_frame.pack(fill=X, pady=4)

        self.global_xmin = tk.DoubleVar(value=0.0)
        self.global_xmax = tk.DoubleVar(value=20.0)
        self.global_ymin = tk.DoubleVar(value=0.0)
        self.global_ymax = tk.DoubleVar(value=100.0)
        self.global_subplot_spacing = tk.DoubleVar(value=0.6)

        tb.Label(axes_frame, text="X-lim min:").pack(side=LEFT)
        tb.Entry(axes_frame, textvariable=self.global_xmin, width=8).pack(side=LEFT, padx=4)
        tb.Label(axes_frame, text="X-lim max:").pack(side=LEFT)
        tb.Entry(axes_frame, textvariable=self.global_xmax, width=8).pack(side=LEFT, padx=4)

        tb.Label(axes_frame, text="Y-lim min:").pack(side=LEFT, padx=(12, 0))
        self.global_ymin_entry = tb.Entry(axes_frame, textvariable=self.global_ymin, width=8)
        self.global_ymin_entry.pack(side=LEFT, padx=4)
        tb.Label(axes_frame, text="Y-lim max:").pack(side=LEFT)
        self.global_ymax_entry = tb.Entry(axes_frame, textvariable=self.global_ymax, width=8)
        self.global_ymax_entry.pack(side=LEFT, padx=4)

        tb.Label(axes_frame, text="Subplot spacing:").pack(side=LEFT, padx=(12, 0))
        tb.Entry(axes_frame, textvariable=self.global_subplot_spacing, width=8).pack(side=LEFT, padx=4)

        style_frame = tb.Frame(lf, bootstyle="secondary")
        style_frame.pack(fill=X, pady=6)
        self.x_label = tk.StringVar()
        self.y_label = tk.StringVar()
        self.line_width = tk.DoubleVar(value=2.0)

        tb.Label(style_frame, text="X label:").pack(side=LEFT)
        tb.Entry(style_frame, textvariable=self.x_label, width=22).pack(side=LEFT, padx=6)
        tb.Label(style_frame, text="Y label:").pack(side=LEFT)
        tb.Entry(style_frame, textvariable=self.y_label, width=22).pack(side=LEFT, padx=6)
        tb.Label(style_frame, text="Line width:").pack(side=LEFT, padx=(10, 0))
        tb.Entry(style_frame, textvariable=self.line_width, width=6).pack(side=LEFT, padx=6)

        toggles = tb.Frame(lf, bootstyle="secondary")
        toggles.pack(fill=X, pady=4)
        self.show_plot = tk.BooleanVar(value=True)
        self.transparent_export = tk.BooleanVar(value=False)
        tb.Checkbutton(
            toggles,
            text="Show Figure in external window (in addition to preview)",
            variable=self.show_plot,
        ).pack(anchor=W, pady=2)
        tb.Checkbutton(
            toggles,
            text="Transparent export background",
            variable=self.transparent_export,
        ).pack(anchor=W, pady=2)

    def _bind_global_y_clamping(self):
        """
        Bind focus-out clamping for global Y-lim entries, using the chrom
        normalization mode as the reference (0–100% vs 0–1). This ensures
        that out-of-range values are corrected automatically.
        """
        if not hasattr(self, "global_ymin_entry") or not hasattr(self, "global_ymax_entry"):
            return
        if not hasattr(self, "chrom_norm_mode"):
            return

        def _on_ymin_focus_out(event=None):
            self._clamp_y_entry(
                entry_var=self.global_ymin,
                mode_var=self.chrom_norm_mode,
                percent_range=(0.0, 100.0),
                unit_range=(0.0, 1.0),
            )

        def _on_ymax_focus_out(event=None):
            self._clamp_y_entry(
                entry_var=self.global_ymax,
                mode_var=self.chrom_norm_mode,
                percent_range=(0.0, 100.0),
                unit_range=(0.0, 1.0),
            )

        self.global_ymin_entry.bind("<FocusOut>", _on_ymin_focus_out)
        self.global_ymax_entry.bind("<FocusOut>", _on_ymax_focus_out)

    # --------------------------- Tab Builders ---------------------------

    def _build_chrom_tab(self, nb: tb.Notebook):
        self.chrom_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.chrom_tab, text="Chromatograms")

        self.chrom_mode = tk.StringVar(value="single")
        self._build_tab_mode_row(
            parent=self.chrom_tab,
            label="Plot Mode:",
            variable=self.chrom_mode,
            options=[("Single Chromatogram", "single"), ("Stacked Chromatograms", "stacked")],
        )

        self._build_labels_colors_section(
            parent=self.chrom_tab,
            title="Chromatogram Labels & Colors",
            text_widget_attr="chrom_labels",
            color_var_attr="chrom_colors",
            default_colors=", ".join(DEFAULT_CHROM_COLORS),
            palette_title="Chromatogram Colors",
        )

        lf_axes = tb.Labelframe(self.chrom_tab, text="Chromatogram Normalization & Standards", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)

        self.chrom_normalize = tk.BooleanVar(value=True)
        self.chrom_norm_mode = tk.StringVar(value="percent")
        self.chrom_show_standards = tk.BooleanVar(value=True)
        self.chrom_standards = tk.StringVar(value=DEFAULT_CHROM_STANDARDS)

        tb.Checkbutton(lf_axes, text="Normalize traces to max", variable=self.chrom_normalize).pack(anchor=W, pady=4)

        self._build_norm_mode_toggle(
            parent=lf_axes,
            label="Normalization scale (for Chromatograms):",
            mode_var=self.chrom_norm_mode,
            ymin_var=self.global_ymin,
            ymax_var=self.global_ymax,
            percent_range=(0.0, 100.0),
            unit_range=(0.0, 1.0),
        )

        std_row = tb.Frame(lf_axes, bootstyle="secondary")
        std_row.pack(fill=X, pady=4)
        tb.Checkbutton(
            std_row,
            text="Show standard volumes (top ticks)",
            variable=self.chrom_show_standards,
        ).pack(side=LEFT)
        tb.Label(std_row, text="Standard Volumes (mL):").pack(side=LEFT, padx=8)
        tb.Entry(std_row, textvariable=self.chrom_standards, width=28).pack(side=LEFT, padx=6)

        self._build_replicates_panel(
            parent=self.chrom_tab,
            title="Replicates (optional)",
            enable_attr="chrom_rep_enable",
            groups_attr="chrom_rep_groups",
            show_mean_attr="chrom_rep_show_mean",
            show_band_attr="chrom_rep_show_band",
            lw_attr="chrom_rep_lw",
            mean_lw_attr="chrom_rep_mean_lw",
        )

    def _build_cdms_tab(self, nb: tb.Notebook):
        self.cdms_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.cdms_tab, text="CD-MS Mass")

        self.cdms_mode = tk.StringVar(value="single")
        self._build_tab_mode_row(
            parent=self.cdms_tab,
            label="Plot Mode:",
            variable=self.cdms_mode,
            options=[("Single CD-MS Mass", "single"), ("Stacked CD-MS Mass", "stacked")],
        )

        self.cdms_two_legends = tk.BooleanVar(value=False)
        self._build_two_legend_labels_section(
            parent=self.cdms_tab,
            title="CD-MS Mass Labels & Colors",
            labels_attr="cdms_labels",
            two_legends_var=self.cdms_two_legends,
            second_labels_attr="cdms_second_labels",
            color_var_attr="cdms_colors",
            default_colors=", ".join(DEFAULT_CDMS_COLORS),
            palette_title="CD-MS Colors",
            toggle_callback=lambda: self._toggle_text_state(
                self.cdms_second_labels,
                self.cdms_two_legends,
            ),
        )

        lf_axes = tb.Labelframe(self.cdms_tab, text="CD-MS Normalization", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)

        self.cdms_normalize = tk.BooleanVar(value=True)
        self.cdms_norm_mode = tk.StringVar(value="percent")

        tb.Checkbutton(lf_axes, text="Normalize to max", variable=self.cdms_normalize).pack(anchor=W, pady=4)

        self._build_norm_mode_toggle(
            parent=lf_axes,
            label="Normalization scale (for CD-MS):",
            mode_var=self.cdms_norm_mode,
            ymin_var=self.global_ymin,
            ymax_var=self.global_ymax,
            percent_range=(0.0, 100.0),
            unit_range=(0.0, 1.0),
        )

        self._build_replicates_panel(
            parent=self.cdms_tab,
            title="Replicates (optional)",
            enable_attr="cdms_rep_enable",
            groups_attr="cdms_rep_groups",
            show_mean_attr="cdms_rep_show_mean",
            show_band_attr="cdms_rep_show_band",
            lw_attr="cdms_rep_lw",
            mean_lw_attr="cdms_rep_mean_lw",
        )

    def _build_mwd_tab(self, nb: tb.Notebook):
        self.mwd_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.mwd_tab, text="MWD Chromatograms")

        lf_files = tb.Labelframe(
            self.mwd_tab, text="MWD Files (drag & drop or browse)", padding=8, bootstyle="secondary"
        )
        lf_files.pack(fill=X, pady=4)

        def file_row(parent, label, var_attr, browse_cmd):
            row = tb.Frame(parent, bootstyle="secondary")
            row.pack(fill=X, pady=2)
            tb.Label(row, text=label).pack(side=LEFT)
            var = tk.StringVar()
            setattr(self, var_attr, var)
            entry = tb.Entry(row, textvariable=var)
            entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
            tb.Button(row, text="Browse", bootstyle="secondary", command=browse_cmd).pack(side=LEFT)
            if DND_AVAILABLE:
                entry.drop_target_register(DND_FILES)
                entry.dnd_bind("<<Drop>>", lambda e, tv=var: self._on_single_file_drop(e, tv))
            return entry

        self.mwd_mau_entry = file_row(lf_files, "mAU file:", "mwd_mau_path", self._browse_mwd_mau)
        self.mwd_mv_entry = file_row(lf_files, "mV file:", "mwd_mv_path", self._browse_mwd_mv)

        lf_labels = tb.Labelframe(self.mwd_tab, text="Labels & Colors", padding=8, bootstyle="secondary")
        lf_labels.pack(fill=X, pady=4)
        row3 = tb.Frame(lf_labels, bootstyle="secondary")
        row3.pack(fill=X, pady=2)

        self.mwd_mau_label = tk.StringVar(value="mAU")
        self.mwd_mv_label = tk.StringVar(value="mV")
        tb.Label(row3, text="mAU label:").pack(side=LEFT)
        tb.Entry(row3, textvariable=self.mwd_mau_label, width=24).pack(side=LEFT, padx=6)
        tb.Label(row3, text="mV label:").pack(side=LEFT)
        tb.Entry(row3, textvariable=self.mwd_mv_label, width=24).pack(side=LEFT, padx=6)

        row4 = tb.Frame(lf_labels, bootstyle="secondary")
        row4.pack(fill=X, pady=2)
        tb.Label(row4, text="Colors (mAU, mV hex):").pack(side=LEFT)
        self.mwd_colors = tk.StringVar(value=", ".join(DEFAULT_MWD_COLORS))
        tb.Entry(row4, textvariable=self.mwd_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row4, text="Pick mAU", bootstyle="info", command=lambda: self._pick_mwd_color(0)).pack(
            side=LEFT, padx=4
        )
        tb.Button(row4, text="Pick mV", bootstyle="info", command=lambda: self._pick_mwd_color(1)).pack(
            side=LEFT, padx=4
        )

        lf_axes = tb.Labelframe(self.mwd_tab, text="MWD Normalization", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=4)

        self.mwd_normalize = tk.BooleanVar(value=True)
        self.mwd_norm_mode = tk.StringVar(value="percent")

        tb.Checkbutton(
            lf_axes, text="Normalize each trace to max", variable=self.mwd_normalize
        ).pack(anchor=W, pady=4)

        self._build_norm_mode_toggle(
            parent=lf_axes,
            label="Normalization scale (for MWD):",
            mode_var=self.mwd_norm_mode,
            ymin_var=self.global_ymin,
            ymax_var=self.global_ymax,
            percent_range=(0.0, 100.0),
            unit_range=(0.0, 1.0),
        )

    def _build_mz_tab(self, nb: tb.Notebook):
        self.mz_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.mz_tab, text="m/z Spectra")

        self.mz_mode = tk.StringVar(value="single")
        self._build_tab_mode_row(
            parent=self.mz_tab,
            label="Plot Mode:",
            variable=self.mz_mode,
            options=[("Single m/z Spectrum", "single"), ("Stacked m/z Spectra", "stacked")],
        )

        self.mz_two_legends = tk.BooleanVar(value=False)
        self._build_two_legend_labels_section(
            parent=self.mz_tab,
            title="m/z Labels & Colors",
            labels_attr="mz_labels",
            two_legends_var=self.mz_two_legends,
            second_labels_attr="mz_second_labels",
            color_var_attr="mz_colors",
            default_colors=", ".join(DEFAULT_MZ_COLORS),
            palette_title="m/z Colors",
            toggle_callback=lambda: self._toggle_text_state(
                self.mz_second_labels,
                self.mz_two_legends,
            ),
        )

        lf_axes = tb.Labelframe(self.mz_tab, text="m/z Normalization", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)

        self.mz_normalize = tk.BooleanVar(value=True)
        self.mz_norm_mode = tk.StringVar(value="percent")

        tb.Checkbutton(lf_axes, text="Normalize to max", variable=self.mz_normalize).pack(anchor=W, pady=4)

        self._build_norm_mode_toggle(
            parent=lf_axes,
            label="Normalization scale (for m/z):",
            mode_var=self.mz_norm_mode,
            ymin_var=self.global_ymin,
            ymax_var=self.global_ymax,
            percent_range=(0.0, 100.0),
            unit_range=(0.0, 1.0),
        )

        self._build_replicates_panel(
            parent=self.mz_tab,
            title="Replicates (optional)",
            enable_attr="mz_rep_enable",
            groups_attr="mz_rep_groups",
            show_mean_attr="mz_rep_show_mean",
            show_band_attr="mz_rep_show_band",
            lw_attr="mz_rep_lw",
            mean_lw_attr="mz_rep_mean_lw",
        )

    def _build_raw_tab(self, nb: tb.Notebook):
        self.raw_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.raw_tab, text="Thermo RAW")

        lf_file = tb.Labelframe(self.raw_tab, text="RAW/mzML File", padding=8, bootstyle="secondary")
        lf_file.pack(fill=X, pady=4)
        row = tb.Frame(lf_file, bootstyle="secondary")
        row.pack(fill=X, pady=2)
        tb.Label(row, text="Path:").pack(side=LEFT)
        self.raw_path = tk.StringVar()
        self.raw_entry = tb.Entry(row, textvariable=self.raw_path)
        self.raw_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row, text="Browse", bootstyle="secondary", command=self._browse_raw_file).pack(side=LEFT)
        if DND_AVAILABLE:
            self.raw_entry.drop_target_register(DND_FILES)
            self.raw_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.raw_path))

        lf_opts = tb.Labelframe(self.raw_tab, text="Options", padding=8, bootstyle="secondary")
        lf_opts.pack(fill=X, pady=4)
        grid = tb.Frame(lf_opts, bootstyle="secondary")
        grid.pack(fill=X)

        self.raw_ms_level = tk.IntVar(value=1)
        self.raw_max_scans = tk.IntVar(value=200)
        self.raw_mode = tk.StringVar(value="apex")
        self.raw_index = tk.IntVar(value=0)
        self.raw_norm = tk.BooleanVar(value=True)
        self.raw_norm_mode = tk.StringVar(value="percent")
        self.raw_ymax_pct = tk.DoubleVar(value=100.0)

        tb.Label(grid, text="MS level:").grid(row=0, column=0, sticky=W)
        tb.Entry(grid, textvariable=self.raw_ms_level, width=6).grid(row=0, column=1, padx=6)
        tb.Label(grid, text="Max scans:").grid(row=0, column=2, sticky=W)
        tb.Entry(grid, textvariable=self.raw_max_scans, width=8).grid(row=0, column=3, padx=6)

        tb.Label(grid, text="Spectrum:").grid(row=1, column=0, sticky=W)
        tb.Radiobutton(grid, text="Apex", variable=self.raw_mode, value="apex").grid(row=1, column=1, sticky=W)
        tb.Radiobutton(grid, text="Index", variable=self.raw_mode, value="index").grid(row=1, column=2, sticky=W)
        tb.Label(grid, text="Index:").grid(row=1, column=3, sticky=W)
        self.raw_index_entry = tb.Entry(grid, textvariable=self.raw_index, width=6)
        self.raw_index_entry.grid(row=1, column=4, padx=6)

        def _toggle_raw_index(*_):
            state = tk.NORMAL if self.raw_mode.get() == "index" else tk.DISABLED
            self.raw_index_entry.config(state=state)

        self.raw_mode.trace_add("write", _toggle_raw_index)
        _toggle_raw_index()

        tb.Checkbutton(lf_opts, text="Normalize spectrum to max", variable=self.raw_norm).pack(anchor=W, pady=4)

        raw_norm_frame = tb.Frame(lf_opts, bootstyle="secondary")
        raw_norm_frame.pack(fill=X, pady=2)
        tb.Label(raw_norm_frame, text="Spectrum normalization scale:").pack(side=LEFT)
        tb.Radiobutton(raw_norm_frame, text="0–100%", variable=self.raw_norm_mode, value="percent").pack(
            side=LEFT, padx=4
        )
        tb.Radiobutton(raw_norm_frame, text="0–1", variable=self.raw_norm_mode, value="unit").pack(
            side=LEFT, padx=4
        )

        self.raw_ymax_entry = tb.Entry(grid, textvariable=self.raw_ymax_pct, width=8)
        tb.Label(grid, text="Spec Y max (normalized):").grid(row=2, column=0, sticky=W, pady=(6, 0))
        self.raw_ymax_entry.grid(row=2, column=1, padx=6, pady=(6, 0))

        self._track_norm_mode(
            mode_var=self.raw_norm_mode,
            min_var=None,
            max_var=self.raw_ymax_pct,
            percent_range=(0.0, 100.0),
            unit_range=(0.0, 1.0),
        )

        def _clamp_raw_ymax(event=None):
            try:
                mode = self.raw_norm_mode.get()
                ymax = float(self.raw_ymax_pct.get())
            except Exception:
                return

            if mode == "percent":
                ymax = max(0.0, min(100.0, ymax))
            else:
                ymax = max(0.0, min(1.0, ymax))

            self.raw_ymax_pct.set(self._format_number_for_entry(ymax))

        self.raw_ymax_entry.bind("<FocusOut>", _clamp_raw_ymax)

        lf_ms = tb.Labelframe(self.raw_tab, text="msconvert (ProteoWizard)", padding=8, bootstyle="secondary")
        lf_ms.pack(fill=X, pady=4)
        row2 = tb.Frame(lf_ms, bootstyle="secondary")
        row2.pack(fill=X, pady=2)
        tb.Label(row2, text="msconvert path (optional):").pack(side=LEFT)
        self.msconvert_path = tk.StringVar()
        self.msconvert_entry = tb.Entry(row2, textvariable=self.msconvert_path)
        self.msconvert_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row2, text="Browse", bootstyle="secondary", command=self._browse_msconvert).pack(side=LEFT)

        self.msconvert_peakpick = tk.BooleanVar(value=True)
        self.msconvert_compress = tk.BooleanVar(value=True)
        self.auto_convert_raw = tk.BooleanVar(value=True)
        self.raw_prefer_native = tk.BooleanVar(value=True)

        tb.Checkbutton(
            lf_ms, text="Centroid (peak-pick) on convert", variable=self.msconvert_peakpick
        ).pack(anchor=W, pady=2)
        tb.Checkbutton(lf_ms, text="Compress mzML (zlib)", variable=self.msconvert_compress).pack(anchor=W, pady=2)
        tb.Checkbutton(
            lf_ms, text="Auto-convert RAW → mzML on open", variable=self.auto_convert_raw
        ).pack(anchor=W, pady=2)
        tb.Checkbutton(
            lf_ms,
            text="Prefer native Thermo MSFileReader on Windows (fallback to msconvert)",
            variable=self.raw_prefer_native,
        ).pack(anchor=W, pady=4)

    # --------------------------- Shared UI Helpers ---------------------------

    def _build_tab_mode_row(
        self,
        parent: tk.Widget,
        label: str,
        variable: tk.StringVar,
        options: Sequence[Tuple[str, str]],
    ):
        frame = tb.Frame(parent, bootstyle="secondary")
        frame.pack(fill=X, pady=4)
        tb.Label(frame, text=label, font=("Arial", 10, "bold")).pack(side=LEFT)
        for text, value in options:
            tb.Radiobutton(frame, text=text, variable=variable, value=value).pack(side=LEFT, padx=8)

    def _build_labels_colors_section(
        self,
        parent: tk.Widget,
        title: str,
        text_widget_attr: str,
        color_var_attr: str,
        default_colors: str,
        palette_title: str,
    ):
        lf = tb.Labelframe(parent, text=title, padding=8, bootstyle="secondary")
        lf.pack(fill=X, pady=4)
        tb.Label(lf, text="Sample Labels (comma-separated; blank = filenames)").pack(anchor=W)
        text = tk.Text(lf, height=2)
        text.pack(fill=BOTH, pady=4)
        setattr(self, text_widget_attr, text)

        tb.Label(lf, text="Colors (comma-separated hex)").pack(anchor=W)
        colors_var = tk.StringVar(value=default_colors)
        setattr(self, color_var_attr, colors_var)

        color_row = tb.Frame(lf, bootstyle="secondary")
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=colors_var).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(
            color_row,
            text="Pick Colors…",
            bootstyle="info",
            command=lambda v=colors_var, t=palette_title: self._open_palette_editor(v, t),
        ).pack(side=LEFT)
        return lf

    def _build_two_legend_labels_section(
        self,
        parent: tk.Widget,
        title: str,
        labels_attr: str,
        two_legends_var: tk.BooleanVar,
        second_labels_attr: str,
        color_var_attr: str,
        default_colors: str,
        palette_title: str,
        toggle_callback: Callable[[], None],
    ):
        lf = tb.Labelframe(parent, text=title, padding=8, bootstyle="secondary")
        lf.pack(fill=X, pady=4)
        tb.Label(lf, text="Labels (comma-separated; blank = filenames)").pack(anchor=W)
        labels_text = tk.Text(lf, height=2)
        labels_text.pack(fill=BOTH, pady=4)
        setattr(self, labels_attr, labels_text)

        tb.Checkbutton(
            lf,
            text="Include second legend near peak",
            variable=two_legends_var,
        ).pack(anchor=W, pady=2)
        tb.Label(
            lf,
            text="Second Labels (one per file; used when second legend is enabled)",
        ).pack(anchor=W)
        second_text = tk.Text(lf, height=2)
        second_text.pack(fill=BOTH, pady=4)
        second_text.config(state=tk.DISABLED)
        setattr(self, second_labels_attr, second_text)
        two_legends_var.trace_add("write", lambda *_: toggle_callback())

        tb.Label(lf, text="Colors (comma-separated hex)").pack(anchor=W)
        colors_var = tk.StringVar(value=default_colors)
        setattr(self, color_var_attr, colors_var)
        color_row = tb.Frame(lf, bootstyle="secondary")
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=colors_var).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(
            color_row,
            text="Pick Colors…",
            bootstyle="info",
            command=lambda v=colors_var, t=palette_title: self._open_palette_editor(v, t),
        ).pack(side=LEFT)
        return lf

    def _build_replicates_panel(
        self,
        parent: tk.Widget,
        title: str,
        enable_attr: str,
        groups_attr: str,
        show_mean_attr: str,
        show_band_attr: str,
        lw_attr: str,
        mean_lw_attr: str,
    ):
        lf = tb.Labelframe(parent, text=title, padding=8, bootstyle="secondary")
        lf.pack(fill=X, pady=6)
        row1 = tb.Frame(lf, bootstyle="secondary")
        row1.pack(fill=X, pady=2)

        enable_var = tk.BooleanVar(value=False)
        setattr(self, enable_attr, enable_var)
        tb.Checkbutton(row1, text="Enable replicates overlay", variable=enable_var).pack(side=LEFT)
        tb.Label(
            row1,
            text="Groups: one line per axis; use file indices or paths (comma-separated)",
        ).pack(side=LEFT, padx=8)

        groups_text = tk.Text(lf, height=3)
        groups_text.pack(fill=X, pady=4)
        setattr(self, groups_attr, groups_text)

        row2 = tb.Frame(lf, bootstyle="secondary")
        row2.pack(fill=X, pady=2)

        show_mean_var = tk.BooleanVar(value=True)
        show_band_var = tk.BooleanVar(value=False)
        rep_lw_var = tk.DoubleVar(value=1.2)
        mean_lw_var = tk.DoubleVar(value=2.2)

        setattr(self, show_mean_attr, show_mean_var)
        setattr(self, show_band_attr, show_band_var)
        setattr(self, lw_attr, rep_lw_var)
        setattr(self, mean_lw_attr, mean_lw_var)

        tb.Checkbutton(row2, text="Show mean line", variable=show_mean_var).pack(side=LEFT)
        tb.Checkbutton(row2, text="Show SEM band", variable=show_band_var).pack(side=LEFT, padx=8)
        tb.Label(row2, text="Replicate LW:").pack(side=LEFT, padx=(12, 2))
        tb.Entry(row2, textvariable=rep_lw_var, width=6).pack(side=LEFT)
        tb.Label(row2, text="Mean LW:").pack(side=LEFT, padx=(12, 2))
        tb.Entry(row2, textvariable=mean_lw_var, width=6).pack(side=LEFT)

    @staticmethod
    def _format_number_for_entry(value: float) -> str:
        if abs(value - round(value)) < 1e-9:
            return str(int(round(value)))
        return f"{value:.6f}".rstrip("0").rstrip(".")

    def _clamp_y_entry(self, entry_var: tk.DoubleVar, mode_var: tk.StringVar,
                       percent_range=(0.0, 100.0), unit_range=(0.0, 1.0)):
        """
        Clamp a Y-limit entry according to the current normalization mode.

        percent_range: bounds when mode == 'percent'
        unit_range:    bounds when mode == 'unit'
        """
        try:
            mode = mode_var.get()
            val = float(entry_var.get())
        except Exception:
            return

        if mode == "percent":
            lo, hi = percent_range
        else:
            lo, hi = unit_range

        val = max(lo, min(hi, val))
        entry_var.set(self._format_number_for_entry(val))

    def _build_norm_mode_toggle(
        self,
        parent: tk.Widget,
        label: str,
        mode_var: tk.StringVar,
        ymin_var: Optional[tk.DoubleVar],
        ymax_var: tk.DoubleVar,
        percent_range: Tuple[float, float],
        unit_range: Tuple[float, float],
    ):
        frame = tb.Frame(parent, bootstyle="secondary")
        frame.pack(fill=X, pady=2)
        tb.Label(frame, text=label).pack(side=LEFT)
        tb.Radiobutton(frame, text="0–100%", variable=mode_var, value="percent").pack(side=LEFT, padx=4)
        tb.Radiobutton(frame, text="0–1", variable=mode_var, value="unit").pack(side=LEFT, padx=4)
        self._track_norm_mode(
            mode_var=mode_var,
            min_var=ymin_var,
            max_var=ymax_var,
            percent_range=percent_range,
            unit_range=unit_range,
        )

    def _track_norm_mode(
        self,
        mode_var: tk.StringVar,
        min_var: Optional[tk.DoubleVar],
        max_var: tk.DoubleVar,
        percent_range: Tuple[float, float],
        unit_range: Tuple[float, float],
    ):
        def _update_on_mode_change(*_):
            try:
                mode = mode_var.get()
                vmax = float(max_var.get())
                vmin = float(min_var.get()) if min_var is not None else None
            except Exception:
                return

            p_lo, p_hi = percent_range
            u_lo, u_hi = unit_range

            if mode == "percent":
                if vmax <= 1.0:
                    vmax *= 100.0
                if vmin is not None and vmin <= 1.0:
                    vmin *= 100.0
                vmax = max(p_lo, min(p_hi, vmax))
                if vmin is not None:
                    vmin = max(p_lo, min(p_hi, vmin))
            else:
                if vmax > 1.0:
                    vmax /= 100.0
                if vmin is not None and vmin > 1.0:
                    vmin /= 100.0
                vmax = max(u_lo, min(u_hi, vmax))
                if vmin is not None:
                    vmin = max(u_lo, min(u_hi, vmin))

            if vmin is not None:
                vmin, vmax2 = clamp_range(
                    vmin, vmax,
                    *(percent_range if mode == "percent" else unit_range)
                )
                min_var.set(self._format_number_for_entry(vmin))
                max_var.set(self._format_number_for_entry(vmax2))
            else:
                max_var.set(self._format_number_for_entry(vmax))

            # If this mode_var is the chrom normalization mode, also
            # clamp the global Y-limits to 0–100 or 0–1.
            if hasattr(self, "chrom_norm_mode") and mode_var is self.chrom_norm_mode:
                if hasattr(self, "global_ymin") and hasattr(self, "global_ymax"):
                    self._clamp_y_entry(
                        entry_var=self.global_ymin,
                        mode_var=self.chrom_norm_mode,
                        percent_range=(0.0, 100.0),
                        unit_range=(0.0, 1.0),
                    )
                    self._clamp_y_entry(
                        entry_var=self.global_ymax,
                        mode_var=self.chrom_norm_mode,
                        percent_range=(0.0, 100.0),
                        unit_range=(0.0, 1.0),
                    )

        mode_var.trace_add("write", _update_on_mode_change)

    def _toggle_text_state(self, text_widget: tk.Text, flag_var: tk.BooleanVar):
        state = tk.NORMAL if flag_var.get() else tk.DISABLED
        text_widget.config(state=state)

    # --------------------------- Files Panel ---------------------------

    def _build_file_panel(self, parent: tk.Widget):
        lf = tb.Labelframe(
            parent,
            text="Files (Chromatograms, CD-MS Mass, MWD, m/z)",
            padding=8,
            bootstyle="secondary",
        )
        lf.pack(fill=BOTH, expand=YES, pady=6)
        btns = tb.Frame(lf, bootstyle="secondary")
        btns.pack(fill=X, pady=4)
        tb.Button(btns, text="Add Files...", bootstyle="success", command=self._add_files).pack(side=LEFT)
        tb.Button(btns, text="Remove Selected", bootstyle="warning", command=self._remove_selected).pack(
            side=LEFT, padx=6
        )
        tb.Button(btns, text="Clear", bootstyle="danger", command=self._clear_files).pack(side=LEFT)

        self.file_tree = tb.Treeview(
            lf, columns=("idx", "path"), show="headings", height=8, bootstyle="info"
        )
        self.file_tree.heading("idx", text="#")
        self.file_tree.column("idx", width=40, anchor=CENTER)
        self.file_tree.heading("path", text="File Path")
        self.file_tree.column("path", width=700, anchor=W)
        self.file_tree.pack(fill=BOTH, expand=YES, pady=6)

        if DND_AVAILABLE:
            self.file_tree.drop_target_register(DND_FILES)
            self.file_tree.dnd_bind("<<Drop>>", self._on_files_drop)

    # --------------------------- Utility / Common Actions ---------------------------

    def _set_status(self, text: str):
        self.status_label.configure(text=f"Status: {text}")
        self.root.update_idletasks()

    def _browse_save(self):
        path = filedialog.asksaveasfilename(
            title="Save figure as...",
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")],
        )
        if path:
            self.save_var.set(path)

    def _refresh_files(self):
        for i in self.file_tree.get_children():
            self.file_tree.delete(i)
        for idx, p in enumerate(self.files):
            self.file_tree.insert("", "end", values=(idx, p))

    def _add_files(self):
        tab = self._active_tab()
        if tab == "Chromatograms":
            paths = filedialog.askopenfilenames(
                title="Select Chromatogram Excel files",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
            )
        elif tab == "CD-MS Mass":
            paths = filedialog.askopenfilenames(
                title="Select CD-MS mass text files",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            )
        elif tab == "m/z Spectra":
            paths = filedialog.askopenfilenames(
                title="Select m/z spectrum text files",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            )
        else:  # MWD
            paths = filedialog.askopenfilenames(
                title="Select mAU/mV Excel files (first=mAU, second=mV)",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
            )
            if paths:
                self.files = list(paths)[:2]
                self._refresh_files()
                return
        if paths:
            self.files.extend(paths)
            self._refresh_files()

    def _remove_selected(self):
        selected = self.file_tree.selection()
        if not selected:
            return
        idxs = sorted(
            (int(self.file_tree.item(i, "values")[0]) for i in selected),
            reverse=True,
        )
        for idx in idxs:
            if 0 <= idx < len(self.files):
                self.files.pop(idx)
        self._refresh_files()

    def _clear_files(self):
        self.files = []
        self._refresh_files()

    def _active_tab(self) -> str:
        if self.chrom_tab.winfo_ismapped():
            return "Chromatograms"
        if self.cdms_tab.winfo_ismapped():
            return "CD-MS Mass"
        if self.mwd_tab.winfo_ismapped():
            return "MWD Chromatograms"
        if self.mz_tab.winfo_ismapped():
            return "m/z Spectra"
        return "Thermo RAW"

    def _on_tab_changed(self, event=None):
        """
        Per-tab default ranges for the global X-lim.

        Chromatograms: 0–20 mL
        CD-MS Mass:    0–500 kDa
        MWD:           0–20 mL
        m/z Spectra:   0–20000 m/z
        Thermo RAW:    leave global X-lim unchanged (RAW uses RT/m/z from data).
        """
        try:
            tab_text = self.nb.tab(self.nb.select(), "text")
        except Exception:
            tab_text = ""

        # Defaults per tab
        if tab_text == "Chromatograms":
            xmin, xmax = 0.0, 20.0
        elif tab_text == "CD-MS Mass":
            xmin, xmax = 0.0, 500.0
        elif tab_text == "MWD Chromatograms":
            xmin, xmax = 0.0, 20.0
        elif tab_text == "m/z Spectra":
            xmin, xmax = 0.0, 20000.0
        else:
            # Thermo RAW or unknown: don't change global X-lim
            return

        # Set global X limits to the tab defaults
        self.global_xmin.set(xmin)
        self.global_xmax.set(xmax)

        # Optionally, you can also reset global Y-lim when switching tabs.
        # Commented out here; uncomment if desired:
        # self.global_ymin.set(0.0)
        # self.global_ymax.set(100.0)

    def _on_files_drop(self, event):
        try:
            paths = list(self.root.splitlist(event.data))
        except Exception:
            paths = event.data.split()
        if not paths:
            return
        tab = self._active_tab()
        if tab == "MWD Chromatograms":
            self.files = paths[:2]
            self._refresh_files()
            self._set_status(f"Set MWD files via drag-and-drop ({len(self.files)} file(s)).")
        else:
            self.files.extend(paths)
            self._refresh_files()
            self._set_status(f"Added {len(paths)} file(s) via drag-and-drop.")

    def _on_single_file_drop(self, event, target_var: tk.StringVar):
        try:
            paths = list(self.root.splitlist(event.data))
        except Exception:
            paths = event.data.split()
        if not paths:
            return
        p = paths[0]
        ext = os.path.splitext(p)[1].lower()
        if ext == ".raw" and bool(self.auto_convert_raw.get()):
            def worker():
                result = self.convert_raw_to_mzml(p)
                if isinstance(result, str) and result.startswith("__CONVERT_ERROR__:"):
                    exc_text = result.split(":", 1)[1]
                    self.root.after(
                        0,
                        lambda: Messagebox.show_error(
                            f"RAW→mzML conversion failed:\n{exc_text}",
                            title="Conversion error",
                        ),
                    )
                    self.root.after(0, lambda: target_var.set(p))
                elif result:
                    self.root.after(0, lambda: target_var.set(result))
                else:
                    self.root.after(0, lambda: target_var.set(p))

            threading.Thread(target=worker, daemon=True).start()
        else:
            target_var.set(p)
            self._set_status("Set file via drag-and-drop.")

    # --------------------------- Color Pickers ---------------------------

    def _pick_color(self, initial: str = "#FFFFFF") -> Optional[str]:
        result = askcolor(color=initial, title="Pick a color")
        return result[1] if result and result[1] else None

    def _open_palette_editor(self, target_var: tk.StringVar, title: str):
        top = tk.Toplevel(self.root)
        top.title(title)
        top.geometry("420x300")
        try:
            top.configure(bg=self.style.colors.bg)
        except Exception:
            pass
        colors = parse_str_list(target_var.get())
        lb_frame = tb.Frame(top, padding=6, bootstyle="secondary")
        lb_frame.pack(fill=BOTH, expand=YES)
        lb = tk.Listbox(lb_frame, height=10)
        lb.pack(fill=BOTH, expand=YES)

        def refresh_list():
            lb.delete(0, tk.END)
            for c in colors:
                lb.insert(tk.END, c)

        refresh_list()

        ctrl = tb.Frame(top, padding=6, bootstyle="secondary")
        ctrl.pack(fill=X)

        tb.Button(
            ctrl,
            text="Add…",
            bootstyle="info",
            command=lambda: self._palette_add(colors, refresh_list),
        ).pack(side=LEFT, padx=4)
        tb.Button(
            ctrl,
            text="Remove",
            bootstyle="warning",
            command=lambda: self._palette_remove(colors, lb, refresh_list),
        ).pack(side=LEFT, padx=4)
        tb.Button(
            ctrl,
            text="Clear",
            bootstyle="danger",
            command=lambda: self._palette_clear(colors, refresh_list),
        ).pack(side=LEFT, padx=4)
        tb.Button(
            ctrl,
            text="Apply",
            bootstyle="success",
            command=lambda: self._palette_apply(colors, target_var, top),
        ).pack(side=RIGHT, padx=4)
        tb.Button(ctrl, text="Cancel", bootstyle="secondary", command=top.destroy).pack(side=RIGHT, padx=4)

    def _palette_add(self, colors: List[str], refresh_list: Callable[[], None]):
        start_color = colors[-1] if colors else "#00A9A5"
        c = self._pick_color(initial=start_color)
        if c:
            colors.append(c)
            refresh_list()

    def _palette_remove(self, colors: List[str], lb: tk.Listbox, refresh_list: Callable[[], None]):
        selection = lb.curselection()
        if selection:
            colors.pop(selection[0])
            refresh_list()

    def _palette_clear(self, colors: List[str], refresh_list: Callable[[], None]):
        colors.clear()
        refresh_list()

    def _palette_apply(self, colors: List[str], target_var: tk.StringVar, top: tk.Toplevel):
        target_var.set(", ".join(colors))
        top.destroy()

    def _pick_mwd_color(self, index: int):
        colors = parse_str_list(self.mwd_colors.get()) or DEFAULT_MWD_COLORS.copy()
        while len(colors) < 2:
            colors.append("#000000")
        start_color = colors[index]
        c = self._pick_color(initial=start_color)
        if c:
            colors[index] = c
            self.mwd_colors.set(", ".join(colors))

    # --------------------------- Browse Helpers ---------------------------

    def _browse_msconvert(self):
        path = filedialog.askopenfilename(
            title="Select msconvert.exe",
            filetypes=[("Executable", "*.exe;*"), ("All files", "*.*")],
        )
        if path:
            self.msconvert_path.set(path)
            self._set_status("msconvert path set.")

    def _browse_raw_file(self):
        path = filedialog.askopenfilename(
            title="Select Thermo RAW or mzML",
            filetypes=[("RAW/mzML", "*.raw *.mzML"), ("All files", "*.*")],
        )
        if not path:
            return
        ext = os.path.splitext(path)[1].lower()
        if ext == ".raw" and bool(self.auto_convert_raw.get()):
            def worker():
                result = self.convert_raw_to_mzml(path)
                if isinstance(result, str) and result.startswith("__CONVERT_ERROR__:"):
                    exc_text = result.split(":", 1)[1]
                    self.root.after(
                        0,
                        lambda: Messagebox.show_error(
                            f"RAW→mzML conversion failed:\n{exc_text}",
                            title="Conversion error",
                        ),
                    )
                    self.root.after(0, lambda: self.raw_path.set(path))
                elif result:
                    self.root.after(0, lambda: self.raw_path.set(result))
                else:
                    self.root.after(0, lambda: self.raw_path.set(path))

            threading.Thread(target=worker, daemon=True).start()
        else:
            self.raw_path.set(path)
            self._set_status("RAW/mzML file selected.")

    def _browse_mwd_mau(self):
        path = filedialog.askopenfilename(
            title="Select mAU Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if path:
            self.mwd_mau_path.set(path)
            self._set_status("mAU file selected.")

    def _browse_mwd_mv(self):
        path = filedialog.askopenfilename(
            title="Select mV Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if path:
            self.mwd_mv_path.set(path)
            self._set_status("mV file selected.")

    # --------------------------- Replicates Helpers ---------------------------

    def _parse_rep_groups(self, text: str) -> List[List[str]]:
        groups: List[List[str]] = []
        for line in text.strip().splitlines():
            items = [s.strip() for s in line.replace(";", ",").split(",") if s.strip()]
            if not items:
                continue
            group: List[str] = []
            for it in items:
                try:
                    idx = int(it)
                    if 0 <= idx < len(self.files):
                        group.append(self.files[idx])
                        continue
                except Exception:
                    pass
                if os.path.isfile(it):
                    group.append(it)
            if group:
                groups.append(group)
        return groups

    @staticmethod
    def _interpolate_to_grid(x: np.ndarray, y: np.ndarray, xgrid: np.ndarray) -> np.ndarray:
        return np.interp(xgrid, x, y, left=np.nan, right=np.nan)

    def _overlay_replicates_generic(
        self,
        axes,
        groups,
        colors,
        normalize: bool,
        xgrid: np.ndarray,
        x_selector: Callable[[str], Tuple[np.ndarray, np.ndarray]],
        norm_mode: str,
        rep_lw: float,
        mean_lw: float,
        show_mean: bool,
        show_band: bool,
    ):
        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
        for i, ax in enumerate(axes_list):
            if i >= len(groups):
                break
            gfiles = groups[i]
            traces = []
            for fp in gfiles:
                try:
                    x, y = x_selector(fp)
                    y = normalize_y_values(y, normalize, norm_mode)
                    ygrid = self._interpolate_to_grid(x, y, xgrid)
                    traces.append(ygrid)
                    ax.plot(x, y, color=colors[i % len(colors)], alpha=0.35, lw=rep_lw)
                except Exception:
                    continue
            if traces and show_mean:
                arr = np.array(traces)
                mean = np.nanmean(arr, axis=0)
                if show_band and arr.shape[0] > 1:
                    sem = np.nanstd(arr, axis=0) / np.sqrt(arr.shape[0])
                    ax.fill_between(
                        xgrid,
                        mean - sem,
                        mean + sem,
                        color=colors[i % len(colors)],
                        alpha=0.15,
                        linewidth=0,
                    )
                ax.plot(xgrid, mean, color=colors[i % len(colors)], lw=mean_lw)

    def _overlay_chrom_replicates(self, axes, groups, colors, normalize: bool, xlim):
        xgrid = np.linspace(xlim[0], xlim[1], 1200)
        self._overlay_replicates_generic(
            axes=axes,
            groups=groups,
            colors=colors,
            normalize=normalize,
            xgrid=xgrid,
            x_selector=lambda fp: (
                read_chrom_file(fp)["ml"].astype(float).values,
                read_chrom_file(fp)["mAU"].astype(float).values,
            ),
            norm_mode=self.chrom_norm_mode.get(),
            rep_lw=float(self.chrom_rep_lw.get()),
            mean_lw=float(self.chrom_rep_mean_lw.get()),
            show_mean=bool(self.chrom_rep_show_mean.get()),
            show_band=bool(self.chrom_rep_show_band.get()),
        )

    def _overlay_cdms_replicates(self, axes, groups, colors, normalize: bool, x_range):
        def selector(fp: str):
            df = _load_mass(fp)
            df = df[(df["mass_kDa"] >= x_range[0]) & (df["mass_kDa"] <= x_range[1])].copy()
            return (
                df["mass_kDa"].astype(float).values,
                df["intensity"].astype(float).values,
            )

        xgrid = np.linspace(x_range[0], x_range[1], 2000)
        self._overlay_replicates_generic(
            axes=axes,
            groups=groups,
            colors=colors,
            normalize=normalize,
            xgrid=xgrid,
            x_selector=selector,
            norm_mode=self.cdms_norm_mode.get(),
            rep_lw=float(self.cdms_rep_lw.get()),
            mean_lw=float(self.cdms_rep_mean_lw.get()),
            show_mean=bool(self.cdms_rep_show_mean.get()),
            show_band=bool(self.cdms_rep_show_band.get()),
        )

    def _overlay_mz_replicates(self, axes, groups, colors, normalize: bool, x_range):
        def selector(fp: str):
            df = read_mz_spectrum_file(fp)
            df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
            return df["mz"].astype(float).values, df["intensity"].astype(float).values

        xgrid = np.linspace(x_range[0], x_range[1], 4000)
        self._overlay_replicates_generic(
            axes=axes,
            groups=groups,
            colors=colors,
            normalize=normalize,
            xgrid=xgrid,
            x_selector=selector,
            norm_mode=self.mz_norm_mode.get(),
            rep_lw=float(self.mz_rep_lw.get()),
            mean_lw=float(self.mz_rep_mean_lw.get()),
            show_mean=bool(self.mz_rep_show_mean.get()),
            show_band=bool(self.mz_rep_show_band.get()),
        )

    # --------------------------- Plot Dispatcher ---------------------------

    def _plot(self):
        cfg = CommonConfig(
            save=self.save_var.get() or None,
            font_name=self.font_name.get(),
            font_size=int(self.font_size.get()),
            dpi=int(self.dpi.get()),
            figsize=(float(self.fig_w.get()), float(self.fig_h.get())),
            transparent=bool(self.transparent_export.get()),
            show=bool(self.show_plot.get()),
        )
        try:
            self.progress.start(10)
            self._set_status("Plotting...")
            tab = self._active_tab()

            if tab == "Chromatograms":
                self._plot_chrom(cfg)
            elif tab == "CD-MS Mass":
                self._plot_cdms(cfg)
            elif tab == "MWD Chromatograms":
                self._plot_mwd(cfg)
            elif tab == "m/z Spectra":
                self._plot_mz(cfg)
            else:
                self._plot_raw(cfg)

        except Exception as e:
            self._set_status("Error")
            Messagebox.show_error(f"Plotting failed:\n{e}", title="Error")
        finally:
            self.progress.stop()

    # ---------------- Per-tab plot methods ----------------

    def _plot_chrom(self, cfg: CommonConfig):
        raw_labels = self.chrom_labels.get("1.0", "end").strip()
        labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
        colors = parse_str_list(self.chrom_colors.get()) or DEFAULT_CHROM_COLORS
        std_vols = parse_float_list(self.chrom_standards.get()) if self.chrom_show_standards.get() else None
        normalize = bool(self.chrom_normalize.get())
        norm_mode = self.chrom_norm_mode.get()

        xlim = (float(self.global_xmin.get()), float(self.global_xmax.get()))
        ylim = (float(self.global_ymin.get()), float(self.global_ymax.get()))
        try:
            hspace = float(self.global_subplot_spacing.get())
        except Exception:
            hspace = 0.55

        lw = float(self.line_width.get())
        x_lab = self.x_label.get().strip() or None
        y_lab = self.y_label.get().strip() or None

        mode = self.chrom_mode.get()
        if mode == "single":
            if len(self.files) != 1:
                Messagebox.show_error(
                    "Select exactly one Excel file for Single Chromatogram.", title="Selection error"
                )
                return
            fig, axes = plot_chrom_stack(
                files=self.files[:1],
                sample_labels=labels[:1],
                colors=colors[:1],
                standard_volumes=std_vols,
                normalize=normalize,
                norm_mode=norm_mode,
                font_name=cfg.font_name,
                font_size=cfg.font_size,
                xlim=xlim,
                ylim=ylim,
                figsize=cfg.figsize,
                dpi=DISPLAY_DPI,
                line_width=lw,
                x_label=x_lab,
                y_label=y_lab,
                save=None,
                show=cfg.show,
            )
        else:
            if len(labels) != len(self.files):
                Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
                return
            fig, axes = plot_chrom_stack(
                files=self.files,
                sample_labels=labels,
                colors=colors,
                standard_volumes=std_vols,
                normalize=normalize,
                norm_mode=norm_mode,
                font_name=cfg.font_name,
                font_size=cfg.font_size,
                xlim=xlim,
                ylim=ylim,
                figsize=cfg.figsize,
                dpi=DISPLAY_DPI,
                hspace=hspace,
                line_width=lw,
                x_label=x_lab,
                y_label=y_lab,
                save=None,
                show=cfg.show,
            )

        if bool(self.chrom_rep_enable.get()):
            groups = self._parse_rep_groups(self.chrom_rep_groups.get("1.0", "end"))
            axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
            if mode == "single":
                groups = groups[:1]
            elif len(groups) != len(axes_list):
                Messagebox.show_error(
                    "Replicate groups must have one line per axis (stacked).",
                    title="Replicates mismatch",
                )
            else:
                self._overlay_chrom_replicates(axes_list, groups, colors, normalize, xlim)

        self._render_figure(fig)
        self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
        self._set_status("Chromatogram plot complete.")

    def _plot_cdms(self, cfg: CommonConfig):
        raw_labels = self.cdms_labels.get("1.0", "end").strip()
        labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
        use_second = bool(self.cdms_two_legends.get())
        second_raw = self.cdms_second_labels.get("1.0", "end").strip()
        second_labels = parse_str_list(second_raw) if (use_second and second_raw) else [None] * len(self.files)
        colors = parse_str_list(self.cdms_colors.get()) or DEFAULT_CDMS_COLORS

        # Global X/Y ranges (normalized Y range)
        x_range = (float(self.global_xmin.get()), float(self.global_xmax.get()))
        normalize = bool(self.cdms_normalize.get())
        norm_mode = self.cdms_norm_mode.get()
        try:
            hspace = float(self.global_subplot_spacing.get())
        except Exception:
            hspace = 0.6
        mode = self.cdms_mode.get()

        lw = float(self.line_width.get())
        x_lab = self.x_label.get().strip() or None
        y_lab = self.y_label.get().strip() or None

        y_min_norm = float(self.global_ymin.get()) if normalize else None
        y_max_norm = float(self.global_ymax.get()) if normalize else None

        if mode == "single":
            if len(self.files) != 1:
                Messagebox.show_error(
                    "Select exactly one text file for Single CD-MS Mass.", title="Selection error"
                )
                return
            second_label_single = second_labels[0] if use_second else None
            fig, ax = plot_single_cdms_mass(
                dataset=self.files[0],
                label=labels[0],
                second_label=second_label_single,
                x_range=x_range,
                color=colors[0],
                normalize=normalize,
                norm_mode=norm_mode,
                y_min=y_min_norm,
                y_max=y_max_norm,
                font_name=cfg.font_name,
                font_size=cfg.font_size,
                figsize=cfg.figsize,
                dpi=DISPLAY_DPI,
                line_width=lw,
                x_label=x_lab,
                y_label=y_lab,
                save=None,
                show=cfg.show,
            )
            if bool(self.cdms_rep_enable.get()):
                groups = self._parse_rep_groups(self.cdms_rep_groups.get("1.0", "end"))[:1]
                self._overlay_cdms_replicates([ax], groups, [colors[0]], normalize, x_range)
            self._render_figure(fig)
            self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
            self._set_status("CD-MS mass (single) plot complete.")
            return

        # stacked
        if len(labels) != len(self.files):
            Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
            return
        if use_second and len(second_labels) != len(self.files):
            Messagebox.show_error(
                "Second labels must match number of files when second legend is enabled.",
                title="Second label mismatch",
            )
            return

        fig, axes = plot_cdms_mass_stack(
            datasets=self.files,
            labels=labels,
            second_labels=second_labels if use_second else [None] * len(self.files),
            x_range=x_range,
            colors=colors,
            y_min=y_min_norm,
            y_max=y_max_norm,
            normalize=normalize,
            norm_mode=norm_mode,
            font_name=cfg.font_name,
            font_size=cfg.font_size,
            figsize=cfg.figsize,
            dpi=DISPLAY_DPI,
            hspace=hspace,
            line_width=lw,
            x_label=x_lab,
            y_label=y_lab,
            save=None,
            show=cfg.show,
        )

        if bool(self.cdms_rep_enable.get()):
            groups = self._parse_rep_groups(self.cdms_rep_groups.get("1.0", "end"))
            axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
            if len(groups) != len(axes_list):
                Messagebox.show_error(
                    "Replicate groups must have one line per axis (stacked).",
                    title="Replicates mismatch",
                )
            else:
                self._overlay_cdms_replicates(axes_list, groups, colors, normalize, x_range)
        self._render_figure(fig)
        self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
        self._set_status("CD-MS mass (stacked) plot complete.")

    def _plot_mwd(self, cfg: CommonConfig):
        path_mau = self.mwd_mau_path.get().strip()
        path_mv = self.mwd_mv_path.get().strip()
        if not path_mau or not path_mv:
            Messagebox.show_error("Please select both mAU and mV files.", title="MWD files missing")
            return
        colors = parse_str_list(self.mwd_colors.get()) or DEFAULT_MWD_COLORS
        if len(colors) < 2:
            Messagebox.show_error("Please provide two colors (mAU, mV).", title="Colors required")
            return

        xlim = (float(self.global_xmin.get()), float(self.global_xmax.get()))
        normalize = bool(self.mwd_normalize.get())
        norm_mode = self.mwd_norm_mode.get()

        y_min_norm = float(self.global_ymin.get()) if normalize else None
        y_max_norm = float(self.global_ymax.get()) if normalize else None

        lw = float(self.line_width.get())
        x_lab = self.x_label.get().strip() or None
        y_lab = self.y_label.get().strip() or None

        fig, (ax1, ax2) = plot_mwd_chromatograms(
            mau_file=path_mau,
            mv_file=path_mv,
            mau_label=self.mwd_mau_label.get(),
            mv_label=self.mwd_mv_label.get(),
            colors=colors[:2],
            xlim=xlim,
            normalize=normalize,
            norm_mode=norm_mode,
            y_min=y_min_norm,
            y_max=y_max_norm,
            font_name=cfg.font_name,
            font_size=cfg.font_size,
            figsize=cfg.figsize,
            dpi=DISPLAY_DPI,
            line_width=lw,
            x_label=x_lab,
            y_label_left=y_lab,
            y_label_right=y_lab,
            save=None,
            show=cfg.show,
        )
        self._render_figure(fig)
        self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
        self._set_status("MWD Chromatograms plot complete.")

    def _plot_mz(self, cfg: CommonConfig):
        raw_labels = self.mz_labels.get("1.0", "end").strip()
        labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
        use_second = bool(self.mz_two_legends.get())
        second_raw = self.mz_second_labels.get("1.0", "end").strip()
        second_labels = parse_str_list(second_raw) if (use_second and second_raw) else [None] * len(self.files)
        colors = parse_str_list(self.mz_colors.get()) or DEFAULT_MZ_COLORS

        x_range = (float(self.global_xmin.get()), float(self.global_xmax.get()))
        normalize = bool(self.mz_normalize.get())
        norm_mode = self.mz_norm_mode.get()
        mode = self.mz_mode.get()
        try:
            hspace = float(self.global_subplot_spacing.get())
        except Exception:
            hspace = 0.6

        lw = float(self.line_width.get())
        x_lab = self.x_label.get().strip() or None
        y_lab = self.y_label.get().strip() or None
        y_min_norm = float(self.global_ymin.get()) if normalize else None
        y_max_norm = float(self.global_ymax.get()) if normalize else None

        if mode == "single":
            if len(self.files) != 1:
                Messagebox.show_error(
                    "Select exactly one text file for Single m/z Spectrum.",
                    title="Selection error",
                )
                return
            second_label_single = second_labels[0] if use_second else None
            fig, ax = plot_single_mz_spectrum(
                dataset=self.files[0],
                label=labels[0],
                second_label=second_label_single,
                x_range=x_range,
                color=colors[0],
                normalize=normalize,
                norm_mode=norm_mode,
                y_min=y_min_norm,
                y_max=y_max_norm,
                font_name=cfg.font_name,
                font_size=cfg.font_size,
                figsize=cfg.figsize,
                dpi=DISPLAY_DPI,
                line_width=lw,
                x_label=x_lab,
                y_label=y_lab,
                save=None,
                show=cfg.show,
            )
            if bool(self.mz_rep_enable.get()):
                groups = self._parse_rep_groups(self.mz_rep_groups.get("1.0", "end"))[:1]
                self._overlay_mz_replicates([ax], groups, [colors[0]], normalize, x_range)
            self._render_figure(fig)
            self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
            self._set_status("m/z spectrum (single) plot complete.")
            return

        # stacked
        if len(labels) != len(self.files):
            Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
            return
        if use_second and len(second_labels) != len(self.files):
            Messagebox.show_error(
                "Second labels must match number of files when second legend is enabled.",
                title="Second label mismatch",
            )
            return

        fig, axes = plot_mz_spectrum_stack(
            datasets=self.files,
            labels=labels,
            second_labels=second_labels if use_second else [None] * len(self.files),
            x_range=x_range,
            colors=colors,
            y_min=y_min_norm,
            y_max=y_max_norm,
            normalize=normalize,
            norm_mode=norm_mode,
            font_name=cfg.font_name,
            font_size=cfg.font_size,
            figsize=cfg.figsize,
            dpi=DISPLAY_DPI,
            hspace=hspace,
            line_width=lw,
            x_label=x_lab,
            y_label=y_lab,
            save=None,
            show=cfg.show,
        )

        if bool(self.mz_rep_enable.get()):
            groups = self._parse_rep_groups(self.mz_rep_groups.get("1.0", "end"))
            axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
            if len(groups) != len(axes_list):
                Messagebox.show_error(
                    "Replicate groups must have one line per axis (stacked).",
                    title="Replicates mismatch",
                )
            else:
                self._overlay_mz_replicates(axes_list, groups, colors, normalize, x_range)
        self._render_figure(fig)
        self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
        self._set_status("m/z spectra (stacked) plot complete.")

    def _plot_raw(self, cfg: CommonConfig):
        raw_path = self.raw_path.get().strip()
        if not raw_path:
            Messagebox.show_error("Please pick a RAW or mzML file.", title="RAW file missing")
            return

        if raw_path.lower().endswith(".raw") and bool(self.auto_convert_raw.get()):
            mz = self.convert_raw_to_mzml(raw_path)
            if mz:
                raw_path = mz
            else:
                return

        norm_mode = self.raw_norm_mode.get()

        lw = float(self.line_width.get())
        x_lab = self.x_label.get().strip() or None
        y_lab = self.y_label.get().strip() or None

        spec_ymin = float(self.global_ymin.get()) if bool(self.raw_norm.get()) else None
        spec_ymax = float(self.global_ymax.get()) if bool(self.raw_norm.get()) else None

        fig, (ax1, ax2) = plot_raw_chrom_and_spectrum(
            raw_path=raw_path,
            ms_level=int(self.raw_ms_level.get()),
            max_scans=int(self.raw_max_scans.get()),
            spectrum_mode=self.raw_mode.get(),
            spectrum_index=(int(self.raw_index.get()) if self.raw_mode.get() == "index" else None),
            normalize_spectrum=bool(self.raw_norm.get()),
            norm_mode=norm_mode,
            spec_y_min=spec_ymin,
            spec_y_max=spec_ymax,
            figsize=cfg.figsize,
            dpi=DISPLAY_DPI,
            font_name=cfg.font_name,
            font_size=cfg.font_size,
            tick_length=5.0,
            tick_width=1.0,
            line_width_tic=lw,
            line_width_spec=lw,
            x_label_tic=x_lab,
            y_label_tic=y_lab,
            x_label_spec=x_lab,
            y_label_spec=y_lab,
            save=None,
            show=cfg.show,
            prefer_native=bool(self.raw_prefer_native.get()),
            msconvert_path=(self.msconvert_path.get().strip() or None),
            peakpick_on_convert=bool(self.msconvert_peakpick.get()),
            compress_mzml=bool(self.msconvert_compress.get()),
        )
        self._render_figure(fig)
        self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
        self._set_status("Thermo RAW plot (TIC + m/z) complete.")

    # --------------------------- Figure Embedding & Export ---------------------------

    def _render_figure(self, fig):
        for child in self.plot_area.winfo_children():
            child.destroy()
        self.plot_area.configure(style="secondary.TFrame")
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_area)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=YES)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_area)
        self.toolbar.update()
        self.toolbar.pack(fill=X)

    def _export_figure(self, fig, save_path: Optional[str], dpi: int, transparent: bool):
        if not save_path:
            return
        try:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=transparent)
            self._set_status(f"Saved to {save_path} (transparent={transparent})")
        except Exception as e:
            Messagebox.show_error(f"Failed to save figure:\n{e}", title="Export error")

    # --------------------------- RAW → mzML Conversion ---------------------------

    def convert_raw_to_mzml(self, raw_path: str) -> Optional[str]:
        try:
            self.root.after(0, lambda: self.plot_btn.configure(state="disabled"))
        except Exception:
            pass

        try:
            self._set_status("Converting RAW → mzML (msconvert)...")
            mzml = _convert_raw_to_mzml(
                raw_path,
                msconvert_path=(self.msconvert_path.get().strip() or None),
                ms_level_filter=str(int(self.raw_ms_level.get())) if hasattr(self, "raw_ms_level") else None,
                peakpick=bool(self.msconvert_peakpick.get()),
                compress=bool(self.msconvert_compress.get()),
            )
            self._set_status(f"Conversion complete: {os.path.basename(mzml)}")
            return mzml
        except Exception as exc:
            self._set_status("Conversion failed")
            return f"__CONVERT_ERROR__:{exc}"
        finally:
            try:
                self.root.after(0, lambda: self.plot_btn.configure(state="normal"))
            except Exception:
                pass


# ------------------------------ Main ------------------------------

def main():
    root = TkinterDnD.Tk() if DND_AVAILABLE else tk.Tk()
    try:
        root.tk.call("package", "require", "msgcat")
    except tk.TclError:
        pass
    app = PlotApp(root, themename="darkly")
    root.mainloop()


if __name__ == "__main__":
    main()
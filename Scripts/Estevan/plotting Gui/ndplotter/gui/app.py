from __future__ import annotations
import os
import tkinter as tk
from tkinter import filedialog
from tkinter.colorchooser import askcolor

import ttkbootstrap as tb
from ttkbootstrap.constants import *
from ttkbootstrap.dialogs import Messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Import core plotting functions with new names/locations
from ndplotter.core.chromatogram import plot_chrom_stack
from ndplotter.core.cdms_mass import plot_cdms_mass_offset_stack, plot_single_cdms_mass
from ndplotter.core.mz_spectrum import plot_single_mz_spectrum, plot_mz_spectrum_offset_stack
from ndplotter.core.mwd_chromatogram import plot_mwd_chromatograms
from ndplotter.core.thermo_raw import plot_raw_chrom_and_spectrum

# Optional drag-and-drop
try:
    from tkinterdnd2 import TkinterDnD, DND_FILES

    DND_AVAILABLE = True
except Exception:
    DND_AVAILABLE = False

# Defaults
DEFAULT_CHROM_COLORS = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CDMS_COLORS = ["#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CHROM_STANDARDS = "13.81, 16.82, 18.09, 19.05, 20.79"
DEFAULT_MWD_COLORS = ["#0F0C08", "#3CAEA3"]
DEFAULT_MZ_COLORS = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442"]
THEMES = ["flatly", "cosmo", "lumen", "minty", "sandstone", "superhero", "cyborg", "darkly", "morph"]


def parse_float_list(text: str) -> list[float]:
    text = text.strip()
    if not text:
        return []
    out = []
    for part in text.replace(";", ",").split(","):
        part = part.strip()
        if part:
            try:
                out.append(float(part))
            except ValueError:
                pass
    return out


def parse_str_list(text: str) -> list[str]:
    text = text.strip()
    if not text:
        return []
    return [s.strip() for s in text.replace(";", ",").split(",") if s.strip()]


def filenames_to_labels(paths: list[str]) -> list[str]:
    return [os.path.splitext(os.path.basename(p))[0] for p in paths]


class PlotApp:
    def __init__(self, root: tk.Tk, themename: str = "darkly"):
        self.root = root
        # Ensure msgcat
        try:
            self.root.tk.call('package', 'require', 'msgcat')
        except tk.TclError:
            pass
        # Style
        try:
            self.style = tb.Style(master=self.root, theme=themename)
        except TypeError:
            self.style = tb.Style(themename)
        try:
            self.style.theme_use(themename)
        except Exception:
            pass

        # Build UI containers
        self.header = tb.Frame(self.root, bootstyle="secondary");
        self.header.pack(fill=X)
        self.centerer = tb.Frame(self.root, bootstyle="secondary");
        self.centerer.pack(fill=BOTH, expand=YES)
        self.footer = tb.Frame(self.root, bootstyle="secondary");
        self.footer.pack(fill=X, side=BOTTOM)

        # Centerer grid for centered content
        self.centerer.grid_columnconfigure(0, weight=1)
        self.centerer.grid_columnconfigure(1, weight=0)
        self.centerer.grid_columnconfigure(2, weight=1)
        self.centerer.grid_rowconfigure(0, weight=1)
        self.body = tb.Frame(self.centerer, bootstyle="secondary");
        self.body.grid(row=0, column=1, sticky="n")

        # Build sections
        self._build_header(self.header)
        self._build_body(self.body)
        self._build_footer(self.footer)

        self.files: list[str] = []
        self._set_status("Ready. Drag files into the list or use the pickers, then click Plot.")

    def _build_header(self, parent):
        header = tb.Frame(parent, padding=10, bootstyle="secondary");
        header.pack(fill=X)
        tb.Label(header, text="Plotting GUI: Chromatograms, CD-MS Mass, MWD, m/z Spectra, and Thermo RAW",
                 font=("Arial", 16, "bold")).pack(side=LEFT)

        theme_frame = tb.Frame(header, bootstyle="secondary");
        theme_frame.pack(side=RIGHT)
        tb.Label(theme_frame, text="Theme:", font=("Arial", 10, "bold")).pack(side=LEFT, padx=(0, 6))
        self.theme_var = tk.StringVar(value=self.style.theme_use())
        combo = tb.Combobox(theme_frame, values=THEMES, textvariable=self.theme_var, state="readonly", width=10)
        combo.pack(side=LEFT)
        combo.bind("<<ComboboxSelected>>", self._on_theme_change)

        bar = tb.Frame(parent, bootstyle="secondary");
        bar.pack(fill=X)
        self.status_label = tb.Label(bar, text="Status: Ready", anchor=W);
        self.status_label.pack(side=LEFT, padx=8, pady=4)

    def _on_theme_change(self, event=None):
        self.style.theme_use(self.theme_var.get())

    def _build_body(self, parent):
        paned = tb.Panedwindow(parent, orient=HORIZONTAL, bootstyle="secondary")
        paned.pack(fill=BOTH, expand=YES, padx=10, pady=6)
        left = tb.Frame(paned, padding=10, bootstyle="secondary");
        paned.add(left, weight=1)
        right = tb.Frame(paned, padding=10, bootstyle="secondary");
        paned.add(right, weight=3)
        self.plot_area = right

        nb = tb.Notebook(left, bootstyle="primary");
        nb.pack(fill=BOTH, expand=YES)
        self._build_common_options(left)

        self._build_chrom_tab(nb)
        self._build_cdms_tab(nb)
        self._build_mwd_tab(nb)
        self._build_mz_tab(nb)
        self._build_raw_tab(nb)

        self._build_file_panel(left)

    def _build_footer(self, parent):
        frame = tb.Frame(parent, bootstyle="secondary");
        frame.pack(fill=X, pady=6)
        self.progress = tb.Progressbar(frame, mode="indeterminate", bootstyle="info")
        self.progress.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(frame, text="Plot", bootstyle="primary", command=self._plot).pack(side=RIGHT, padx=6)

    def _build_common_options(self, parent):
        lf = tb.Labelframe(parent, text="Common Options", padding=10, bootstyle="secondary");
        lf.pack(fill=X, pady=6)
        sp = tb.Frame(lf, bootstyle="secondary");
        sp.pack(fill=X, pady=4)
        tb.Label(sp, text="Save as:", width=12).pack(side=LEFT)
        self.save_var = tk.StringVar()
        tb.Entry(sp, textvariable=self.save_var).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(sp, text="Browse", bootstyle="secondary", command=self._browse_save).pack(side=LEFT)

        grid = tb.Frame(lf, bootstyle="secondary");
        grid.pack(fill=X, pady=4)
        tb.Label(grid, text="Font:", width=8).grid(row=0, column=0, sticky=W)
        self.font_name = tk.StringVar(value="Arial");
        tb.Entry(grid, textvariable=self.font_name, width=14).grid(row=0, column=1, sticky=W, padx=4)
        tb.Label(grid, text="Size:", width=6).grid(row=0, column=2, sticky=W)
        self.font_size = tk.IntVar(value=10);
        tb.Entry(grid, textvariable=self.font_size, width=6).grid(row=0, column=3, sticky=W, padx=4)
        tb.Label(grid, text="DPI:", width=6).grid(row=0, column=4, sticky=W)
        self.dpi = tk.IntVar(value=300);
        tb.Entry(grid, textvariable=self.dpi, width=8).grid(row=0, column=5, sticky=W, padx=4)
        tb.Label(grid, text="Figsize (W,H):", width=14).grid(row=0, column=6, sticky=W)
        self.fig_w = tk.DoubleVar(value=10.0);
        tb.Entry(grid, textvariable=self.fig_w, width=6).grid(row=0, column=7, sticky=W, padx=2)
        self.fig_h = tk.DoubleVar(value=6.0);
        tb.Entry(grid, textvariable=self.fig_h, width=6).grid(row=0, column=8, sticky=W, padx=2)
        toggles = tb.Frame(lf, bootstyle="secondary");
        toggles.pack(fill=X, pady=4)
        self.show_plot = tk.BooleanVar(value=True)
        tb.Checkbutton(toggles, text="Show Figure in external window (in addition to preview)",
                       variable=self.show_plot).pack(anchor=W, pady=2)
        self.transparent_export = tk.BooleanVar(value=False)
        tb.Checkbutton(toggles, text="Transparent export background", variable=self.transparent_export).pack(anchor=W,
                                                                                                             pady=2)

    def _build_chrom_tab(self, nb):
        self.chrom_tab = tb.Frame(nb, padding=8, bootstyle="secondary");
        nb.add(self.chrom_tab, text="Chromatograms")
        mode = tb.Frame(self.chrom_tab, bootstyle="secondary");
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.chrom_mode = tk.StringVar(value="single");
        tb.Radiobutton(mode, text="Single Chromatogram", variable=self.chrom_mode, value="single").pack(side=LEFT,
                                                                                                        padx=8)
        tb.Radiobutton(mode, text="Stacked Chromatograms", variable=self.chrom_mode, value="stacked").pack(side=LEFT)

        lf_labels = tb.Labelframe(self.chrom_tab, text="Chromatogram Labels & Colors", padding=8,
                                  bootstyle="secondary");
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Sample Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.chrom_labels = tk.Text(lf_labels, height=2);
        self.chrom_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.chrom_colors = tk.StringVar(value=", ".join(DEFAULT_CHROM_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary");
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.chrom_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info",
                  command=lambda: self._open_palette_editor(self.chrom_colors, "Chromatogram Colors")).pack(side=LEFT)

        lf_axes = tb.Labelframe(self.chrom_tab, text="Chromatogram Axes & Standards", padding=8, bootstyle="secondary");
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary");
        grid.pack(fill=X)
        tb.Label(grid, text="X-lim (mL):").grid(row=0, column=0, sticky=W)
        self.chrom_xmin = tk.DoubleVar(value=10.0);
        tb.Entry(grid, textvariable=self.chrom_xmin, width=8).grid(row=0, column=1, padx=4)
        self.chrom_xmax = tk.DoubleVar(value=20.0);
        tb.Entry(grid, textvariable=self.chrom_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-lim:").grid(row=0, column=3, sticky=W)
        self.chrom_ymin = tk.DoubleVar(value=0.0);
        tb.Entry(grid, textvariable=self.chrom_ymin, width=8).grid(row=0, column=4, padx=4)
        self.chrom_ymax = tk.DoubleVar(value=1.0);
        tb.Entry(grid, textvariable=self.chrom_ymax, width=8).grid(row=0, column=5, padx=4)
        self.chrom_normalize = tk.BooleanVar(value=True);
        tb.Checkbutton(lf_axes, text="Normalize traces to max", variable=self.chrom_normalize).pack(anchor=W, pady=4)

        std_row = tb.Frame(lf_axes, bootstyle="secondary");
        std_row.pack(fill=X, pady=4)
        self.chrom_show_standards = tk.BooleanVar(value=True)
        tb.Checkbutton(std_row, text="Show standard volumes (top ticks)", variable=self.chrom_show_standards).pack(
            side=LEFT)
        tb.Label(std_row, text="Standard Volumes (mL):").pack(side=LEFT, padx=8)
        self.chrom_standards = tk.StringVar(value=DEFAULT_CHROM_STANDARDS)
        tb.Entry(std_row, textvariable=self.chrom_standards, width=28).pack(side=LEFT, padx=6)

    def _build_cdms_tab(self, nb):
        self.cdms_tab = tb.Frame(nb, padding=8, bootstyle="secondary");
        nb.add(self.cdms_tab, text="CD-MS Mass")
        mode = tb.Frame(self.cdms_tab, bootstyle="secondary");
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.cdms_mode = tk.StringVar(value="single");
        tb.Radiobutton(mode, text="Single CD-MS Mass", variable=self.cdms_mode, value="single").pack(side=LEFT, padx=8)
        tb.Radiobutton(mode, text="Stacked CD-MS Mass", variable=self.cdms_mode, value="stacked").pack(side=LEFT)

        lf_labels = tb.Labelframe(self.cdms_tab, text="CD-MS Mass Labels & Colors", padding=8, bootstyle="secondary");
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.cdms_labels = tk.Text(lf_labels, height=2);
        self.cdms_labels.pack(fill=X, pady=4)
        self.cdms_two_legends = tk.BooleanVar(value=False)
        tb.Checkbutton(lf_labels, text="Include second legend near peak", variable=self.cdms_two_legends).pack(anchor=W,
                                                                                                               pady=2)
        tb.Label(lf_labels, text="Second Labels (one per file; used when second legend is enabled)").pack(anchor=W)
        self.cdms_second_labels = tk.Text(lf_labels, height=2);
        self.cdms_second_labels.pack(fill=X, pady=4)
        self.cdms_second_labels.config(state=tk.DISABLED)
        self.cdms_two_legends.trace_add("write", lambda *args: self._toggle_second_labels_input())
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.cdms_colors = tk.StringVar(value=", ".join(DEFAULT_CDMS_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary");
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.cdms_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info",
                  command=lambda: self._open_palette_editor(self.cdms_colors, "CD-MS Colors")).pack(side=LEFT)

        lf_axes = tb.Labelframe(self.cdms_tab, text="CD-MS Mass Axes & Layout", padding=8, bootstyle="secondary");
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary");
        grid.pack(fill=X)
        tb.Label(grid, text="X-range (kDa):").grid(row=0, column=0, sticky=W)
        self.cdms_xmin = tk.DoubleVar(value=0.0);
        tb.Entry(grid, textvariable=self.cdms_xmin, width=8).grid(row=0, column=1, padx=4)
        self.cdms_xmax = tk.DoubleVar(value=500.0);
        tb.Entry(grid, textvariable=self.cdms_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-offset step (stacked):").grid(row=0, column=3, sticky=W)
        self.cdms_y_offset = tk.DoubleVar(value=120.0);
        tb.Entry(grid, textvariable=self.cdms_y_offset, width=8).grid(row=0, column=4, padx=4)
        self.cdms_normalize = tk.BooleanVar(value=True);
        tb.Checkbutton(lf_axes, text="Normalize to max (percent)", variable=self.cdms_normalize).pack(anchor=W, pady=4)

    def _build_mwd_tab(self, nb):
        self.mwd_tab = tb.Frame(nb, padding=8, bootstyle="secondary");
        nb.add(self.mwd_tab, text="MWD Chromatograms")
        lf_files = tb.Labelframe(self.mwd_tab, text="MWD Files (drag & drop or browse)", padding=8,
                                 bootstyle="secondary");
        lf_files.pack(fill=X, pady=4)
        row1 = tb.Frame(lf_files, bootstyle="secondary");
        row1.pack(fill=X, pady=2)
        tb.Label(row1, text="mAU file:").pack(side=LEFT)
        self.mwd_mau_path = tk.StringVar()
        self.mwd_mau_entry = tb.Entry(row1, textvariable=self.mwd_mau_path)
        self.mwd_mau_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row1, text="Browse", bootstyle="secondary", command=self._browse_mwd_mau).pack(side=LEFT)
        if DND_AVAILABLE:
            self.mwd_mau_entry.drop_target_register(DND_FILES)
            self.mwd_mau_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.mwd_mau_path))
        row2 = tb.Frame(lf_files, bootstyle="secondary");
        row2.pack(fill=X, pady=2)
        tb.Label(row2, text="mV file:").pack(side=LEFT)
        self.mwd_mv_path = tk.StringVar()
        self.mwd_mv_entry = tb.Entry(row2, textvariable=self.mwd_mv_path)
        self.mwd_mv_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row2, text="Browse", bootstyle="secondary", command=self._browse_mwd_mv).pack(side=LEFT)
        if DND_AVAILABLE:
            self.mwd_mv_entry.drop_target_register(DND_FILES)
            self.mwd_mv_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.mwd_mv_path))
        lf_labels = tb.Labelframe(self.mwd_tab, text="Labels & Colors", padding=8, bootstyle="secondary");
        lf_labels.pack(fill=X, pady=4)
        row3 = tb.Frame(lf_labels, bootstyle="secondary");
        row3.pack(fill=X, pady=2)
        tb.Label(row3, text="mAU label:").pack(side=LEFT)
        self.mwd_mau_label = tk.StringVar(value="AqpZ (mAU)")
        tb.Entry(row3, textvariable=self.mwd_mau_label, width=24).pack(side=LEFT, padx=6)
        tb.Label(row3, text="mV label:").pack(side=LEFT)
        self.mwd_mv_label = tk.StringVar(value="AqpZ (mV)")
        tb.Entry(row3, textvariable=self.mwd_mv_label, width=24).pack(side=LEFT, padx=6)
        row4 = tb.Frame(lf_labels, bootstyle="secondary");
        row4.pack(fill=X, pady=2)
        tb.Label(row4, text="Colors (mAU, mV hex):").pack(side=LEFT)
        self.mwd_colors = tk.StringVar(value=", ".join(DEFAULT_MWD_COLORS))
        tb.Entry(row4, textvariable=self.mwd_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row4, text="Pick mAU", bootstyle="info", command=lambda: self._pick_mwd_color(index=0)).pack(
            side=LEFT, padx=4)
        tb.Button(row4, text="Pick mV", bootstyle="info", command=lambda: self._pick_mwd_color(index=1)).pack(side=LEFT,
                                                                                                              padx=4)
        lf_axes = tb.Labelframe(self.mwd_tab, text="Axes & Options", padding=8, bootstyle="secondary");
        lf_axes.pack(fill=X, pady=4)
        grid = tb.Frame(lf_axes, bootstyle="secondary");
        grid.pack(fill=X)
        tb.Label(grid, text="X-lim (mL):").grid(row=0, column=0, sticky=W)
        self.mwd_xmin = tk.DoubleVar(value=0.0);
        tb.Entry(grid, textvariable=self.mwd_xmin, width=8).grid(row=0, column=1, padx=4)
        self.mwd_xmax = tk.DoubleVar(value=26.0);
        tb.Entry(grid, textvariable=self.mwd_xmax, width=8).grid(row=0, column=2, padx=4)
        self.mwd_normalize = tk.BooleanVar(value=True);
        tb.Checkbutton(lf_axes, text="Normalize each trace to max", variable=self.mwd_normalize).pack(anchor=W, pady=4)

    def _build_mz_tab(self, nb):
        self.mz_tab = tb.Frame(nb, padding=8, bootstyle="secondary");
        nb.add(self.mz_tab, text="m/z Spectra")
        mode = tb.Frame(self.mz_tab, bootstyle="secondary");
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.mz_mode = tk.StringVar(value="single");
        tb.Radiobutton(mode, text="Single m/z Spectrum", variable=self.mz_mode, value="single").pack(side=LEFT, padx=8)
        tb.Radiobutton(mode, text="Stacked m/z Spectra", variable=self.mz_mode, value="stacked").pack(side=LEFT)
        lf_labels = tb.Labelframe(self.mz_tab, text="m/z Labels & Colors", padding=8, bootstyle="secondary");
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.mz_labels = tk.Text(lf_labels, height=2);
        self.mz_labels.pack(fill=X, pady=4)
        self.mz_two_legends = tk.BooleanVar(value=False)
        tb.Checkbutton(lf_labels, text="Include second legend near peak", variable=self.mz_two_legends).pack(anchor=W,
                                                                                                             pady=2)
        tb.Label(lf_labels, text="Second Labels (one per file; used when second legend is enabled)").pack(anchor=W)
        self.mz_second_labels = tk.Text(lf_labels, height=2);
        self.mz_second_labels.pack(fill=X, pady=4)
        self.mz_second_labels.config(state=tk.DISABLED)
        self.mz_two_legends.trace_add("write", lambda *args: self._toggle_mz_second_labels())
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.mz_colors = tk.StringVar(value=", ".join(DEFAULT_MZ_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary");
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.mz_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info",
                  command=lambda: self._open_palette_editor(self.mz_colors, "m/z Colors")).pack(side=LEFT)
        lf_axes = tb.Labelframe(self.mz_tab, text="m/z Axes & Layout", padding=8, bootstyle="secondary");
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary");
        grid.pack(fill=X)
        tb.Label(grid, text="X-range (m/z):").grid(row=0, column=0, sticky=W)
        self.mz_xmin = tk.DoubleVar(value=100.0);
        tb.Entry(grid, textvariable=self.mz_xmin, width=8).grid(row=0, column=1, padx=4)
        self.mz_xmax = tk.DoubleVar(value=2000.0);
        tb.Entry(grid, textvariable=self.mz_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-offset step (stacked):").grid(row=0, column=3, sticky=W)
        self.mz_y_offset = tk.DoubleVar(value=120.0);
        tb.Entry(grid, textvariable=self.mz_y_offset, width=8).grid(row=0, column=4, padx=4)
        self.mz_normalize = tk.BooleanVar(value=True);
        tb.Checkbutton(lf_axes, text="Normalize to max (percent)", variable=self.mz_normalize).pack(anchor=W, pady=4)

    def _build_raw_tab(self, nb):
        self.raw_tab = tb.Frame(nb, padding=8, bootstyle="secondary");
        nb.add(self.raw_tab, text="Thermo RAW")
        lf_file = tb.Labelframe(self.raw_tab, text="RAW/mzML File", padding=8, bootstyle="secondary");
        lf_file.pack(fill=X, pady=4)
        row = tb.Frame(lf_file, bootstyle="secondary");
        row.pack(fill=X, pady=2)
        tb.Label(row, text="Path:").pack(side=LEFT)
        self.raw_path = tk.StringVar()
        self.raw_entry = tb.Entry(row, textvariable=self.raw_path)
        self.raw_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row, text="Browse", bootstyle="secondary", command=self._browse_raw_file).pack(side=LEFT)
        if DND_AVAILABLE:
            self.raw_entry.drop_target_register(DND_FILES)
            self.raw_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.raw_path))
        lf_opts = tb.Labelframe(self.raw_tab, text="Options", padding=8, bootstyle="secondary");
        lf_opts.pack(fill=X, pady=4)
        grid = tb.Frame(lf_opts, bootstyle="secondary");
        grid.pack(fill=X)
        tb.Label(grid, text="MS level:").grid(row=0, column=0, sticky=W)
        self.raw_ms_level = tk.IntVar(value=1);
        tb.Entry(grid, textvariable=self.raw_ms_level, width=6).grid(row=0, column=1, padx=6)
        tb.Label(grid, text="Max scans:").grid(row=0, column=2, sticky=W)
        self.raw_max_scans = tk.IntVar(value=200);
        tb.Entry(grid, textvariable=self.raw_max_scans, width=8).grid(row=0, column=3, padx=6)
        tb.Label(grid, text="Spectrum:").grid(row=1, column=0, sticky=W)
        self.raw_mode = tk.StringVar(value="apex")
        tb.Radiobutton(grid, text="Apex", variable=self.raw_mode, value="apex").grid(row=1, column=1, sticky=W)
        tb.Radiobutton(grid, text="Index", variable=self.raw_mode, value="index").grid(row=1, column=2, sticky=W)
        tb.Label(grid, text="Index:").grid(row=1, column=3, sticky=W)
        self.raw_index = tk.IntVar(value=0);
        self.raw_index_entry = tb.Entry(grid, textvariable=self.raw_index, width=6)
        self.raw_index_entry.grid(row=1, column=4, padx=6)

        def _toggle_raw_index(*args):
            self.raw_index_entry.config(state=(tk.NORMAL if self.raw_mode.get() == "index" else tk.DISABLED))

        self.raw_mode.trace_add("write", _toggle_raw_index);
        _toggle_raw_index()
        self.raw_norm = tk.BooleanVar(value=True);
        tb.Checkbutton(lf_opts, text="Normalize spectrum to max (percent)", variable=self.raw_norm).pack(anchor=W,
                                                                                                         pady=4)
        lf_ms = tb.Labelframe(self.raw_tab, text="msconvert (ProteoWizard)", padding=8, bootstyle="secondary");
        lf_ms.pack(fill=X, pady=4)
        row2 = tb.Frame(lf_ms, bootstyle="secondary");
        row2.pack(fill=X, pady=2)
        tb.Label(row2, text="msconvert path (optional):").pack(side=LEFT)
        self.msconvert_path = tk.StringVar()
        self.msconvert_entry = tb.Entry(row2, textvariable=self.msconvert_path)
        self.msconvert_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row2, text="Browse", bootstyle="secondary", command=self._browse_msconvert).pack(side=LEFT)
        self.raw_prefer_native = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_ms, text="Prefer native Thermo MSFileReader on Windows (fallback to msconvert)",
                       variable=self.raw_prefer_native).pack(anchor=W, pady=4)

    def _build_file_panel(self, parent):
        lf = tb.Labelframe(parent, text="Files (Chromatograms, CD-MS Mass, MWD, m/z)", padding=8,
                           bootstyle="secondary");
        lf.pack(fill=BOTH, expand=YES, pady=6)
        btns = tb.Frame(lf, bootstyle="secondary");
        btns.pack(fill=X, pady=4)
        tb.Button(btns, text="Add Files...", bootstyle="success", command=self._add_files).pack(side=LEFT)
        tb.Button(btns, text="Remove Selected", bootstyle="warning", command=self._remove_selected).pack(side=LEFT,
                                                                                                         padx=6)
        tb.Button(btns, text="Clear", bootstyle="danger", command=self._clear_files).pack(side=LEFT)
        self.file_tree = tb.Treeview(lf, columns=("idx", "path"), show="headings", height=8, bootstyle="info")
        self.file_tree.heading("idx", text="#");
        self.file_tree.column("idx", width=40, anchor=CENTER)
        self.file_tree.heading("path", text="File Path");
        self.file_tree.column("path", width=700, anchor=W)
        self.file_tree.pack(fill=BOTH, expand=YES, pady=6)
        if DND_AVAILABLE:
            self.file_tree.drop_target_register(DND_FILES)
            self.file_tree.dnd_bind("<<Drop>>", self._on_files_drop)

    def _browse_save(self):
        path = filedialog.asksaveasfilename(title="Save figure as...", defaultextension=".png",
                                            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if path:
            self.save_var.set(path)

    def _add_files(self):
        tab = self._active_tab()
        if tab == "Chromatograms":
            paths = filedialog.askopenfilenames(title="Select Chromatogram Excel files",
                                                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
        elif tab == "CD-MS Mass":
            paths = filedialog.askopenfilenames(title="Select CD-MS mass text files",
                                                filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        elif tab == "m/z Spectra":
            paths = filedialog.askopenfilenames(title="Select m/z spectrum text files",
                                                filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        else:
            paths = filedialog.askopenfilenames(title="Select mAU/mV Excel files (first=mAU, second=mV)",
                                                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
            if paths:
                self.files = list(paths)[:2];
                self._refresh_files();
                return
        if paths:
            self.files.extend(paths);
            self._refresh_files()

    def _remove_selected(self):
        selected = self.file_tree.selection()
        if not selected: return
        idxs = sorted([int(self.file_tree.item(i, "values")[0]) for i in selected], reverse=True)
        for idx in idxs:
            if 0 <= idx < len(self.files):
                self.files.pop(idx)
        self._refresh_files()

    def _clear_files(self):
        self.files = [];
        self._refresh_files()

    def _refresh_files(self):
        for i in self.file_tree.get_children():
            self.file_tree.delete(i)
        for idx, p in enumerate(self.files):
            self.file_tree.insert("", "end", values=(idx, p))

    def _active_tab(self) -> str:
        if self.chrom_tab.winfo_ismapped(): return "Chromatograms"
        if self.cdms_tab.winfo_ismapped():  return "CD-MS Mass"
        if self.mwd_tab.winfo_ismapped():   return "MWD Chromatograms"
        if self.mz_tab.winfo_ismapped():    return "m/z Spectra"
        return "Thermo RAW"

    def _on_files_drop(self, event):
        try:
            paths = list(self.root.splitlist(event.data))
        except Exception:
            paths = event.data.split()
        if not paths: return
        tab = self._active_tab()
        if tab == "MWD Chromatograms":
            self.files = paths[:2];
            self._refresh_files()
            self._set_status(f"Set MWD files via drag-and-drop ({len(self.files)} file(s)).")
        else:
            self.files.extend(paths);
            self._refresh_files()
            self._set_status(f"Added {len(paths)} file(s) via drag-and-drop.")

    def _on_single_file_drop(self, event, target_var: tk.StringVar):
        try:
            paths = list(self.root.splitlist(event.data))
        except Exception:
            paths = event.data.split()
        if paths:
            target_var.set(paths[0]);
            self._set_status("Set file via drag-and-drop.")

    def _set_status(self, text: str):
        self.status_label.configure(text=f"Status: {text}")
        self.root.update_idletasks()

    def _open_palette_editor(self, target_var: tk.StringVar, title: str):
        top = tk.Toplevel(self.root);
        top.title(title);
        top.geometry("420x300")
        colors = parse_str_list(target_var.get())
        lb_frame = tb.Frame(top, padding=6, bootstyle="secondary");
        lb_frame.pack(fill=BOTH, expand=YES)
        lb = tk.Listbox(lb_frame, height=10);
        lb.pack(fill=BOTH, expand=YES)

        def refresh_list():
            lb.delete(0, tk.END)
            for c in colors: lb.insert(tk.END, c)

        refresh_list()
        ctrl = tb.Frame(top, padding=6, bootstyle="secondary");
        ctrl.pack(fill=X)
        tb.Button(ctrl, text="Add…", bootstyle="info", command=lambda: self._palette_add(colors, refresh_list)).pack(
            side=LEFT, padx=4)
        tb.Button(ctrl, text="Remove", bootstyle="warning",
                  command=lambda: self._palette_remove(colors, lb, refresh_list)).pack(side=LEFT, padx=4)
        tb.Button(ctrl, text="Clear", bootstyle="danger",
                  command=lambda: self._palette_clear(colors, refresh_list)).pack(side=LEFT, padx=4)
        tb.Button(ctrl, text="Apply", bootstyle="success",
                  command=lambda: self._palette_apply(colors, target_var, top)).pack(side=RIGHT, padx=4)
        tb.Button(ctrl, text="Cancel", bootstyle="secondary", command=top.destroy).pack(side=RIGHT, padx=4)

    def _palette_add(self, colors: list[str], refresh_list):
        start_color = colors[-1] if colors else "#00A9A5"
        c = askcolor(color=start_color, title="Pick a color")[1]
        if c: colors.append(c); refresh_list()

    def _palette_remove(self, colors: list[str], lb: tk.Listbox, refresh_list):
        sel = lb.curselection()
        if sel:
            idx = sel[0];
            colors.pop(idx);
            refresh_list()

    def _palette_clear(self, colors: list[str], refresh_list):
        colors.clear();
        refresh_list()

    def _palette_apply(self, colors: list[str], target_var: tk.StringVar, top: tk.Toplevel):
        target_var.set(", ".join(colors));
        top.destroy()

    def _pick_mwd_color(self, index: int):
        colors = parse_str_list(self.mwd_colors.get()) or DEFAULT_MWD_COLORS.copy()
        while len(colors) < 2: colors.append("#000000")
        start_color = colors[index]
        c = askcolor(color=start_color, title="Pick a color")[1]
        if c: colors[index] = c; self.mwd_colors.set(", ".join(colors))

    def _browse_msconvert(self):
        path = filedialog.askopenfilename(title="Select msconvert.exe",
                                          filetypes=[("Executable", "*.exe;*"), ("All files", "*.*")])
        if path:
            self.msconvert_path.set(path)
            self._set_status("msconvert path set.")

    def _browse_raw_file(self):
        path = filedialog.askopenfilename(title="Select Thermo RAW or mzML",
                                          filetypes=[("RAW/mzML", "*.raw *.mzML"), ("All files", "*.*")])
        if path:
            self.raw_path.set(path)
            self._set_status("RAW/mzML file selected.")

    def _browse_mwd_mau(self):
        path = filedialog.askopenfilename(title="Select mAU Excel file",
                                          filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
        if path: self.mwd_mau_path.set(path); self._set_status("mAU file selected.")

    def _browse_mwd_mv(self):
        path = filedialog.askopenfilename(title="Select mV Excel file",
                                          filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
        if path: self.mwd_mv_path.set(path); self._set_status("mV file selected.")

    def _plot(self):
        try:
            self.progress.start(10);
            self._set_status("Plotting...")
            tab = self._active_tab()
            # Common config
            font_name = self.font_name.get();
            font_size = int(self.font_size.get())
            dpi = int(self.dpi.get());
            figsize = (float(self.fig_w.get()), float(self.fig_h.get()))
            save = self.save_var.get() or None;
            show = bool(self.show_plot.get())

            if tab == "Chromatograms":
                raw_labels = self.chrom_labels.get("1.0", "end").strip()
                labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
                colors = parse_str_list(self.chrom_colors.get()) or DEFAULT_CHROM_COLORS
                std_vols = parse_float_list(self.chrom_standards.get()) if bool(
                    self.chrom_show_standards.get()) else None
                normalize = bool(self.chrom_normalize.get())
                if self.chrom_mode.get() == "single":
                    if len(self.files) != 1:
                        Messagebox.show_error("Select exactly one Excel file for Single Chromatogram.",
                                              title="Selection error");
                        return
                    fig, axes = plot_chrom_stack(
                        files=self.files[:1], sample_labels=labels[:1], colors=colors[:1],
                        standard_volumes=std_vols, normalize=normalize,
                        font_name=font_name, font_size=font_size,
                        xlim=(float(self.chrom_xmin.get()), float(self.chrom_xmax.get())),
                        ylim=(float(self.chrom_ymin.get()), float(self.chrom_ymax.get())),
                        figsize=figsize, dpi=dpi, save=None, show=show,
                    )
                else:
                    if len(labels) != len(self.files):
                        Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch");
                        return
                    fig, axes = plot_chrom_stack(
                        files=self.files, sample_labels=labels, colors=colors,
                        standard_volumes=std_vols, normalize=normalize,
                        font_name=font_name, font_size=font_size,
                        xlim=(float(self.chrom_xmin.get()), float(self.chrom_xmax.get())),
                        ylim=(float(self.chrom_ymin.get()), float(self.chrom_ymax.get())),
                        figsize=figsize, dpi=dpi, save=None, show=show,
                    )
                self._render_figure(fig);
                if save: self._export_figure(fig, save, dpi, bool(self.transparent_export.get()))
                self._set_status("Chromatogram plot complete.")

            elif tab == "CD-MS Mass":
                raw_labels = self.cdms_labels.get("1.0", "end").strip()
                labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
                use_second = bool(self.cdms_two_legends.get())
                second_raw = self.cdms_second_labels.get("1.0", "end").strip()
                second_labels = parse_str_list(second_raw) if (use_second and second_raw) else [None] * len(self.files)
                colors = parse_str_list(self.cdms_colors.get()) or DEFAULT_CDMS_COLORS
                x_range = (float(self.cdms_xmin.get()), float(self.cdms_xmax.get()))
                normalize = bool(self.cdms_normalize.get())
                if self.cdms_mode.get() == "single":
                    if len(self.files) != 1:
                        Messagebox.show_error("Select exactly one text file for Single CD-MS Mass.",
                                              title="Selection error");
                        return
                    second_label_single = second_labels[0] if use_second else None
                    fig, ax = plot_single_cdms_mass(
                        dataset=self.files[0], label=labels[0], second_label=second_label_single,
                        x_range=x_range, color=colors[0], normalize=normalize,
                        font_name=font_name, font_size=font_size,
                        figsize=figsize, dpi=dpi, save=None, show=show,
                    )
                    self._render_figure(fig);
                    if save: self._export_figure(fig, save, dpi, bool(self.transparent_export.get()))
                    self._set_status("CD-MS mass (single) plot complete.")
                else:
                    if len(labels) != len(self.files):
                        Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch");
                        return
                    if use_second and len(second_labels) != len(self.files):
                        Messagebox.show_error("Second labels must match number of files when second legend is enabled.",
                                              title="Second label mismatch");
                        return
                    y_ticks = ("0", "%", "100") if normalize else ("", "", "")
                    fig, axes = plot_cdms_mass_offset_stack(
                        datasets=self.files, labels=labels,
                        second_labels=second_labels if use_second else [None] * len(self.files),
                        x_range=x_range, colors=colors, y_offset_step=float(self.cdms_y_offset.get()),
                        y_tick_labels=y_ticks, normalize=normalize,
                        font_name=font_name, font_size=font_size,
                        figsize=figsize, dpi=dpi, save=None, show=show,
                    )
                    self._render_figure(fig);
                    if save: self._export_figure(fig, save, dpi, bool(self.transparent_export.get()))
                    self._set_status("CD-MS mass (stacked) plot complete.")

                    fig, (ax1, ax2) = plot_mwd_chromatograms(

                        mau_file=path_mau,

                        mv_file=path_mv,

                        mau_label=self.mwd_mau_label.get(),

                        mv_label=self.mwd_mv_label.get(),

                        colors=colors[:2],

                        xlim=xlim,

                        normalize=normalize,

                        font_name=font_name,

                        font_size=font_size,

                        figsize=figsize,

                        dpi=dpi,

                        save=None,

                        show=show,

                    )

                    self._render_figure(fig)

                    if save:
                        self._export_figure(fig, save, dpi, bool(self.transparent_export.get()))

                    self._set_status("MWD Chromatograms plot complete.")

                elif tab == "m/z Spectra":

                raw_labels = self.mz_labels.get("1.0", "end").strip()

                labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)

                use_second = bool(self.mz_two_legends.get())

                second_raw = self.mz_second_labels.get("1.0", "end").strip()

                second_labels = parse_str_list(second_raw) if (use_second and second_raw) else [None] * len(self.files)

                colors = parse_str_list(self.mz_colors.get()) or DEFAULT_MZ_COLORS

                x_range = (float(self.mz_xmin.get()), float(self.mz
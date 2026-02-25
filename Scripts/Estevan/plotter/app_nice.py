# app_nice.py (Chunk 1 of 2)
from __future__ import annotations
import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from tkinter.colorchooser import askcolor
from typing import List, Sequence, Tuple, Optional
from dataclasses import dataclass

# Drag-and-drop root (optional)
try:
    from tkinterdnd2 import TkinterDnD, DND_FILES
    DND_AVAILABLE = True
except Exception:
    DND_AVAILABLE = False

import ttkbootstrap as tb
from ttkbootstrap.constants import *
from ttkbootstrap.dialogs import Messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Preview DPI (moderate size on-screen); export uses CommonConfig.dpi
DISPLAY_DPI = 110

# Plotting utilities (ensure these files exist in the same folder)
from chromatogram_plot_utils import plot_chrom_stack, read_chrom_file
from cdms_mass_plot_utils import (
    plot_cdms_mass_offset_stack, plot_single_cdms_mass, read_cdms_mass_file
)
from mwd_chromatogram_utils import plot_mwd_chromatograms
from mz_spectrum_plot_utils import (
    plot_single_mz_spectrum, plot_mz_spectrum_offset_stack, read_mz_spectrum_file
)
from thermo_raw_utils import plot_raw_chrom_and_spectrum

# Defaults
DEFAULT_CHROM_COLORS = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CDMS_COLORS  = ["#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CHROM_STANDARDS = "13.81, 16.82, 18.09, 19.05, 20.79"
DEFAULT_MWD_COLORS   = ["#0F0C08", "#3CAEA3"]
DEFAULT_MZ_COLORS    = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442"]

import subprocess
import configparser

# Popular themes used for fallback selection
POPULAR_THEMES_RANKED = [
    # Light
    "litera", "flatly", "cosmo", "minty", "yeti", "sandstone", "journal", "united",
    "cerulean", "lux", "quartz", "zephyr",
    # Dark
    "darkly", "morph", "vapor", "cyborg", "superhero",
]

def is_system_dark_mode() -> bool:
    """
    Best-effort cross-platform detection of user dark mode preference.
    Returns True for dark mode, False for light/unknown.
    """
    # Windows 10/11: HKCU\Software\Microsoft\Windows\CurrentVersion\Themes\Personalize\AppsUseLightTheme
    try:
        import winreg
        with winreg.OpenKey(
            winreg.HKEY_CURRENT_USER,
            r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize"
        ) as key:
            val, _ = winreg.QueryValueEx(key, "AppsUseLightTheme")
            return int(val) == 0  # 0 = dark, 1 = light
    except Exception:
        pass

    # macOS: 'defaults read -g AppleInterfaceStyle' returns 'Dark' when dark mode is active
    try:
        proc = subprocess.run(
            ["defaults", "read", "-g", "AppleInterfaceStyle"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if proc.returncode == 0:
            return proc.stdout.strip().lower() == "dark" or proc.stdout.strip() != ""
    except Exception:
        pass

    # GNOME: org.gnome.desktop.interface color-scheme: 'prefer-dark' or 'prefer-light'
    try:
        proc = subprocess.run(
            ["gsettings", "get", "org.gnome.desktop.interface", "color-scheme"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if proc.returncode == 0:
            val = proc.stdout.strip().strip("'").lower()
            if "prefer-dark" in val:
                return True
            if "prefer-light" in val:
                return False
    except Exception:
        pass

    # GNOME fallback: check gtk-theme contains 'dark'
    try:
        proc = subprocess.run(
            ["gsettings", "get", "org.gnome.desktop.interface", "gtk-theme"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if proc.returncode == 0:
            theme = proc.stdout.strip().strip("'").lower()
            if "dark" in theme:
                return True
    except Exception:
        pass

    # KDE: ~/.config/kdeglobals [General] ColorScheme, often contains 'Dark'
    try:
        import os
        kdeglobals = os.path.expanduser("~/.config/kdeglobals")
        if os.path.isfile(kdeglobals):
            cp = configparser.ConfigParser()
            cp.read(kdeglobals)
            scheme = cp.get("General", "ColorScheme", fallback="")
            if "dark" in scheme.lower():
                return True
    except Exception:
        pass

    return False

def choose_default_theme() -> str:
    """
    Choose litera for light or darkly for dark; if unavailable, pick a close modern alternative.
    Ensures the theme exists in this ttkbootstrap installation.
    """
    try:
        tmp_style = tb.Style()
        available = set(tmp_style.theme_names())
    except Exception:
        # Fallback: assume all popular themes may be available
        available = set(POPULAR_THEMES_RANKED)

    system_dark = is_system_dark_mode()
    preferred = "darkly" if system_dark else "litera"

    if preferred not in available:
        fallbacks = (
            ["morph", "vapor", "cyborg", "superhero", "flatly"] if system_dark
            else ["flatly", "cosmo", "minty", "yeti", "sandstone", "journal"]
        )
        for t in fallbacks:
            if t in available:
                preferred = t
                break
        else:
            # Last-resort: pick first available known theme or the style's default
            preferred = next(iter(available)) if available else "flatly"

    return preferred

def parse_float_list(text: str) -> List[float]:
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

def parse_str_list(text: str) -> List[str]:
    text = text.strip()
    if not text:
        return []
    return [s.strip() for s in text.replace(";", ",").split(",") if s.strip()]

def filenames_to_labels(paths: Sequence[str]) -> List[str]:
    return [os.path.splitext(os.path.basename(p))[0] for p in paths]

@dataclass
class CommonConfig:
    save: Optional[str]
    font_name: str
    font_size: int
    dpi: int
    figsize: Tuple[float, float]
    transparent: bool
    show: bool

class PlotApp:
    def __init__(self, root: tk.Tk, themename: str = "flatly"):
        self.root = root
        try:
            self.root.tk.call('package', 'require', 'msgcat')
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

        self.body = tb.Frame(self.root, bootstyle="secondary")
        self.body.pack(fill=BOTH, expand=YES)

        self.footer = tb.Frame(self.root, bootstyle="secondary")
        self.footer.pack(fill=X, side=BOTTOM)

        self._build_header(self.header)
        self._build_body(self.body)
        self._build_footer(self.footer)

        self._center_on_screen()

        self.files: List[str] = []
        self._set_status("Ready. Drag files into the list or use the pickers, then click Plot.")

        # Keyboard shortcuts
        self.root.bind("<F5>", lambda e: self._plot())
        self.root.bind("<Control-Return>", lambda e: self._plot())
        self.root.bind("<Control-KP_Enter>", lambda e: self._plot())

    # ---------- Header ----------
    def _build_header(self, parent):
        header = tb.Frame(parent, padding=10, bootstyle="secondary")
        header.pack(fill=X)
        tb.Label(header, text="Plotting GUI: Chromatograms, CD-MS Mass, MWD, m/z Spectra, and Thermo RAW",
                 font=("Arial", 16, "bold")).pack(side=LEFT)

        theme_frame = tb.Frame(header, bootstyle="secondary")
        theme_frame.pack(side=RIGHT)
        tb.Label(theme_frame, text="Theme:", font=("Arial", 10, "bold")).pack(side=LEFT, padx=(0, 6))

        # Detect available themes and keep only the popular ones (preserving our popularity order)
        try:
            available_set = set(self.style.theme_names())
        except Exception:
            available_set = set(POPULAR_THEMES_RANKED)
        themes_available = [t for t in POPULAR_THEMES_RANKED if t in available_set]

        # Default intelligently: keep current, otherwise best popular theme available
        current_theme = self.style.theme_use()
        if current_theme not in themes_available and themes_available:
            # Prefer a modern dark default if present, else the first light theme
            preferred_defaults = ["darkly", "morph", "flatly", "litera"]
            for pd in preferred_defaults:
                if pd in themes_available:
                    current_theme = pd
                    break
            else:
                current_theme = themes_available[0]

        self.theme_var = tk.StringVar(value=current_theme)
        combo = tb.Combobox(
            theme_frame,
            values=themes_available,
            textvariable=self.theme_var,
            state="readonly",
            width=14
        )
        combo.pack(side=LEFT)
        combo.bind("<<ComboboxSelected>>", self._on_theme_change)

        bar = tb.Frame(parent, bootstyle="secondary")
        bar.pack(fill=X)
        self.status_label = tb.Label(bar, text="Status: Ready", anchor=W)
        self.status_label.pack(side=LEFT, padx=8, pady=4)

        # Plot button: create, then apply theme-aware sizing
        self.plot_btn = tb.Button(bar, text="Plot", bootstyle="success", command=self._plot)
        self._apply_plot_button_sizing()
        self.plot_btn.pack(side=RIGHT, padx=8, pady=4)

    def _on_theme_change(self, event=None):
        self.style.theme_use(self.theme_var.get())
        try:
            self.root.configure(bg=self.style.colors.bg)
        except Exception:
            pass
        # Re-apply plot button sizing because a theme change can adjust widget metrics
        self._apply_plot_button_sizing()

    def _apply_plot_button_sizing(self):
        """
        Make the Plot button bigger in a theme-aware, DPI-aware way without
        affecting other buttons or losing the 'success' bootstyle colors.
        """
        if not hasattr(self, "plot_btn") or self.plot_btn is None:
            return

        # Tk scaling accounts for DPI / OS scaling; ensures consistent physical size
        try:
            scale = float(self.root.tk.call("tk", "scaling"))
        except Exception:
            scale = 1.0

        # Character width ensures a consistent minimum width across themes
        # (ttk Button 'width' is in characters, not pixels)
        self.plot_btn.configure(width=12)

        # Padding enlarges the clickable area; scale it to match theme density
        # Adjust these bases if you’d like it even bigger/smaller
        base_pad_x = 14
        base_pad_y = 10

        pad_x = int(base_pad_x * scale)
        pad_y = int(base_pad_y * scale)

        # Per-widget padding avoids modifying global styles or other success buttons
        try:
            self.plot_btn.configure(padding=(pad_x, pad_y))
        except Exception:
            # Very old Tk builds may not accept widget padding; gracefully ignore
            pass

    # ---------- Scrollable left helper ----------
    def _build_scrollable_left(self, paned):
        """
        Create a scrollable left pane and return (container_frame, inner_content_frame).
        Build all left-side controls inside 'inner_content_frame'.
        """
        container = tb.Frame(paned, bootstyle="secondary")

        # Canvas + vertical scrollbar
        canvas = tk.Canvas(container, borderwidth=0, highlightthickness=0)
        vsb = tb.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side=RIGHT, fill=Y)
        canvas.pack(side=LEFT, fill=BOTH, expand=YES)

        # Inner frame holds your actual UI
        inner = tb.Frame(canvas, padding=10, bootstyle="secondary")
        window_id = canvas.create_window((0, 0), window=inner, anchor="nw")

        # Helpers to keep scrollregion tight to content and width synced to viewport
        def _update_scrollregion():
            # Requested content size
            content_w = inner.winfo_reqwidth()
            content_h = inner.winfo_reqheight()
            # Viewport size
            viewport_w = canvas.winfo_width()
            viewport_h = canvas.winfo_height()

            # Make inner window width track the canvas width (no horizontal scrolling)
            canvas.itemconfigure(window_id, width=viewport_w)

            # Clamp scrollregion to the exact content height (no phantom top/bottom)
            canvas.configure(scrollregion=(0, 0, content_w, content_h))

            # If content fits vertically, reset view to top and effectively disable scrolling
            if content_h <= viewport_h:
                try:
                    canvas.yview_moveto(0.0)
                except Exception:
                    pass

        # Keep scrollregion synced to inner contents
        inner.bind("<Configure>", lambda e: _update_scrollregion())

        # Make inner frame width match the container width and update scrollregion on viewport resize
        def _on_canvas_configure(event):
            _update_scrollregion()

        canvas.bind("<Configure>", _on_canvas_configure)

        # Mouse wheel scrolling only when pointer is over the canvas
        def _on_mousewheel(event):
            # If content fits, ignore scroll
            if inner.winfo_reqheight() <= canvas.winfo_height():
                return
            delta = event.delta if event.delta else 0
            if delta:
                # Windows/macOS standard scaling
                canvas.yview_scroll(-int(delta / 120) if abs(delta) > 1 else -delta, "units")

        def _bind_mousewheel(_):
            # Bind when entering canvas
            canvas.bind_all("<MouseWheel>", _on_mousewheel)

        def _unbind_mousewheel(_):
            # Unbind when leaving canvas
            canvas.unbind_all("<MouseWheel>")

        canvas.bind("<Enter>", _bind_mousewheel)
        canvas.bind("<Leave>", _unbind_mousewheel)

        # Linux (X11) support: bind only to canvas, and only scroll if content exceeds viewport
        def _on_linux_scroll_up(event):
            if inner.winfo_reqheight() > canvas.winfo_height():
                canvas.yview_scroll(-3, "units")

        def _on_linux_scroll_down(event):
            if inner.winfo_reqheight() > canvas.winfo_height():
                canvas.yview_scroll(+3, "units")

        canvas.bind("<Button-4>", _on_linux_scroll_up)
        canvas.bind("<Button-5>", _on_linux_scroll_down)

        return container, inner

    # ---------- Body ----------
    def _build_body(self, parent):
        # Split: left control panel (scrollable) + right plot area
        paned = tb.Panedwindow(parent, orient=HORIZONTAL, bootstyle="secondary")
        # Before: paned.pack(fill=BOTH, expand=YES, padx=10, pady=6)
        paned.pack(fill=BOTH, expand=YES, padx=10, pady=0)  # remove vertical gap

        # Scrollable left
        left_container, left = self._build_scrollable_left(paned)
        paned.add(left_container, weight=1)

        # Right plot area
        right = tb.Frame(paned, padding=10, bootstyle="secondary")
        paned.add(right, weight=3)
        self.plot_area = right

        # Build left content inside 'left'
        nb = tb.Notebook(left, bootstyle="primary")
        nb.pack(fill=BOTH, expand=YES)

        # Common options
        self._build_common_options(left)

        # Tabs
        self._build_chrom_tab(nb)
        self._build_cdms_tab(nb)
        self._build_mwd_tab(nb)
        self._build_mz_tab(nb)
        self._build_raw_tab(nb)  # Thermo RAW

        # Files panel
        self._build_file_panel(left)

    # ---------- Persistent footer ----------
    def _build_footer(self, parent):
        frame = tb.Frame(parent, bootstyle="secondary")
        # Before: frame.pack(fill=X, pady=6)
        frame.pack(fill=X)  # remove extra gap below body
        self.progress = tb.Progressbar(frame, mode="indeterminate", bootstyle="info")
        self.progress.pack(side=LEFT, fill=X, expand=YES, padx=6)

    def _center_on_screen(self):
        self.root.update_idletasks()
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        x = (sw - self.win_w) // 2
        y = (sh - self.win_h) // 2
        self.root.geometry(f"{self.win_w}x{self.win_h}")
        try:
            self.root.state('zoomed')  # maximize on launch (Windows)
        except Exception:
            pass

    # ---------- Common options ----------
    def _build_common_options(self, parent):
        lf = tb.Labelframe(parent, text="Common Options", padding=10, bootstyle="secondary")
        lf.pack(fill=X, pady=6)

        sp = tb.Frame(lf, bootstyle="secondary")
        sp.pack(fill=X, pady=4)
        tb.Label(sp, text="Save as:", width=12).pack(side=LEFT)
        self.save_var = tk.StringVar()
        tb.Entry(sp, textvariable=self.save_var).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(sp, text="Browse", bootstyle="secondary", command=self._browse_save).pack(side=LEFT)

        grid = tb.Frame(lf, bootstyle="secondary")
        grid.pack(fill=X, pady=4)
        tb.Label(grid, text="Font:", width=8).grid(row=0, column=0, sticky=W)
        self.font_name = tk.StringVar(value="Arial")
        tb.Entry(grid, textvariable=self.font_name, width=14).grid(row=0, column=1, sticky=W, padx=4)

        tb.Label(grid, text="Size:", width=6).grid(row=0, column=2, sticky=W)
        self.font_size = tk.IntVar(value=10)
        tb.Entry(grid, textvariable=self.font_size, width=6).grid(row=0, column=3, sticky=W, padx=4)

        tb.Label(grid, text="DPI:", width=6).grid(row=0, column=4, sticky=W)
        self.dpi = tk.IntVar(value=300)
        tb.Entry(grid, textvariable=self.dpi, width=8).grid(row=0, column=5, sticky=W, padx=4)

        tb.Label(grid, text="Figsize (W,H):", width=14).grid(row=0, column=6, sticky=W)
        self.fig_w = tk.DoubleVar(value=10.0)
        tb.Entry(grid, textvariable=self.fig_w, width=6).grid(row=0, column=7, sticky=W, padx=2)
        self.fig_h = tk.DoubleVar(value=6.0)
        tb.Entry(grid, textvariable=self.fig_h, width=6).grid(row=0, column=8, sticky=W, padx=2)

        style_frame = tb.Frame(lf, bootstyle="secondary")
        style_frame.pack(fill=X, pady=6)
        tb.Label(style_frame, text="X label override:").pack(side=LEFT)
        self.x_label_override = tk.StringVar()
        tb.Entry(style_frame, textvariable=self.x_label_override, width=22).pack(side=LEFT, padx=6)
        tb.Label(style_frame, text="Y label override:").pack(side=LEFT)
        self.y_label_override = tk.StringVar()
        tb.Entry(style_frame, textvariable=self.y_label_override, width=22).pack(side=LEFT, padx=6)
        tb.Label(style_frame, text="Line width:").pack(side=LEFT, padx=(10, 0))
        self.line_width = tk.DoubleVar(value=2.0)
        tb.Entry(style_frame, textvariable=self.line_width, width=6).pack(side=LEFT, padx=6)

        toggles = tb.Frame(lf, bootstyle="secondary")
        toggles.pack(fill=X, pady=4)
        self.show_plot = tk.BooleanVar(value=True)
        tb.Checkbutton(toggles, text="Show Figure in external window (in addition to preview)", variable=self.show_plot).pack(anchor=W, pady=2)
        self.transparent_export = tk.BooleanVar(value=False)
        tb.Checkbutton(toggles, text="Transparent export background", variable=self.transparent_export).pack(anchor=W, pady=2)

    # ---------- Chromatograms tab ----------
    def _build_chrom_tab(self, nb):
        self.chrom_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.chrom_tab, text="Chromatograms")

        mode = tb.Frame(self.chrom_tab, bootstyle="secondary")
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.chrom_mode = tk.StringVar(value="single")
        tb.Radiobutton(mode, text="Single Chromatogram", variable=self.chrom_mode, value="single").pack(side=LEFT, padx=8)
        tb.Radiobutton(mode, text="Stacked Chromatograms", variable=self.chrom_mode, value="stacked").pack(side=LEFT)

        lf_labels = tb.Labelframe(self.chrom_tab, text="Chromatogram Labels & Colors", padding=8, bootstyle="secondary")
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Sample Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.chrom_labels = tk.Text(lf_labels, height=2)
        self.chrom_labels.pack(fill=BOTH, pady=4)
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.chrom_colors = tk.StringVar(value=", ".join(DEFAULT_CHROM_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary")
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.chrom_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info", command=lambda: self._open_palette_editor(self.chrom_colors, "Chromatogram Colors")).pack(side=LEFT)

        lf_axes = tb.Labelframe(self.chrom_tab, text="Chromatogram Axes & Standards", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary")
        grid.pack(fill=X)
        tb.Label(grid, text="X-lim (mL):").grid(row=0, column=0, sticky=W)
        self.chrom_xmin = tk.DoubleVar(value=10.0)
        tb.Entry(grid, textvariable=self.chrom_xmin, width=8).grid(row=0, column=1, padx=4)
        self.chrom_xmax = tk.DoubleVar(value=20.0)
        tb.Entry(grid, textvariable=self.chrom_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-lim:").grid(row=0, column=3, sticky=W)
        self.chrom_ymin = tk.DoubleVar(value=0.0)
        tb.Entry(grid, textvariable=self.chrom_ymin, width=8).grid(row=0, column=4, padx=4)
        self.chrom_ymax = tk.DoubleVar(value=1.0)
        tb.Entry(grid, textvariable=self.chrom_ymax, width=8).grid(row=0, column=5, padx=4)
        self.chrom_normalize = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_axes, text="Normalize traces to max", variable=self.chrom_normalize).pack(anchor=W, pady=4)

        std_row = tb.Frame(lf_axes, bootstyle="secondary")
        std_row.pack(fill=X, pady=4)
        self.chrom_show_standards = tk.BooleanVar(value=True)
        tb.Checkbutton(std_row, text="Show standard volumes (top ticks)", variable=self.chrom_show_standards).pack(side=LEFT)
        tb.Label(std_row, text="Standard Volumes (mL):").pack(side=LEFT, padx=8)
        self.chrom_standards = tk.StringVar(value=DEFAULT_CHROM_STANDARDS)
        tb.Entry(std_row, textvariable=self.chrom_standards, width=28).pack(side=LEFT, padx=6)

        # Replicates panel
        lf_rep = tb.Labelframe(self.chrom_tab, text="Replicates (optional)", padding=8, bootstyle="secondary")
        lf_rep.pack(fill=X, pady=6)
        row_rep1 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep1.pack(fill=X, pady=2)
        self.chrom_rep_enable = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep1, text="Enable replicates overlay", variable=self.chrom_rep_enable).pack(side=LEFT)
        tb.Label(row_rep1, text="Groups: one line per axis; use file indices or paths (comma-separated)").pack(side=LEFT, padx=8)
        self.chrom_rep_groups = tk.Text(lf_rep, height=3)
        self.chrom_rep_groups.pack(fill=X, pady=4)
        row_rep2 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep2.pack(fill=X, pady=2)
        self.chrom_rep_show_mean = tk.BooleanVar(value=True)
        tb.Checkbutton(row_rep2, text="Show mean line", variable=self.chrom_rep_show_mean).pack(side=LEFT)
        self.chrom_rep_show_band = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep2, text="Show SEM band", variable=self.chrom_rep_show_band).pack(side=LEFT, padx=8)
        tb.Label(row_rep2, text="Replicate LW:").pack(side=LEFT, padx=(12, 2))
        self.chrom_rep_lw = tk.DoubleVar(value=1.2)
        tb.Entry(row_rep2, textvariable=self.chrom_rep_lw, width=6).pack(side=LEFT)
        tb.Label(row_rep2, text="Mean LW:").pack(side=LEFT, padx=(12, 2))
        self.chrom_rep_mean_lw = tk.DoubleVar(value=2.2)
        tb.Entry(row_rep2, textvariable=self.chrom_rep_mean_lw, width=6).pack(side=LEFT)
    # ---------- CD-MS tab ----------
    def _build_cdms_tab(self, nb):
        self.cdms_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.cdms_tab, text="CD-MS Mass")

        mode = tb.Frame(self.cdms_tab, bootstyle="secondary")
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.cdms_mode = tk.StringVar(value="single")
        tb.Radiobutton(mode, text="Single CD-MS Mass", variable=self.cdms_mode, value="single").pack(side=LEFT, padx=8)
        tb.Radiobutton(mode, text="Stacked CD-MS Mass", variable=self.cdms_mode, value="stacked").pack(side=LEFT)

        lf_labels = tb.Labelframe(self.cdms_tab, text="CD-MS Mass Labels & Colors", padding=8, bootstyle="secondary")
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.cdms_labels = tk.Text(lf_labels, height=2)
        self.cdms_labels.pack(fill=BOTH, pady=4)
        self.cdms_two_legends = tk.BooleanVar(value=False)
        tb.Checkbutton(lf_labels, text="Include second legend near peak", variable=self.cdms_two_legends).pack(anchor=W, pady=2)
        tb.Label(lf_labels, text="Second Labels (one per file; used when second legend is enabled)").pack(anchor=W)
        self.cdms_second_labels = tk.Text(lf_labels, height=2)
        self.cdms_second_labels.pack(fill=BOTH, pady=4)
        self.cdms_second_labels.config(state=tk.DISABLED)
        self.cdms_two_legends.trace_add("write", lambda *args: self._toggle_second_labels_input())
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.cdms_colors = tk.StringVar(value=", ".join(DEFAULT_CDMS_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary")
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.cdms_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info",
                  command=lambda: self._open_palette_editor(self.cdms_colors, "CD-MS Colors")).pack(side=LEFT)

        lf_axes = tb.Labelframe(self.cdms_tab, text="CD-MS Mass Axes & Layout", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary")
        grid.pack(fill=X)
        tb.Label(grid, text="X-range (kDa):").grid(row=0, column=0, sticky=W)
        self.cdms_xmin = tk.DoubleVar(value=0.0)
        tb.Entry(grid, textvariable=self.cdms_xmin, width=8).grid(row=0, column=1, padx=4)
        self.cdms_xmax = tk.DoubleVar(value=500.0)
        tb.Entry(grid, textvariable=self.cdms_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-max (%):").grid(row=0, column=3, sticky=W)
        self.cdms_ymax_pct = tk.DoubleVar(value=100.0)
        tb.Entry(grid, textvariable=self.cdms_ymax_pct, width=8).grid(row=0, column=4, padx=4)
        self.cdms_normalize = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_axes, text="Normalize to max (percent)", variable=self.cdms_normalize).pack(anchor=W, pady=4)

        # Replicates panel
        lf_rep = tb.Labelframe(self.cdms_tab, text="Replicates (optional)", padding=8, bootstyle="secondary")
        lf_rep.pack(fill=X, pady=6)
        row_rep1 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep1.pack(fill=X, pady=2)
        self.cdms_rep_enable = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep1, text="Enable replicates overlay", variable=self.cdms_rep_enable).pack(side=LEFT)
        tb.Label(row_rep1, text="Groups: one line per axis; use file indices or paths (comma-separated)").pack(side=LEFT, padx=8)
        self.cdms_rep_groups = tk.Text(lf_rep, height=3)
        self.cdms_rep_groups.pack(fill=X, pady=4)
        row_rep2 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep2.pack(fill=X, pady=2)
        self.cdms_rep_show_mean = tk.BooleanVar(value=True)
        tb.Checkbutton(row_rep2, text="Show mean line", variable=self.cdms_rep_show_mean).pack(side=LEFT)
        self.cdms_rep_show_band = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep2, text="Show SEM band", variable=self.cdms_rep_show_band).pack(side=LEFT, padx=8)
        tb.Label(row_rep2, text="Replicate LW:").pack(side=LEFT, padx=(12, 2))
        self.cdms_rep_lw = tk.DoubleVar(value=1.2)
        tb.Entry(row_rep2, textvariable=self.cdms_rep_lw, width=6).pack(side=LEFT)
        tb.Label(row_rep2, text="Mean LW:").pack(side=LEFT, padx=(12, 2))
        self.cdms_rep_mean_lw = tk.DoubleVar(value=2.2)
        tb.Entry(row_rep2, textvariable=self.cdms_rep_mean_lw, width=6).pack(side=LEFT)

    # ---------- MWD tab ----------
    def _build_mwd_tab(self, nb):
        self.mwd_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.mwd_tab, text="MWD Chromatograms")

        lf_files = tb.Labelframe(self.mwd_tab, text="MWD Files (drag & drop or browse)", padding=8, bootstyle="secondary")
        lf_files.pack(fill=X, pady=4)
        row1 = tb.Frame(lf_files, bootstyle="secondary")
        row1.pack(fill=X, pady=2)
        tb.Label(row1, text="mAU file:").pack(side=LEFT)
        self.mwd_mau_path = tk.StringVar()
        self.mwd_mau_entry = tb.Entry(row1, textvariable=self.mwd_mau_path)
        self.mwd_mau_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row1, text="Browse", bootstyle="secondary", command=self._browse_mwd_mau).pack(side=LEFT)
        if DND_AVAILABLE:
            self.mwd_mau_entry.drop_target_register(DND_FILES)
            self.mwd_mau_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.mwd_mau_path))

        row2 = tb.Frame(lf_files, bootstyle="secondary")
        row2.pack(fill=X, pady=2)
        tb.Label(row2, text="mV file:").pack(side=LEFT)
        self.mwd_mv_path = tk.StringVar()
        self.mwd_mv_entry = tb.Entry(row2, textvariable=self.mwd_mv_path)
        self.mwd_mv_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row2, text="Browse", bootstyle="secondary", command=self._browse_mwd_mv).pack(side=LEFT)
        if DND_AVAILABLE:
            self.mwd_mv_entry.drop_target_register(DND_FILES)
            self.mwd_mv_entry.dnd_bind("<<Drop>>", lambda e: self._on_single_file_drop(e, self.mwd_mv_path))

        lf_labels = tb.Labelframe(self.mwd_tab, text="Labels & Colors", padding=8, bootstyle="secondary")
        lf_labels.pack(fill=X, pady=4)
        row3 = tb.Frame(lf_labels, bootstyle="secondary")
        row3.pack(fill=X, pady=2)
        tb.Label(row3, text="mAU label:").pack(side=LEFT)
        self.mwd_mau_label = tk.StringVar(value="AqpZ (mAU)")
        tb.Entry(row3, textvariable=self.mwd_mau_label, width=24).pack(side=LEFT, padx=6)
        tb.Label(row3, text="mV label:").pack(side=LEFT)
        self.mwd_mv_label = tk.StringVar(value="AqpZ (mV)")
        tb.Entry(row3, textvariable=self.mwd_mv_label, width=24).pack(side=LEFT, padx=6)

        row4 = tb.Frame(lf_labels, bootstyle="secondary")
        row4.pack(fill=X, pady=2)
        tb.Label(row4, text="Colors (mAU, mV hex):").pack(side=LEFT)
        self.mwd_colors = tk.StringVar(value=", ".join(DEFAULT_MWD_COLORS))
        tb.Entry(row4, textvariable=self.mwd_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row4, text="Pick mAU", bootstyle="info", command=lambda: self._pick_mwd_color(index=0)).pack(side=LEFT, padx=4)
        tb.Button(row4, text="Pick mV", bootstyle="info", command=lambda: self._pick_mwd_color(index=1)).pack(side=LEFT, padx=4)

        lf_axes = tb.Labelframe(self.mwd_tab, text="Axes & Options", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=4)
        grid = tb.Frame(lf_axes, bootstyle="secondary")
        grid.pack(fill=X)
        tb.Label(grid, text="X-lim (mL):").grid(row=0, column=0, sticky=W)
        self.mwd_xmin = tk.DoubleVar(value=0.0)
        tb.Entry(grid, textvariable=self.mwd_xmin, width=8).grid(row=0, column=1, padx=4)
        self.mwd_xmax = tk.DoubleVar(value=26.0)
        tb.Entry(grid, textvariable=self.mwd_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-max (%):").grid(row=0, column=3, sticky=W)
        self.mwd_ymax_pct = tk.DoubleVar(value=100.0)
        tb.Entry(grid, textvariable=self.mwd_ymax_pct, width=8).grid(row=0, column=4, padx=4)
        self.mwd_normalize = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_axes, text="Normalize each trace to max", variable=self.mwd_normalize).pack(anchor=W, pady=4)

    # ---------- m/z tab ----------
    def _build_mz_tab(self, nb):
        self.mz_tab = tb.Frame(nb, padding=8, bootstyle="secondary")
        nb.add(self.mz_tab, text="m/z Spectra")

        mode = tb.Frame(self.mz_tab, bootstyle="secondary")
        mode.pack(fill=X, pady=4)
        tb.Label(mode, text="Plot Mode:", font=("Arial", 10, "bold")).pack(side=LEFT)
        self.mz_mode = tk.StringVar(value="single")
        tb.Radiobutton(mode, text="Single m/z Spectrum", variable=self.mz_mode, value="single").pack(side=LEFT, padx=8)
        tb.Radiobutton(mode, text="Stacked m/z Spectra", variable=self.mz_mode, value="stacked").pack(side=LEFT)

        lf_labels = tb.Labelframe(self.mz_tab, text="m/z Labels & Colors", padding=8, bootstyle="secondary")
        lf_labels.pack(fill=X, pady=4)
        tb.Label(lf_labels, text="Labels (comma-separated; blank = filenames)").pack(anchor=W)
        self.mz_labels = tk.Text(lf_labels, height=2)
        self.mz_labels.pack(fill=BOTH, pady=4)
        self.mz_two_legends = tk.BooleanVar(value=False)
        tb.Checkbutton(lf_labels, text="Include second legend near peak", variable=self.mz_two_legends).pack(anchor=W, pady=2)
        tb.Label(lf_labels, text="Second Labels (one per file; used when second legend is enabled)").pack(anchor=W)
        self.mz_second_labels = tk.Text(lf_labels, height=2)
        self.mz_second_labels.pack(fill=BOTH, pady=4)
        self.mz_second_labels.config(state=tk.DISABLED)
        self.mz_two_legends.trace_add("write", lambda *args: self._toggle_mz_second_labels())
        tb.Label(lf_labels, text="Colors (comma-separated hex)").pack(anchor=W)
        self.mz_colors = tk.StringVar(value=", ".join(DEFAULT_MZ_COLORS))
        color_row = tb.Frame(lf_labels, bootstyle="secondary")
        color_row.pack(fill=X, pady=2)
        tb.Entry(color_row, textvariable=self.mz_colors).pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(color_row, text="Pick Colors…", bootstyle="info",
                  command=lambda: self._open_palette_editor(self.mz_colors, "m/z Colors")).pack(side=LEFT)

        lf_axes = tb.Labelframe(self.mz_tab, text="m/z Axes & Layout", padding=8, bootstyle="secondary")
        lf_axes.pack(fill=X, pady=6)
        grid = tb.Frame(lf_axes, bootstyle="secondary")
        grid.pack(fill=X)
        tb.Label(grid, text="X-range (m/z):").grid(row=0, column=0, sticky=W)
        self.mz_xmin = tk.DoubleVar(value=100.0)
        tb.Entry(grid, textvariable=self.mz_xmin, width=8).grid(row=0, column=1, padx=4)
        self.mz_xmax = tk.DoubleVar(value=2000.0)
        tb.Entry(grid, textvariable=self.mz_xmax, width=8).grid(row=0, column=2, padx=4)
        tb.Label(grid, text="Y-max (%):").grid(row=0, column=3, sticky=W)
        self.mz_ymax_pct = tk.DoubleVar(value=100.0)
        tb.Entry(grid, textvariable=self.mz_ymax_pct, width=8).grid(row=0, column=4, padx=4)
        self.mz_normalize = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_axes, text="Normalize to max (percent)", variable=self.mz_normalize).pack(anchor=W, pady=4)

        # Replicates panel
        lf_rep = tb.Labelframe(self.mz_tab, text="Replicates (optional)", padding=8, bootstyle="secondary")
        lf_rep.pack(fill=X, pady=6)
        row_rep1 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep1.pack(fill=X, pady=2)
        self.mz_rep_enable = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep1, text="Enable replicates overlay", variable=self.mz_rep_enable).pack(side=LEFT)
        tb.Label(row_rep1, text="Groups: one line per axis; use file indices or paths (comma-separated)").pack(side=LEFT, padx=8)
        self.mz_rep_groups = tk.Text(lf_rep, height=3)
        self.mz_rep_groups.pack(fill=X, pady=4)
        row_rep2 = tb.Frame(lf_rep, bootstyle="secondary")
        row_rep2.pack(fill=X, pady=2)
        self.mz_rep_show_mean = tk.BooleanVar(value=True)
        tb.Checkbutton(row_rep2, text="Show mean line", variable=self.mz_rep_show_mean).pack(side=LEFT)
        self.mz_rep_show_band = tk.BooleanVar(value=False)
        tb.Checkbutton(row_rep2, text="Show SEM band", variable=self.mz_rep_show_band).pack(side=LEFT, padx=8)
        tb.Label(row_rep2, text="Replicate LW:").pack(side=LEFT, padx=(12, 2))
        self.mz_rep_lw = tk.DoubleVar(value=1.2)
        tb.Entry(row_rep2, textvariable=self.mz_rep_lw, width=6).pack(side=LEFT)
        tb.Label(row_rep2, text="Mean LW:").pack(side=LEFT, padx=(12, 2))
        self.mz_rep_mean_lw = tk.DoubleVar(value=2.2)
        tb.Entry(row_rep2, textvariable=self.mz_rep_mean_lw, width=6).pack(side=LEFT)

    # ---------- Thermo RAW tab ----------
    def _build_raw_tab(self, nb):
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
        tb.Label(grid, text="MS level:").grid(row=0, column=0, sticky=W)
        self.raw_ms_level = tk.IntVar(value=1)
        tb.Entry(grid, textvariable=self.raw_ms_level, width=6).grid(row=0, column=1, padx=6)
        tb.Label(grid, text="Max scans:").grid(row=0, column=2, sticky=W)
        self.raw_max_scans = tk.IntVar(value=200)
        tb.Entry(grid, textvariable=self.raw_max_scans, width=8).grid(row=0, column=3, padx=6)

        tb.Label(grid, text="Spectrum:").grid(row=1, column=0, sticky=W)
        self.raw_mode = tk.StringVar(value="apex")
        tb.Radiobutton(grid, text="Apex", variable=self.raw_mode, value="apex").grid(row=1, column=1, sticky=W)
        tb.Radiobutton(grid, text="Index", variable=self.raw_mode, value="index").grid(row=1, column=2, sticky=W)
        tb.Label(grid, text="Index:").grid(row=1, column=3, sticky=W)
        self.raw_index = tk.IntVar(value=0)
        self.raw_index_entry = tb.Entry(grid, textvariable=self.raw_index, width=6)
        self.raw_index_entry.grid(row=1, column=4, padx=6)

        def _toggle_raw_index(*args):
            state = tk.NORMAL if self.raw_mode.get() == "index" else tk.DISABLED
            self.raw_index_entry.config(state=state)

        self.raw_mode.trace_add("write", _toggle_raw_index)
        _toggle_raw_index()

        self.raw_norm = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_opts, text="Normalize spectrum to max (percent)", variable=self.raw_norm).pack(anchor=W, pady=4)

        tb.Label(grid, text="Spec Y-max (%):").grid(row=2, column=0, sticky=W, pady=(6, 0))
        self.raw_ymax_pct = tk.DoubleVar(value=100.0)
        tb.Entry(grid, textvariable=self.raw_ymax_pct, width=8).grid(row=2, column=1, padx=6, pady=(6, 0))

        lf_ms = tb.Labelframe(self.raw_tab, text="msconvert (ProteoWizard)", padding=8, bootstyle="secondary")
        lf_ms.pack(fill=X, pady=4)
        row2 = tb.Frame(lf_ms, bootstyle="secondary")
        row2.pack(fill=X, pady=2)
        tb.Label(row2, text="msconvert path (optional):").pack(side=LEFT)
        self.msconvert_path = tk.StringVar()
        self.msconvert_entry = tb.Entry(row2, textvariable=self.msconvert_path)
        self.msconvert_entry.pack(side=LEFT, fill=X, expand=YES, padx=6)
        tb.Button(row2, text="Browse", bootstyle="secondary", command=self._browse_msconvert).pack(side=LEFT)

        self.raw_prefer_native = tk.BooleanVar(value=True)
        tb.Checkbutton(lf_ms, text="Prefer native Thermo MSFileReader on Windows (fallback to msconvert)",
                       variable=self.raw_prefer_native).pack(anchor=W, pady=4)

    # ---------- Files panel ----------
    def _build_file_panel(self, parent):
        lf = tb.Labelframe(parent, text="Files (Chromatograms, CD-MS Mass, MWD, m/z)", padding=8, bootstyle="secondary")
        lf.pack(fill=BOTH, expand=YES, pady=6)
        btns = tb.Frame(lf, bootstyle="secondary")
        btns.pack(fill=X, pady=4)
        tb.Button(btns, text="Add Files...", bootstyle="success", command=self._add_files).pack(side=LEFT)
        tb.Button(btns, text="Remove Selected", bootstyle="warning", command=self._remove_selected).pack(side=LEFT, padx=6)
        tb.Button(btns, text="Clear", bootstyle="danger", command=self._clear_files).pack(side=LEFT)

        self.file_tree = tb.Treeview(lf, columns=("idx", "path"), show="headings", height=8, bootstyle="info")
        self.file_tree.heading("idx", text="#")
        self.file_tree.column("idx", width=40, anchor=CENTER)
        self.file_tree.heading("path", text="File Path")
        self.file_tree.column("path", width=700, anchor=W)
        self.file_tree.pack(fill=BOTH, expand=YES, pady=6)

        if DND_AVAILABLE:
            self.file_tree.drop_target_register(DND_FILES)
            self.file_tree.dnd_bind("<<Drop>>", self._on_files_drop)

    # ---------- Utility ----------
    def _set_status(self, text: str):
        self.status_label.configure(text=f"Status: {text}")
        self.root.update_idletasks()

    def _browse_save(self):
        path = filedialog.asksaveasfilename(title="Save figure as...", defaultextension=".png",
                                            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
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
            paths = filedialog.askopenfilenames(title="Select Chromatogram Excel files",
                                                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
        elif tab == "CD-MS Mass":
            paths = filedialog.askopenfilenames(title="Select CD-MS mass text files",
                                                filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        elif tab == "m/z Spectra":
            paths = filedialog.askopenfilenames(title="Select m/z spectrum text files",
                                                filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        else:  # MWD
            paths = filedialog.askopenfilenames(title="Select mAU/mV Excel files (first=mAU, second=mV)",
                                                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
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
        idxs = sorted([int(self.file_tree.item(i, "values")[0]) for i in selected], reverse=True)
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
        if paths:
            target_var.set(paths[0])
            self._set_status("Set file via drag-and-drop.")

    # ---------- Color pickers ----------
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
        tb.Button(ctrl, text="Add…", bootstyle="info",
                  command=lambda: self._palette_add(colors, refresh_list)).pack(side=LEFT, padx=4)
        tb.Button(ctrl, text="Remove", bootstyle="warning",
                  command=lambda: self._palette_remove(colors, lb, refresh_list)).pack(side=LEFT, padx=4)
        tb.Button(ctrl, text="Clear", bootstyle="danger",
                  command=lambda: self._palette_clear(colors, refresh_list)).pack(side=LEFT, padx=4)
        tb.Button(ctrl, text="Apply", bootstyle="success",
                  command=lambda: self._palette_apply(colors, target_var, top)).pack(side=RIGHT, padx=4)
        tb.Button(ctrl, text="Cancel", bootstyle="secondary", command=top.destroy).pack(side=RIGHT, padx=4)

    def _palette_add(self, colors: List[str], refresh_list):
        start_color = colors[-1] if colors else "#00A9A5"
        c = self._pick_color(initial=start_color)
        if c:
            colors.append(c)
            refresh_list()

    def _palette_remove(self, colors: List[str], lb: tk.Listbox, refresh_list):
        selection = lb.curselection()
        if selection:
            idx = selection[0]
            colors.pop(idx)
            refresh_list()

    def _palette_clear(self, colors: List[str], refresh_list):
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

    # ---------- Toggle helpers ----------
    def _toggle_second_labels_input(self):
        state = tk.NORMAL if self.cdms_two_legends.get() else tk.DISABLED
        self.cdms_second_labels.config(state=state)

    def _toggle_mz_second_labels(self):
        state = tk.NORMAL if self.mz_two_legends.get() else tk.DISABLED
        self.mz_second_labels.config(state=state)

    # ---------- Browse helpers ----------
    def _browse_msconvert(self):
        path = filedialog.askopenfilename(
            title="Select msconvert.exe",
            filetypes=[("Executable", "*.exe;*"), ("All files", "*.*")]
        )
        if path:
            self.msconvert_path.set(path)
            self._set_status("msconvert path set.")

    def _browse_raw_file(self):
        path = filedialog.askopenfilename(
            title="Select Thermo RAW or mzML",
            filetypes=[("RAW/mzML", "*.raw *.mzML"), ("All files", "*.*")]
        )
        if path:
            self.raw_path.set(path)
            self._set_status("RAW/mzML file selected.")

    def _browse_mwd_mau(self):
        path = filedialog.askopenfilename(
            title="Select mAU Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.mwd_mau_path.set(path)
            self._set_status("mAU file selected.")

    def _browse_mwd_mv(self):
        path = filedialog.askopenfilename(
            title="Select mV Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.mwd_mv_path.set(path)
            self._set_status("mV file selected.")

    # ---------- Replicates helpers ----------
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

    def _interpolate_to_grid(self, x: np.ndarray, y: np.ndarray, xgrid: np.ndarray) -> np.ndarray:
        return np.interp(xgrid, x, y, left=np.nan, right=np.nan)

    def _overlay_chrom_replicates(self, axes, groups, colors, normalize: bool, xlim):
        xgrid = np.linspace(xlim[0], xlim[1], 1200)
        rep_lw = float(self.chrom_rep_lw.get())
        mean_lw = float(self.chrom_rep_mean_lw.get())
        show_mean = bool(self.chrom_rep_show_mean.get())
        show_band = bool(self.chrom_rep_show_band.get())
        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]

        for i, ax in enumerate(axes_list):
            if i >= len(groups): break
            gfiles = groups[i]
            traces = []
            for fp in gfiles:
                try:
                    df = read_chrom_file(fp)
                    x = df["ml"].astype(float).values
                    y = df["mAU"].astype(float).values
                    if normalize:
                        m = float(np.nanmax(y))
                        if m > 0: y = y / m
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
                    ax.fill_between(xgrid, mean - sem, mean + sem, color=colors[i % len(colors)], alpha=0.15, linewidth=0)
                ax.plot(xgrid, mean, color=colors[i % len(colors)], lw=mean_lw)

    def _overlay_cdms_replicates(self, axes, groups, colors, normalize: bool, x_range):
        xgrid = np.linspace(x_range[0], x_range[1], 2000)
        rep_lw = float(self.cdms_rep_lw.get())
        mean_lw = float(self.cdms_rep_mean_lw.get())
        show_mean = bool(self.cdms_rep_show_mean.get())
        show_band = bool(self.cdms_rep_show_band.get())
        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]

        for i, ax in enumerate(axes_list):
            if i >= len(groups): break
            gfiles = groups[i]
            traces = []
            for fp in gfiles:
                try:
                    df = read_cdms_mass_file(fp)
                    df = df[(df["mass_kDa"] >= x_range[0]) & (df["mass_kDa"] <= x_range[1])].copy()
                    x = df["mass_kDa"].astype(float).values
                    y = df["intensity"].astype(float).values
                    if normalize:
                        m = float(np.nanmax(y))
                        if m > 0: y = (y / m) * 100.0
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
                    ax.fill_between(xgrid, mean - sem, mean + sem, color=colors[i % len(colors)], alpha=0.15, linewidth=0)
                ax.plot(xgrid, mean, color=colors[i % len(colors)], lw=mean_lw)

    def _overlay_mz_replicates(self, axes, groups, colors, normalize: bool, x_range):
        xgrid = np.linspace(x_range[0], x_range[1], 4000)
        rep_lw = float(self.mz_rep_lw.get())
        mean_lw = float(self.mz_rep_mean_lw.get())
        show_mean = bool(self.mz_rep_show_mean.get())
        show_band = bool(self.mz_rep_show_band.get())
        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]

        for i, ax in enumerate(axes_list):
            if i >= len(groups): break
            gfiles = groups[i]
            traces = []
            for fp in gfiles:
                try:
                    df = read_mz_spectrum_file(fp)
                    df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
                    x = df["mz"].astype(float).values
                    y = df["intensity"].astype(float).values
                    if normalize:
                        m = float(np.nanmax(y))
                        if m > 0: y = (y / m) * 100.0
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
                    ax.fill_between(xgrid, mean - sem, mean + sem, color=colors[i % len(colors)], alpha=0.15, linewidth=0)
                ax.plot(xgrid, mean, color=colors[i % len(colors)], lw=mean_lw)

    # ---------- Post-process overrides ----------
    def _apply_overrides(self, fig, axes_obj, tab_name: str, normalize: bool):
        if axes_obj is None:
            ax_list = []
        elif isinstance(axes_obj, (list, tuple)):
            ax_list = list(axes_obj)
        else:
            ax_list = [axes_obj]

        # Line width override
        try:
            lw = float(self.line_width.get())
            for ax in ax_list:
                for ln in ax.get_lines():
                    ln.set_linewidth(lw)
        except Exception:
            pass

        # Axis label overrides
        xlab = self.x_label_override.get().strip()
        ylab = self.y_label_override.get().strip()

        if xlab:
            applied = False
            try:
                fig.supxlabel(xlab)
                applied = True
            except Exception:
                applied = False
            if not applied:
                try:
                    (ax_list[-1] if ax_list else fig).set_xlabel(xlab)
                except Exception:
                    fig.text(0.5, 0.02, xlab, ha="center")

        if ylab:
            applied = False
            try:
                fig.supylabel(ylab)
                applied = True
            except Exception:
                applied = False
            if not applied:
                try:
                    if tab_name in ("MWD Chromatograms", "Thermo RAW") and len(ax_list) >= 2:
                        ax_list[0].set_ylabel(ylab)
                        ax_list[1].set_ylabel(ylab)
                    elif ax_list:
                        ax_list[0].set_ylabel(ylab)
                    else:
                        fig.text(0.02, 0.5, ylab, va="center", rotation="vertical")
                except Exception:
                    pass

        # Clamp normalized Y
        try:
            if tab_name == "CD-MS Mass" and normalize:
                ymax = min(float(self.cdms_ymax_pct.get()), 100.0)
                for ax in ax_list: ax.set_ylim(0, ymax)
            elif tab_name == "m/z Spectra" and normalize:
                ymax = min(float(self.mz_ymax_pct.get()), 100.0)
                for ax in ax_list: ax.set_ylim(0, ymax)
            elif tab_name == "MWD Chromatograms" and normalize and len(ax_list) == 2:
                ymax = min(float(self.mwd_ymax_pct.get()), 100.0)
                ax_list[0].set_ylim(0, ymax)
                ax_list[1].set_ylim(0, ymax)
            elif tab_name == "Chromatograms":
                ymin = float(self.chrom_ymin.get())
                ymax = float(self.chrom_ymax.get())
                if normalize: ymax = min(ymax, 1.0)
                for ax in ax_list: ax.set_ylim(ymin, ymax)
            elif tab_name == "Thermo RAW" and normalize and len(ax_list) == 2:
                ymax = min(float(self.raw_ymax_pct.get()), 100.0)
                ax_list[1].set_ylim(0, ymax)
        except Exception:
            pass

        try:
            fig.canvas.draw_idle()
        except Exception:
            pass

    # ---------- Plot ----------
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

            # Chromatograms
            if tab == "Chromatograms":
                raw_labels = self.chrom_labels.get("1.0", "end").strip()
                labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
                colors = parse_str_list(self.chrom_colors.get()) or DEFAULT_CHROM_COLORS
                std_vols = parse_float_list(self.chrom_standards.get()) if self.chrom_show_standards.get() else None
                normalize = bool(self.chrom_normalize.get())
                if self.chrom_mode.get() == "single":
                    if len(self.files) != 1:
                        Messagebox.show_error("Select exactly one Excel file for Single Chromatogram.", title="Selection error")
                        return
                    fig, axes = plot_chrom_stack(
                        files=self.files[:1], sample_labels=labels[:1], colors=colors[:1],
                        standard_volumes=std_vols, normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        xlim=(float(self.chrom_xmin.get()), float(self.chrom_xmax.get())),
                        ylim=(float(self.chrom_ymin.get()), float(self.chrom_ymax.get())),
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                else:
                    if len(labels) != len(self.files):
                        Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
                        return
                    fig, axes = plot_chrom_stack(
                        files=self.files, sample_labels=labels, colors=colors,
                        standard_volumes=std_vols, normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        xlim=(float(self.chrom_xmin.get()), float(self.chrom_xmax.get())),
                        ylim=(float(self.chrom_ymin.get()), float(self.chrom_ymax.get())),
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                # Replicates
                if bool(self.chrom_rep_enable.get()):
                    groups = self._parse_rep_groups(self.chrom_rep_groups.get("1.0", "end"))
                    axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
                    if self.chrom_mode.get() == "single":
                        groups = groups[:1]
                    elif len(groups) != len(axes_list):
                        Messagebox.show_error("Replicate groups must have one line per axis (stacked).", title="Replicates mismatch")
                    else:
                        self._overlay_chrom_replicates(axes_list, groups, colors, normalize,
                                                       (float(self.chrom_xmin.get()), float(self.chrom_xmax.get())))
                self._apply_overrides(fig, axes, tab_name="Chromatograms", normalize=normalize)
                self._render_figure(fig)
                self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                self._set_status("Chromatogram plot complete.")

            # CD-MS Mass
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
                        Messagebox.show_error("Select exactly one text file for Single CD-MS Mass.", title="Selection error")
                        return
                    second_label_single = second_labels[0] if use_second else None
                    fig, ax = plot_single_cdms_mass(
                        dataset=self.files[0], label=labels[0], second_label=second_label_single,
                        x_range=x_range, color=colors[0], normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                    if bool(self.cdms_rep_enable.get()):
                        groups = self._parse_rep_groups(self.cdms_rep_groups.get("1.0", "end"))[:1]
                        self._overlay_cdms_replicates([ax], groups, [colors[0]], normalize, x_range)
                    self._apply_overrides(fig, ax, tab_name="CD-MS Mass", normalize=normalize)
                    self._render_figure(fig)
                    self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                    self._set_status("CD-MS mass (single) plot complete.")

                else:
                    if len(labels) != len(self.files):
                        Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
                        return
                    if use_second and len(second_labels) != len(self.files):
                        Messagebox.show_error(
                            "Second labels must match number of files when second legend is enabled.",
                            title="Second label mismatch"
                        )
                        return
                    y_ticks = ("0", "%", "100") if normalize else ("", "", "")
                    fig, axes = plot_cdms_mass_offset_stack(
                        datasets=self.files, labels=labels,
                        second_labels=second_labels if use_second else [None] * len(self.files),
                        x_range=x_range, colors=colors, y_offset_step=float(self.cdms_y_offset.get()),
                        y_tick_labels=y_ticks, normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                    if bool(self.cdms_rep_enable.get()):
                        groups = self._parse_rep_groups(self.cdms_rep_groups.get("1.0", "end"))
                        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
                        if len(groups) != len(axes_list):
                            Messagebox.show_error("Replicate groups must have one line per axis (stacked).",
                                                  title="Replicates mismatch")
                        else:
                            self._overlay_cdms_replicates(axes_list, groups, colors, normalize, x_range)
                    self._apply_overrides(fig, axes, tab_name="CD-MS Mass", normalize=normalize)
                    self._render_figure(fig)
                    self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                    self._set_status("CD-MS mass (stacked) plot complete.")

            # MWD Chromatograms
            elif tab == "MWD Chromatograms":
                path_mau = self.mwd_mau_path.get().strip()
                path_mv = self.mwd_mv_path.get().strip()
                if not path_mau or not path_mv:
                    Messagebox.show_error("Please select both mAU and mV files.", title="MWD files missing")
                    return
                colors = parse_str_list(self.mwd_colors.get()) or DEFAULT_MWD_COLORS
                if len(colors) < 2:
                    Messagebox.show_error("Please provide two colors (mAU, mV).", title="Colors required")
                    return
                xlim = (float(self.mwd_xmin.get()), float(self.mwd_xmax.get()))
                normalize = bool(self.mwd_normalize.get())
                fig, (ax1, ax2) = plot_mwd_chromatograms(
                    mau_file=path_mau, mv_file=path_mv,
                    mau_label=self.mwd_mau_label.get(), mv_label=self.mwd_mv_label.get(),
                    colors=colors[:2], xlim=xlim, normalize=normalize,
                    font_name=cfg.font_name, font_size=cfg.font_size,
                    figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                )
                self._apply_overrides(fig, (ax1, ax2), tab_name="MWD Chromatograms", normalize=normalize)
                self._render_figure(fig)
                self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                self._set_status("MWD Chromatograms plot complete.")

            # m/z Spectra
            elif tab == "m/z Spectra":
                raw_labels = self.mz_labels.get("1.0", "end").strip()
                labels = parse_str_list(raw_labels) if raw_labels else filenames_to_labels(self.files)
                use_second = bool(self.mz_two_legends.get())
                second_raw = self.mz_second_labels.get("1.0", "end").strip()
                second_labels = parse_str_list(second_raw) if (use_second and second_raw) else [None] * len(self.files)
                colors = parse_str_list(self.mz_colors.get()) or DEFAULT_MZ_COLORS
                x_range = (float(self.mz_xmin.get()), float(self.mz_xmax.get()))
                normalize = bool(self.mz_normalize.get())

                if self.mz_mode.get() == "single":
                    if len(self.files) != 1:
                        Messagebox.show_error("Select exactly one text file for Single m/z Spectrum.", title="Selection error")
                        return
                    second_label_single = second_labels[0] if use_second else None
                    fig, ax = plot_single_mz_spectrum(
                        dataset=self.files[0], label=labels[0], second_label=second_label_single,
                        x_range=x_range, color=colors[0], normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                    if bool(self.mz_rep_enable.get()):
                        groups = self._parse_rep_groups(self.mz_rep_groups.get("1.0", "end"))[:1]
                        self._overlay_mz_replicates([ax], groups, [colors[0]], normalize, x_range)
                    self._apply_overrides(fig, ax, tab_name="m/z Spectra", normalize=normalize)
                    self._render_figure(fig)
                    self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                    self._set_status("m/z spectrum (single) plot complete.")

                else:
                    if len(labels) != len(self.files):
                        Messagebox.show_error("Number of labels must match number of files.", title="Label mismatch")
                        return
                    if use_second and len(second_labels) != len(self.files):
                        Messagebox.show_error(
                            "Second labels must match number of files when second legend is enabled.",
                            title="Second label mismatch"
                        )
                        return
                    y_ticks = ("0", "%", "100") if normalize else ("", "", "")
                    fig, axes = plot_mz_spectrum_offset_stack(
                        datasets=self.files, labels=labels,
                        second_labels=second_labels if use_second else [None] * len(self.files),
                        x_range=x_range, colors=colors, y_offset_step=float(self.mz_y_offset.get()),
                        y_tick_labels=y_ticks, normalize=normalize,
                        font_name=cfg.font_name, font_size=cfg.font_size,
                        figsize=cfg.figsize, dpi=DISPLAY_DPI, save=None, show=cfg.show,
                    )
                    if bool(self.mz_rep_enable.get()):
                        groups = self._parse_rep_groups(self.mz_rep_groups.get("1.0", "end"))
                        axes_list = axes if isinstance(axes, (list, tuple)) else [axes]
                        if len(groups) != len(axes_list):
                            Messagebox.show_error("Replicate groups must have one line per axis (stacked).",
                                                  title="Replicates mismatch")
                        else:
                            self._overlay_mz_replicates(axes_list, groups, colors, normalize, x_range)
                    self._apply_overrides(fig, axes, tab_name="m/z Spectra", normalize=normalize)
                    self._render_figure(fig)
                    self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                    self._set_status("m/z spectra (stacked) plot complete.")

            # Thermo RAW (TIC + m/z)
            else:
                raw_path = self.raw_path.get().strip()
                if not raw_path:
                    Messagebox.show_error("Please pick a RAW or mzML file.", title="RAW file missing")
                    return
                fig, (ax1, ax2) = plot_raw_chrom_and_spectrum(
                    raw_path=raw_path,
                    ms_level=int(self.raw_ms_level.get()),
                    max_scans=int(self.raw_max_scans.get()),
                    spectrum_mode=self.raw_mode.get(),
                    spectrum_index=(int(self.raw_index.get()) if self.raw_mode.get() == "index" else None),
                    normalize_spectrum=bool(self.raw_norm.get()),
                    figsize=cfg.figsize, dpi=DISPLAY_DPI,
                    font_name=cfg.font_name, font_size=cfg.font_size,
                    save=None, show=cfg.show,
                    prefer_native=bool(self.raw_prefer_native.get()),
                    msconvert_path=(self.msconvert_path.get().strip() or None),
                )
                self._apply_overrides(fig, (ax1, ax2), tab_name="Thermo RAW", normalize=bool(self.raw_norm.get()))
                self._render_figure(fig)
                self._export_figure(fig, cfg.save, cfg.dpi, cfg.transparent)
                self._set_status("Thermo RAW plot (TIC + m/z) complete.")

        except Exception as e:
            self._set_status("Error")
            Messagebox.show_error(f"Plotting failed:\n{e}", title="Error")
        finally:
            self.progress.stop()

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

def main():
    root = TkinterDnD.Tk() if DND_AVAILABLE else tk.Tk()
    try:
        root.tk.call('package', 'require', 'msgcat')
    except tk.TclError:
        pass

    # Detect system light/dark and choose a suitable default theme
    default_theme = choose_default_theme()
    app = PlotApp(root, themename=default_theme)

    root.mainloop()


if __name__ == "__main__":
    main()
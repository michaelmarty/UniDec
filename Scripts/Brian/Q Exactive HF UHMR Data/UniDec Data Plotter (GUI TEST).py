"""Modern GUI for plotting UniDec two-column text exports."""

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Optional
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter


FIGURE_SIZE = (16, 9)
FONT_FAMILY = "Arial"
FONT_SIZE = 28
TICK_SIZE = 12
TICK_WIDTH = 2
STAGGER_STEP_NORMALIZED = 0.1
STAGGER_STEP_FRACTION = 0.1

THEMES = {
    "Light": {
        "app_bg": "#f4f6f8",
        "panel_bg": "#ffffff",
        "text_primary": "#172033",
        "text_secondary": "#586174",
        "accent": "#2563eb",
        "accent_active": "#1d4ed8",
        "border": "#d7dce5",
        "input_bg": "#ffffff",
        "button_bg": "#eef2f7",
        "button_active": "#e1e7f0",
        "disabled_text": "#8d96a8",
        "plot_bg": "#ffffff",
        "empty_plot_text": "#6b7280",
    },
    "Dark": {
        "app_bg": "#111827",
        "panel_bg": "#1f2937",
        "text_primary": "#f8fafc",
        "text_secondary": "#cbd5e1",
        "accent": "#3b82f6",
        "accent_active": "#2563eb",
        "border": "#374151",
        "input_bg": "#111827",
        "button_bg": "#374151",
        "button_active": "#4b5563",
        "disabled_text": "#64748b",
        "plot_bg": "#ffffff",
        "empty_plot_text": "#6b7280",
    },
}


@dataclass(frozen=True)
class PlotConfig:
    """Settings that differ between the three supported plot types."""

    name: str
    x_label: str
    output_stem: str
    single_output_stem: str
    zoom_start_label: str
    zoom_end_label: str
    zoom_value_error: str
    zoom_order_error: str
    x_scale: float = 1.0


@dataclass(frozen=True)
class Dataset:
    path: Path
    values: np.ndarray
    legend_label: Optional[str] = None


@dataclass
class DatasetInput:
    """Editable dataset values captured by the GUI before files are loaded."""

    path: str = ""
    legend_label: str = ""


@dataclass(frozen=True)
class PlotOptions:
    normalize: bool
    stagger: bool
    use_zoom: bool
    zoom_start: Optional[float] = None
    zoom_end: Optional[float] = None


@dataclass(frozen=True)
class SaveOptions:
    output_dir: Path
    output_stem: str
    save_png: bool
    save_pdf: bool


class UserInputError(ValueError):
    """Raised when GUI input cannot be turned into a valid plot request."""


PLOT_CONFIGS = (
    PlotConfig(
        name="m/z Spectra",
        x_label=r"$m/z$",
        output_stem="mz_spectra_comparison",
        single_output_stem="mass_spectrum",
        zoom_start_label="Start m/z",
        zoom_end_label="End m/z",
        zoom_value_error="Enter numeric values for the start and end m/z values.",
        zoom_order_error="End m/z must be greater than start m/z.",
    ),
    PlotConfig(
        name="Mass Distribution",
        x_label=r"Mass (kDa)",
        output_stem="mass_distribution_comparison",
        single_output_stem="mass_distribution",
        zoom_start_label="Start mass (kDa)",
        zoom_end_label="End mass (kDa)",
        zoom_value_error="Enter numeric values for the start and end masses.",
        zoom_order_error="End mass must be greater than start mass.",
        x_scale=1 / 1000,
    ),
    PlotConfig(
        name="Total Ion Chromatogram",
        x_label=r"Time (minutes)",
        output_stem="total_ion_chromatogram_comparison",
        single_output_stem="total_ion_chromatogram",
        zoom_start_label="Start time (min)",
        zoom_end_label="End time (min)",
        zoom_value_error="Enter numeric values for the start and end times.",
        zoom_order_error="End time must be greater than start time.",
    ),
)
PLOT_CONFIG_BY_NAME = {config.name: config for config in PLOT_CONFIGS}


def clean_file_path(raw_path):
    return Path(raw_path.strip().strip('"').strip("'"))


def clean_output_stem(raw_stem):
    output_stem = raw_stem.strip()
    if not output_stem:
        return ""

    output_path = Path(output_stem)
    if output_path.suffix.lower() in (".png", ".pdf"):
        output_stem = output_path.stem

    invalid_characters = '<>:"/\\|?*'
    if any(character in output_stem for character in invalid_characters):
        raise UserInputError(f"Output file name cannot contain: {invalid_characters}")

    return output_stem


def compact_display_name(raw_path, max_length=36):
    name = Path(raw_path).name
    if len(name) <= max_length:
        return name

    head_length = (max_length - 3) // 2
    tail_length = max_length - 3 - head_length
    return f"{name[:head_length]}...{name[-tail_length:]}"


def detect_system_theme():
    if sys.platform.startswith("win"):
        try:
            import winreg

            key_path = r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize"
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, key_path) as key:
                value, _value_type = winreg.QueryValueEx(key, "AppsUseLightTheme")
            return "Light" if value else "Dark"
        except OSError:
            pass

    return "Light"


def load_dataset(path, dataset_number):
    """Load a UniDec-style text export with x values in column 1 and intensity in column 2."""
    data = np.loadtxt(path)

    if data.ndim == 1:
        if data.size == 2:
            data = data.reshape(1, 2)
        else:
            raise ValueError(
                f"Dataset {dataset_number} must have two columns: x-axis values and intensity."
            )

    if data.ndim != 2:
        raise ValueError(
            f"Dataset {dataset_number} must have two columns: x-axis values and intensity."
        )

    if data.shape[1] < 2:
        raise ValueError(
            f"Dataset {dataset_number} must have at least two columns: x-axis values and intensity."
        )

    return data


def configure_matplotlib():
    """Apply the shared visual style from the TUI plotter."""
    mpl.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": FONT_SIZE,
            "xtick.major.size": TICK_SIZE,
            "xtick.major.width": TICK_WIDTH,
            "ytick.major.size": TICK_SIZE,
            "ytick.major.width": TICK_WIDTH,
            "xtick.major.pad": 8,
            "ytick.major.pad": 8,
            "axes.labelpad": 18,
            "axes.titlepad": 18,
            "axes.linewidth": TICK_WIDTH,
            "pdf.fonttype": 42,
        }
    )


def get_x_values(dataset, config):
    # Mass distribution exports are in Da; plot them in kDa.
    return dataset.values[:, 0] * config.x_scale


def get_intensity_values(dataset):
    return dataset.values[:, 1].astype(float, copy=True)


def normalize_intensity(intensity, dataset_number):
    max_intensity = np.max(intensity)
    if max_intensity <= 0:
        raise ValueError(f"Cannot normalize dataset {dataset_number}; maximum intensity is {max_intensity}.")
    return intensity / max_intensity


def get_stagger_offset_step(all_intensities, normalize):
    if normalize:
        return STAGGER_STEP_NORMALIZED
    return STAGGER_STEP_FRACTION * np.max(all_intensities)


def is_single_dataset(datasets):
    return len(datasets) == 1


def plot_datasets(ax, config, datasets, options):
    all_raw_intensities = np.concatenate([get_intensity_values(dataset) for dataset in datasets])
    stagger_offset_step = get_stagger_offset_step(all_raw_intensities, options.normalize)
    plotted_intensities = []
    single_dataset = is_single_dataset(datasets)

    for dataset_number, dataset in enumerate(datasets, start=1):
        intensity = get_intensity_values(dataset)

        if options.normalize:
            intensity = normalize_intensity(intensity, dataset_number)
        if options.stagger:
            intensity = intensity + ((dataset_number - 1) * stagger_offset_step)

        plotted_intensities.append(intensity)
        if single_dataset:
            ax.plot(get_x_values(dataset, config), intensity, color="black")
        else:
            ax.plot(get_x_values(dataset, config), intensity, label=dataset.legend_label)

    return all_raw_intensities, plotted_intensities


def format_axes(ax, config, datasets, options, all_raw_intensities, plotted_intensities):
    ax.set_xlabel(config.x_label)

    y_label = "Normalized Intensity" if options.normalize else "Intensity"
    if options.stagger:
        y_label += " + Offset"
    ax.set_ylabel(y_label)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    set_y_limits(ax, options, all_raw_intensities, plotted_intensities)
    set_x_limits(ax, config, datasets, options)
    format_x_axis_without_scientific_notation(ax)


def set_y_limits(ax, options, all_raw_intensities, plotted_intensities):
    if options.stagger:
        y_max = np.max(np.concatenate(plotted_intensities))
    elif options.normalize:
        y_max = 1.0
    else:
        y_max = np.max(all_raw_intensities)

    ax.set_ylim(0, y_max)


def set_x_limits(ax, config, datasets, options):
    if options.use_zoom:
        ax.set_xlim(options.zoom_start, options.zoom_end)
        return

    all_x_values = np.concatenate([get_x_values(dataset, config) for dataset in datasets])
    ax.set_xlim(np.min(all_x_values), np.max(all_x_values))


def format_x_axis_without_scientific_notation(ax):
    formatter = ScalarFormatter(useOffset=False, useMathText=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)


def build_output_stem(config, datasets, options):
    if is_single_dataset(datasets):
        output_stem = config.single_output_stem
    else:
        output_stem = config.output_stem

    if options.normalize:
        output_stem += "_normalized"
    if options.stagger:
        output_stem += "_staggered"
    if options.use_zoom:
        output_stem += "_zoomed"
    return output_stem


def save_figure(figure, save_options):
    saved_paths = []

    if save_options.save_png:
        output_path_png = save_options.output_dir / f"{save_options.output_stem}.png"
        figure.savefig(output_path_png, dpi=300, transparent=True)
        saved_paths.append(output_path_png)

    if save_options.save_pdf:
        output_path_pdf = save_options.output_dir / f"{save_options.output_stem}.pdf"
        figure.savefig(output_path_pdf, dpi=300)
        saved_paths.append(output_path_pdf)

    return saved_paths


def render_plot_on_figure(figure, config, datasets, options):
    configure_matplotlib()
    figure.clear()

    ax = figure.add_subplot(111)
    all_raw_intensities, plotted_intensities = plot_datasets(ax, config, datasets, options)
    format_axes(ax, config, datasets, options, all_raw_intensities, plotted_intensities)

    if not is_single_dataset(datasets):
        ax.legend()

    figure.tight_layout()
    return ax


def build_fixed_output_figure(config, datasets, options):
    """Render a fixed-size figure so saved files do not depend on GUI preview size."""
    figure = Figure(figsize=FIGURE_SIZE, dpi=100, facecolor=THEMES["Light"]["plot_bg"])
    render_plot_on_figure(figure, config, datasets, options)
    return figure


class DatasetRow:
    def __init__(self, app, parent, dataset_input):
        self.app = app
        self.dataset_input = dataset_input
        self.frame = ttk.Frame(parent, style="Dataset.TFrame", padding=(0, 2))
        self.path_var = tk.StringVar(value=dataset_input.path)
        self.label_var = tk.StringVar(value=dataset_input.legend_label)

        self.index_label = ttk.Label(self.frame, width=3, anchor="center", style="Index.TLabel")
        self.path_entry = ttk.Entry(self.frame, textvariable=self.path_var)
        self.browse_button = ttk.Button(self.frame, text="Browse", style="Dataset.TButton", command=self.browse)
        self.label_entry = ttk.Entry(self.frame, textvariable=self.label_var)
        self.remove_button = ttk.Button(self.frame, text="Remove", style="Dataset.TButton", command=self.remove)

        self.index_label.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.path_entry.grid(row=0, column=1, sticky="ew", padx=(0, 8))
        self.browse_button.grid(row=0, column=2, sticky="ew", padx=(0, 8))
        self.label_entry.grid(row=0, column=3, sticky="ew", padx=(0, 8))
        self.remove_button.grid(row=0, column=4, sticky="ew")

        self.frame.columnconfigure(1, weight=1)
        self.frame.columnconfigure(3, weight=1)

        self.path_var.trace_add("write", self.on_field_changed)
        self.label_var.trace_add("write", self.on_field_changed)

    def on_field_changed(self, *_args):
        self.sync_to_input()
        self.app.update_dataset_summary()

    def browse(self):
        self.app.browse_for_row(self)

    def clear(self):
        self.path_var.set("")
        self.label_var.set("")
        self.sync_to_input()
        self.app.update_dataset_summary()

    def sync_to_input(self):
        self.dataset_input.path = self.path_var.get()
        self.dataset_input.legend_label = self.label_var.get()

    def grid(self, row):
        self.frame.grid(row=row, column=0, sticky="ew", pady=(0, 4))

    def remove(self):
        self.app.remove_dataset_row(self)

    def destroy(self):
        self.frame.destroy()

    def update_index(self, index, labels_enabled):
        self.index_label.configure(text=str(index))
        self.label_entry.configure(state="normal" if labels_enabled else "disabled")
        self.remove_button.configure(text="Clear" if len(self.app.dataset_rows) == 1 else "Remove")


class UniDecDataPlotterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("UniDec Data Plotter")
        self.root.geometry("1280x780")
        self.root.minsize(1060, 680)

        self.style = ttk.Style(self.root)
        self.style.theme_use("clam")
        self.dataset_inputs = [DatasetInput()]
        self.dataset_rows = []
        self.dataset_editor_window = None
        self.plot_type_var = tk.StringVar(value=PLOT_CONFIGS[0].name)
        self.theme_var = tk.StringVar(value=detect_system_theme())
        self.normalize_var = tk.BooleanVar(value=True)
        self.stagger_var = tk.BooleanVar(value=False)
        self.zoom_var = tk.BooleanVar(value=False)
        self.zoom_start_var = tk.StringVar()
        self.zoom_end_var = tk.StringVar()
        self.zoom_start_label_var = tk.StringVar()
        self.zoom_end_label_var = tk.StringVar()
        self.output_dir_var = tk.StringVar()
        self.output_stem_var = tk.StringVar()
        self.save_png_var = tk.BooleanVar(value=True)
        self.save_pdf_var = tk.BooleanVar(value=True)
        self.dataset_summary_var = tk.StringVar()
        self.dataset_detail_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Ready.")

        self.configure_ttk_style()
        self.build_layout()
        self.update_dataset_summary()
        self.update_zoom_labels()
        self.update_zoom_state()
        self.draw_empty_preview()

    def configure_ttk_style(self):
        colors = THEMES[self.theme_var.get()]
        self.root.configure(bg=colors["app_bg"])

        default_font = ("Segoe UI", 10)
        heading_font = ("Segoe UI Semibold", 11)
        title_font = ("Segoe UI Semibold", 18)

        self.style.configure(
            ".",
            font=default_font,
            background=colors["app_bg"],
            foreground=colors["text_primary"],
        )
        self.style.configure("App.TFrame", background=colors["app_bg"])
        self.style.configure(
            "Panel.TFrame",
            background=colors["panel_bg"],
            bordercolor=colors["border"],
            borderwidth=1,
            relief="solid",
        )
        self.style.configure("Surface.TFrame", background=colors["panel_bg"])
        self.style.configure("Dataset.TFrame", background=colors["panel_bg"])
        self.style.configure(
            "Header.TLabel",
            background=colors["panel_bg"],
            foreground=colors["text_primary"],
            font=title_font,
        )
        self.style.configure(
            "Section.TLabel",
            background=colors["panel_bg"],
            foreground=colors["text_primary"],
            font=heading_font,
        )
        self.style.configure(
            "Field.TLabel",
            background=colors["panel_bg"],
            foreground=colors["text_secondary"],
        )
        self.style.configure(
            "Index.TLabel",
            background=colors["panel_bg"],
            foreground=colors["text_secondary"],
        )
        self.style.configure(
            "Status.TLabel",
            background=colors["app_bg"],
            foreground=colors["text_secondary"],
        )
        self.style.configure(
            "TCheckbutton",
            background=colors["panel_bg"],
            foreground=colors["text_primary"],
        )
        self.style.map(
            "TCheckbutton",
            foreground=[("disabled", colors["disabled_text"])],
            background=[("active", colors["panel_bg"])],
        )
        self.style.configure(
            "TEntry",
            fieldbackground=colors["input_bg"],
            foreground=colors["text_primary"],
            bordercolor=colors["border"],
            lightcolor=colors["border"],
            darkcolor=colors["border"],
            insertcolor=colors["text_primary"],
        )
        self.style.map(
            "TEntry",
            fieldbackground=[("disabled", colors["app_bg"])],
            foreground=[("disabled", colors["disabled_text"])],
        )
        self.style.configure(
            "TCombobox",
            fieldbackground=colors["input_bg"],
            foreground=colors["text_primary"],
            background=colors["input_bg"],
            bordercolor=colors["border"],
            lightcolor=colors["border"],
            darkcolor=colors["border"],
            arrowcolor=colors["text_primary"],
        )
        self.style.map(
            "TCombobox",
            fieldbackground=[("readonly", colors["input_bg"])],
            foreground=[("readonly", colors["text_primary"])],
            selectbackground=[("readonly", colors["input_bg"])],
            selectforeground=[("readonly", colors["text_primary"])],
        )
        self.style.configure(
            "TButton",
            padding=(10, 6),
            background=colors["button_bg"],
            foreground=colors["text_primary"],
            bordercolor=colors["border"],
        )
        self.style.configure(
            "Dataset.TButton",
            padding=(8, 3),
            background=colors["button_bg"],
            foreground=colors["text_primary"],
            bordercolor=colors["border"],
        )
        self.style.map(
            "TButton",
            background=[("active", colors["button_active"]), ("disabled", colors["app_bg"])],
            foreground=[("disabled", colors["disabled_text"])],
        )
        self.style.map(
            "Dataset.TButton",
            background=[("active", colors["button_active"]), ("disabled", colors["app_bg"])],
            foreground=[("disabled", colors["disabled_text"])],
        )
        self.style.configure(
            "Primary.TButton",
            padding=(14, 8),
            background=colors["accent"],
            foreground="#ffffff",
            bordercolor=colors["accent"],
        )
        self.style.map(
            "Primary.TButton",
            background=[("active", colors["accent_active"]), ("disabled", colors["button_bg"])],
            foreground=[("disabled", colors["disabled_text"])],
        )

    def build_layout(self):
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        main = ttk.Frame(self.root, style="App.TFrame", padding=18)
        main.grid(row=0, column=0, sticky="nsew")
        main.columnconfigure(0, minsize=520)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        controls = ttk.Frame(main, style="Panel.TFrame", padding=18)
        controls.grid(row=0, column=0, sticky="nsew", padx=(0, 16))
        controls.columnconfigure(0, weight=1)

        preview = ttk.Frame(main, style="Panel.TFrame", padding=12)
        preview.grid(row=0, column=1, sticky="nsew")
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(1, weight=1)

        self.build_controls(controls)
        self.build_preview(preview)

        status = ttk.Label(main, textvariable=self.status_var, style="Status.TLabel", anchor="w")
        status.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(12, 0))

    def build_controls(self, parent):
        header = ttk.Frame(parent, style="Surface.TFrame")
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(0, weight=1)

        ttk.Label(header, text="UniDec Data Plotter", style="Header.TLabel").grid(
            row=0, column=0, sticky="w"
        )
        ttk.Label(header, text="Theme", style="Field.TLabel").grid(row=0, column=1, sticky="e", padx=(12, 8))
        self.theme_combo = ttk.Combobox(
            header,
            textvariable=self.theme_var,
            values=list(THEMES.keys()),
            state="readonly",
            width=8,
        )
        self.theme_combo.grid(row=0, column=2, sticky="e")
        self.theme_combo.bind("<<ComboboxSelected>>", self.on_theme_changed)

        plot_type_frame = ttk.Frame(parent, style="Surface.TFrame")
        plot_type_frame.grid(row=1, column=0, sticky="ew", pady=(18, 0))
        plot_type_frame.columnconfigure(0, weight=1)

        ttk.Label(plot_type_frame, text="Plot type", style="Field.TLabel").grid(row=0, column=0, sticky="w")
        self.plot_type_combo = ttk.Combobox(
            plot_type_frame,
            textvariable=self.plot_type_var,
            values=[config.name for config in PLOT_CONFIGS],
            state="readonly",
        )
        self.plot_type_combo.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        self.plot_type_combo.bind("<<ComboboxSelected>>", self.on_plot_type_changed)

        datasets = ttk.Frame(parent, style="Surface.TFrame")
        datasets.grid(row=2, column=0, sticky="ew", pady=(18, 0))
        datasets.columnconfigure(0, weight=1)
        datasets.columnconfigure(1, weight=0)

        ttk.Label(datasets, text="Datasets", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Button(datasets, text="Manage Datasets", command=self.open_dataset_editor).grid(
            row=0,
            column=1,
            padx=(8, 0),
        )
        ttk.Label(datasets, textvariable=self.dataset_summary_var, style="Field.TLabel").grid(
            row=1,
            column=0,
            columnspan=2,
            sticky="w",
            pady=(8, 0),
        )
        ttk.Label(datasets, textvariable=self.dataset_detail_var, style="Field.TLabel", wraplength=460).grid(
            row=2,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(3, 0),
        )

        self.options_frame = ttk.Frame(parent, style="Surface.TFrame")
        options = self.options_frame
        options.grid(row=3, column=0, sticky="ew", pady=(16, 0))
        options.columnconfigure(0, weight=1)
        options.columnconfigure(1, weight=1)

        ttk.Label(options, text="Options", style="Section.TLabel").grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 8)
        )
        ttk.Checkbutton(options, text="Normalize intensity", variable=self.normalize_var).grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(0, 6)
        )
        self.stagger_check = ttk.Checkbutton(options, text="Stagger overlay traces", variable=self.stagger_var)
        self.stagger_check.grid(row=2, column=0, columnspan=2, sticky="w", pady=(0, 6))
        ttk.Checkbutton(
            options,
            text="Use zoomed x-axis range",
            variable=self.zoom_var,
            command=self.update_zoom_state,
        ).grid(row=3, column=0, columnspan=2, sticky="w", pady=(0, 10))

        self.zoom_start_label = ttk.Label(options, textvariable=self.zoom_start_label_var, style="Field.TLabel")
        self.zoom_start_label.grid(row=4, column=0, sticky="w", padx=(0, 8))
        self.zoom_end_label = ttk.Label(options, textvariable=self.zoom_end_label_var, style="Field.TLabel")
        self.zoom_end_label.grid(row=4, column=1, sticky="w")

        self.zoom_start_entry = ttk.Entry(options, textvariable=self.zoom_start_var)
        self.zoom_start_entry.grid(row=5, column=0, sticky="ew", padx=(0, 8), pady=(5, 0))
        self.zoom_end_entry = ttk.Entry(options, textvariable=self.zoom_end_var)
        self.zoom_end_entry.grid(row=5, column=1, sticky="ew", pady=(5, 0))

        self.output_frame = ttk.Frame(parent, style="Surface.TFrame")
        output = self.output_frame
        output.grid(row=4, column=0, sticky="ew", pady=(14, 0))
        output.columnconfigure(0, weight=1)

        ttk.Label(output, text="Output", style="Section.TLabel").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Label(output, text="Output folder", style="Field.TLabel").grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(8, 0)
        )
        self.output_dir_entry = ttk.Entry(output, textvariable=self.output_dir_var)
        self.output_dir_entry.grid(row=2, column=0, columnspan=2, sticky="ew", padx=(0, 8), pady=(5, 0))
        ttk.Button(output, text="Browse", command=self.browse_output_folder).grid(
            row=2, column=2, sticky="ew", pady=(5, 0)
        )

        ttk.Label(output, text="Base filename (blank uses default)", style="Field.TLabel").grid(
            row=3, column=0, columnspan=3, sticky="w", pady=(10, 0)
        )
        self.output_stem_entry = ttk.Entry(output, textvariable=self.output_stem_var)
        self.output_stem_entry.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(5, 0))

        ttk.Checkbutton(output, text="PNG", variable=self.save_png_var).grid(
            row=5, column=0, sticky="w", pady=(10, 0)
        )
        ttk.Checkbutton(output, text="PDF", variable=self.save_pdf_var).grid(
            row=5, column=1, sticky="w", pady=(10, 0)
        )

        self.actions_frame = ttk.Frame(parent, style="Surface.TFrame")
        actions = self.actions_frame
        actions.grid(row=5, column=0, sticky="ew", pady=(16, 0))
        actions.columnconfigure(0, weight=1)
        actions.columnconfigure(1, weight=1)

        ttk.Button(actions, text="Preview Plot", command=self.preview_plot).grid(
            row=0, column=0, sticky="ew", padx=(0, 8)
        )
        ttk.Button(actions, text="Save Plot", style="Primary.TButton", command=self.save_plot).grid(
            row=0, column=1, sticky="ew"
        )

    def build_preview(self, parent):
        ttk.Label(parent, text="Preview", style="Section.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 8))

        colors = THEMES[self.theme_var.get()]
        self.figure = Figure(figsize=FIGURE_SIZE, dpi=100, facecolor=colors["plot_bg"])
        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.configure(bg=colors["panel_bg"], highlightbackground=colors["panel_bg"])
        self.canvas_widget.grid(row=1, column=0, sticky="nsew")

        toolbar_frame = ttk.Frame(parent, style="Surface.TFrame")
        toolbar_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.LEFT)

    def draw_empty_preview(self):
        colors = THEMES[self.theme_var.get()]
        self.figure.clear()
        self.figure.set_facecolor(colors["plot_bg"])
        ax = self.figure.add_subplot(111)
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "Add a dataset and preview the plot.",
            ha="center",
            va="center",
            color=colors["empty_plot_text"],
            fontsize=16,
            transform=ax.transAxes,
        )
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def get_selected_config(self):
        return PLOT_CONFIG_BY_NAME[self.plot_type_var.get()]

    def on_plot_type_changed(self, _event=None):
        self.update_zoom_labels()
        self.status_var.set(f"Plot type set to {self.plot_type_var.get()}.")

    def on_theme_changed(self, _event=None):
        self.configure_ttk_style()
        self.update_canvas_theme()
        self.status_var.set(f"{self.theme_var.get()} theme applied.")

    def update_canvas_theme(self):
        colors = THEMES[self.theme_var.get()]
        if self.dataset_editor_is_open():
            self.dataset_editor_window.configure(bg=colors["app_bg"])
            self.dataset_rows_canvas.configure(
                bg=colors["panel_bg"],
                highlightbackground=colors["border"],
            )
        if hasattr(self, "canvas_widget"):
            self.canvas_widget.configure(bg=colors["panel_bg"], highlightbackground=colors["panel_bg"])

        if hasattr(self, "figure") and self.figure.axes and not self.figure.axes[0].lines:
            self.draw_empty_preview()

    def dataset_editor_is_open(self):
        return self.dataset_editor_window is not None and self.dataset_editor_window.winfo_exists()

    def open_dataset_editor(self):
        if self.dataset_editor_is_open():
            self.sync_dataset_editor_rows()
            self.refresh_dataset_rows()
            self.dataset_editor_window.lift()
            self.dataset_editor_window.focus_force()
            return

        self.sync_dataset_editor_rows()
        colors = THEMES[self.theme_var.get()]
        self.dataset_editor_window = tk.Toplevel(self.root)
        self.dataset_editor_window.title("Manage Datasets")
        self.dataset_editor_window.geometry("920x420")
        self.dataset_editor_window.minsize(760, 320)
        self.dataset_editor_window.configure(bg=colors["app_bg"])
        self.dataset_editor_window.transient(self.root)
        self.dataset_editor_window.protocol("WM_DELETE_WINDOW", self.close_dataset_editor)
        self.dataset_editor_window.columnconfigure(0, weight=1)
        self.dataset_editor_window.rowconfigure(1, weight=1)

        header = ttk.Frame(self.dataset_editor_window, style="Panel.TFrame", padding=(14, 12))
        header.grid(row=0, column=0, sticky="ew", padx=14, pady=(14, 10))
        header.columnconfigure(0, weight=1)
        ttk.Label(header, text="Datasets", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Button(header, text="Add Files", command=self.add_files).grid(row=0, column=1, padx=(8, 0))
        ttk.Button(header, text="Add Row", command=self.add_dataset_row).grid(row=0, column=2, padx=(8, 0))
        ttk.Button(header, text="Close", command=self.close_dataset_editor).grid(row=0, column=3, padx=(8, 0))

        list_panel = ttk.Frame(self.dataset_editor_window, style="Panel.TFrame", padding=14)
        list_panel.grid(row=1, column=0, sticky="nsew", padx=14, pady=(0, 14))
        list_panel.columnconfigure(0, weight=1)
        list_panel.rowconfigure(1, weight=1)

        column_labels = ttk.Frame(list_panel, style="Surface.TFrame")
        column_labels.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        column_labels.columnconfigure(1, weight=1)
        column_labels.columnconfigure(3, weight=1)
        ttk.Label(column_labels, text="#", width=3, anchor="center", style="Field.TLabel").grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 8),
        )
        ttk.Label(column_labels, text="File path", style="Field.TLabel").grid(
            row=0,
            column=1,
            sticky="w",
            padx=(0, 8),
        )
        ttk.Label(column_labels, text="", style="Field.TLabel").grid(row=0, column=2, padx=(0, 8))
        ttk.Label(column_labels, text="Legend label", style="Field.TLabel").grid(
            row=0,
            column=3,
            sticky="w",
            padx=(0, 8),
        )

        dataset_list = ttk.Frame(list_panel, style="Surface.TFrame")
        dataset_list.grid(row=1, column=0, sticky="nsew")
        dataset_list.columnconfigure(0, weight=1)
        dataset_list.rowconfigure(0, weight=1)

        self.dataset_rows_canvas = tk.Canvas(
            dataset_list,
            bg=colors["panel_bg"],
            borderwidth=0,
            highlightthickness=1,
            highlightbackground=colors["border"],
        )
        self.dataset_rows_scrollbar = ttk.Scrollbar(
            dataset_list,
            orient="vertical",
            command=self.dataset_rows_canvas.yview,
        )
        self.dataset_rows_frame = ttk.Frame(self.dataset_rows_canvas, style="Surface.TFrame")
        self.dataset_rows_window = self.dataset_rows_canvas.create_window(
            (0, 0),
            window=self.dataset_rows_frame,
            anchor="nw",
        )

        self.dataset_rows_canvas.configure(yscrollcommand=self.dataset_rows_scrollbar.set)
        self.dataset_rows_canvas.grid(row=0, column=0, sticky="nsew")
        self.dataset_rows_scrollbar.grid(row=0, column=1, sticky="ns", padx=(8, 0))
        self.dataset_rows_frame.columnconfigure(0, weight=1)

        self.dataset_rows_frame.bind("<Configure>", lambda _event: self.refresh_dataset_viewport())
        self.dataset_rows_canvas.bind(
            "<Configure>",
            lambda event: self.dataset_rows_canvas.itemconfigure(self.dataset_rows_window, width=event.width),
        )
        self.bind_mousewheel(self.dataset_rows_canvas)
        self.refresh_dataset_rows()

    def bind_mousewheel(self, canvas):
        canvas.bind(
            "<Enter>",
            lambda _event: self.root.bind_all(
                "<MouseWheel>",
                lambda event, target=canvas: self.scroll_canvas(event, target),
            ),
        )
        canvas.bind("<Leave>", lambda _event: self.root.unbind_all("<MouseWheel>"))

    def scroll_canvas(self, event, canvas):
        bbox = canvas.bbox("all")
        if not bbox or bbox[3] - bbox[1] <= canvas.winfo_height():
            return "break"

        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        return "break"

    def refresh_dataset_viewport(self):
        if not self.dataset_editor_is_open() or not hasattr(self, "dataset_rows_canvas"):
            return

        bbox = self.dataset_rows_canvas.bbox("all")
        content_height = bbox[3] - bbox[1] if bbox else 0
        viewport_height = self.dataset_rows_canvas.winfo_height()
        self.dataset_rows_canvas.configure(scrollregion=bbox)

        if content_height > viewport_height:
            self.dataset_rows_scrollbar.grid()
        else:
            self.dataset_rows_scrollbar.grid_remove()
            self.dataset_rows_canvas.yview_moveto(0.0)

    def update_zoom_labels(self):
        config = self.get_selected_config()
        self.zoom_start_label_var.set(config.zoom_start_label)
        self.zoom_end_label_var.set(config.zoom_end_label)

    def update_zoom_state(self):
        state = "normal" if self.zoom_var.get() else "disabled"
        self.zoom_start_entry.configure(state=state)
        self.zoom_end_entry.configure(state=state)

    def get_initial_output_dir(self):
        output_dir = self.output_dir_var.get().strip()
        if output_dir:
            path = clean_file_path(output_dir)
            if path.exists():
                return path

        self.sync_dataset_editor_rows()
        for dataset_input in self.dataset_inputs:
            file_path = dataset_input.path.strip()
            if file_path:
                path = clean_file_path(file_path)
                if path.parent.exists():
                    return path.parent

        return Path.home()

    def browse_for_row(self, row):
        path = filedialog.askopenfilename(
            title="Select UniDec data file",
            filetypes=(
                ("Text data files", "*.txt *.dat"),
                ("All files", "*.*"),
            ),
        )
        if not path:
            return

        selected_path = Path(path)
        row.path_var.set(str(selected_path))
        if not row.label_var.get().strip():
            row.label_var.set(selected_path.stem)

        row.sync_to_input()
        self.update_dataset_summary()
        self.status_var.set(f"Selected {selected_path.name}.")

    def browse_output_folder(self):
        output_dir = filedialog.askdirectory(
            title="Select output folder",
            initialdir=str(self.get_initial_output_dir()),
        )
        if not output_dir:
            return

        selected_dir = Path(output_dir)
        self.output_dir_var.set(str(selected_dir))
        self.status_var.set(f"Output folder set to {selected_dir}.")

    def add_files(self):
        self.sync_dataset_editor_rows()
        paths = filedialog.askopenfilenames(
            title="Select UniDec data files",
            filetypes=(
                ("Text data files", "*.txt *.dat"),
                ("All files", "*.*"),
            ),
        )
        if not paths:
            return

        blank_inputs = [
            dataset_input
            for dataset_input in self.dataset_inputs
            if not dataset_input.path.strip()
        ]
        for raw_path in paths:
            path = Path(raw_path)
            if blank_inputs:
                dataset_input = blank_inputs.pop(0)
            else:
                dataset_input = DatasetInput()
                self.dataset_inputs.append(dataset_input)

            dataset_input.path = str(path)
            if not dataset_input.legend_label.strip():
                dataset_input.legend_label = path.stem

        self.update_dataset_summary()
        self.refresh_dataset_rows()
        if self.dataset_editor_is_open():
            self.dataset_editor_window.after_idle(self.scroll_dataset_rows_to_bottom)
        self.status_var.set(f"Added {len(paths)} dataset file(s).")

    def add_dataset_row(self):
        self.sync_dataset_editor_rows()
        self.dataset_inputs.append(DatasetInput())
        self.refresh_dataset_rows()
        if self.dataset_editor_is_open():
            self.dataset_editor_window.after_idle(self.scroll_dataset_rows_to_bottom)
            return self.dataset_rows[-1]
        return None

    def remove_dataset_row(self, row):
        self.sync_dataset_editor_rows()
        if len(self.dataset_inputs) == 1:
            row.clear()
        else:
            if row.dataset_input in self.dataset_inputs:
                self.dataset_inputs.remove(row.dataset_input)
            row.destroy()

        self.refresh_dataset_rows()
        self.update_dataset_summary()

    def refresh_dataset_rows(self):
        if not self.dataset_inputs:
            self.dataset_inputs.append(DatasetInput())

        if not self.dataset_editor_is_open():
            self.update_dataset_summary()
            return

        for row in self.dataset_rows:
            row.destroy()

        self.dataset_rows = []
        labels_enabled = len(self.dataset_inputs) > 1

        for index, dataset_input in enumerate(self.dataset_inputs, start=1):
            row = DatasetRow(self, self.dataset_rows_frame, dataset_input)
            self.dataset_rows.append(row)
            row.grid(index - 1)
            row.update_index(index, labels_enabled)

        self.update_dataset_summary()
        self.dataset_rows_frame.update_idletasks()
        self.refresh_dataset_viewport()

    def sync_dataset_editor_rows(self):
        for row in self.dataset_rows:
            row.sync_to_input()

    def close_dataset_editor(self):
        self.sync_dataset_editor_rows()
        self.update_dataset_summary()
        if self.dataset_editor_is_open():
            self.dataset_editor_window.destroy()

        self.dataset_editor_window = None
        self.dataset_rows = []

    def update_dataset_summary(self):
        active_inputs = [
            dataset_input
            for dataset_input in self.dataset_inputs
            if dataset_input.path.strip()
        ]
        dataset_count = len(active_inputs)

        if dataset_count == 0:
            self.dataset_summary_var.set("No datasets selected.")
            self.dataset_detail_var.set("Use Manage Datasets to choose UniDec text exports.")
        elif dataset_count == 1:
            self.dataset_summary_var.set("1 dataset selected.")
            self.dataset_detail_var.set(compact_display_name(active_inputs[0].path, max_length=54))
        else:
            first_name = compact_display_name(active_inputs[0].path, max_length=38)
            remaining_count = dataset_count - 1
            self.dataset_summary_var.set(f"{dataset_count} datasets selected.")
            self.dataset_detail_var.set(f"{first_name}, +{remaining_count} more")

        if not hasattr(self, "stagger_check"):
            return

        if dataset_count > 1:
            self.stagger_check.configure(state="normal")
        else:
            self.stagger_var.set(False)
            self.stagger_check.configure(state="disabled")

    def scroll_dataset_rows_to_bottom(self):
        if not self.dataset_editor_is_open():
            return

        self.refresh_dataset_viewport()
        bbox = self.dataset_rows_canvas.bbox("all")
        if bbox and bbox[3] - bbox[1] > self.dataset_rows_canvas.winfo_height():
            self.dataset_rows_canvas.yview_moveto(1.0)

    def collect_inputs(self):
        config = self.get_selected_config()
        datasets = self.collect_datasets()
        options = self.collect_plot_options(config, is_single_dataset(datasets))
        return config, datasets, options

    def collect_datasets(self):
        self.sync_dataset_editor_rows()
        active_inputs = [
            dataset_input
            for dataset_input in self.dataset_inputs
            if dataset_input.path.strip()
        ]
        if not active_inputs:
            raise UserInputError("Add at least one dataset file.")

        multiple_datasets = len(active_inputs) > 1
        datasets = []
        for index, dataset_input in enumerate(active_inputs, start=1):
            path = clean_file_path(dataset_input.path)
            if not path.exists():
                raise UserInputError(f"Dataset {index} file does not exist:\n{path}")

            values = load_dataset(path, index)
            legend_label = None
            if multiple_datasets:
                legend_label = dataset_input.legend_label.strip() or path.stem

            datasets.append(Dataset(path=path, values=values, legend_label=legend_label))

        return datasets

    def collect_plot_options(self, config, single_dataset):
        use_zoom = self.zoom_var.get()
        zoom_start = None
        zoom_end = None

        if use_zoom:
            try:
                zoom_start = float(self.zoom_start_var.get())
                zoom_end = float(self.zoom_end_var.get())
            except ValueError as exc:
                raise UserInputError(config.zoom_value_error) from exc

            if zoom_end <= zoom_start:
                raise UserInputError(config.zoom_order_error)

        return PlotOptions(
            normalize=self.normalize_var.get(),
            stagger=False if single_dataset else self.stagger_var.get(),
            use_zoom=use_zoom,
            zoom_start=zoom_start,
            zoom_end=zoom_end,
        )

    def collect_save_options(self, config, datasets, options):
        save_png = self.save_png_var.get()
        save_pdf = self.save_pdf_var.get()
        if not save_png and not save_pdf:
            raise UserInputError("Select PNG, PDF, or both before saving.")

        output_dir_text = self.output_dir_var.get().strip()
        if not output_dir_text:
            raise UserInputError("Select an output folder before saving.")

        output_dir = clean_file_path(output_dir_text)
        if not output_dir.exists():
            raise UserInputError(f"Output folder does not exist:\n{output_dir}")
        if not output_dir.is_dir():
            raise UserInputError(f"Output path is not a folder:\n{output_dir}")

        output_stem = clean_output_stem(self.output_stem_var.get())
        if not output_stem:
            output_stem = build_output_stem(config, datasets, options)

        return SaveOptions(
            output_dir=output_dir,
            output_stem=output_stem,
            save_png=save_png,
            save_pdf=save_pdf,
        )

    def preview_plot(self):
        try:
            config, datasets, options = self.collect_inputs()
            render_plot_on_figure(self.figure, config, datasets, options)
            self.canvas.draw_idle()
            self.status_var.set("Preview updated.")
        except Exception as exc:
            self.show_error("Cannot create plot", exc)

    def save_plot(self):
        try:
            config, datasets, options = self.collect_inputs()
            save_options = self.collect_save_options(config, datasets, options)
            output_figure = build_fixed_output_figure(config, datasets, options)
            saved_paths = save_figure(output_figure, save_options)
            render_plot_on_figure(self.figure, config, datasets, options)
            self.canvas.draw_idle()
            saved_names = " and ".join(path.name for path in saved_paths)
            saved_message = "\n".join(str(path) for path in saved_paths)
            self.status_var.set(f"Saved: {saved_names}")
            messagebox.showinfo(
                "Plot saved",
                f"Saved plot file(s):\n{saved_message}",
            )
        except Exception as exc:
            self.show_error("Cannot save plot", exc)

    def show_error(self, title, error):
        self.status_var.set(str(error))
        messagebox.showerror(title, str(error))


def main():
    root = tk.Tk()
    UniDecDataPlotterApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()

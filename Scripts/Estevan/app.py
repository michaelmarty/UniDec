# app.py
from __future__ import annotations
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Sequence, List
import os

# Import your plotting utilities
# Ensure these modules are in the same folder or on PYTHONPATH
try:
    from chromatogram_plot_utils import plot_sec_stack
except ImportError:
    plot_sec_stack = None

try:
    from cdms_mass_plot_utils import plot_mass_offset_stack, plot_single_mass
except ImportError:
    plot_mass_offset_stack = None
    plot_single_mass = None

DEFAULT_COLORS = [
    "#333333", "#00A9A5", "#C76BA0", "#C76BA0",
    "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5",
]

DEFAULT_MASS_COLORS = ["#CC79A7", "#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"]

PLOT_TYPES = [
    "Single SEC",
    "Stacked SEC",
    "Single Mass",
    "Stacked Mass",
]

def parse_float_list(text: str) -> List[float]:
    text = text.strip()
    if not text:
        return []
    out = []
    for part in text.replace(";", ",").split(","):
        try:
            out.append(float(part.strip()))
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

class PlotGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SEC & Mass Plotter")
        self.geometry("920x720")

        # State
        self.files: List[str] = []
        self.plot_type = tk.StringVar(value=PLOT_TYPES[0])

        # Build UI
        self._build_header()
        self._build_file_section()
        self._build_options_section()
        self._build_sec_options()
        self._build_mass_options()
        self._build_action_buttons()
        self._update_visibility()

    def _build_header(self):
        header = ttk.Label(self, text="Plotting GUI: SEC chromatograms and mass distributions", font=("Arial", 14, "bold"))
        header.pack(pady=10)

        # Plot type selection
        frame = ttk.Frame(self)
        frame.pack(fill="x", padx=12, pady=6)
        ttk.Label(frame, text="Plot Type:", font=("Arial", 11, "bold")).pack(side="left")
        combo = ttk.Combobox(frame, values=PLOT_TYPES, textvariable=self.plot_type, state="readonly")
        combo.pack(side="left", padx=10)
        combo.bind("<<ComboboxSelected>>", lambda e: self._update_visibility())

    def _build_file_section(self):
        frame = ttk.LabelFrame(self, text="Files")
        frame.pack(fill="x", padx=12, pady=6)

        # Buttons for selecting files
        btns = ttk.Frame(frame)
        btns.pack(fill="x", pady=4)
        ttk.Button(btns, text="Add Files...", command=self._add_files).pack(side="left")
        ttk.Button(btns, text="Clear Files", command=self._clear_files).pack(side="left", padx=6)

        # File listbox
        self.file_list = tk.Listbox(frame, height=8)
        self.file_list.pack(fill="both", padx=6, pady=6)

    def _build_options_section(self):
        frame = ttk.LabelFrame(self, text="Common Options")
        frame.pack(fill="x", padx=12, pady=6)

        # Save path
        path_frame = ttk.Frame(frame)
        path_frame.pack(fill="x", pady=4)
        ttk.Label(path_frame, text="Save as (optional):").pack(side="left")
        self.save_path = tk.StringVar()
        ttk.Entry(path_frame, textvariable=self.save_path, width=40).pack(side="left", padx=6)
        ttk.Button(path_frame, text="Browse...", command=self._browse_save).pack(side="left")

        # Show checkbox
        self.show_plot = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, text="Show Figure", variable=self.show_plot).pack(anchor="w", pady=2)

        # Font settings
        font_frame = ttk.Frame(frame)
        font_frame.pack(fill="x", pady=4)
        ttk.Label(font_frame, text="Font Name:").pack(side="left")
        self.font_name = tk.StringVar(value="Arial")
        ttk.Entry(font_frame, textvariable=self.font_name, width=16).pack(side="left", padx=6)

        ttk.Label(font_frame, text="Font Size:").pack(side="left")
        self.font_size = tk.IntVar(value=10)
        ttk.Entry(font_frame, textvariable=self.font_size, width=6).pack(side="left", padx=6)

        # DPI and figsize
        fig_frame = ttk.Frame(frame)
        fig_frame.pack(fill="x", pady=4)
        ttk.Label(fig_frame, text="Figure Size (W,H):").pack(side="left")
        self.fig_w = tk.DoubleVar(value=9.0)
        self.fig_h = tk.DoubleVar(value=6.0)
        ttk.Entry(fig_frame, textvariable=self.fig_w, width=6).pack(side="left", padx=3)
        ttk.Entry(fig_frame, textvariable=self.fig_h, width=6).pack(side="left", padx=3)
        ttk.Label(fig_frame, text="DPI:").pack(side="left", padx=6)
        self.dpi = tk.IntVar(value=300)
        ttk.Entry(fig_frame, textvariable=self.dpi, width=8).pack(side="left")

    def _build_sec_options(self):
        # SEC options frame
        self.sec_frame = ttk.LabelFrame(self, text="SEC Options")
        self.sec_frame.pack(fill="x", padx=12, pady=6)

        # Labels
        lab_frame = ttk.Frame(self.sec_frame)
        lab_frame.pack(fill="x", pady=4)
        ttk.Label(lab_frame, text="Sample Labels (comma-separated; leave blank to derive from filenames):").pack(anchor="w")
        self.sec_labels_text = tk.Text(lab_frame, height=2)
        self.sec_labels_text.pack(fill="x", padx=6, pady=4)

        # Colors
        color_frame = ttk.Frame(self.sec_frame)
        color_frame.pack(fill="x", pady=4)
        ttk.Label(color_frame, text="Colors (comma-separated hex):").pack(anchor="w")
        self.sec_colors = tk.StringVar(value=", ".join(DEFAULT_COLORS))
        ttk.Entry(color_frame, textvariable=self.sec_colors).pack(fill="x", padx=6, pady=4)

        # Limits
        limits_frame = ttk.Frame(self.sec_frame)
        limits_frame.pack(fill="x", pady=4)
        ttk.Label(limits_frame, text="X-lim (mL):").pack(side="left")
        self.sec_xmin = tk.DoubleVar(value=10.0)
        self.sec_xmax = tk.DoubleVar(value=20.0)
        ttk.Entry(limits_frame, textvariable=self.sec_xmin, width=8).pack(side="left", padx=3)
        ttk.Entry(limits_frame, textvariable=self.sec_xmax, width=8).pack(side="left", padx=3)
        ttk.Label(limits_frame, text="Normalize traces to max:").pack(side="left", padx=12)
        self.sec_normalize = tk.BooleanVar(value=True)
        ttk.Checkbutton(limits_frame, variable=self.sec_normalize).pack(side="left")

        # Y-lim
        y_frame = ttk.Frame(self.sec_frame)
        y_frame.pack(fill="x", pady=4)
        ttk.Label(y_frame, text="Y-lim:").pack(side="left")
        self.sec_ymin = tk.DoubleVar(value=0.0)
        self.sec_ymax = tk.DoubleVar(value=1.0)
        ttk.Entry(y_frame, textvariable=self.sec_ymin, width=8).pack(side="left", padx=3)
        ttk.Entry(y_frame, textvariable=self.sec_ymax, width=8).pack(side="left", padx=3)

        # Standards
        std_frame = ttk.Frame(self.sec_frame)
        std_frame.pack(fill="x", pady=4)
        ttk.Label(std_frame, text="Standard Volumes (mL, comma-separated):").pack(anchor="w")
        self.sec_standards = tk.StringVar(value="13.81, 16.82, 18.09, 19.05, 20.79")
        ttk.Entry(std_frame, textvariable=self.sec_standards).pack(fill="x", padx=6, pady=4)

    def _build_mass_options(self):
        # Mass options frame
        self.mass_frame = ttk.LabelFrame(self, text="Mass Distribution Options")
        self.mass_frame.pack(fill="x", padx=12, pady=6)

        # Labels
        lab_frame = ttk.Frame(self.mass_frame)
        lab_frame.pack(fill="x", pady=4)
        ttk.Label(lab_frame, text="Labels (comma-separated; leave blank to derive from filenames):").pack(anchor="w")
        self.mass_labels_text = tk.Text(lab_frame, height=2)
        self.mass_labels_text.pack(fill="x", padx=6, pady=4)

        # Second labels (stacked mass)
        sec_lab_frame = ttk.Frame(self.mass_frame)
        sec_lab_frame.pack(fill="x", pady=4)
        ttk.Label(sec_lab_frame, text="Second Labels (stacked mass; comma-separated):").pack(anchor="w")
        self.mass_second_labels_text = tk.Text(sec_lab_frame, height=2)
        self.mass_second_labels_text.pack(fill="x", padx=6, pady=4)

        # Colors
        color_frame = ttk.Frame(self.mass_frame)
        color_frame.pack(fill="x", pady=4)
        ttk.Label(color_frame, text="Colors (comma-separated hex):").pack(anchor="w")
        self.mass_colors = tk.StringVar(value=", ".join(DEFAULT_MASS_COLORS))
        ttk.Entry(color_frame, textvariable=self.mass_colors).pack(fill="x", padx=6, pady=4)

        # X-range and offset
        x_frame = ttk.Frame(self.mass_frame)
        x_frame.pack(fill="x", pady=4)
        ttk.Label(x_frame, text="X-range (kDa):").pack(side="left")
        self.mass_xmin = tk.DoubleVar(value=0.0)
        self.mass_xmax = tk.DoubleVar(value=500.0)
        ttk.Entry(x_frame, textvariable=self.mass_xmin, width=8).pack(side="left", padx=3)
        ttk.Entry(x_frame, textvariable=self.mass_xmax, width=8).pack(side="left", padx=3)
        ttk.Label(x_frame, text="Y-offset step (stacked):").pack(side="left", padx=12)
        self.mass_y_offset = tk.DoubleVar(value=120.0)
        ttk.Entry(x_frame, textvariable=self.mass_y_offset, width=8).pack(side="left")

    def _build_action_buttons(self):
        frame = ttk.Frame(self)
        frame.pack(fill="x", padx=12, pady=12)
        ttk.Button(frame, text="Plot", command=self._plot).pack(side="left")
        ttk.Button(frame, text="Quit", command=self.destroy).pack(side="right")

    def _update_visibility(self):
        pt = self.plot_type.get()
        # Show/hide frames based on plot type
        if pt in ("Single SEC", "Stacked SEC"):
            self.sec_frame.pack_configure(fill="x", padx=12, pady=6)
            self.mass_frame.pack_forget()
        else:
            self.mass_frame.pack_configure(fill="x", padx=12, pady=6)
            self.sec_frame.pack_forget()

    def _add_files(self):
        pt = self.plot_type.get()
        if pt in ("Single SEC", "Stacked SEC"):
            paths = filedialog.askopenfilenames(
                title="Select SEC Excel files",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
            )
        else:
            paths = filedialog.askopenfilenames(
                title="Select mass distribution text files",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
            )
        if paths:
            self.files.extend(paths)
            self._refresh_file_list()

    def _clear_files(self):
        self.files = []
        self._refresh_file_list()

    def _refresh_file_list(self):
        self.file_list.delete(0, tk.END)
        for p in self.files:
            self.file_list.insert(tk.END, p)

    def _browse_save(self):
        path = filedialog.asksaveasfilename(
            title="Save figure as...",
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")]
        )
        if path:
            self.save_path.set(path)

    def _plot(self):
        pt = self.plot_type.get()
        files = list(self.files)
        if not files:
            messagebox.showerror("No files", "Please add at least one file.")
            return

        # Common settings
        font_name = self.font_name.get()
        font_size = int(self.font_size.get())
        figsize = (float(self.fig_w.get()), float(self.fig_h.get()))
        dpi = int(self.dpi.get())
        save = self.save_path.get() or None
        show = bool(self.show_plot.get())

        try:
            if pt == "Single SEC":
                if plot_sec_stack is None:
                    raise ImportError("sec_plot_utils.plot_sec_stack not found.")
                # Labels: single; derive if blank
                labels_text = self.sec_labels_text.get("1.0", "end").strip()
                labels = parse_str_list(labels_text) if labels_text else filenames_to_labels(files)
                if len(files) != 1:
                    messagebox.showerror("Selection error", "Select exactly one SEC Excel file for 'Single SEC'.")
                    return
                # Use stacked function with one file for consistency
                colors = parse_str_list(self.sec_colors.get()) or DEFAULT_COLORS
                std_vols = parse_float_list(self.sec_standards.get())
                fig, axes = plot_sec_stack(
                    files=files,
                    sample_labels=labels[:1],
                    colors=colors[:1],
                    standard_volumes=std_vols,
                    normalize=bool(self.sec_normalize.get()),
                    font_name=font_name,
                    font_size=font_size,
                    xlim=(float(self.sec_xmin.get()), float(self.sec_xmax.get())),
                    ylim=(float(self.sec_ymin.get()), float(self.sec_ymax.get())),
                    figsize=figsize,
                    dpi=dpi,
                    save=save,
                    show=show
                )

            elif pt == "Stacked SEC":
                if plot_sec_stack is None:
                    raise ImportError("sec_plot_utils.plot_sec_stack not found.")
                labels_text = self.sec_labels_text.get("1.0", "end").strip()
                labels = parse_str_list(labels_text) if labels_text else filenames_to_labels(files)
                if len(labels) != len(files):
                    messagebox.showerror("Label mismatch", "Number of labels must match number of files.")
                    return
                colors = parse_str_list(self.sec_colors.get()) or DEFAULT_COLORS
                std_vols = parse_float_list(self.sec_standards.get())
                fig, axes = plot_sec_stack(
                    files=files,
                    sample_labels=labels,
                    colors=colors,
                    standard_volumes=std_vols,
                    normalize=bool(self.sec_normalize.get()),
                    font_name=font_name,
                    font_size=font_size,
                    xlim=(float(self.sec_xmin.get()), float(self.sec_xmax.get())),
                    ylim=(float(self.sec_ymin.get()), float(self.sec_ymax.get())),
                    figsize=figsize,
                    dpi=dpi,
                    save=save,
                    show=show
                )

            elif pt == "Single Mass":
                # Use single mass function if available; otherwise emulate via stack with one dataset
                labels_text = self.mass_labels_text.get("1.0", "end").strip()
                labels = parse_str_list(labels_text) if labels_text else filenames_to_labels(files)
                if len(files) != 1:
                    messagebox.showerror("Selection error", "Select exactly one mass file for 'Single Mass'.")
                    return
                color_list = parse_str_list(self.mass_colors.get()) or DEFAULT_MASS_COLORS
                x_range = (float(self.mass_xmin.get()), float(self.mass_xmax.get()))

                if plot_single_mass is not None:
                    fig, ax = plot_single_mass(
                        dataset=files[0],
                        label=labels[0],
                        second_label=None,  # Optional
                        x_range=x_range,
                        color=color_list[0],
                        font_name=font_name,
                        font_size=font_size,
                        figsize=figsize,
                        dpi=dpi,
                        save=save,
                        show=show
                    )
                else:
                    if plot_mass_offset_stack is None:
                        raise ImportError("mass_plot_utils not found.")
                    fig, axes = plot_mass_offset_stack(
                        datasets=files[:1],
                        labels=labels[:1],
                        second_labels=[None],
                        x_range=x_range,
                        colors=color_list[:1],
                        y_offset_step=float(self.mass_y_offset.get()),
                        font_name=font_name,
                        font_size=font_size,
                        figsize=figsize,
                        dpi=dpi,
                        save=save,
                        show=show
                    )

            elif pt == "Stacked Mass":
                if plot_mass_offset_stack is None:
                    raise ImportError("mass_plot_utils.plot_mass_offset_stack not found.")
                labels_text = self.mass_labels_text.get("1.0", "end").strip()
                labels = parse_str_list(labels_text) if labels_text else filenames_to_labels(files)
                if len(labels) != len(files):
                    messagebox.showerror("Label mismatch", "Number of labels must match number of files.")
                    return

                second_text = self.mass_second_labels_text.get("1.0", "end").strip()
                second_labels = parse_str_list(second_text) if second_text else [None] * len(files)
                if len(second_labels) != len(files):
                    messagebox.showerror("Second label mismatch", "Second labels must match number of files or be empty.")
                    return

                colors = parse_str_list(self.mass_colors.get()) or DEFAULT_MASS_COLORS
                x_range = (float(self.mass_xmin.get()), float(self.mass_xmax.get()))
                fig, axes = plot_mass_offset_stack(
                    datasets=files,
                    labels=labels,
                    second_labels=second_labels,
                    x_range=x_range,
                    colors=colors,
                    y_offset_step=float(self.mass_y_offset.get()),
                    font_name=font_name,
                    font_size=font_size,
                    figsize=figsize,
                    dpi=dpi,
                    save=save,
                    show=show
                )

            else:
                messagebox.showerror("Unknown plot type", f"Unsupported plot type: {pt}")
                return

        except Exception as e:
            messagebox.showerror("Error", f"Plotting failed:\n{e}")

def main():
    app = PlotGUI()
    app.mainloop()

if __name__ == "__main__":
    main()
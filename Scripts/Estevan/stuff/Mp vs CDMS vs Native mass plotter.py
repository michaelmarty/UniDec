import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# USER SETTINGS
# =====================================================
font_name = "Arial"
font_size = 10
x_range = (50, 250)  # x-axis limits in kDa

# Gaussian bandwidth used to convert raw MP events into a smooth curve.
# Smaller values give sharper peaks; larger values give smoother peaks.
mp_gaussian_bandwidth_kda = 3.0

# Number of points used for the smooth MP curve.
mp_curve_points = 2000

# Set True to fill under MP curves.
fill_mp_curves = False

# Optional output file. Set to None to avoid saving.
output_file = r"G:\EC\Mp data\MP_and_mass_spec_overlay.pdf"

# =====================================================
# FILES TO PLOT
# =====================================================
files = [
    r"Z:\Group Share\Estevan\Q Exactive Data\20260330\20260330_EC_AW_MSP1D1_1_to_80_DMPC_unidecfiles\20260330_EC_AW_MSP1D1_1_to_80_DMPC_mass.txt",
    r"Z:\Group Share\Estevan\Q Exactive Data\20260330\20260330_EC_AW_MSP1D1_sfCherryS_1_to_80_DMPC_CDMS_unidecfiles\20260330_EC_AW_MSP1D1_sfCherryS_1_to_80_DMPC_CDMS_mass.txt",
    r"G:\EC\Mp data\MP TRY1.txt",
]

# =====================================================
# FILE TYPES
# Use "mass_spec" for UniDec two-column mass files.
# Use "mp" for raw mass-photometry files containing calibrated_values.
# =====================================================
file_types = [
    "mass_spec",
    "mass_spec",
    "mp",
]

# =====================================================
# LEGEND LABELS
# One label per file
# =====================================================
labels = [
    "MSP1D1-sfGFP DMPC + AM2",
    "MSP1D1 DMPC CD-MS",
    "MSP1D1 DMPC MP",
]

# =====================================================
# CUSTOM COLORS
# One color per file
# =====================================================
colors = [
    "#68BB59",
    "#080808",
    "#FF0000",
]

# =====================================================
# MASS PHOTOMETRY SETTINGS
# =====================================================
mp_mass_column = "calibrated_values"

# Change to "Da" if calibrated_values are stored in daltons.
# Keep as "kDa" if they are already stored in kilodaltons.
mp_input_mass_unit = "kDa"

# Set True to use only events marked as selected.
mp_selected_only = False


# =====================================================
# HELPER FUNCTIONS
# =====================================================
def detect_delimiter(file_path, encodings=("utf-8", "utf-8-sig", "latin1", "utf-16")):
    """Detect a likely delimiter in a text file."""
    for encoding in encodings:
        try:
            with open(file_path, "r", encoding=encoding, errors="strict") as handle:
                sample = handle.read(4096)

            if "\t" in sample:
                return "\t"
            if "," in sample:
                return ","
            if ";" in sample:
                return ";"
            return None
        except Exception:
            continue

    return None


def read_mass_spec_file(file_path):
    """
    Read a two-column UniDec-style mass file.

    Expected columns:
        mass in Da, intensity
    """
    df = pd.read_csv(
        file_path,
        sep=r"\s+",
        header=None,
        comment="#",
        engine="python",
    )

    if df.shape[1] < 2:
        raise ValueError(
            f"Mass-spec file must contain at least two columns:\n{file_path}"
        )

    df = df.iloc[:, :2].copy()
    df.columns = ["Mass_Da", "Intensity"]

    df["Mass_Da"] = pd.to_numeric(df["Mass_Da"], errors="coerce")
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    if df.empty:
        raise ValueError(f"No numeric mass-spec data found in:\n{file_path}")

    df = df[df["Intensity"] >= 0].copy()
    df["Mass_kDa"] = df["Mass_Da"] / 1000.0

    maximum = df["Intensity"].max()
    if maximum <= 0 or np.isnan(maximum):
        raise ValueError(f"Mass-spec intensity has no positive values in:\n{file_path}")

    df["Percent"] = (df["Intensity"] / maximum) * 100.0
    return df[["Mass_kDa", "Percent"]]


def read_mp_events(file_path):
    """Read raw MP events and return calibrated masses in kDa."""
    extension = os.path.splitext(file_path)[1].lower()

    if extension in (".xls", ".xlsx"):
        df = pd.read_excel(file_path)
    else:
        encodings = ("utf-8", "utf-8-sig", "latin1", "utf-16")
        delimiter = detect_delimiter(file_path, encodings)
        df = None
        last_error = None

        for encoding in encodings:
            try:
                if delimiter is None:
                    df = pd.read_csv(
                        file_path,
                        sep=None,
                        engine="python",
                        encoding=encoding,
                    )
                else:
                    df = pd.read_csv(
                        file_path,
                        sep=delimiter,
                        encoding=encoding,
                    )
                break
            except Exception as exc:
                last_error = exc

        if df is None:
            raise RuntimeError(
                f"Could not read MP file:\n{file_path}\nLast error: {last_error}"
            )

    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.replace("\ufeff", "", regex=False)
    )

    lookup = {column.lower(): column for column in df.columns}
    actual_mass_column = lookup.get(mp_mass_column.lower())

    if actual_mass_column is None:
        available = ", ".join(df.columns)
        raise KeyError(
            f"MP mass column '{mp_mass_column}' was not found in:\n{file_path}\n"
            f"Available columns: {available}"
        )

    if mp_selected_only:
        selection_column = lookup.get("selections")
        if selection_column is None:
            raise KeyError(
                f"mp_selected_only is True, but no selections column exists in:\n{file_path}"
            )

        selected_values = {
            "1", "1.0", "true", "yes", "y", "selected", "keep"
        }
        selection_text = (
            df[selection_column]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
        )
        df = df.loc[selection_text.isin(selected_values)].copy()

    masses = pd.to_numeric(df[actual_mass_column], errors="coerce")
    masses = masses.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)

    if masses.size == 0:
        raise ValueError(f"No numeric MP masses found in:\n{file_path}")

    if mp_input_mass_unit.lower() == "da":
        masses = masses / 1000.0
    elif mp_input_mass_unit.lower() != "kda":
        raise ValueError("mp_input_mass_unit must be 'Da' or 'kDa'.")

    return masses


def gaussian_mp_curve(masses_kda, x_min, x_max, bandwidth_kda, points=2000):
    """
    Convert individual MP event masses into a Gaussian kernel density curve.

    The output is normalized so that its tallest point equals 100 percent.
    """
    if bandwidth_kda <= 0:
        raise ValueError("mp_gaussian_bandwidth_kda must be greater than zero.")

    masses_kda = np.asarray(masses_kda, dtype=float)
    masses_kda = masses_kda[np.isfinite(masses_kda)]

    if masses_kda.size == 0:
        raise ValueError("No finite MP masses are available for Gaussian smoothing.")

    x_values = np.linspace(x_min, x_max, points)
    intensity = np.zeros_like(x_values)

    # Process in chunks to avoid creating a very large temporary matrix.
    chunk_size = 5000
    denominator = 2.0 * bandwidth_kda ** 2

    for start in range(0, masses_kda.size, chunk_size):
        chunk = masses_kda[start:start + chunk_size]
        differences = x_values[:, None] - chunk[None, :]
        intensity += np.exp(-(differences ** 2) / denominator).sum(axis=1)

    maximum = intensity.max()
    if maximum <= 0 or np.isnan(maximum):
        raise ValueError("Could not generate a positive MP Gaussian curve.")

    percent = (intensity / maximum) * 100.0

    return pd.DataFrame({
        "Mass_kDa": x_values,
        "Percent": percent,
    })


# =====================================================
# VALIDATE SETTINGS
# =====================================================
if not (len(files) == len(file_types) == len(labels)):
    raise ValueError(
        "files, file_types, and labels must contain the same number of items."
    )

if len(colors) == 0:
    raise ValueError("At least one color is required.")

if x_range[1] <= x_range[0]:
    raise ValueError("x_range maximum must be greater than its minimum.")


# =====================================================
# FIGURE SETUP
# =====================================================
plt.rcParams.update({
    "font.size": font_size,
    "font.family": font_name,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=300)

handles = []
legend_labels = []


# =====================================================
# READ AND PLOT FILES
# =====================================================
for i, (file_path, file_type, label) in enumerate(
    zip(files, file_types, labels)
):
    color = colors[i % len(colors)]
    file_type = file_type.strip().lower()

    if file_type == "mass_spec":
        plot_df = read_mass_spec_file(file_path)

    elif file_type == "mp":
        mp_masses = read_mp_events(file_path)

        # Keep events close enough to the requested range to contribute to
        # Gaussian tails at the visible boundaries.
        margin = 4.0 * mp_gaussian_bandwidth_kda
        mp_masses = mp_masses[
            (mp_masses >= x_range[0] - margin)
            & (mp_masses <= x_range[1] + margin)
        ]

        if mp_masses.size == 0:
            raise ValueError(
                f"No MP events fall near the selected x-axis range for:\n{file_path}"
            )

        plot_df = gaussian_mp_curve(
            masses_kda=mp_masses,
            x_min=x_range[0],
            x_max=x_range[1],
            bandwidth_kda=mp_gaussian_bandwidth_kda,
            points=mp_curve_points,
        )

    else:
        raise ValueError(
            f"Unknown file type '{file_type}' for:\n{file_path}\n"
            "Use 'mass_spec' or 'mp'."
        )

    line, = ax.plot(
        plot_df["Mass_kDa"],
        plot_df["Percent"],
        color=color,
        linewidth=1.5,
        label=label,
    )

    if file_type == "mp" and fill_mp_curves:
        ax.fill_between(
            plot_df["Mass_kDa"],
            plot_df["Percent"],
            0,
            color=color,
            alpha=0.15,
            linewidth=0,
        )

    handles.append(line)
    legend_labels.append(label)


# =====================================================
# FINAL AXIS STYLING
# =====================================================
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(x_range)
ax.set_ylim(0, 100)

ax.set_yticks([0, 50, 100])
ax.set_yticklabels(["0", "%", "100"])

ax.set_xlabel(
    "Mass (kDa)",
    fontname=font_name,
    fontsize=font_size,
)

ax.set_ylabel(
    "Relative Intensity (%)",
    fontname=font_name,
    fontsize=font_size,
)

ax.legend(
    handles,
    legend_labels,
    loc="upper left",
    frameon=False,
    fontsize=font_size,
)

plt.tight_layout()


# =====================================================
# SAVE AND SHOW
# =====================================================
if output_file is not None:
    output_directory = os.path.dirname(output_file)
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    plt.savefig(
        output_file,
        bbox_inches="tight",
    )

    print(f"Figure saved to:\n{output_file}")

plt.show()
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "public" / "UniDec"))
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# Make the repository importable when this script is run directly.
repository_root = Path(__file__).resolve().parents[3]
if str(repository_root) not in sys.path:
    sys.path.insert(0, str(repository_root))

from PublicScripts.Lipids.LipidFunctions import network_analysis  # noqa: E402


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

input_file = "Combined_Fragment_Assignment_Colored.xlsx"
removed_file = "Removed_Identifications.xlsx"

df = pd.read_excel(input_file, sheet_name="Sheet1")


def normalized_text(series):
    return (
        series
        .fillna("")
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
        .str.casefold()
    )


# Apex is outside the dataset comparison used in the manuscript figures.
if "Dataset" not in df.columns:
    raise ValueError("The input workbook is missing required column: Dataset")
apex_mask = normalized_text(df["Dataset"]).eq("apex")
apex_count = int(apex_mask.sum())
df = df.loc[~apex_mask].copy()
print(f"Excluded {apex_count} Apex IDs from analysis.")


# Exclude every exact entry recorded in the manual-removal workbook. The
# category rules below also prevent equivalent rows from reappearing if a new
# export assigns them different alignment IDs.
exact_manual_removal = pd.Series(False, index=df.index)
if os.path.exists(removed_file):
    removed_df = pd.read_excel(removed_file, sheet_name="Sheet1")
    removal_key_columns = [
        "Dataset",
        "Spectrum reference file name",
        "Alignment ID",
        "Metabolite name",
        "Adduct type",
    ]
    if all(column in df.columns for column in removal_key_columns) and all(
        column in removed_df.columns for column in removal_key_columns
    ):
        def make_removal_keys(frame):
            key_frame = pd.DataFrame(index=frame.index)
            for column in removal_key_columns:
                if column == "Alignment ID":
                    key_frame[column] = (
                        pd.to_numeric(frame[column], errors="coerce")
                        .astype("Int64")
                        .astype(str)
                    )
                else:
                    key_frame[column] = normalized_text(frame[column])
            return key_frame.agg("||".join, axis=1)

        removed_keys = set(make_removal_keys(removed_df))
        exact_manual_removal = make_removal_keys(df).isin(removed_keys)

class_text = normalized_text(df["Ontology"])
adduct_text = normalized_text(df["Adduct type"]).str.replace(
    " ", "", regex=False
)
lipid_name_text = normalized_text(df["Metabolite name"])

# Rules documented in Removed_Identifications.xlsx:
#   - DG sodium adduct excluded
#   - Formate adduct excluded
#   - FAHFA excluded
dg_sodium_removal = (
    class_text.eq("dg")
    & adduct_text.str.contains(r"m\+na", regex=True, na=False)
)
formate_removal = adduct_text.str.contains(
    r"hcoo|cho2|formate",
    regex=True,
    na=False,
)
fahfa_removal = (
    class_text.str.contains("fahfa", regex=False, na=False)
    | lipid_name_text.str.contains("fahfa", regex=False, na=False)
)

manual_removal_mask = (
    exact_manual_removal
    | dg_sodium_removal
    | formate_removal
    | fahfa_removal
)
manual_removal_count = int(manual_removal_mask.sum())
df = df.loc[~manual_removal_mask].copy().reset_index(drop=True)
print(f"Excluded {manual_removal_count} manually removed IDs from analysis.")

required_columns = {
    "Dataset",
    "Metabolite name",
    "Ontology",
    "Polarity",
    "Adduct type",
    "Average Mz",
    "MS/MS spectrum",
    "Overall Match Status",
    "Reverse dot product",
    "Fill %",
}
missing_columns = required_columns.difference(df.columns)
if missing_columns:
    raise ValueError(
        "The input workbook is missing required columns: "
        + ", ".join(sorted(missing_columns))
    )

true_values = [True, "TRUE", "True", "true", 1, "1"]
df["ID Status"] = np.where(
    df["Overall Match Status"].isin(true_values),
    "True",
    "False",
)

status_colors = {
    "True": "#0072B2",
    "False": "#D55E00",
}


# Calculate each ID's class strength in a spectral-similarity network containing
# every candidate from the same lipid class and polarity. Raw strength is useful
# as a predictor, while within-group percentiles make different-sized classes
# comparable in the all-class plots below.
network_tolerance = 0.3
minimum_network_size = 3


def calculate_all_class_strengths(frame):
    frame = frame.copy()
    frame["_Figure Row"] = np.arange(len(frame))
    frame["_Network Class"] = normalized_text(frame["Ontology"])
    frame["_Network Polarity"] = normalized_text(frame["Polarity"])
    frame["Class Strength"] = np.nan

    analyzed_groups = 0
    for _, network_input in frame.groupby(
        ["_Network Class", "_Network Polarity"],
        sort=True,
    ):
        if len(network_input) < minimum_network_size:
            continue

        local_input = network_input.reset_index(drop=True)
        network_df, _ = network_analysis(
            local_input,
            mz_tol=network_tolerance,
            header="Class_",
            add_all_cols=True,
        )
        frame.loc[
            network_df["_Figure Row"].astype(int).to_numpy(),
            "Class Strength",
        ] = network_df["Class_strength"].to_numpy()
        analyzed_groups += 1

    frame["Class Strength Percentile"] = (
        frame.groupby(
            ["_Network Class", "_Network Polarity"],
        )["Class Strength"]
        .rank(method="average", pct=True)
        * 100
    )
    frame["Network Strength Outlier"] = (
        frame["Class Strength Percentile"] <= 5
    )
    frame = frame.drop(
        columns=["_Figure Row", "_Network Class", "_Network Polarity"]
    )

    print(
        f"Calculated class strength for {analyzed_groups} "
        f"class/polarity networks."
    )
    return frame


df = calculate_all_class_strengths(df)
df["Square Root Class Strength"] = np.sqrt(
    pd.to_numeric(df["Class Strength"], errors="coerce").clip(lower=0)
)

class_strength_valid = df.dropna(
    subset=["Square Root Class Strength"]
).copy()
class_strength_auc = roc_auc_score(
    class_strength_valid["ID Status"].eq("True"),
    class_strength_valid["Square Root Class Strength"],
)

# PLS regression with binary acceptance status as the supervised response.
# This is commonly visualized as PLS discriminant analysis (PLS-DA).
# Restore the three predictors used in the original figure.
pls_features = [
    "Square Root Class Strength",
    "Reverse dot product",
    "Fill %",
]

for column in pls_features:
    df[column] = pd.to_numeric(df[column], errors="coerce")

pls_df = df.dropna(subset=pls_features).copy()
pls_df["Accepted Response"] = (
    pls_df["ID Status"].eq("True").astype(float)
)

feature_scaler = StandardScaler()
scaled_features = feature_scaler.fit_transform(pls_df[pls_features])
accepted_response = pls_df["Accepted Response"].to_numpy()

pls = PLSRegression(n_components=2, scale=False)
pls.fit(scaled_features, accepted_response)
pls_scores = pls.x_scores_.copy()

# Orient component 1 so positive scores point toward accepted IDs. PLS signs
# are otherwise arbitrary and can flip between equivalent fits.
accepted_component_1 = pls_scores[accepted_response == 1, 0].mean()
excluded_component_1 = pls_scores[accepted_response == 0, 0].mean()
if accepted_component_1 < excluded_component_1:
    pls_scores[:, 0] *= -1

pls_df["PLS1"] = pls_scores[:, 0]
pls_df["PLS2"] = pls_scores[:, 1]

# Cross-validated predictions provide a less optimistic summary than scoring
# the same rows used to fit the displayed PLS model.
cross_validation = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42,
)
cross_validated_predictions = cross_val_predict(
    make_pipeline(
        StandardScaler(),
        PLSRegression(n_components=2, scale=False),
    ),
    pls_df[pls_features],
    accepted_response,
    cv=cross_validation,
).ravel()
cross_validated_auc = roc_auc_score(
    accepted_response,
    cross_validated_predictions,
)
cross_validated_accuracy = accuracy_score(
    accepted_response,
    cross_validated_predictions >= 0.5,
)

pls_fig = plt.figure(figsize=(6.5, 4.77), layout="constrained")
figure_grid = pls_fig.add_gridspec(
    nrows=2,
    ncols=2,
    height_ratios=[1, 1.05],
    width_ratios=[1.35, 1],
)
score_ax = pls_fig.add_subplot(figure_grid[0, 0])
coefficient_ax = pls_fig.add_subplot(figure_grid[0, 1])
strength_distribution_ax = pls_fig.add_subplot(figure_grid[1, 0])
class_strength_ax = pls_fig.add_subplot(figure_grid[1, 1])

status_markers = {
    "False": ("X", "Excluded", 40),
    "True": ("o", "Accepted", 28),
}

for status, (marker, label, marker_size) in status_markers.items():
    status_df = pls_df[pls_df["ID Status"] == status]
    score_ax.scatter(
        status_df["PLS1"],
        status_df["PLS2"],
        color=status_colors[status],
        marker=marker,
        s=marker_size,
        alpha=0.55,
        linewidths=0.25,
        label=label,
    )

score_ax.set_xlabel("PLS Component 1")
score_ax.set_ylabel("PLS Component 2")
score_ax.set_title("Top-3-Feature Identification-Quality PLS-DA")
score_ax.spines["top"].set_visible(False)
score_ax.spines["right"].set_visible(False)
score_ax.legend(
    frameon=True,
    facecolor="white",
    edgecolor="none",
    framealpha=0.90,
    loc="upper left",
)
score_ax.text(
    0.02,
    0.02,
    (
        f"5-fold CV AUC = {cross_validated_auc:.2f}\n"
        f"5-fold CV accuracy = {100 * cross_validated_accuracy:.1f}%"
    ),
    transform=score_ax.transAxes,
    fontsize=8,
    ha="left",
    va="bottom",
    bbox={
        "facecolor": "white",
        "edgecolor": "none",
        "alpha": 0.90,
        "pad": 3,
    },
)

feature_labels = {
    "Reverse dot product": "Reverse dot product",
    "Fill %": "Fill %",
    "Square Root Class Strength": "Square-root class strength",
}

# Because all predictors were standardized, coefficient magnitudes can be
# compared directly. Positive values predict acceptance; negative values
# predict exclusion after accounting for the other features.
pls_coefficients = np.asarray(pls.coef_).reshape(-1)
coefficient_positions = np.arange(len(pls_features))
coefficient_colors = np.where(
    pls_coefficients >= 0,
    status_colors["True"],
    status_colors["False"],
)

coefficient_ax.barh(
    coefficient_positions,
    pls_coefficients,
    color=coefficient_colors,
    edgecolor="black",
    linewidth=0.6,
)
coefficient_ax.set_yticks(
    coefficient_positions,
    labels=[feature_labels[feature] for feature in pls_features],
)
coefficient_ax.invert_yaxis()
coefficient_ax.axvline(0, color="#777777", linewidth=0.9)
coefficient_ax.set_xlabel("Standardized PLS Coefficient")
coefficient_ax.set_ylabel("Quality Feature")
coefficient_ax.set_title("Top 3 ID-Status Predictors")
coefficient_ax.spines["top"].set_visible(False)
coefficient_ax.spines["right"].set_visible(False)
coefficient_ax.text(
    0.98,
    0.04,
    "Blue: Accepted\nOrange: Excluded",
    transform=coefficient_ax.transAxes,
    fontsize=8,
    ha="right",
    va="bottom",
)


# C: a square-root transformation retains zero-strength IDs and controls the
# long tail without using a logarithm or adding a pseudocount.
strength_plot_df = df.dropna(
    subset=["Square Root Class Strength", "Class Strength Percentile"]
).copy()
violin_values = [
    strength_plot_df.loc[
        strength_plot_df["ID Status"].eq(status),
        "Square Root Class Strength",
    ].to_numpy()
    for status in ["False", "True"]
]
violins = strength_distribution_ax.violinplot(
    violin_values,
    positions=[0, 1],
    widths=0.76,
    showmeans=False,
    showmedians=False,
    showextrema=False,
)
for violin_body, status in zip(violins["bodies"], ["False", "True"]):
    violin_body.set_facecolor(status_colors[status])
    violin_body.set_edgecolor("black")
    violin_body.set_linewidth(0.8)
    violin_body.set_alpha(0.48)

strength_distribution_ax.boxplot(
    violin_values,
    positions=[0, 1],
    widths=0.18,
    patch_artist=True,
    showfliers=False,
    boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 1.0},
    medianprops={"color": "black", "linewidth": 1.5},
    whiskerprops={"color": "black", "linewidth": 1.0},
    capprops={"color": "black", "linewidth": 1.0},
)

random_generator = np.random.default_rng(42)
for position, status in enumerate(["False", "True"]):
    status_df = strength_plot_df[
        strength_plot_df["ID Status"].eq(status)
    ]
    jitter = random_generator.normal(0, 0.065, len(status_df))
    strength_distribution_ax.scatter(
        position + jitter,
        status_df["Square Root Class Strength"],
        s=8,
        color=status_colors[status],
        alpha=0.14,
        linewidths=0,
        rasterized=True,
    )
    outlier_mask = status_df["Network Strength Outlier"].to_numpy()
    strength_distribution_ax.scatter(
        position + jitter[outlier_mask],
        status_df.loc[outlier_mask, "Square Root Class Strength"],
        s=19,
        facecolors="none",
        edgecolors="black",
        linewidths=0.55,
        alpha=0.80,
        rasterized=True,
    )

strength_distribution_ax.set_xticks(
    [0, 1],
    labels=["Excluded", "Accepted"],
)
strength_distribution_ax.set_ylabel("Square-Root Class Strength")
strength_distribution_ax.set_title(
    "Class-Network Strength by ID Status",
    pad=10,
)
strength_distribution_ax.set_ylim(
    0,
    strength_plot_df["Square Root Class Strength"].max() * 1.28,
)
strength_distribution_ax.spines["top"].set_visible(False)
strength_distribution_ax.spines["right"].set_visible(False)
strength_distribution_ax.text(
    0.50,
    0.98,
    f"Class strength ROC AUC = {class_strength_auc:.2f}",
    transform=strength_distribution_ax.transAxes,
    fontsize=8,
    ha="center",
    va="top",
)
strength_distribution_ax.legend(
    handles=[
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=0.8,
            markersize=7,
            label="Lowest 5% within class/polarity",
        )
    ],
    loc="upper center",
    bbox_to_anchor=(0.5, 0.92),
    frameon=False,
    fontsize=8,
)

# D: show whether the strength shift is consistent across sufficiently sampled
# lipid classes. Percentiles normalize strength within class and polarity.
class_strength_summary = (
    strength_plot_df
    .groupby(["Ontology", "ID Status"])
    .agg(
        Count=("Class Strength Percentile", "size"),
        Median_Percentile=("Class Strength Percentile", "median"),
    )
    .reset_index()
)
class_counts = class_strength_summary.pivot(
    index="Ontology",
    columns="ID Status",
    values="Count",
)
class_medians = class_strength_summary.pivot(
    index="Ontology",
    columns="ID Status",
    values="Median_Percentile",
)
eligible_classes = class_counts.index[
    class_counts.get("False", pd.Series(index=class_counts.index, dtype=float)).ge(10)
    & class_counts.get("True", pd.Series(index=class_counts.index, dtype=float)).ge(10)
]
class_medians = class_medians.loc[eligible_classes].dropna(
    subset=["False", "True"]
)
class_medians["Strength Shift"] = (
    class_medians["True"] - class_medians["False"]
)
class_medians = class_medians.sort_values("Strength Shift", ascending=False)

class_positions = np.arange(len(class_medians))
class_strength_ax.hlines(
    class_positions,
    class_medians["False"],
    class_medians["True"],
    color="#A5A5A5",
    linewidth=1.5,
    zorder=1,
)
class_strength_ax.scatter(
    class_medians["False"],
    class_positions,
    color=status_colors["False"],
    edgecolor="black",
    linewidth=0.5,
    s=45,
    label="Excluded",
    zorder=2,
)
class_strength_ax.scatter(
    class_medians["True"],
    class_positions,
    color=status_colors["True"],
    edgecolor="black",
    linewidth=0.5,
    s=45,
    label="Accepted",
    zorder=2,
)
class_strength_ax.set_yticks(class_positions, labels=class_medians.index)
class_strength_ax.invert_yaxis()
class_strength_ax.set_xlim(0, 100)
class_strength_ax.set_xlabel(
    "Median Strength Percentile\nWithin Class/Polarity"
)
class_strength_ax.set_ylabel("Lipid Class")
class_strength_ax.set_title("Class-Specific Network Strength")
class_strength_ax.spines["top"].set_visible(False)
class_strength_ax.spines["right"].set_visible(False)
# Keep the status key outside panel D so its markers cannot be mistaken for data points.
class_strength_ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    ncol=1,
    frameon=False,
    fontsize=8,
    borderaxespad=0,
)


for panel_label, ax in zip(
    ["A", "B", "C", "D"],
    [score_ax, coefficient_ax, strength_distribution_ax, class_strength_ax],
):
    ax.text(
        0.02,
        1.18,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

pls_fig.savefig(
    "Figure5.png",
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.08,
)
pls_fig.savefig(
    "Figure5.pdf",
    bbox_inches="tight",
    pad_inches=0.08,
)

print("\nIDs included in PLS-DA:")
print(pls_df.groupby("ID Status").size().to_string())
print("PLS-DA features: " + ", ".join(pls_features))
print(
    f"5-fold cross-validated AUC: {cross_validated_auc:.3f}\n"
    f"5-fold cross-validated accuracy: {cross_validated_accuracy:.3f}"
)
print(
    f"\nAll-class network strength "
    f"(m/z tolerance = {network_tolerance:.2f} Da):"
)
print(f"Class-strength ROC AUC: {class_strength_auc:.3f}")
print("Lowest-5% network-strength outliers by ID status:")
print(
    strength_plot_df
    .groupby("ID Status")["Network Strength Outlier"]
    .agg(["sum", "count", "mean"])
    .to_string()
)

plt.show(block=True)

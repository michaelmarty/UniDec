import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

# ------------------------------
# Load data
# ------------------------------

folder = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\20251216_MoreLipidDiscs_FTs\DDM Comparisons"
file_path = f"{folder}\\1%vs2%_DDM_Comparison_Stats.xlsx"

df = pd.read_excel(file_path)

# ------------------------------
# Clean column names
# ------------------------------
df.columns = (
    df.columns
      .str.strip()
      .str.replace(r'\s+', '_', regex=True)
      .str.replace('-', '_')
      .str.replace('+', '_')
)

# Remove rows with NaNs in key columns
df_clean = df.dropna(subset=['Avg_log2fc', 'p', 'BH_p']).copy()

# ------------------------------
# Assign specific colors per class
# ------------------------------
class_colors = {
    'CL': '#F25C54',
    'PG': '#356935',
    'PE': '#6017A9',
}


def volcano_plot(df, y_col, title, fc_thresh=0.263, p_thresh=0.05):
    plt.figure(figsize=(10, 6))

    for cls in df['Class'].unique():
        subset = df[df['Class'] == cls]
        color = class_colors.get(cls, '#7f7f7f')
        plt.scatter(
            subset['Avg_log2fc'],
            -np.log10(subset[y_col]),
            color=color,
            label=cls,
            alpha=0.7
        )

    # Filter for significant lipids
    sig_df = df[(df[y_col] < p_thresh) & (abs(df['Avg_log2fc']) >= fc_thresh)]

    for i, row in sig_df.iterrows():
        plt.text(
            row['Avg_log2fc'],
            -np.log10(row[y_col]) + 0.1,  # small offset
            row['Lipid'],
            fontsize=8,
            ha='center',
            rotation=30
        )

    # Horizontal significance line
    plt.axhline(-np.log10(p_thresh), color='grey', linestyle='--', linewidth=0.8)

    # Vertical fold-change thresholds
    plt.axvline(fc_thresh, color='grey', linestyle='--', linewidth=0.8)
    plt.axvline(-fc_thresh, color='grey', linestyle='--', linewidth=0.8)
    plt.xlim(-3,3)
    plt.xlabel('Average log2 Fold Change')
    plt.ylabel(f'-log10({y_col})')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

# ------------------------------
# Volcano plots
# ------------------------------
# volcano_plot(df_clean, 'p', 'Volcano Plot (Raw p-values)')
volcano_plot(df_clean, 'BH_p', 'Volcano Plot (BH-adjusted FDR)')

# ------------------------------
# Fold-change plot per class
# ------------------------------
plt.figure(figsize=(8, 6))

class_means = df_clean.groupby('Class')['Avg_log2fc'].mean()
classes = class_means.index
avg_fc = class_means.values
colors = [class_colors.get(c, '#7f7f7f') for c in classes]

bars = plt.bar(classes, avg_fc, color=colors)

# Significance labels (BH-adjusted)
for i, cls in enumerate(classes):
    sig = (df_clean[df_clean['Class'] == cls]['BH_p'] < 0.05).any()
    label = '*' if sig else 'ns'
    plt.text(
        i,
        avg_fc[i] + 0.05,
        label,
        ha='center',
        va='bottom',
        fontsize=12,
        color='red' if sig else 'black'
    )

plt.axhline(0, color='black', linewidth=0.8)
plt.ylabel('Average log2 Fold Change')
plt.xlabel('Lipid Class')
plt.title('Average Fold Change per Class')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
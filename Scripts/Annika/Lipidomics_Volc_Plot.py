import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# ---- Load CSV ----
file_path = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\Analysis\NewStats.csv"
df = pd.read_csv(file_path)


# ---- Columns ----
log2fc_col = "Avg log2fc"
raw_p_col = "p value"
class_col = "Class"
name_col = "Lipid"
# Convert to numeric
df[raw_p_col] = pd.to_numeric(df[raw_p_col], errors='coerce')
# ---- Ensure numeric ----
df[log2fc_col] = pd.to_numeric(df[log2fc_col], errors='coerce')
df[raw_p_col] = pd.to_numeric(df[raw_p_col], errors='coerce')
df = df.dropna(subset=[log2fc_col, raw_p_col, class_col])

# ---- Replace zero p-values ----
df[raw_p_col] = df[raw_p_col].replace(0, 1e-300)

# ---- Compute -log10(p-value) ----
df['NegLog10P'] = -np.log10(df[raw_p_col])

# ---- Volcano parameters ----
pval_thresh = 0.05
fc_threshold = 1

# ---- Class colors ----
class_colors = {'CL': '#F25C54', 'PE': '#6017A9', 'PG':'#356935'}
df['color'] = df[class_col].map(class_colors).fillna('gray')

# ---- SIGNIFICANT LIPIDS (CIRCLE THESE) ----
sig_points = df[df[raw_p_col] < pval_thresh]

# ---- Volcano plot ----
plt.figure(figsize=(8, 6))

# Scatter all points
plt.scatter(df[log2fc_col], df['NegLog10P'],
            color=df['color'], alpha=0.7, s=100)

# Circle significant points
plt.scatter(sig_points[log2fc_col], sig_points['NegLog10P'],
            facecolors='none', edgecolors='black',
            linewidths=1.8, s=200)

# Threshold lines
plt.axvline(fc_threshold, color='grey', linestyle='--', linewidth=1)
plt.axvline(-fc_threshold, color='grey', linestyle='--', linewidth=1)
neglog10_cutoff = -np.log10(pval_thresh)
plt.axhline(neglog10_cutoff, color='grey', linestyle=':', linewidth=1.5)

# ---- Label significant lipids ----
for i, row in sig_points.iterrows():
    plt.text(row[log2fc_col], row['NegLog10P'] + 0.1,  # slightly above the point
             row[name_col],
             fontsize=9,
             ha='center',
             va='bottom',
             rotation=45)  # optional: rotate for readability

plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
plt.title('Volcano Plot (p < 0.05 circled & labeled)')
plt.tight_layout()
plt.show()

from scipy.stats import combine_pvalues

# ================================
#      FOLD CHANGE BAR PLOT
# ================================

# ---- Standardize class names ----
df[class_col] = df[class_col].astype(str).str.replace(r'\s+', '', regex=True).str.upper()
classes = ['PE', 'PG', 'CL']

# ---- Class colors ----
class_colors = {'CL': '#F25C54', 'PE': '#6017A9', 'PG': '#356935'}

# ---- Prepare lists to store class stats ----
mean_fc = []
yerr_lower = []
yerr_upper = []
overall_p = []

# ---- Loop through each class ----
for cls in classes:
    class_data = df[df[class_col] == cls]

    if class_data.empty:
        print(f"Warning: No data for class {cls}")
        mean_fc.append(0)
        yerr_lower.append(0)
        yerr_upper.append(0)
        overall_p.append(1.0)
        continue

    # Mean log2 fold change
    fc_vals = class_data[log2fc_col].astype(float)
    mean_val = fc_vals.mean()
    mean_fc.append(mean_val)

    # Confidence intervals from Excel
    ci_lower = class_data['Lower CI'].astype(float).mean()  # mean if multiple rows per class
    ci_upper = class_data['Upper CI'].astype(float).mean()

    # Compute distances from mean for yerr
    lower_err = max(mean_val - ci_lower, 0)
    upper_err = max(ci_upper - mean_val, 0)
    yerr_lower.append(lower_err)
    yerr_upper.append(upper_err)

    # Overall p-value using Fisher's method
    pvals = class_data[raw_p_col].astype(float).values
    if len(pvals) > 0:
        stat, p_combined = combine_pvalues(pvals)
        overall_p.append(p_combined)
    else:
        overall_p.append(1.0)


# ---- Convert p-values to significance stars ----
def p_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


stars = [p_to_stars(p) for p in overall_p]

# ---- Convert error bars to numpy array for asymmetric bars ----
yerr = [np.array(yerr_lower), np.array(yerr_upper)]

# ---- Plot bar chart ----
plt.figure(figsize=(6, 6))
bars = plt.bar(classes, mean_fc, yerr=yerr,
               color=[class_colors[c] for c in classes],
               capsize=5)
plt.ylabel('Mean Log2 Fold Change')
plt.title('Class-level Fold Change (with CIs)')
plt.axhline(0, color='black', linewidth=0.8)

# ---- Add stars above/below bars ----
for i, (bar, star) in enumerate(zip(bars, stars)):
    height = bar.get_height()
    if height >= 0:
        y = height + yerr_upper[i] + 0.05
        va = 'bottom'
    else:
        y = height - yerr_lower[i] - 0.05
        va = 'top'
    plt.text(bar.get_x() + bar.get_width() / 2, y, star,
             ha='center', va=va, fontsize=16, fontweight='bold')

plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
# ---- Load CSV ----
file_path = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\Analysis\THISONESamples1,2,3,4,5,6_ANALYSIS.csv"
df = pd.read_csv(file_path)

# ---- Columns ----
log2fc_col = "Log2FC"
raw_p_col = "P-value"
class_col = "Class"
name_col = "Lipid"

# ---- Ensure numeric ----
df[log2fc_col] = pd.to_numeric(df[log2fc_col], errors='coerce')
df[raw_p_col] = pd.to_numeric(df[raw_p_col], errors='coerce')
df = df.dropna(subset=[log2fc_col, raw_p_col, class_col])

# ---- Replace zero p-values ----
df[raw_p_col] = df[raw_p_col].replace(0, 1e-300)

# ---- Compute -log10(p-value) ----
df['NegLog10P'] = -np.log10(df[raw_p_col])

# ---- Parameters ----
fc_threshold = 1
pval_thresh = 0.05
top_n = 10

# ---- Color mapping by class ----
class_colors = {'CL': '#F25C54', 'PE': '#356935', 'PG': '#6017A9'}
df['color'] = df[class_col].map(class_colors).fillna('gray')

# ---- Identify top N lipids by significance ----
df['rank_p'] = df['NegLog10P'].rank(method='min', ascending=False)
top_points = df[df['rank_p'] <= top_n]

# ---- Volcano plot ----
plt.figure(figsize=(8,6))
plt.scatter(df[log2fc_col], df['NegLog10P'], color=df['color'], alpha=0.7, s=100)

# Highlight top points
plt.scatter(top_points[log2fc_col], top_points['NegLog10P'],
            facecolors='none', edgecolors='black', s=120, linewidths=1.5)

# Fold-change threshold lines
plt.axvline(fc_threshold, color='grey', linestyle='--', linewidth=1)
plt.axvline(-fc_threshold, color='grey', linestyle='--', linewidth=1)

# P-value threshold line
neglog10_cutoff = -np.log10(pval_thresh)
plt.axhline(neglog10_cutoff, color='grey', linestyle=':', linewidth=1.5)
plt.text(plt.xlim()[1] - 1, neglog10_cutoff + 0.1, f'p = {pval_thresh}', color='grey', fontsize=10)

# Label top points
x_offset = 0.15
y_offset = 0.05
for idx, (i, row) in enumerate(top_points.iterrows()):
    x = row[log2fc_col]
    y = row['NegLog10P']
    if idx % 2 == 0:
        plt.text(x + x_offset, y + y_offset, row[name_col], ha='left', va='bottom', fontsize=9)
    else:
        plt.text(x - x_offset, y + y_offset, row[name_col], ha='right', va='bottom', fontsize=9)

plt.xlabel('Log2 Fold Change')


plt.ylabel('-Log10(p-value)')
plt.title('Volcano Plot Colored by Lipid Class')
plt.tight_layout()
plt.show()
extreme = df[df[log2fc_col].abs() > 5]
print(extreme[[name_col, log2fc_col, raw_p_col, class_col]])



# ---- Define thresholds ----
fc_thresh = 1  # optional: consider only significant fold-changes
pval_thresh = 0.05  # optional: consider only significant lipids

# ---- Compute mean Log2 fold change per class ----
classes = ['PE', 'PG', 'CL']
mean_fc = []
std_fc = []
p_vals = []

for cls in classes:
    class_data = df[df['Class'] == cls]

    # Optionally, only use significant lipids:
    sig = (abs(class_data['Log2FC']) >= fc_thresh) & (class_data['NegLog10P'] >= pval_thresh)
    fc_values = class_data.loc[sig, 'Log2FC']
    p_values = 10 ** (-class_data.loc[sig, 'NegLog10P'])  # Convert -log10(p) → p
    # If no significant lipids, take all
    if fc_values.empty:
        fc_values = class_data['Log2FC']
        p_values = 10 ** (-class_data['NegLog10P'])
    mean_fc.append(fc_values.mean())
    std_fc.append(fc_values.std())
    p_vals.append(np.nanmean(p_values))  # average p-value per class

    # ---- Determine significance stars ----
    def p_to_stars(p):
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return 'ns'


    stars = [p_to_stars(p) for p in p_vals]
# ---- Plot bar chart ----
plt.figure(figsize=(6, 6))

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 14
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
bars = plt.bar(classes, mean_fc, yerr=std_fc, color=['royalblue', 'darkslateblue', 'Magenta'], capsize=5)
plt.ylabel('Mean Log2 Fold Change')
plt.title('Class-level Fold Change')
plt.axhline(0, color='black', linewidth=0.8)  # zero line
#plt.grid(axis='y', linestyle='--', alpha=0.7)
# ---- Add significance stars above each bar ----
for i, (bar, star) in enumerate(zip(bars, stars)):
    height = bar.get_height()
    err = std_fc[i] if not np.isnan(std_fc[i]) else 0  # handle NaNs safely
    y = height + err + 0.15  # space above error bar
    plt.text(bar.get_x() + bar.get_width()/2, y, star,
             ha='center', va='bottom', fontsize=20, fontweight='bold')
plt.tight_layout()
plt.show()


import pandas as pd
import matplotlib.pyplot as plt

# ---- Step 1: Load your processed data ----
# Your CSV should have columns: 'Class', 'Log2FC', 'NegLog10P'
file_path = r"Z:\Group Share\Annika\Stellar\EC\Replicates\MembDisc_vs_FlowThruDisc.csv"
data = pd.read_csv(file_path)
# Ensure numeric columns are read correctly
data['Log2FC'] = pd.to_numeric(data['Log2FC'], errors='coerce')
data['NegLog10P'] = pd.to_numeric(data['NegLog10P'], errors='coerce')

# ---- Step 2: Volcano plot by lipid class ----
plt.figure(figsize=(10,6))

# Define colors for each class
colors = {'PE':'royalblue', 'CL':'magenta', 'PG':'darkslateblue'}

# Plot points by class
for lipid_class in data['Class'].unique():
    class_data = data[data['Class'] == lipid_class]
    plt.scatter(class_data['Log2FC'], class_data['NegLog10P'],
                color=colors.get(lipid_class, 'gray'), label=lipid_class, s=50, alpha=0.7)

# Highlight significant points
sig = (abs(data['Log2FC']) >= 1) & (data['NegLog10P'] >= 1.3)
plt.scatter(data.loc[sig, 'Log2FC'], data.loc[sig, 'NegLog10P'],
            facecolors='none', edgecolors='none', s=80, linewidths=1.5)

# Threshold lines
plt.axhline(y=1.3, color='black', linestyle='--', linewidth=1)
plt.axvline(x=1, color='black', linestyle='--', linewidth=1)
plt.axvline(x=-1, color='black', linestyle='--', linewidth=1)

# Add labels for significant points
sig = (abs(data['Log2FC']) >= 1) & (data['NegLog10P'] >= 1.3)
plt.scatter(data.loc[sig, 'Log2FC'], data.loc[sig, 'NegLog10P'],
            color=[colors.get(c,'gray') for c in data.loc[sig, 'Class']], s=80, alpha=0.9)

# Label only the top N significant points by -log10(p-value)
top_n = 5
top_sig = data.loc[sig].nlargest(top_n, 'NegLog10P')
for i, row in top_sig.iterrows():
    plt.text(row['Log2FC'] + 0.05, row['NegLog10P'] + 0.05, row['Lipid'], fontsize=9)

# Labels and legend
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
plt.title('Volcano Plot by Lipid Class')
plt.legend(title='Class')
plt.grid(True)
plt.tight_layout()
plt.show()



# ---- Define thresholds ----
fc_thresh = 1  # optional: consider only significant fold-changes
pval_thresh = 1.3  # optional: consider only significant lipids

# ---- Compute mean Log2 fold change per class ----
classes = ['PE', 'PG', 'CL']
mean_fc = []
std_fc = []

for cls in classes:
    class_data = data[data['Class'] == cls]

    # Optionally, only use significant lipids:
    sig = (abs(class_data['Log2FC']) >= fc_thresh) & (class_data['NegLog10P'] >= pval_thresh)
    fc_values = class_data.loc[sig, 'Log2FC']

    # If no significant lipids, take all
    if fc_values.empty:
        fc_values = class_data['Log2FC']

    mean_fc.append(fc_values.mean())
    std_fc.append(fc_values.std())

# ---- Plot bar chart ----
plt.figure(figsize=(6, 6))

plt.bar(classes, mean_fc, yerr=std_fc, color=['royalblue', 'darkslateblue', 'Magenta'], capsize=5)
plt.ylabel('Mean Log2 Fold Change')
plt.title('Class-level Fold Change')
plt.axhline(0, color='black', linewidth=0.8)  # zero line
#plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
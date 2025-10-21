import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Step 1: Load CSV
# -----------------------------
folder = r"Z:\\Group Share\\Annika\\Stellar\\EC\\Replicates"
file_path = f"{folder}\\MembDisc_vs_FlowThruDisc.csv"

data = pd.read_csv(file_path, header=0)


# -----------------------------
# Step 2: Define sample columns
# -----------------------------
# Make sure these match your CSV exactly
samples = ['MD1','MD2','MD3','FD1','FD2','FD3']

# Convert to numeric, coerce errors, fill NaNs
for col in samples:
    data[col] = pd.to_numeric(data[col], errors='coerce')
data[samples] = data[samples].fillna(0)

# -----------------------------
# Step 3: Normalize per sample (mTIC)
# -----------------------------
data_norm = data[samples].div(data[samples].sum(axis=0), axis=1)

# -----------------------------
# Step 4: Prepare matrix for PCA
# -----------------------------
X = data_norm.T.values  # samples x lipids
X_centered = X - X.mean(axis=0)  # mean-center each lipid

# -----------------------------
# Step 5: Covariance and PCA
# -----------------------------
cov_matrix = np.cov(X_centered, rowvar=False)
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

# Sort eigenvectors by descending eigenvalues
sorted_idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[sorted_idx]
eigenvectors = eigenvectors[:, sorted_idx]

# Project onto first 2 principal components
pc1_pc2 = X_centered @ eigenvectors[:, :2]

# -----------------------------
# Step 6: Prepare DataFrame for plotting
# -----------------------------
groups = ['MD','MD','MD','FD','FD','FD']  # matches 'samples' order
pca_df = pd.DataFrame(pc1_pc2, columns=['PC1','PC2'])
pca_df['Sample'] = samples
pca_df['Group'] = groups

# -----------------------------
# Step 7: Plot PCA
# -----------------------------
plt.figure(figsize=(8,6))
colors = {'MD':'blue','FD':'red'}

for i, row in pca_df.iterrows():
    plt.scatter(row['PC1'], row['PC2'], color=colors[row['Group']], s=200)
    plt.text(row['PC1'] + 0.05, row['PC2'] + 0.05, row['Sample'], fontsize=12)
plt.scatter(pc1_pc2[:,0], pc1_pc2[:,1], s=200)
plt.xlim(pc1_pc2[:,0].min() - 1, pc1_pc2[:,0].max() + 1)
plt.ylim(pc1_pc2[:,1].min() - 1, pc1_pc2[:,1].max() + 1)

plt.xlabel(f"PC1 ({100*eigenvalues[0]/eigenvalues.sum():.1f}%)")
plt.ylabel(f"PC2 ({100*eigenvalues[1]/eigenvalues.sum():.1f}%)")
plt.title("PCA of Normalized Lipidomics Data")
plt.grid(True)
plt.tight_layout()
plt.show()

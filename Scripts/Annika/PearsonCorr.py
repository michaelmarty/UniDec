import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your data
folder = r"Z:\\Group Share\\Annika\\Stellar\\EC\\Replicates"
file_path = f"{folder}\\2AqpZ_Memb_Ext_Disc_Merged.csv"

df = pd.read_csv(file_path)  # row 4 is header

# Define replicates
sampleA = ["ME1", "ME2", "ME3"]
sampleB = ["D1", "D2", "D3"]

# Average replicates for each sample
df["SampleA_mean"] = df[sampleA].mean(axis=1)
df["SampleB_mean"] = df[sampleB].mean(axis=1)

# Pearson correlation between the two mean profiles
corr_value = df[["SampleA_mean", "SampleB_mean"]].corr(method="pearson").iloc[0,1]

print("Pearson correlation between Sample A and B (averaged replicates):", corr_value)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your data
file_path = f"{folder}\\2AqpZ_Memb_Ext_Disc_Merged.csv"
data = pd.read_csv(file_path, index_col=0)  # lipids in rows, samples in columns

# Select only the ME and D replicates
samples = ['ME1', 'ME2', 'ME3', 'D1', 'D2', 'D3']
data_subset = data[samples]

# Transpose so samples are rows and features are columns
X = data_subset.T.values

# Standardize data (mean=0, std=1)
X_mean = X.mean(axis=0)
X_std = X.std(axis=0, ddof=1)
X_scaled = (X - X_mean) / X_std

# Compute covariance matrix
cov_matrix = np.cov(X_scaled, rowvar=False)

# Eigen decomposition
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

# Sort eigenvalues and eigenvectors in descending order
sorted_idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[sorted_idx]
eigenvectors = eigenvectors[:, sorted_idx]

# Project data onto first two principal components
pc1_pc2 = X_scaled @ eigenvectors[:, :2]

# Create DataFrame for plotting
pca_df = pd.DataFrame(pc1_pc2, columns=['PC1', 'PC2'])
pca_df['Sample'] = samples
pca_df['Group'] = ['ME']*3 + ['D']*3

# Plot PCA
plt.figure(figsize=(8,6))
colors = {'ME':'orange', 'D':'purple'}
for group in pca_df['Group'].unique():
    subset = pca_df[pca_df['Group'] == group]
    plt.scatter(subset['PC1'], subset['PC2'], label=group, color=colors[group], s=100)

plt.xlabel(f"PC1 ({100*eigenvalues[0]/eigenvalues.sum():.1f}%)")
plt.ylabel(f"PC2 ({100*eigenvalues[1]/eigenvalues.sum():.1f}%)")
plt.title("PCA of ME vs D Replicates")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

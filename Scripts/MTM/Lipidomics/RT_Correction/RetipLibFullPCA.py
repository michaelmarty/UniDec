import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

path = r"C:\Python\pyRetip\notebooks\lipidomics\SplitCSVs\Pos\PS_data_preprocessed_full_training.csv"
path = r"C:\Python\pyRetip\notebooks\lipidomics\SplitCSVs\Pos_preprocessed_full_training_merged.csv"

df = pd.read_csv(path)

# Select only columns 5 and onwards (assuming 0-based indexing)
data = df.iloc[:, 5:].values

# Drop columns with any NaN values
data = data[:, ~np.isnan(data).any(axis=0)]

# Drop VR1_A and VR2_A columns if they exist
drop_cols = ['VR1_A', 'VR2_A', 'VE2_A', 'VE3_A']
for col in drop_cols:
    if col in df.columns:
        col_index = df.columns.get_loc(col) - 5  # Adjust index for data array
        if 0 <= col_index < data.shape[1]:
            data = np.delete(data, col_index, axis=1)

# Standardize the data (mean=0, variance=1)
# data_meaned = data - np.mean(data, axis=0)
# data_std = data_meaned / np.std(data_meaned, axis=0)

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(data)

# Create a DataFrame with PCA results
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_df['CompoundClass'] = df['CompoundClass']

print(pca_df.head())

# plot
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 8))
plt.subplot(1, 2, 1)
classes = pca_df['CompoundClass'].unique()
colors = plt.cm.get_cmap('tab10', len(classes))

for i, class_name in enumerate(classes):
    subset = pca_df[pca_df['CompoundClass'] == class_name]
    plt.scatter(subset['PC1'], subset['PC2'], label=class_name, color=colors(i))
plt.title('PCA of Lipidomics Data')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend()
plt.grid()

# Top 10 components contributing to PC1
loading_scores = pca.components_[0]
abs_loading_scores = np.abs(loading_scores)
sorted_indices = np.argsort(abs_loading_scores)[::-1]
top_10_indices = sorted_indices[:10]
top_10_features = df.columns[5:][top_10_indices]
top_10_scores = loading_scores[top_10_indices]
plt.subplot(1, 2, 2)
plt.barh(top_10_features, top_10_scores)
plt.title('Top 10 Features Contributing to PC1')
plt.xlabel('Loading Score')
plt.ylabel('Feature')
plt.tight_layout()

plt.show()
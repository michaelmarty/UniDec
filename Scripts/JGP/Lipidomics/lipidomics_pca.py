import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

samples = [1,2,3,4,5,6,7,8,9]

data = pd.DataFrame({
    "Tail Length" : [33.6,33.97,34.34,35.945,35.095,38.905,34.325,41.1,37.945],
    "Tail Unsaturation" : [1.26,1.373,1.26,1.23,1.249,1.17,1.09,1.23,1.14],
    "Percent CL" : [5.41,0.37,2.01,4.82,3.83,11.85,1.71,15.33,9.94],
    "Percent PE" : [79.04,37.39,72.65,62.19,77.78,61.07,86.51,66.00,63.23],
    "Percent PG" : [15.55,62.24,25.34,32.99,18.39,27.08,11.78,18.67,26.83]
})

scaler = StandardScaler()
scaled_data = scaler.fit_transform(data)

pca = PCA(n_components=2)
principal_components = pca.fit_transform(scaled_data)

pca_df = pd.DataFrame(data=principal_components, columns=["PC1", "PC2"])
print(pca_df)


plt.scatter(pca_df["PC1"], pca_df["PC2"])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()

print(f'Explained Variance Ratio: {pca.explained_variance_ratio_}')
import pandas as pd
import numpy as np
import peptidy as pep
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import time

amino_acid_letters = "ACDEFGHIKLMNPQRSTVWY"

def monte_carlo_peptides(num, length):
    peptides = []
    for _ in range(num):
        peptide = ''.join(np.random.choice(list(amino_acid_letters), length))
        peptides.append(peptide)
    return peptides

def generate_descriptors(peptides):
    descriptors = []
    for peptide in peptides:
        desc = pep.descriptors.compute_descriptors(peptide)
        descriptors.append(desc)
    df = pd.DataFrame(descriptors)
    # Drop columns with freq_ in name
    df = df.loc[:, ~df.columns.str.contains('freq_')]
    print(df.columns)
    return df

def pca_of_df(descriptors_df, plot=True):
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(descriptors_df)
    pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    if plot:
        plt.figure(figsize=(8,6))
        plt.subplot(121)
        plt.scatter(pc_df['PC1'], pc_df['PC2'], alpha=0.1)
        plt.title('PCA of Peptide Descriptors')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.grid()

        # Plot feature loadings in bar chart
        plt.subplot(122)
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        # Sort by magnitude of first principal component
        loadings = loadings[np.argsort(-np.abs(loadings[:, 0]))]
        labels = descriptors_df.columns[np.argsort(-np.abs(loadings[:, 0]))]
        # If too many features, limit to first 5
        loadings = loadings[:5, :]
        labels = labels[:5]
        num_features = loadings.shape[0]
        x = np.arange(num_features)
        plt.bar(x - 0.2, loadings[:, 0], width=0.4, label='PC1')
        plt.bar(x + 0.2, loadings[:, 1], width=0.4, label='PC2')
        plt.xticks(x, labels, rotation=90)
        plt.title('PCA Feature Loadings')
        plt.xlabel('Features')
        plt.ylabel('Loading Score')
        plt.legend()
        plt.tight_layout()

        plt.show()

    return pc_df

if __name__ == "__main__":
    starttime = time.perf_counter()
    peptides = monte_carlo_peptides(10000, 7)
    descriptors_df = generate_descriptors(peptides)
    print(descriptors_df)
    endtime = time.perf_counter()
    print("Descriptor generation time:", endtime - starttime)
    pc_df = pca_of_df(descriptors_df)
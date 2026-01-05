import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class_color_map = {
    "PC": "#788BFF",
    "PC O": "#BFD7FF",
    "LPC-O/P": "#BFD7FF",
    "LPC": "#3019E3",
    "PE": "#6017A9",
    "PE-O": "#A276E5",
    "LPE": "#8F5CD2",
    "PG": "#356935",
    "LPG": "#B1D8B1",
    "LPS": "#FDF32F",
    "PS": "#C89504",
    "PA": "#95F4D4",
    "LPA": "#42EBB3",
    "PI": "#8FFBFF",
    "LPI": "#00CAD1",
    "CL": "#FFA200",
    "Cer": "#AB30A1",
    "HexCer": "#F708E3",
    "SM": "#FFDB73",
    "CoQ": "#F25C54",
    "MAG": "#F7B267",
    "DG": "#F4845F",
    "DAG": "#F4845F",
    "TG": "#8F1458",
    "TAG": "#8F1458",
    "FA": "#FFA1DB",
    "CE": "#D10081",
    "NAE": "#777777",
    "ST": "#B40808",
    "Standard": "#0D0D0D"
}

def sum_transitions(df, mode="Products", drop_IS=True, normalize_IS=True, normalize_TIC=True):
    molecules = df["Molecule"].unique()
    replicates = df["Replicate Name"].unique()

    newrows = []
    for m in molecules:
        newrow = {"Molecule": m, "Molecule List Name": df[df["Molecule"] == m]["Molecule List Name"].values[0]}
        for r in replicates:
            mask = (df["Molecule"] == m) & (df["Replicate Name"] == r)
            subset = df[mask]
            if mode == "Products":
                # Drop any with "precursor" in "Fragment Ion"
                subset = subset[~subset["Fragment Ion"].str.contains("precursor", case=False, na=False)]
            elif mode == "Precursors":
                # Keep only those with "precursor" in "Fragment Ion"
                subset = subset[subset["Fragment Ion"].str.contains("precursor", case=False, na=False)]
            elif mode == "Tails":
                # Keep only those with "T" at the start of "Fragment Ion"
                subset = subset[subset["Fragment Ion"].str.startswith("T", na=False)]
            elif mode == "Heads":
                # Keep only those with "H" at the start of "Fragment Ion"
                subset = subset[subset["Fragment Ion"].str.startswith("H", na=False)]
            elif mode == "All":
                pass

            total_area = subset["Area"].sum()
            newrow[r] = total_area
        newrows.append(newrow)

    outdf = pd.DataFrame(newrows)
    if normalize_IS:
        outdf = normalize_is(outdf, replicates)

    if drop_IS:
        # Drop any rows where Molecule starts with "IS_"
        outdf = outdf[~outdf["Molecule"].str.startswith("IS")]

    if normalize_TIC:
        outdf = normalize_tic(outdf, replicates)
    return outdf, replicates


def normalize_tic(df, replicates):
    normdf = df.copy()
    for r in replicates:
        total = normdf[r].sum()
        normdf[r] = normdf[r] / total
    return normdf


def normalize_is(df, replicates):
    normdf = df.copy()
    classes = normdf["Molecule List Name"].unique()
    for r in replicates:
        normdf[r] = normdf[r].astype(float)

    for c in classes:
        mask = normdf["Molecule List Name"] == c
        subdf = normdf[normdf["Molecule List Name"] == c]
        # Find IS as molecule that starts with "IS_"
        isrow = subdf[subdf["Molecule"].str.startswith("IS")]
        for r in replicates:
            isval = float(isrow[r].values[0])
            if isval == 0:
                print(
                    f"Warning: IS value is zero for class {c} in replicate {r}. Skipping normalization for this class and replicate.")
                continue
            normdf.loc[mask, r] = normdf.loc[mask, r].astype(float) / isval
    return normdf


def stats_calcs(df, replicates_set1, replicates_set2, paired=False, bh_correction=True):
    if paired:
        if len(replicates_set1) != len(replicates_set2):
            raise ValueError("For paired t-test, the number of replicates in both sets must be equal.")
        set1_values = np.array(df[replicates_set1], dtype=np.float64)
        set2_values = np.array(df[replicates_set2], dtype=np.float64)
        ratios = set2_values / set1_values
        # Put ratios in their own columns
        for i in range(ratios.shape[1]):
            df[f"Ratio_{i + 1}"] = ratios[:, i]
        ratiomean = ratios.mean(axis=1)
        df["FC"] = ratiomean
        df["Log2FC"] = np.log2(np.array(ratiomean) + 1e-10)
        df["StdDevFC"] = np.std(ratios, axis=1, ddof=1)

    else:
        set1_means = df[replicates_set1].mean(axis=1)
        set2_means = df[replicates_set2].mean(axis=1)
        df["Set1_Mean"] = set1_means
        df["Set2_Mean"] = set2_means
        df["FC"] = set2_means / set1_means
        df["Log2FC"] = np.log2(set2_means + 1e-10) - np.log2(set1_means + 1e-10)
        # Propagate standard deviation
        set1_std = df[replicates_set1].std(axis=1, ddof=1)
        set2_std = df[replicates_set2].std(axis=1, ddof=1)
        df["StdDevFC"] = df["FC"] * np.sqrt((set1_std / set1_means) ** 2 + (set2_std / set2_means) ** 2)

    p_values = []
    for index, row in df.iterrows():
        set1_values = row[replicates_set1].values
        set2_values = row[replicates_set2].values

        set1_values = np.array(set1_values, dtype=np.float64)
        set2_values = np.array(set2_values, dtype=np.float64)

        if paired:
            _, p = stats.ttest_rel(set1_values, set2_values)
        else:
            _, p = stats.ttest_ind(set1_values, set2_values)
        p_values.append(p)
    df["p-value"] = p_values

    if bh_correction:
        corrected = multipletests(df["p-value"], method='fdr_bh')
        df["p-value"] = corrected[1]

    df["Log10p"] = -np.log10(df["p-value"] + 1e-10)

    return df


def class_calcs(df, replicate_set1, replicate_set2, class_col="Molecule List Name", normalize=True, paired=True,
                bh_correction=True):
    classes = df[class_col].unique()
    replicates = replicate_set1 + replicate_set2

    results = []
    for cls in classes:
        newrow = {class_col: cls}
        for r in replicates:
            mask = df[class_col] == cls
            class_sum = df[mask][r].sum()
            newrow[r] = class_sum
        results.append(newrow)
    outdf = pd.DataFrame(results)

    if normalize:
        outdf = normalize_tic(outdf, replicates)

    outdf = stats_calcs(outdf, replicate_set1, replicate_set2, paired=paired, bh_correction=bh_correction)
    print(outdf)
    return outdf


def make_volcano_plot(df, log2fc_col="Log2FC", neglogp_col="Log10p", title="Volcano Plot", ax=None):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)
    # Color by Molecule List Name
    classes = df["Molecule List Name"].unique()
    colors = plt.cm.get_cmap('tab10', len(classes))
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    plt.scatter(df[log2fc_col], df[neglogp_col], alpha=0.7, c=df["Molecule List Name"].map(class_color_map))
    plt.title(title)
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 p-value")
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    plt.axvline(x=1, color='g', linestyle='--')
    plt.axvline(x=-1, color='g', linestyle='--')


def make_class_plots(classdf, ax=None):
    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    classes = classdf["Molecule List Name"].unique()
    colors = plt.cm.get_cmap('tab10', len(classes))
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    colors = [class_color_map[cls] for cls in classdf["Molecule List Name"]]

    # Labels are class name plus p-value
    labels = [f"{cls}\n(p={pval:.2f})" for cls, pval in zip(classdf["Molecule List Name"], classdf["p-value"])]

    plt.bar(classdf["Molecule List Name"], classdf["FC"], yerr=classdf["StdDevFC"], color=colors,
            label=labels, capsize=5)
    # line at y=1
    plt.axhline(y=1, color='r', linestyle='--')
    plt.ylabel("Fold Change")
    plt.title("Class Fold Changes")
    plt.legend()


def pca_plot(df, replicates_set1, replicates_set2, title="PCA Plot", ax=None):
    replicates = replicates_set1 + replicates_set2
    data = df[replicates].values
    data = StandardScaler().fit_transform(data.T)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(data)
    pc_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    pc_df['Group'] = ['Set1'] * len(replicates_set1) + ['Set2'] * len(replicates_set2)

    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    colors = {'Set1': 'b', 'Set2': 'r'}
    for group in pc_df['Group'].unique():
        indicesToKeep = pc_df['Group'] == group
        plt.scatter(pc_df.loc[indicesToKeep, 'PC1'], pc_df.loc[indicesToKeep, 'PC2'], c=colors[group], s=50)
    plt.title(title)
    plt.xlabel(f'PC1 - {pca.explained_variance_ratio_[0] * 100:.2f}%')
    plt.ylabel(f'PC2 - {pca.explained_variance_ratio_[1] * 100:.2f}%')
    plt.legend(pc_df['Group'].unique())


def full_pipeline(filepath, set1, set2, mode="Products", drop_IS=True, normalize_IS=True,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True):

    df = pd.read_csv(filepath)
    normdf, replicates = sum_transitions(df, mode=mode, drop_IS=drop_IS, normalize_IS=normalize_IS,
                                         normalize_TIC=normalize_TIC)
    classdf = class_calcs(normdf, set1, set2, paired=paired)
    normdf = stats_calcs(normdf, set1, set2, paired=paired, bh_correction=bh_correction)

    # Write classdf to CSV
    classdf.to_csv(os.path.splitext(filepath)[0] + "_class_stats.csv", index=False)
    # Write normdf to CSV
    normdf.to_csv(os.path.splitext(filepath)[0] + "_molecule_stats.csv", index=False)

    if plot_results:
        plt.figure(figsize=(12, 5))
        ax1 = plt.subplot(1, 3, 1)
        make_class_plots(classdf, ax=ax1)

        ax2 = plt.subplot(1, 3, 2)
        make_volcano_plot(normdf, title="Volcano Plot", ax=ax2)

        ax3 = plt.subplot(1, 3, 3)
        pca_plot(normdf, set1, set2, title="PCA Plot", ax=ax3)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\DDM 05 AqpZ Molecule Transition Results2.csv"

    set1 = ["2", "3", "4"]
    set2 = ["5_D1_EC_2", "7_D2_2", "9_D3_EC_2"]

    full_pipeline(file, set1, set2, mode="All", drop_IS=True, normalize_IS=False,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True)

import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl
# Import s

mpl.use("WxAgg")

def get_row_dict(id, mdf):
    row = mdf[mdf["Alignment ID"] == id]
    if len(row) == 0:
        return None
    return row.iloc[0].to_dict()

def translate_to_skyline(vdf, mdf):
    """
    Translates a DataFrame from MS-DIAL format to Skyline format. Takes dfs from PeakValue and PeakMaster files.
    """
    rows = []
    for i, row in vdf.iterrows():
        peakid = row["ID"]
        rowdict = get_row_dict(peakid, mdf)

        name = rowdict["Metabolite name"] if rowdict is not None else f"Unknown_{peakid}"
        listname = rowdict["Ontology"] if rowdict is not None else "Unknown"
        fragmentname = "None"
        area = row["Area"]
        replicate = row["File"]

        rows.append({"Molecule List Name": listname,
                          "Molecule": name,
                          "Fragment Name": fragmentname,
                          "Area": area,
                          "Replicate Name": replicate})

    # Create a new DataFrame for Skyline
    sdf = pd.DataFrame(rows)
    return sdf


class_color_map = {
    "PC": "#788BFF",
    "PC O": "#BFD7FF",
    "EtherPC": "#BFD7FF",
    "LPC-O/P": "#27D3F5",
    "LPC": "#3019E3",
    "EtherLPC": "#27D3F5",
    "PE": "#6017A9",
    "PE-O": "#A276E5",
    "EtherPE": "#A276E5",
    "LPE": "#8F5CD2",
    "EtherLPE": "#8F5CD2",
    "PEtOH" : "#D291EF",
    "PMeOH": "#D291EF",
    "PG": "#356935",
    "LPG": "#B1D8B1",
    "LPS": "#FDF32F",
    "PS": "#C89504",
    "PA": "#95F4D4",
    "LPA": "#42EBB3",
    "PI": "#8FFBFF",
    "EtherPI": "#8FFBFF",
    "LPI": "#00CAD1",
    "CL": "#FFA200",
    "MLCL": "#FFDAA5",
    "DLCL": "#FFE6C2",
    "Cer": "#AB30A1",
    "Cer_NS": "#AB30A1",
    "Cer_NDS": "#AB30A1",
    "Cer_EODS": "#AB30A1",
    "Cer_EOS": "#AB30A1",
    "Cer_EBDS": "#AB30A1",
    "Cer_AP": "#AB30A1",
    "Cer_AS": "#AB30A1",
    "CAR": "#AB30A1",
    "HexCer": "#F708E3",
    "HexCer_NS": "#F708E3",
    "HexCer_NDS": "#F708E3",
    "HexCer_HDS": "#F708E3",
    "HexCer_AP": "#F708E3",
    "HexCer_HS": "#F708E3",
    "Hex2Cer": "#FF73FF",
    "SM": "#FFDB73",
    "ASM": "#FFDB73",
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
    "Standard": "#0D0D0D",
    "GM3": "#FF73FF",
    "HBMP": "#FFA200",
    "SHexCer": "#F708E3",
    "OxPI": "#8FFBFF",
    "PI_Cer": "#8FFBFF",
    "PE_Cer": "#6017A9",
    "AHexBRS": "#FF73FF",
    "AHexCAS": "#FF73FF",
    "EtherPG": "#356935",
    "OxPE": "#6017A9",
    "AHexCer": "#FF73FF",
    "EtherLPG": "#B1D8B1",
    "EtherPG": "#356935",
    "EtherTG": "#8F1458",
    "BileAcid": "#FF8800",
    "SL": "#FFDB73",
    "LNAPE": "#8F5CD2",
    "LNAPS": "#FDF32F",
    "EtherDGDG": "#FFA200",
    "EtherMGDG": "#FFA200",
    "DGGA": "#FFA200",
    "SQDG": "#FFA200",
    "MGDG": "#FFA200",
    "OxPS": "#C89504",
}

# python
default_color = "#0D0D0D"
lower_color_map = {k.strip().lower(): v for k, v in class_color_map.items()}

def get_color_for_class(cls):
    key = str(cls).strip().lower()
    color = lower_color_map.get(key)
    if color is None:
        print(f"WARNING: No color found for class '{cls}'. Using default color.")
        return default_color
    return color


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
        if len(outdf) == 0:
            raise ValueError("No data left after IS normalization. Check if IS compounds are present.")

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
        if len(isrow) == 0:
            print(f"WARNING: No IS found for class {c}. Skipping normalization for this class.")
            continue
        for r in replicates:
            isval = float(isrow[r].values[0])
            if isval == 0:
                print(
                    f"WARNING: IS value is zero for class {c} in replicate {r}. Skipping normalization for this class and replicate.")
                continue
            normdf.loc[mask, r] = normdf.loc[mask, r].astype(float) / isval
    return normdf


def stats_calcs(df, replicates_set1, replicates_set2, paired=False, bh_correction=True):

    df["Mean"] = df[replicates_set1 + replicates_set2].mean(axis=1)
    set1_means = df[replicates_set1].mean(axis=1)
    set2_means = df[replicates_set2].mean(axis=1)
    df["Set1_Mean"] = set1_means
    df["Set2_Mean"] = set2_means

    summarycols = []
    if paired:
        if len(replicates_set1) != len(replicates_set2):
            raise ValueError("For paired t-test, the number of replicates in both sets must be equal.")
        set1_values = np.array(df[replicates_set1], dtype=np.float64)
        set2_values = np.array(df[replicates_set2], dtype=np.float64)
        ratios = set2_values / set1_values
        # Put ratios in their own columns
        for i in range(ratios.shape[1]):
            df[f"Ratio_{i + 1}"] = ratios[:, i]
            summarycols.append(f"Ratio_{i + 1}")
        ratiomean = ratios.mean(axis=1)
        df["FC"] = ratiomean
        df["Log2FC"] = np.log2(np.array(ratiomean) + 1e-10)
        df["StdDevFC"] = np.std(ratios, axis=1, ddof=1)
        df["CI95_FC"] = stats.t.ppf(0.975, df=len(replicates_set1) - 1) * (df["StdDevFC"] / np.sqrt(len(replicates_set1)))

    else:
        summarycols.append("Set1_Mean")
        summarycols.append("Set2_Mean")
        df["FC"] = set2_means / set1_means
        df["Log2FC"] = np.log2(set2_means + 1e-10) - np.log2(set1_means + 1e-10)
        # Propagate standard deviation
        set1_std = df[replicates_set1].std(axis=1, ddof=1)
        set2_std = df[replicates_set2].std(axis=1, ddof=1)
        df["StdDevFC"] = df["FC"] * np.sqrt((set1_std / set1_means) ** 2 + (set2_std / set2_means) ** 2)
        df["CI95_FC"] = stats.t.ppf(0.975, df=len(replicates_set1) + len(replicates_set2) - 2) * (
                    df["StdDevFC"] / np.sqrt(len(replicates_set1) + len(replicates_set2)))

    p_values = []
    for index, row in df.iterrows():
        set1_values = row[replicates_set1].values
        set2_values = row[replicates_set2].values

        set1_values = np.array(set1_values, dtype=np.float64)
        set2_values = np.array(set2_values, dtype=np.float64)

        if index == 0:
            print("N Set1:", len(set1_values), "N Set2:", len(set2_values))

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

    return df, summarycols


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

    outdf, summarycols = stats_calcs(outdf, replicate_set1, replicate_set2, paired=paired, bh_correction=bh_correction)
    print(outdf)
    return outdf, summarycols


def make_volcano_plot(df, log2fc_col="Log2FC", neglogp_col="Log10p", title="Volcano Plot", ax=None, xlims=[-12, 12]):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)
    # Color by Molecule List Name
    classes = df["Molecule List Name"].unique()
    # colors = plt.cm.get_cmap('tab10', len(classes))
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    sizes = df["Mean"]
    # Scale sizes to between 20 and 50
    sizes = 20 + (sizes - sizes.min()) / (sizes.max() - sizes.min()) * (50 - 20)
    plt.scatter(df[log2fc_col], df[neglogp_col], alpha=0.7, c=df["Molecule List Name"].map(get_color_for_class), s=sizes)
    plt.title(title)
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 p-value")
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    plt.axvline(x=1, color='g', linestyle='--')
    plt.axvline(x=-1, color='g', linestyle='--')
    plt.xlim(xlims)

    # Add cursor hover to show molecule names
    annot = plt.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                         bbox=dict(boxstyle="round", fc="w"),
                         arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "\n".join([df["Molecule"].iloc[n] for n in ind["ind"]])
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor("lightyellow")
        annot.get_bbox_patch().set_alpha(0.9)
    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                plt.draw()
            else:
                if vis:
                    annot.set_visible(False)
                    plt.draw()

    sc = plt.scatter(df[log2fc_col], df[neglogp_col], alpha=0)  # Invisible scatter for hover
    plt.gcf().canvas.mpl_connect("motion_notify_event", hover)




def make_class_plots(classdf, ax=None, fclimit=10):
    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    # classes = classdf["Molecule List Name"].unique()
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    colors = [get_color_for_class(cls) for cls in classdf["Molecule List Name"]]

    # Labels are class name plus p-value
    labels = [f"{cls}\n(p={pval:.2f})" for cls, pval in zip(classdf["Molecule List Name"], classdf["p-value"])]

    if len(labels) > 10:
        rotation=90
    else:
        rotation=0

    plt.bar(classdf["Molecule List Name"], classdf["FC"], yerr=classdf["StdDevFC"], color=colors,
            label=labels, capsize=5)
    plt.xticks(rotation=rotation)
    # line at y=1
    plt.axhline(y=1, color='r', linestyle='--')
    plt.ylabel("Fold Change")
    plt.title("Class Fold Changes")
    plt.ylim(0, fclimit)
    if len(labels) <= 10:
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

def heatmap_plot(df, replicates, title="Heatmap", ax=None, lindex=-3, colormap='viridis', normtype="minmax", scale="linear"):
    data = df[replicates].values
    vmax = 1-1e-12
    vmin = 1e-12
    if normtype=="minmax":
        # Normalize data to 0-1 for better color mapping
        data = (data - np.min(data)) / (np.max(data) - np.min(data))
    elif normtype=="zscore":
        data = StandardScaler().fit_transform(data)
    elif normtype=="zerocentered":
        data = (data / (2 * np.amax(np.abs(data)))) + 0.5
    elif normtype=="onecentered":
        data = ((data - 1) / (2 * np.amax(np.abs(data - 1)))) + 0.5
    elif normtype=="zerocentered_clip":
        data = data + 0.5
        data[data < vmin] = vmin
        data[data > vmax] = vmax
    elif normtype=="onecentered_clip":
        data = data - 0.5
        data[data < vmin] = vmin
        data[data > vmax] = vmax

    print("Min max of heatmap data:", np.min(data), np.max(data))
    if ax is None:
        plt.figure(figsize=(10, 8))
    else:
        plt.sca(ax)

    plt.imshow(data, aspect='auto', cmap=colormap, norm=scale, vmax=vmax, vmin=vmin)
    # plt.colorbar(label='Normalized Intensity')
    plt.yticks(ticks=np.arange(len(df)), labels=df['Molecule'], fontsize=5)
    # Labels are last 3 characters of replicates
    labels = [r[lindex:] for r in replicates]
    plt.xticks(ticks=np.arange(len(replicates)), labels=labels, rotation=0, fontsize=8)

    # Draw colored rectangles around each class
    classes = df["Molecule List Name"].values
    class_changes = [0]
    for i in range(1, len(classes)):
        if classes[i] != classes[i - 1]:
            class_changes.append(i)
    class_changes.append(len(classes))
    for i in range(len(class_changes) - 1):
        start = class_changes[i] - 0.5
        end = class_changes[i + 1] - 0.5
        rect = mpl.patches.Rectangle(( -0.5, start), len(replicates), end - start,
                                     linewidth=1, edgecolor=get_color_for_class(classes[class_changes[i]]),
                                     facecolor='none')
        plt.gca().add_patch(rect)
        # print(start, end)
        # # Add a thick colored line to the right of the rectangle
        # line = mpl.lines.Line2D([len(replicates) - 0.5, len(replicates) - 0.5], [start, end],
        #                         linewidth=3, color=class_color_map[classes[class_changes[i]]])
        # plt.gca().add_line(line)

    plt.title(title)

def pie_chart(df, set, ax=None, title="Class Distribution", otherthresh=0.035):
    class_sums = df.groupby("Molecule List Name")[set].sum()
    if ax is None:
        plt.figure(figsize=(8, 8))
    else:
        plt.sca(ax)

    # For sums that are less than 1% of total, group into "Other"
    total = class_sums.sum()
    b1 = class_sums / total >= otherthresh
    class_sums_filtered = class_sums[b1]
    colors = [get_color_for_class(cls) for cls in class_sums_filtered.index]
    other_sum = total - class_sums_filtered.sum()
    if other_sum > 0:
        class_sums_filtered["Other"] = other_sum
        colors.append("#CCCCCC")


    plt.pie(class_sums_filtered, labels=class_sums_filtered.index, colors=colors, autopct='%1.0f%%', startangle=140)
    plt.title(title)


def full_pipeline(filepath, set1, set2, mode="Products", drop_IS=True, normalize_IS=True,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True, heatplots=True,
                  piecharts=True, figsavepath=None):

    df = pd.read_csv(filepath)
    normdf, replicates = sum_transitions(df, mode=mode, drop_IS=drop_IS, normalize_IS=normalize_IS,
                                         normalize_TIC=normalize_TIC)
    classdf, classsummarycols = class_calcs(normdf, set1, set2, paired=paired)
    normdf, summarycols = stats_calcs(normdf, set1, set2, paired=paired, bh_correction=bh_correction)

    # Write to two sheets in an Excel file
    with pd.ExcelWriter(os.path.splitext(filepath)[0] + "_stats.xlsx") as writer:
        classdf.to_excel(writer, sheet_name="Class Stats", index=False)
        normdf.to_excel(writer, sheet_name="Molecule Stats", index=False)
        df.to_excel(writer, sheet_name="Raw Data", index=False)

    # # Write classdf to CSV
    # classdf.to_csv(os.path.splitext(filepath)[0] + "_class_stats.csv", index=False)
    # # Write normdf to CSV
    # normdf.to_csv(os.path.splitext(filepath)[0] + "_molecule_stats.csv", index=False)

    if plot_results:
        plt.figure(figsize=(18, 10))
        nrow = 1
        if heatplots:
            nrow +=1
        if piecharts:
            nrow +=1
        ncol = 3
        ax1 = plt.subplot(nrow, ncol, 1)
        make_class_plots(classdf, ax=ax1)

        ax2 = plt.subplot(nrow, ncol, 2)
        make_volcano_plot(normdf, title="Volcano Plot", ax=ax2)

        ax3 = plt.subplot(nrow, ncol, 3)
        pca_plot(normdf, set1, set2, title="PCA Plot", ax=ax3)

        i=3
        if heatplots:
            ax4 = plt.subplot(nrow, ncol, i+1)
            # Sort normdf by class and then by Mean
            normdfsorted = normdf.sort_values(by=["Molecule List Name", "Mean"], ascending=[True, False])
            heatmap_plot(normdfsorted, replicates, title="Heatmap", ax=ax4, normtype="minmax", scale="linear")

            ax5 = plt.subplot(nrow, ncol, i+2)
            if paired:
                titleval = "Paired Ratio Heatmap"
                normtype = "onecentered"
                cmap = "RdBu_r"
            else:
                titleval = "Set Means Heatmap"
                normtype = "minmax"
                cmap = "viridis"
            heatmap_plot(normdfsorted, summarycols, title=titleval, ax=ax5, lindex=0, colormap=cmap, normtype=normtype)

            ax6 = plt.subplot(nrow, ncol, i+3)
            cols = ["Log2FC", "StdDevFC"]
            heatmap_plot(normdfsorted, cols, title="Fold Change Heatmap", ax=ax6, lindex=0, colormap="RdBu_r", normtype="zerocentered")
            i += 3

        if piecharts:
            ax7 = plt.subplot(nrow, ncol, i+1)
            pie_chart(classdf, "Set1_Mean", ax=ax7, title="Class Distribution Set 1")
            ax8 = plt.subplot(nrow, ncol, i+2)
            pie_chart(classdf, "Set2_Mean", ax=ax8, title="Class Distribution Set 2")

        plt.tight_layout()

        if figsavepath is not None:
            plt.savefig(figsavepath)

        plt.show()


if __name__ == "__main__":
    file = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\DDM 05 AqpZ Molecule Transition Results2.csv"

    set1 = ["2", "3", "4"]
    set2 = ["5_D1_EC_2", "7_D2_2", "9_D3_EC_2"]

    full_pipeline(file, set1, set2, mode="All", drop_IS=True, normalize_IS=False,
                  normalize_TIC=True, paired=True, bh_correction=True, plot_results=True)

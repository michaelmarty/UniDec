import os
from PublicScripts.Lipids.LipidFunctions import *
import PublicScripts.Lipids.LipidOutlierAnalysis as ot
from venn import venn
import matplotlib.pyplot as plt
import matplotlib as mpl

def prep_datasets(topdir, datasets, min_count=None, drop_d7=True, write_output=True,
                  dtype="msdial", rttol=0.2):
    os.chdir(topdir)
    dfs = []
    # Data sets
    for d in datasets:
        if dtype=="msdial":
            mergeddf, posdf, negdf = cleanup_and_merge_msdial(f"{d}_pos.csv", f"{d}_neg.csv", f"{d}_pos.xml",
                                                              f"{d}_neg.xml", rttol=rttol)
        elif dtype=="mzmine":
            mergeddf, posdf, negdf = cleanup_and_merge_mzmine(f"{d}_Pos_features_data_msdial.csv", f"{d}_Neg_features_data_msdial.csv", rttol=rttol)
        # Write out merged dataframes for reference
        # mergeddf.to_csv(f"{d}_Merged.csv", index=False)
        mergeddf["Dataset"] = d
        dfs.append(mergeddf)

    # Combine all datasets
    combineddf = pd.concat(dfs)
    # Sort by class, then by name
    combineddf = combineddf.sort_values(by=["Ontology", "Metabolite name"])

    if drop_d7:
        # Drop any with (d7) in name
        combineddf = combineddf[~combineddf["Metabolite name"].str.contains(r"\(d7\)", na=False)]

    # Drop any where the class has less than min_count total entries across all datasets, to focus on more common classes
    if min_count is not None:
        class_counts = combineddf["Ontology"].value_counts()
        classes_to_keep = class_counts[class_counts >= 2].index
        combineddf = combineddf[combineddf["Ontology"].isin(classes_to_keep)]

    if write_output:
        # Write out combined dataframe for reference
        combineddf.to_csv("Combined_DDA_Full.csv", index=False)
    return combineddf

def full_name_assignment(df, rttol=0.1, rtcol="Average Rt(min)"):
    newrows = []
    for i, row in df.iterrows():
        name = row["Metabolite name"]
        if "|" in name:
            simple_name = name.split("|")[0]
        else:
            simple_name = name

        candidates = df[df["Metabolite name"].str.startswith(simple_name)]
        # Take only candidates not equal to the current row
        candidates = candidates[candidates["Metabolite name"] != name]
        if len(candidates) > 0:
            rt = row[rtcol]
            rtdiffs = abs(candidates[rtcol] - rt)
            if any(rtdiffs < rttol):
                # Filter candidates below rttol and sort by RTdiff
                candidates = candidates.assign(RTdiff=rtdiffs)
                candidates = candidates[candidates["RTdiff"] < rttol].sort_values(by="RTdiff")
                for j, candidate in candidates.iterrows():
                    # For cls and other partially named candidates, check if the candidate name is longer
                    if len(candidate["Metabolite name"]) > len(name):
                        # print(f"Assigning {candidate['Metabolite name']} to {row['Metabolite name']} based on RT match within {rttol} min")
                        row["Metabolite name"] = candidate["Metabolite name"]
                        break

        newrows.append(row)

    df = pd.DataFrame(newrows)
    return df

def get_number_of_datasets(df, namecol="Metabolite name"):
    df = df.copy()
    for i, row in df.iterrows():
        name = row[namecol]
        count = len(df[df[namecol] == name]["Dataset"].unique())
        df.at[i, "Dataset Count"] = count
    return df

def split_dfs(df):
    max_num_datasets = int(df["Dataset Count"].max())
    dfs = {}
    # For everything with 2, 3, 4, etc. datasets, split into separate dfs
    for i in range(2, max_num_datasets + 1):
        subdf = df[df["Dataset Count"] == i]
        dfs[i] = subdf
    # For everything with only 1 dataset, split by dataset
    subdf = df[df["Dataset Count"] == 1]
    for dataset in subdf["Dataset"].unique():
        dataset_subdf = subdf[subdf["Dataset"] == dataset]
        dfs[dataset] = dataset_subdf
    return dfs

def write_dfs_to_excel(df, filename="Combined_DDA_Full.xlsx"):
    subdfs = split_dfs(df)

    # Write df to first sheet, then each subdf to separate sheets in the same Excel file
    with pd.ExcelWriter(filename) as writer:
        df.to_excel(writer, sheet_name="Full Data", index=False)
        for key, subdf in subdfs.items():
            if isinstance(key, int):
                sheet_name = f"{key} Datasets"
            else:
                sheet_name = f"{key} Only"
            subdf.to_excel(writer, sheet_name=sheet_name, index=False)
        # Split df by dataset and write each to a separate sheet
        for dataset in df["Dataset"].unique():
            dataset_subdf = df[df["Dataset"] == dataset]
            sheet_name = f"{dataset} Total"
            dataset_subdf.to_excel(writer, sheet_name=sheet_name, index=False)

def dataset_split(df, rttol=0.2, rtcol="Average Rt(min)", namecol="Metabolite name", write_output=True):
    df = full_name_assignment(df, rttol=rttol, rtcol=rtcol)
    df = get_number_of_datasets(df, namecol=namecol)
    if write_output:
        write_dfs_to_excel(df)
    return df

def outlier_setup_dda(df, mztol=0.3, write_output=True, do_heads=True, do_tails=False):
    df = smiles_analysis(df)

    posdf = df[df["Polarity"] == "Positive"]
    negdf = df[df["Polarity"] == "Negative"]

    outdf = ot.outlier_analysis(posdf, negdf, tol=mztol, contribs=False, add_all_cols_heads=True,
                                add_all_cols_tails=False, drop_cols=False, do_tails=do_tails, do_heads=do_heads,
                                auto_tolcol=True)

    if write_output:
        # Write out combined dataframe for reference
        outdf.to_excel("Combined_DDA_OutlierAnalysis.xlsx", index=False)

    return outdf


def venn_diagram(df, sets, namecol="Metabolite name", title="", ax=None):
    vdata = {}
    for s in sets:
        subdf = df[df["Dataset"] == s]
        unique_names = subdf[namecol].unique()
        l = len(unique_names)
        vdata[s + " " + str(l)] = set(unique_names)

    if ax is None:
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
    venn(vdata, ax=ax)
    plt.title(title)

def count_number_of_datasets(df, namecol="Metabolite name", dataset_col="Dataset", class_col="Ontology"):
    # If Dataset Count column already exists, drop it to avoid confusion
    if "Dataset Count" in df.columns:
        df = df.drop(columns=["Dataset Count"])
    # Count the number of unique datasets each metabolite is found in
    dataset_counts = df.groupby(namecol)[dataset_col].nunique()
    # Addback into the dataframe
    df = df.merge(dataset_counts.rename("Dataset Count"), left_on=namecol, right_index=True)
    return df

def get_class_counts(df, name_col="Metabolite name", dataset_col="Dataset", class_col="Ontology"):
    datasets = sorted(df[dataset_col].unique())
    counts = []

    for i in range(2, len(datasets) + 1):
        subdf = df[df["Dataset Count"] == i]
        # Drop repeated names
        subdf = subdf.drop_duplicates(subset=[class_col, name_col])
        class_counts = subdf.groupby([dataset_col, class_col]).size().unstack(fill_value=0)
        # Sum across datasets to get total counts per class for coloring
        class_counts = class_counts.sum(axis=0)
        print(f"Number of unique names in {i} datasets:", subdf[name_col].nunique())
        cdf = pd.DataFrame(class_counts, columns=[str(i) + " Datasets"]).transpose()
        counts.append(cdf)

    counts = counts[::-1]  # Reverse to have higher counts on bottom of stack
    class_counts_1 = df[df["Dataset Count"] == 1].groupby([dataset_col, class_col]).size().unstack(fill_value=0)
    counts.append(class_counts_1)
    # Combine counts for 1, 2, and 3 datasets
    class_counts = pd.concat(counts, axis=1).fillna(0)
    # Merge duplicate classes together by summing across them
    class_counts = class_counts.T.groupby(class_counts.columns).sum().T
    # Rename rows that have a dataset in them to say dataset + only
    newrows = []
    for idx in class_counts.index:
        if any(dataset in idx for dataset in datasets):
            newrows.append(idx + " Only")
        else:
            newrows.append(idx)
    class_counts.index = newrows
    return class_counts

def bar_chart_of_classes(
        df,
        datasets,
        dataset_col="Dataset",
        class_col="Ontology",
        name_col="Metabolite name",
        title="",
        ax=None,
        legend_fontsize=8,
        legend_bbox_to_anchor=(1.05, 1),
        legend_loc="upper left",
        hide_top_and_right=False
):
    class_counts = get_class_counts(df, name_col=name_col, dataset_col=dataset_col, class_col=class_col)
    # Sort order of class counts rows to match datasets order with only at the end, add the other classes at the beginning in the same order they appear in the dataframe
    new_order = []
    for idx in class_counts.index:
        if "Only" not in idx:
            new_order.append(idx)
    for dataset in datasets:
        only_row = dataset + " Only"
        if only_row in class_counts.index:
            new_order.append(only_row)
    class_counts = class_counts.reindex(new_order)

    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    colors = [class_color_map.get(cls, "#333333") for cls in class_counts.columns]
    plot_ax = class_counts.plot(kind="bar", stacked=True, ax=ax, color=colors)
    # rotate x tick labels by 45 degrees
    plot_ax.tick_params(axis="x", rotation=45)
    for label in plot_ax.get_xticklabels():
        label.set_horizontalalignment("right")

    plot_ax.set_title(title)
    plot_ax.set_ylabel("Count")
    # plt.xlabel(dataset_col)
    plot_ax.legend(
        bbox_to_anchor=legend_bbox_to_anchor,
        loc=legend_loc,
        fontsize=legend_fontsize,
        reverse=True
    )

    if hide_top_and_right:
        plot_ax.spines["top"].set_visible(False)
        plot_ax.spines["right"].set_visible(False)


def shared_unique_plot(df, sets, namecol="Metabolite name", title="", ax=None):
    """Plot each dataset total as unique-only plus shared confirmed IDs."""
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7))

    dataset_counts = df.groupby(namecol)["Dataset"].nunique()
    unique_names = dataset_counts[dataset_counts == 1].index

    total_counts = (
        df.groupby("Dataset")[namecol]
        .nunique()
        .reindex(sets, fill_value=0)
    )
    unique_counts = (
        df[df[namecol].isin(unique_names)]
        .groupby("Dataset")[namecol]
        .nunique()
        .reindex(sets, fill_value=0)
    )
    shared_counts = total_counts - unique_counts

    x_positions = list(range(len(sets)))
    unique_bars = ax.bar(
        x_positions,
        unique_counts.values,
        color="#E69F00",
        edgecolor="black",
        linewidth=1,
        label="Unique to Dataset"
    )
    shared_bars = ax.bar(
        x_positions,
        shared_counts.values,
        bottom=unique_counts.values,
        color="#035D99",
        edgecolor="black",
        linewidth=1,
        label="Shared with Another Dataset"
    )

    ax.bar_label(
        shared_bars,
        labels=[f"{value} shared" for value in shared_counts.values],
        label_type="center",
        color="white",
        fontsize=11
    )
    ax.bar_label(
        shared_bars,
        labels=[str(value) for value in total_counts.values],
        padding=4,
        fontsize=14
    )
    for bar, value in zip(unique_bars, unique_counts.values):
        ax.annotate(
            f"{value} unique",
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11,
            color="#7A4C00",
            fontweight="bold"
        )

    ax.set_xticks(x_positions, labels=sets, rotation=45, ha="right")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Confirmed Lipid IDs")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(y=0.12)

    return total_counts, unique_counts, shared_counts


def upset_diagram(df, sets, namecol="Metabolite name", title="", subplot_spec=None):
    """Plot exact set intersections and dataset totals without extra dependencies."""
    if subplot_spec is None:
        fig = plt.figure(figsize=(12, 8))
        subplot_spec = fig.add_gridspec(1, 1)[0]
    else:
        fig = plt.gcf()

    membership_counts = {}
    for _, group in df.groupby(namecol, sort=False):
        observed_sets = set(group["Dataset"])
        membership = tuple(dataset for dataset in sets if dataset in observed_sets)
        if membership:
            membership_counts[membership] = membership_counts.get(membership, 0) + 1

    intersections = sorted(
        membership_counts.items(),
        key=lambda item: (-item[1], -len(item[0]), item[0])
    )

    grid = subplot_spec.subgridspec(
        2,
        2,
        height_ratios=(3, 2),
        width_ratios=(2, 7),
        hspace=0.05,
        wspace=0.05
    )
    label_ax = fig.add_subplot(grid[0, 0])
    intersection_ax = fig.add_subplot(grid[0, 1])
    set_size_ax = fig.add_subplot(grid[1, 0])
    matrix_ax = fig.add_subplot(grid[1, 1], sharex=intersection_ax)
    label_ax.axis("off")

    x_positions = list(range(len(intersections)))
    intersection_values = [count for _, count in intersections]
    intersection_bars = intersection_ax.bar(
        x_positions,
        intersection_values,
        color="#035D99",
        edgecolor="black",
        linewidth=0.8
    )
    intersection_ax.bar_label(
        intersection_bars,
        padding=2,
        fontsize=10,
        rotation=90
    )
    intersection_ax.set_ylabel("Intersection Size")
    intersection_ax.set_title(title)
    intersection_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    intersection_ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    intersection_ax.spines["top"].set_visible(False)
    intersection_ax.spines["right"].set_visible(False)
    intersection_ax.margins(y=0.15)

    row_positions = list(range(len(sets)))
    set_sizes = [df.loc[df["Dataset"] == dataset, namecol].nunique() for dataset in sets]
    set_bars = set_size_ax.barh(
        row_positions,
        set_sizes,
        color="#076E3D",
        edgecolor="black",
        linewidth=0.8
    )
    set_size_ax.bar_label(set_bars, padding=3, fontsize=12)
    set_size_ax.set_yticks(row_positions, labels=sets)
    set_size_ax.tick_params(axis="y", length=0)
    set_size_ax.invert_yaxis()
    set_size_ax.invert_xaxis()
    set_size_ax.set_xlabel("Dataset Size")
    set_size_ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, nbins=4))
    set_size_ax.spines["top"].set_visible(False)
    set_size_ax.spines["left"].set_visible(False)

    for row in row_positions:
        if row % 2 == 0:
            matrix_ax.axhspan(row - 0.5, row + 0.5, color="#F2F2F2", zorder=0)

    for column, (membership, _) in enumerate(intersections):
        active_rows = [sets.index(dataset) for dataset in membership]
        matrix_ax.scatter(
            [column] * len(row_positions),
            row_positions,
            s=28,
            color="#D3D3D3",
            zorder=1
        )
        if len(active_rows) > 1:
            matrix_ax.plot(
                [column, column],
                [min(active_rows), max(active_rows)],
                color="black",
                linewidth=1.5,
                zorder=2
            )
        matrix_ax.scatter(
            [column] * len(active_rows),
            active_rows,
            s=48,
            color="black",
            zorder=3
        )

    matrix_ax.set_yticks(row_positions, labels=[])
    matrix_ax.tick_params(axis="y", left=False, labelleft=False)
    matrix_ax.set_xticks([])
    matrix_ax.set_xlabel("Exact Dataset Combination")
    matrix_ax.set_ylim(len(sets) - 0.5, -0.5)
    matrix_ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

    return membership_counts


def compare_classes_plot(
        df,
        datasets=None,
        use_simple_names=True,
        use_simple_classes=True,
        drop_low_quality=True,
        output_basename="DDA_Overlap_and_Class_Distribution",
        overlap_style="venn",
        class_legend_fontsize=8,
        class_legend_bbox_to_anchor=(1.05, 1),
        class_legend_loc="upper left",
        hide_class_top_and_right=False
):
    df = df.copy()
    if datasets is None:
        datasets = sorted(df["Dataset"].unique())

    if drop_low_quality:
        # Remove low quality ids where comment includes "Low quality"
        df = df[~df["Comment"].str.contains("Low quality", case=False, na=False)]

    if use_simple_names:
        simpnames = df["Metabolite name"].apply(lambda x: x.split("|")[0] if "|" in x else x)
        df["Simple name"] = simpnames
        namecol = "Simple name"
    else:
        namecol = "Metabolite name"

    # Drop duplicate rows based on dataset and class to avoid counting the same metabolite multiple times
    df = df.drop_duplicates(subset=["Dataset", namecol])

    # Simplify classes by applying hg_simplifier dict to Ontology column
    df["Ontology"] = df["Ontology"].replace(hg_simplifier)

    if use_simple_classes:
        # Simplify classes by applying class_simplifier dict to Ontology column
        df = simplify_class_df(df)

    # Count number of datasets for each metabolite name and add as a column to the dataframe
    df = count_number_of_datasets(df, namecol=namecol)
    # Rename datasets to add / to OTOT, OTIT, and ITIT
    df["Dataset"] = df["Dataset"].replace({"OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT"})
    datasets = [d.replace("OTOT", "OT/OT").replace("OTIT", "OT/IT").replace("ITIT", "IT/IT") for d in datasets]

    if overlap_style == "shared_unique":
        fig, axes = plt.subplots(
            ncols=2,
            figsize=(18, 8),
            layout="constrained"
        )
        shared_unique_plot(
            df,
            datasets,
            namecol=namecol,
            title="Confirmed Lipid IDs: Unique vs Shared",
            ax=axes[0]
        )
        bar_chart_of_classes(
            df,
            datasets,
            dataset_col="Dataset",
            class_col="Ontology",
            title="Class Distribution by Overlap Category",
            name_col=namecol,
            ax=axes[1],
            legend_fontsize=class_legend_fontsize,
            legend_bbox_to_anchor=class_legend_bbox_to_anchor,
            legend_loc=class_legend_loc,
            hide_top_and_right=hide_class_top_and_right
        )
    elif overlap_style == "upset":
        fig = plt.figure(figsize=(22, 10), layout="constrained")
        outer_grid = fig.add_gridspec(
            1,
            2,
            width_ratios=(1.5, 1),
            wspace=0.25
        )
        upset_diagram(
            df,
            datasets,
            namecol=namecol,
            title="Overlap of Confirmed Lipid IDs",
            subplot_spec=outer_grid[0]
        )
        class_ax = fig.add_subplot(outer_grid[1])
        bar_chart_of_classes(
            df,
            datasets,
            dataset_col="Dataset",
            class_col="Ontology",
            title="Class Distribution by Overlap Category",
            name_col=namecol,
            ax=class_ax,
            legend_fontsize=class_legend_fontsize,
            legend_bbox_to_anchor=class_legend_bbox_to_anchor,
            legend_loc=class_legend_loc,
            hide_top_and_right=hide_class_top_and_right
        )
    elif overlap_style == "venn":
        fig = plt.figure(figsize=(18, 10))
        plt.subplot(1, 2, 1)
        venn_diagram(df, datasets, namecol=namecol, title="Overlap of Confirmed Lipid IDs", ax=plt.gca())
        plt.subplot(1, 2, 2)
        bar_chart_of_classes(
            df,
            datasets,
            dataset_col="Dataset",
            class_col="Ontology",
            title="Class Distribution by Dataset",
            name_col=namecol,
            ax=plt.gca(),
            legend_fontsize=class_legend_fontsize,
            legend_bbox_to_anchor=class_legend_bbox_to_anchor,
            legend_loc=class_legend_loc,
            hide_top_and_right=hide_class_top_and_right
        )
    else:
        raise ValueError(
            "overlap_style must be 'venn', 'upset', or 'shared_unique'"
        )

    if overlap_style == "venn":
        fig.tight_layout()

    fig.savefig(f"{output_basename}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{output_basename}.pdf", bbox_inches="tight")
    plt.show()

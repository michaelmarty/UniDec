import os
import sys
import pandas as pd
import matplotlib

matplotlib.use('WxAgg')
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

data = [
    ("Brodbelt", 2, 2),
    ("Marty", 2, 2),
    ("Mullins", 2, 2),
    ("Ren", 3, 3),
    ("Samanta", 2, 2),
    ("Callmann", 1, 1),
    ("Hsu", 1, 3),
    ("Lu", 2, 2),
    ("Newberry", 1, 2),
    ("Aubrey", 1, 1),
    ("Carsch", 5, 5),
    ("Humphrey", 2, 2),
    ("Anslyn", 0, 0),
    ("Krische", 4, 4),
    ("Page", 1, 2),
    ("Baiz", 1, 1),
    ("Henkelman", 1, 2),
    ("Lin", 3, 3),
    ("Makarov", 1, 2),
    ("Raccah", 2, 4),
    ("Roberts", 1, 1),
    ("Thirumalai", 0, 1),
    ("Webb", 2, 2),
    ("Li", 1, 1),
    ("Brenna", 1, 1),
    ("Ellington", 0, 0),
    ("Milliron", 0, 0),
    ("Yu", 0, 0),
    ("Rose", 0, 0),
    ("Que", 1, 1),
    ("Hull", 0, 2),
    ("Sessler", 1, 2),
    ("Belardi", 1, 1),
]
facdf = pd.DataFrame(data, columns=["Names", "MIN", "MAX"])
facnames = facdf["Names"].tolist()


# Add one to the MAX to make it inclusive
# facdf["MAX"] = facdf["MAX"] + 1

def split_names(cell):
    """Split a cell (string) by commas into cleaned name tokens.
    Handles NaN and empty strings. Trims whitespace and ignores empty tokens.
    """
    if pd.isna(cell):
        return []
    # If already a list, flatten
    if isinstance(cell, list):
        tokens = cell
    else:
        tokens = str(cell).split(",")
    # Strip whitespace and filter empties
    cleaned = [t.strip() for t in tokens if t is not None and str(t).strip() != ""]
    # If space in any split that
    final_tokens = []
    for token in cleaned:
        sub_tokens = token.split()
        final_tokens.extend([st.strip() for st in sub_tokens if st.strip() != ""])
    # Drop any that are "-"
    final_tokens = [t for t in final_tokens if t != "--"]
    return final_tokens


def compute_counts(df):
    """Compute per-name counts using row-wise presence:
    - names: sorted list of unique names (by row-wise counts desc)
    - counts1: list of row-wise counts (number of rows where name appears in either column)
    - counts2: None (we plot only row-wise counts). For backward compatibility we still return a same-length list of zeros if requested by plotting logic.
    - counts_df: DataFrame with name, row_count, col1_count, col2_count
    """
    c_row = Counter()

    # Iterate rows and count presence per row (deduplicate within the row)
    rows = []
    for _, row in df.iterrows():
        names_in_row_col1 = set(split_names(row["col1"]))
        names_in_row_col2 = set(split_names(row["col2"]))
        # union of names in this row across both columns
        names_in_row = names_in_row_col1 | names_in_row_col2
        # print(names_in_row)
        # increment row-wise counter once per name present in this row
        for name in names_in_row:
            c_row[name] += 1
        row["Names"] = ",".join(names_in_row)
        rows.append(row)

    df = pd.DataFrame(rows)

    # Unique names union
    all_names = set(c_row.keys())

    # Sort by row-wise counts descending
    names_sorted = sorted(all_names, key=lambda n: (-(c_row[n]), n))

    counts1 = [c_row[name] for name in names_sorted]

    return names_sorted, counts1, df


def plot_grouped_counts(names, counts1, title=None):
    """Create a grouped bar chart for counts. If counts2 is None, plot a single-series bar chart.
    Saves the figure to outpath if provided, and also shows it.
    """
    x = range(len(names))
    plt.figure(figsize=(max(8, len(names) * 0.4), 6))
    plt.bar(x, counts1, width=0.6, label="Row-wise count", color="#4C72B0")
    plt.xticks(x, names, rotation=90)
    plt.xlabel("Name")
    plt.ylabel("Count")
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.show()


def check_full(facname, faclist):
    if facname not in facdf["Names"].values:
        return False
    maxtaking = facdf[facdf["Names"] == facname]["MAX"].values[0]
    currentcount = faclist.count(facname)
    return currentcount >= maxtaking


def check_all_full(facinterest, faclist):
    for fac in facinterest:
        if not check_full(fac, faclist):
            return False
    return True


def monte_carlo(df, matriculation_range=[.4, .6], num=1000, limit=True):
    """Perform Monte Carlo simulation of faculty assignments based on matriculation rates."""
    faculty = df["Names"].str.split(",").explode().str.strip().unique()
    faculty = [f for f in faculty if f != ""]
    faculty.append("FAILED")  # To account for students who cannot be assigned
    faccounts = []
    failed_interests = []
    for i in range(num):
        # Uniform:
        matriculation_rate = np.random.uniform(matriculation_range[0], matriculation_range[1])
        # Gaussian
        # matriculation_rate = np.random.normal(loc=(matriculation_range[0] + matriculation_range[1]) / 2,
        #     scale=(matriculation_range[1] - matriculation_range[0]) / 3 )
        # Select students based on the matriculation rate
        sampled_df = df.sample(frac=matriculation_rate, replace=True)
        # Randomly shuffle the sampled dataframe
        sampled_df = sampled_df.sample(frac=1).reset_index(drop=True)
        # For each student in the sampled_df, pick a random faculty
        faculty_list = []
        for _, row in sampled_df.iterrows():
            facinterest = row["Names"].split(",")
            allfull = check_all_full(facinterest, faculty_list)
            if allfull:
                chosen_faculty = "FAILED"
                faculty_list.append(chosen_faculty)
                failed_interests.append(facinterest)
                continue

            if len(facinterest) == 0:
                raise ValueError("No faculty interest found for a student.")
            chosen_faculty = np.random.choice(facinterest)

            labfull = False

            if limit:
                labfull = check_full(chosen_faculty, faculty_list)

            while not allfull and labfull:
                # if labfull:
                #     print("Initial choice:", chosen_faculty, "Lab full:", labfull, "Faculty list:", faculty_list)
                chosen_faculty = np.random.choice(facinterest)
                if limit:
                    labfull = check_full(chosen_faculty, faculty_list)

            if allfull:
                chosen_faculty = "FAILED"
                failed_interests.append(facinterest)
            faculty_list.append(chosen_faculty)

        # For each faculty, count how many times they were chosen
        fac_counter = Counter(faculty_list)
        flist = []
        for f in faculty:
            flist.append(fac_counter[f])
        faccounts.append(flist)

    # Make a bar plot of the mean and std dev of faculty counts
    faccounts = np.array(faccounts)
    fac_means = np.mean(faccounts, axis=0)
    fac_stds = np.std(faccounts, axis=0)
    # Sort by means descending
    sorted_indices = np.argsort(-fac_means)
    faculty = [faculty[i] for i in sorted_indices]
    fac_means = fac_means[sorted_indices]
    fac_stds = fac_stds[sorted_indices]

    # Drop any who are zero mean
    nonzero_indices = [i for i in range(len(fac_means)) if fac_means[i] > 0]
    faculty = [faculty[i] for i in nonzero_indices]
    fac_means = fac_means[nonzero_indices]
    fac_stds = fac_stds[nonzero_indices]

    # Print max and min values for each faculty
    for i, fac in enumerate(faculty):
        print(f"{fac}: Max={np.max(faccounts[:, sorted_indices[i]])}, Min={np.min(faccounts[:, sorted_indices[i]])}")

    x = range(len(faculty))
    plt.figure(figsize=(max(8, len(faculty) * 0.4), 6))
    plt.subplot(2,1,1)
    plt.bar(x, fac_means, yerr=fac_stds, capsize=5, color="#55A868", alpha=0.7)

    # Add lines on each bar for min and max from facdf
    for i, fac in enumerate(faculty):
        if fac == "FAILED" or fac not in facdf["Names"].values:
            continue
        max_allowed = facdf[facdf["Names"] == fac]["MAX"].values[0]
        min_allowed = facdf[facdf["Names"] == fac]["MIN"].values[0]

        # Use fill betewen for min and max
        plt.fill_between([x[i] - 0.2, x[i] + 0.2], [min_allowed, min_allowed], [max_allowed, max_allowed], color='gray',
                         alpha=0.2, label='Allowed Range' if i == 0 else "")
        plt.plot([x[i] - 0.2, x[i] + 0.2], [max_allowed, max_allowed], color='r', linestyle='-',
                 label='Max Allowed' if i == 0 else "")
        plt.plot([x[i] - 0.2, x[i] + 0.2], [min_allowed, min_allowed], color='b', linestyle='-',
                 label='Min Allowed' if i == 0 else "")
    plt.xticks(x, faculty, rotation=90)
    plt.xlabel("Faculty")
    plt.ylabel("Average Assigned Students")
    plt.title("Monte Carlo Simulation of Faculty Assignments")

    plt.subplot(2,1,2)
    # Histogram of failed assignments
    failed_interests_flat = [fac for sublist in failed_interests for fac in sublist]
    failed_counter = Counter(failed_interests_flat)
    failed_names = list(failed_counter.keys())
    failed_counts = [failed_counter[name] for name in failed_names]
    # Sort by counts descending
    sorted_failed_indices = np.argsort(-np.array(failed_counts))
    failed_names = [failed_names[i] for i in sorted_failed_indices]
    failed_counts = [failed_counts[i] for i in sorted_failed_indices]
    x2 = range(len(failed_names))
    plt.bar(x2, failed_counts, color="#C44E52", alpha=0.7)
    plt.xticks(x2, failed_names, rotation=90)
    plt.xlabel("Faculty")
    plt.ylabel("Number of Failed Assignments")
    plt.title("Faculty Interest of Students Who Could Not Be Assigned")


    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    filepath = r"C:\Users\mm96978\Downloads\FacInterest.xlsx"
    df = pd.read_excel(filepath, dtype=str, header=None)
    df.columns = ["col1", "col2"]

    names, counts1, df = compute_counts(df)
    print(df)

    monte_carlo(df)
    exit()

    # Plot and save figure
    plot_grouped_counts(names, counts1, title="Name counts by column")

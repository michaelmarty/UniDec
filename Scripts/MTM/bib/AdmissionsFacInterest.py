import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('WxAgg')
import matplotlib.pyplot as plt
from collections import Counter

# Default input file path (can be overridden by command-line argument)
filepath = r"C:\Users\mm96978\Downloads\FacInterest.xlsx"


def load_two_columns(path):
    """Load the Excel and return the first two columns as strings.
    If the file has headers, pandas will keep them; we select the first two columns by position.
    Returns a DataFrame with columns ['col1', 'col2'].
    """
    df = pd.read_excel(path, dtype=str, header=None)
    df2 = df.iloc[:, 0:2].copy()
    df2.columns = ["col1", "col2"]
    return df2


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
    c1 = Counter()
    c2 = Counter()

    # Iterate rows and count presence per row (deduplicate within the row)
    for _, row in df.iterrows():
        names_in_row_col1 = set(split_names(row["col1"]))
        names_in_row_col2 = set(split_names(row["col2"]))
        # union of names in this row across both columns
        names_in_row = names_in_row_col1 | names_in_row_col2
        print(names_in_row)
        # increment row-wise counter once per name present in this row
        for name in names_in_row:
            c_row[name] += 1

    # Unique names union
    all_names = set(c_row.keys()) | set(c1.keys()) | set(c2.keys())

    # Sort by row-wise counts descending
    names_sorted = sorted(all_names, key=lambda n: (-(c_row[n]), n))

    counts1 = [c_row[name] for name in names_sorted]

    counts_df = pd.DataFrame({
        "name": names_sorted,
        "row_count": counts1,
        "col1_count": [c1[name] for name in names_sorted],
        "col2_count": [c2[name] for name in names_sorted],
    })
    return names_sorted, counts1, counts_df


def plot_grouped_counts(names, counts1, outpath=None, title=None):
    """Create a grouped bar chart for counts. If counts2 is None, plot a single-series bar chart.
    Saves the figure to outpath if provided, and also shows it.
    """
    x = range(len(names))
    width = 0.35

    plt.figure(figsize=(max(8, len(names) * 0.4), 6))

    plt.bar(x, counts1, width=0.6, label="Row-wise count", color="#4C72B0")



    plt.xticks(x, names, rotation=90)
    plt.xlabel("Name")
    plt.ylabel("Count")
    if title:
        plt.title(title)
    # plt.legend()
    plt.tight_layout()

    # if outpath:
    #     plt.savefig(outpath)
    #     print(f"Saved plot to: {outpath}")
    plt.show()


def main(path=None):
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return

    try:
        df = load_two_columns(path)
    except Exception as e:
        print(f"Error reading Excel: {e}")
        return

    names, counts1, counts_df = compute_counts(df)

    if counts_df.empty:
        print("No names found in the first two columns.")
        return

    # # Save counts to CSV next to the input file
    # out_csv = os.path.splitext(path)[0] + "_name_counts.csv"
    # counts_df.to_csv(out_csv, index=False)
    # print(f"Saved counts CSV to: {out_csv}")

    # Plot and save figure
    out_png = os.path.splitext(path)[0] + "_name_counts.png"
    plot_grouped_counts(names, counts1, outpath=out_png, title="Name counts by column")


if __name__ == "__main__":
    # Allow optional path override: python AdmissionsFacInterest.py <path_to_excel>
    arg_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(arg_path or filepath)

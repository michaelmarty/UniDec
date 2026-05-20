import os
import pandas as pd

# -----------------------------
# SET WORKING DIRECTORY
# -----------------------------
topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Brain"
os.chdir(topdir)

# -----------------------------
# INPUT FILES
# -----------------------------
negfile = "20260506_AS_03_negBPE_50uM_UntDDA.csv"
posfile = "20260506_AS_03_posBPE_50uM_UntDDA.csv"

# -----------------------------
# HELPERS
# -----------------------------
def pick_col(df, options):
    for c in options:
        if c in df.columns:
            return c
    return None


def parse_msms(s):
    """
    MS-DIAL format:
    '184.0732:1200;104.107:300'
    """
    if pd.isna(s) or s == "":
        return []
    out = []
    for item in str(s).split(";"):
        try:
            mz, inten = item.split(":")
            out.append((float(mz), float(inten)))
        except:
            continue
    return out


def build_transition_df(df, molecule_list_name):
    mz_col = pick_col(df, ["Average Mz", "m/z", "Precursor m/z"])
    rt_col = pick_col(df, ["RT (min)", "Retention time", "RT"])

    # FIXED: hard map MS-DIAL Name column
    if "Name" not in df.columns:
        raise ValueError("Missing required column: 'Name'")
    name_col = "Name"

    adduct_col = pick_col(df, ["Adduct type"])
    formula_col = pick_col(df, ["Formula", "Precursor formula"])
    spec_col = pick_col(df, ["MS/MS spectrum", "MSMS spectrum", "MS2 spectrum", "Ref Spec"])

    transitions = []

    for _, row in df.iterrows():

        precursor_mz = row[mz_col]
        rt = row[rt_col] if rt_col else np.nan
        name = row[name_col]  # <-- NOW CORRECT
        adduct = row[adduct_col] if adduct_col else ""
        formula = row[formula_col] if formula_col else ""

        spectrum = row[spec_col] if spec_col else ""
        frags = parse_msms(spectrum)

        # retention time window (fixed default ±0.5 min)
        rt_window = 0.5

        # -----------------------------
        # if no MSMS → pseudo transition
        # -----------------------------
        if len(frags) == 0:
            transitions.append([
                molecule_list_name,
                name,
                formula,
                adduct,
                precursor_mz,
                1,
                name,
                formula,
                adduct,
                precursor_mz,
                1,
                rt,
                rt_window
            ])
        else:
            for mz, inten in frags:

                transitions.append([
                    molecule_list_name,
                    name,
                    formula,
                    adduct,
                    precursor_mz,
                    1,
                    name + "_frag",
                    "",                 # product formula unknown (standard)
                    "fragment",
                    mz,
                    1,
                    rt,
                    rt_window
                ])

    return pd.DataFrame(transitions, columns=[
        "Molecule List Name",
        "Molecule Name",
        "PrecursorFormula",
        "Precursor Adduct",
        "Precursor Mz",
        "Precursor Charge",
        "ProductName",
        "ProductFormula",
        "Product Adduct",
        "Product Mz",
        "Product Charge",
        "Explicit Retention Time",
        "Explicit Retention Time Window"
    ])


# -----------------------------
# LOAD DATA
# -----------------------------
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

# -----------------------------
# CLEANING (optional but same as yours)
# -----------------------------
if "Comment" in posdf.columns:
    posdf = posdf[~posdf["Comment"].str.contains("Low Quality", case=False, na=False)]
    negdf = negdf[~negdf["Comment"].str.contains("Low Quality", case=False, na=False)]

if "Ontology" in posdf.columns:
    posdf = posdf[~posdf["Ontology"].str.startswith("CL", na=False)]
    posdf = posdf[~posdf["Ontology"].str.startswith("PS", na=False)]

# -----------------------------
# BUILD
# -----------------------------
postldf = build_transition_df(posdf, "Brain_Lipid_POS")
negtldf = build_transition_df(negdf, "Brain_Lipid_NEG")

tldf = pd.concat([postldf, negtldf], ignore_index=True)

# -----------------------------
# RT WINDOW (optional tweak)
# -----------------------------
tldf["Explicit Retention Time Window"] = 0.5

# -----------------------------
# SAVE
# -----------------------------
tldf.to_csv("Combined_TransitionList.csv", index=False)

# -----------------------------
# STATS
# -----------------------------
print("Total transitions:", len(tldf))
print("Unique molecules:", tldf["Molecule Name"].nunique())
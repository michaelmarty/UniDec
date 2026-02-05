import os
import pandas as pd
from PublicScripts.Lipids.LipidFunctions import add_ref_spec_xml, class_network_analysis, network_plot, tail_network_analysis, merge_pos_neg, detect_outliers_class, get_unique_names
from PublicScripts.Lipids.LipidFragPrediction import assign_df_fragments
import matplotlib.pyplot as plt
import numpy as np

drop_columns = ["MS/MS assigned", "m/z matched", "Manually modified for annotation", "RT matched", "MS/MS matched",
                "Post curation result", "Fill %", "Manually modified for quantification",
                "Isotope tracking parent ID", "Isotope tracking weight number", "Spectrum reference file name",
                "INCHIKEY", "Annotation tag (VS1.0)", "Matched peaks percentage", "Matched peaks count",
                "RT similarity", "m/z similarity", "Formula"]

def comment_parser(df):
    if "Comment" not in df.columns:
        return df
    if df["Comment"].isnull().all():
        return df
    # See if any of the comments contain "Low quality"
    ucomments = df["Comment"].unique()
    if not any("Low quality" in str(c) for c in ucomments):
        return df
    # Create a new column "Low_Quality_Flag"
    df["Low_Quality_Flag"] = df["Comment"].apply(lambda x: True if "Low quality" in str(x) else False)
    return df

def outlier_pipeline(posfile, negfile, posxml="", negxml="", tol=0.05):
    print("Importing Files:", posfile, negfile)
    posdf = pd.read_csv(posfile)
    negdf = pd.read_csv(negfile)

    # Drop unnecessary columns if present
    for col in drop_columns:
        if col in posdf.columns:
            posdf = posdf.drop(columns=[col])
        if col in negdf.columns:
            negdf = negdf.drop(columns=[col])

    print("Adding Reference Spectra from XML if needed...")
    # Add in Reference Spectra if not already present
    if posxml != "" and "Ref Spec" not in posdf.columns:
        posdf = add_ref_spec_xml(posdf, posxml)
    if negxml != "" and "Ref Spec" not in negdf.columns:
        negdf = add_ref_spec_xml(negdf, negxml)

    print("Headgroup and Tail Fragment Checks...")
    # Perform Tail and Headgroup Fragment Checks
    posdf = assign_df_fragments(posdf, tol=tol)
    negdf = assign_df_fragments(negdf, tol=tol)

    print("MSMS Networking Analysis...")
    posdf, posgraphs = class_network_analysis(posdf, tol=tol)
    negdf, neggraphs = class_network_analysis(negdf, tol=tol)
    # headgraphs = posgraphs + neggraphs
    #
    posdf, posgraphs = tail_network_analysis(posdf, min_count=5, mz_tol=tol)
    negdf, neggraphs = tail_network_analysis(negdf, min_count=5, mz_tol=tol)
    # tailgraphs = posgraphs + neggraphs

    # Merge DFs
    rttol = 0.1
    posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

    # Combined Pos and Neg
    combined_df = pd.concat([posdf, negdf], ignore_index=True)

    # Outlier Detection
    rtdiff = combined_df['Average Rt(min)'] - combined_df['Reference RT']
    combined_df['RTdiff'] = rtdiff
    combined_df["LogSN"] = np.log10(combined_df["S/N average"] + 1)

    features = ['RTdiff', 'Total score', 'LogSN', 'Headgroup_Match', 'Tail_Match_Fraction',
                "Class_z_Anomaly_Score", "Tail_z_Anomaly_Score"]
    combined_df = detect_outliers_class(combined_df, features, separate_adducts=True, contribs=True)

    # Drop all columns that start with "T1" or so on
    cols_to_drop = [col for col in combined_df.columns if col.startswith(('T1', 'T2', 'T3', 'T4', 'T5'))]
    combined_df = combined_df.drop(columns=cols_to_drop)

    # Sort by Class and then Metabolite Name
    combined_df = combined_df.sort_values(by=["Ontology", "Metabolite name"])

    # Parse Comments if present
    combined_df = comment_parser(combined_df)

    # Write Outputs
    combined_df.to_excel("Combined_OutlierAnalysis.xlsx", index=False, sheet_name="Outlier_Analysis")

    namedf = get_unique_names(combined_df)
    # Write to separate sheet in the same Excel file
    with pd.ExcelWriter("Combined_OutlierAnalysis.xlsx", engine='openpyxl', mode='a') as writer:
        namedf.to_excel(writer, sheet_name='Unique_Names', index=False)


if __name__ == "__main__":

    topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK\New Libs Tests"
    posfile = "PosIDsRTP2.csv"
    negfile = "NegIDsRTP2.csv"
    posfile = "PosIDsC1_Ref.csv"
    negfile = "NegIDsC1_Ref.csv"
    posxml = ""
    negxml = ""

    os.chdir(topdir)
    outdf = outlier_pipeline(posfile, negfile, posxml, negxml)


import numpy as np
import pandas as pd
import molmass as mm
import unidec.IsoDec.match as match
import os

topdir = "C:\\Data\\Yuri\\"
os.chdir(topdir)

output_file = "C:\\Data\\Yuri\\known_masses.csv"

if os.path.exists(output_file):
    # import it
    df = pd.read_csv(output_file)

else:
    file = "C:\\Data\\Yuri\\CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_Toshare (1).xlsx"

    df = pd.read_excel(file)

    print("Columns:", df.columns)

    masses = []
    for index, row in df.iterrows():
        formula = row["Formula"]
        mass = mm.Formula(formula).monoisotopic_mass

        adduct = row["Adduct"]
        if type(adduct) == str:
            adduct_mass = mm.Formula(adduct).monoisotopic_mass
        else :
            adduct_mass = 0.0

        mass += adduct_mass
        masses.append(mass)


    df["Mass"] = masses
    df["m/z"] = (df["Mass"] + df["Charge"] * 1.007276467) / df["Charge"]

    # Sort df by masses
    df = df.sort_values(by="Mass")

    # Remove duplicates
    # df = df.drop_duplicates(subset=["Mass"])

    # Export to CSV
    df.to_csv(output_file, index=False)

print(len(df["Mass"]))

# Limit df to only those with "Ion" containing "c"
df = df[df["Ion"].str.contains("c", na=False)]

mc = match.df_to_matchedcollection(df, monoiso="Mass", peakmz = "m/z", z="Charge", scan=None, intensity=None, peakmass=None, ion=None)
print("Data MatchedCollection created with", len(mc), "peaks.")

path = "C:\\Data\\Yuri\\"

fname = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_100k_Charge Reduction-yes_adduct_N-term_MP2_pi200_IntGaus_NoP_N3_1_centroided"
fname = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms_1_centroided"
fname = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms_NT3p5_1_centroided"
fname = "CA2_Sim_c_34+_LTQ-FT 21T_Charge Reduction-no_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL0_1_centroided"

datafile = os.path.join(path, fname + "_unidecfiles", fname + ".tsv")

data = pd.read_csv(datafile, sep="\t")

dc = match.df_to_matchedcollection(data, ion=None, checkmultiplemonoisos=True)

shared, u1, u2 = match.compare_matchedcollections(mc, dc, ignorescan=True)

# Save TSV of matched collection for shared
shared_file = os.path.join(topdir, "shared_peaks.tsv")
shared.to_df().to_csv(shared_file, sep="\t", index=False)
print(f"Shared peaks saved to {shared_file}")

# Save TSV of unmatched peaks for u1 and u2
u1_file = os.path.join(topdir, "unmatched_known.tsv")
u2_file = os.path.join(topdir, "unmatched_decon.tsv")
u1.to_df().to_csv(u1_file, sep="\t", index=False)
u2.to_df().to_csv(u2_file, sep="\t", index=False)
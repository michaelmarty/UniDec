import os
import pandas as pd
import sqlite3

topdir = r"Z:\Group Share\Annika\Stellar"
os.chdir(topdir)
rtlibfile = "RT_Library.csv"

###
###
###
# posfile = "poslib.csv"
# negfile = "neglib.csv"
# posdf = pd.read_csv(posfile)
# negdf = pd.read_csv(negfile)
#
# # For posdf, keep only columns "compound_name" and "retention_time"
# posdf = posdf[["compound_name", "retention_time"]]
# posdf = posdf.rename(columns={"compound_name": "Molecule Name", "retention_time": "RT_Pos"})
# negdf = negdf[["compound_name", "retention_time"]]
# negdf = negdf.rename(columns={"compound_name": "Molecule Name", "retention_time": "RT_Neg"})
# # Merge posdf and negdf on "Molecule Name"
# rtdf = pd.merge(posdf, negdf, on="Molecule Name", how="outer")
#
# rtdf = msdialnames_to_LC(rtdf)
# print(rtdf)
# # Write merged RT library
# rtdf.to_csv(rtlibfile, index=False)
###
###
###

rtlib = pd.read_csv(rtlibfile)

infile = "LipidList_Comparison.csv"
outfile = "LipidList_Comparison_RT.csv"

df = pd.read_csv(infile)

df = pd.merge(df, rtlib, on="Molecule Name", how="left")


# Print maximum difference between RT_Pos and RT_Neg
df["RT_Diff"] = abs(df["RT_Pos"] - df["RT_Neg"])
print("Max RT difference:", df["RT_Diff"].max())

# Set RT as either RT_Pos or RT_Neg, whichever is not NaN. If both are not NaN, take the average
df["RT"] = df[["RT_Pos", "RT_Neg"]].mean(axis=1)

# Drop RT_Pos, RT_Neg, RT_Diff
df = df.drop(columns=["RT_Pos", "RT_Neg", "RT_Diff"])

# Drop files with NaN RT
df = df.dropna(subset=["RT"])

# Export to outfile
print(df)
df.to_csv(outfile, index=False)

# Write to iRTDb file, which is a SQLite tab-delimited file with IrtLibrary structure and columns "PeptideModSeq" and "Irt"
irdbfile = "Lipid2.irtdb"
irdb = pd.DataFrame()
irdb["PeptideModSeq"] = df["Molecule Name"]
irdb["Irt"] = df["RT"]
# Create column of "Standard" with all values as 1
irdb["Standard"] = 1
# Create column of "TimeSource" with all values as 1
irdb["TimeSource"] = 1
# Reset index
# irdb = irdb.reset_index()
# Set index column to "Id"
irdb = irdb.reset_index().rename(columns={"index": "Id"})
print(irdb)

# Connect to SQLite database (it will be created if it doesn't exist)
conn = sqlite3.connect(irdbfile)
# Write the DataFrame to a table named 'IrtLibrary'
irdb.to_sql('IrtLibrary', conn, if_exists='replace', index=False)
# Rename index column to 'Id'
# conn.execute('ALTER TABLE IrtLibrary RENAME COLUMN "index" TO "Id";')
# Close the connection
conn.close()




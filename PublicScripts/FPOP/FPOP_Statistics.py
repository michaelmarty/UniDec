import pandas as pd
import scipy.stats as stats
import os
import numpy as np
import pingouin as pg


# The Big Stats Tester: Tests for different types of things.
def stats_test(df, ratio=None, ltype=None, column="+Ox", type="levene", pcutoff=0.05):
    # Raise an error if neither ratio or lipid type are specified or if both are
    if (ratio is None and ltype is None) or (ratio is not None and ltype is not None):
        print("Error: Need to select either lipid type or ratio but not both", ltype, ratio)
        return -1

    # If ratios is specified, break out all lipid types that have that ratio and group them. Also, generate a DataFrame
    if ratio is not None:
        between = "Lipid Type"
        const = ratio
        df2 = df[df["Ratio"] == ratio]
        uvals = np.unique(df2["Lipid Type"])
        data_sets = []
        for l in uvals:
            df3 = df2[df2["Lipid Type"] == l]
            oxdat = df3[column].to_numpy()
            data_sets.append(oxdat)

    # If ltype is specified, break out all ratios that have that lipid type and group them
    elif ltype is not None:
        between = "Ratio"
        const = ltype
        df2 = df[df["Lipid Type"] == ltype]
        uvals = np.unique(df2["Ratio"])
        data_sets = []
        for l in uvals:
            df3 = df2[df2["Ratio"] == l]
            oxdat = df3[column].to_numpy()
            data_sets.append(oxdat)

    # Run the appropriate statistical test
    print("Ratio:", ratio, "\tLipid Type:", ltype)
    # Levene Test
    if type == "levene":
        stat, p = stats.levene(*data_sets, center="mean")
        print("Levene Pvalue:", p)
    # Bartlett Test
    elif type == "bartlett":
        stat, p = stats.bartlett(*data_sets)
        print("Bartlett Pvalue:", p)
    # Tukey HSD Test
    elif type == "tukey":
        results = pg.pairwise_tukey(df2, dv=column, between=between)
        results["Constant"] = const  # Set a constant column for reading the output
        # Filter the results to only statistically significant things
        p = results["p-tukey"]
        b1 = p.to_numpy() < pcutoff
        print("Tukey HSD Test:", results[b1])
        p = results[b1]
    # Games-Howell Test
    elif type == "games":
        results = pg.pairwise_gameshowell(df2, dv=column, between=between)
        results["Constant"] = const  # Set a constant column for reading the output
        # Filter the results to only statistically significant things
        p = results["pval"]
        b1 = p.to_numpy() < pcutoff
        print("Games-Howell Test:", results[b1])
        p = results[b1]
    # Shapiro-Wilk Test for Normality
    elif type == "norm":
        results = []
        for d in data_sets:
            if len(d) > 3:
                results.append(float(pg.normality(d, alpha=pcutoff)["pval"]))
            else:
                results.append(1)
        print("Normality Test:", results)
        p = np.array(results)
    else:
        print("Type Not Recognized:", type)
        return -1
    return p


def stats_flow(df, ratio=None, ltype=None, column="+Ox", pcutoff=0.05):
    # First, test for normality
    p0 = stats_test(df, ratio=ratio, ltype=ltype, column=column, type="norm", pcutoff=pcutoff)

    # Second, test variances
    # Test if any of the p values are significantly non-normal. If so, use the Levene test
    if np.any(p0 < pcutoff):
        p1 = stats_test(df, ratio=ratio, ltype=ltype, column=column, type="levene", pcutoff=pcutoff)
    # If all are normal or <3 data points, use bartlett
    else:
        p1 = stats_test(df, ratio=ratio, ltype=ltype, column=column, type="bartlett", pcutoff=pcutoff)

    # Third, test difference in means
    # If the same variances, use Tukey HSD
    if p1 > pcutoff:
        p2 = stats_test(df, ratio=ratio, ltype=ltype, column=column, type="tukey", pcutoff=pcutoff)
    # If different, use Games-Howell
    else:
        p2 = stats_test(df, ratio=ratio, ltype=ltype, column=column, type="games", pcutoff=pcutoff)

    return p2


# Set the file and directory
topdir = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\For Public Repository\\FPOP Data"
os.chdir(topdir)
file = "RawExtractedPercents.xlsx"

# Read File
df = pd.read_excel(file)

# Set List of Lipids and Ratios to Consider
lipids = ["DLPC", "DMPG", "DMPC", "POPC", "DPPC", "Ecoli", "Free"]
ratios = [3, 9, 18]

# Set the column and cutoff
cutoff = 0.05
column = "+Ox"

# Create Output
outdf = pd.DataFrame()
# Loop through all the ratios
for r in ratios:
    p2 = stats_flow(df, ratio=r, pcutoff=cutoff, column=column)
    outdf = pd.concat([outdf, p2])

# Same as for the ratios but looping over lipid types
for l in lipids:
    p2 = stats_flow(df, ltype=l, pcutoff=cutoff, column=column)
    outdf = pd.concat([outdf, p2])

# Output results
print(outdf)
outdf.to_excel("FPOP_Stats.xlsx")

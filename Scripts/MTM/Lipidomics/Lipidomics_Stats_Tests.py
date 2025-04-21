import pandas as pd
import scipy.stats as stats
import os
import numpy as np
import pingouin as pg


# The Big Stats Tester: Tests for different types of things.
def stats_test(df, type2=None, type1=None, column="Intensity", stats_type="levene", pcutoff=0.05, name1="Temperature",
               name2="Belt", silent=True):
    if (type2 is not None and type1 is not None) and (stats_type == "tukey" or stats_type == "games"):
        print("Error: Cannot run Tukey or Games-Howell tests on both type1 and type2. Please select only one.")
        return -1


    data_sets = []
    between = ""
    const = ""

    if type1 is None and type2 is None:
        uvals1 = np.unique(df[name1])
        uvals2 = np.unique(df[name2])
        for l1 in uvals1:
            df2 = df[df[name1] == l1]
            for l2 in uvals2:
                df3 = df2[df2[name2] == l2]
                data_sets.append(df3[column].to_numpy())
        between = name1
        df2 = df
    # If type2s is specified, break out all lipid types that have that type2 and group them. Also, generate a DataFrame
    elif type2 is not None:
        between = name1
        const = type2
        df2 = df[df[name2] == type2]
        uvals = np.unique(df2[name1])
        for l in uvals:
            df3 = df2[df2[name1] == l]
            oxdat = df3[column].to_numpy()
            data_sets.append(oxdat)

    # If type1 is specified, break out all type2s that have that lipid type and group them
    elif type1 is not None:
        between = name2
        const = type1
        df2 = df[df[name1] == type1]
        uvals = np.unique(df2[name2])
        data_sets = []
        for l in uvals:
            df3 = df2[df2[name2] == l]
            oxdat = df3[column].to_numpy()
            data_sets.append(oxdat)

    # Run the appropriate statistical test
    if not silent:
        print("type2:", type2, "\ttype1:", type1)
    # Levene Test
    if stats_type == "levene":
        stat, p = stats.levene(*data_sets, center="mean")
        if not silent:
            print("Levene Pvalue:", p)
    # Bartlett Test
    elif stats_type == "bartlett":
        stat, p = stats.bartlett(*data_sets)
        if not silent:
            print("Bartlett Pvalue:", p)
    # Tukey HSD Test
    elif stats_type == "tukey":
        results = pg.pairwise_tukey(df2, dv=column, between=between)
        results["Constant"] = const  # Set a constant column for reading the output
        # Rename p-tukey column to pval
        results = results.rename(columns={"p-tukey": "pval"})
        if not silent:
            print("Tukey HSD Test:", results)
        p = results
    # Games-Howell Test
    elif stats_type == "games":
        results = pg.pairwise_gameshowell(df2, dv=column, between=between)
        results["Constant"] = const  # Set a constant column for reading the output
        if not silent:
            print("Games-Howell Test:", results)
        p = results
    # Shapiro-Wilk Test for Normality
    elif stats_type == "norm":
        results = []
        for d in data_sets:
            if len(d) > 3:
                results.append(float(pg.normality(d, alpha=pcutoff)["pval"]))
            else:
                results.append(1)
        if not silent:
            print("Normality Test:", results)
        p = np.array(results)
    else:
        print("Type Not Recognized:", type)
        return -1

    return p


def stats_equal_variance(df, type1=None, type2=None, column="Intensity", pcutoff=0.05, test_type=None,
                         name1="Temperature", name2="Belt"):
    if test_type is None:
        # First, test for normality
        p0 = stats_test(df, type1=type1, type2=type2, column=column, stats_type="norm", pcutoff=pcutoff, name1=name1,
                        name2=name2)

        # Second, test variances
        # Test if any of the p values are significantly non-normal. If so, use the Levene test
        if np.any(p0 < pcutoff):
            test_type = "levene"
        # If all are normal or <=3 data points, use bartlett
        else:
            test_type = "bartlett"
    p1 = stats_test(df, type1=type1, type2=type2, column=column, stats_type=test_type, pcutoff=pcutoff, name1=name1,
                    name2=name2)
    return p1


def post_hoc(df, type1=None, type2=None, column="Intensity", pcutoff=0.05, p1=None, name1="Temperature",
             name2="Belt"):
    if p1 is None:
        p1 = stats_equal_variance(df, type1=type1, type2=type2, column=column, pcutoff=pcutoff, name1=name1,
                                  name2=name2)

    # Third, test difference in means
    # If the same variances, use Tukey HSD
    if p1 > pcutoff:
        t2 = "tukey"
    # If different, use Games-Howell
    else:
        t2 = "games"
    p2 = stats_test(df, type1=type1, type2=type2, column=column, stats_type=t2, pcutoff=pcutoff, name1=name1,
                    name2=name2)

    # print("Stats Completed:", type1, type2, t1, t2)
    return p2

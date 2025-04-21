import os
from Lipidomics_Stats_Tests import *

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_Brain_Classes.xlsx"
path = "NA_Ecoli_Classes.xlsx"

df = pd.read_excel(path, sheet_name="Reformat")

# Strip off the class from the name
classes = df["Class Name"]
# Get unique classes while preserving order
indexes = np.unique(classes, return_index=True)[1]
uclasses = np.array([classes[index] for index in sorted(indexes)])

int_names = ["Intensity", "Tail Length", "Tail Unsaturation"]

for int_name in int_names:
    # int_name = "Tail Unsaturation"

    temps = np.unique(df["Temperature"])
    belts = np.unique(df["Belt"])

    full_results = []

    anova_p_critical = 0.1
    # Belt Comparisons
    for c in uclasses:
        subdf = df[df["Class Name"] == c]

        p_var = stats_equal_variance(subdf, column=int_name, test_type="levene")
        print(c, p_var)
        if np.isnan(p_var):
            continue

        if p_var < anova_p_critical:
            print("Variance is not equal")
            bresults = pg.welch_anova(subdf, dv=int_name, between="Belt")
            tresults = pg.welch_anova(subdf, dv=int_name, between="Temperature")
        else:
            print("Variance is equal")
            bresults = pg.anova(subdf, dv=int_name, between="Belt")
            tresults = pg.anova(subdf, dv=int_name, between="Temperature")
        print(bresults, tresults)
        p_belt = bresults["p-unc"][0]
        p_temp = tresults["p-unc"][0]
        print(p_belt, p_temp)
        if p_belt < anova_p_critical:
            for t in temps:
                results = post_hoc(subdf, type1=t, column=int_name)
                if len(results) > 0:
                    for i, row in results.iterrows():
                        a = row["A"]
                        b = row["B"]
                        pval = row["pval"]
                        print(c, t, a, b, pval)
                        full_results.append([c, t, a, b, pval])

        if p_temp < anova_p_critical:
            for belt in belts:
                subdf = df[df["Class Name"] == c]
                results = post_hoc(subdf, type2=belt, column=int_name)
                if len(results) > 0:
                    for i, row in results.iterrows():
                        a = row["A"]
                        b = row["B"]
                        pval = row["pval"]
                        print(c, belt, a, b, pval)
                        full_results.append([c, belt, a, b, pval])

    full_results = np.array(full_results)
    if len(full_results) == 0:
        print("No Significance Found", int_name)
        full_results = pd.DataFrame()
    else:
        pvals = full_results[:, -1].astype(float)
        b1 = pvals < anova_p_critical
        full_results = full_results[b1]
        # adjusted_pvals = pg.multicomp(pvals[b1], method="fdr_bh")[1]

        full_results = pd.DataFrame(full_results, columns=["Class Name", "Constant", "A", "B", "pval"])
        # full_results["Adjusted pval"] = adjusted_pvals

        print(full_results)

    writer = pd.ExcelWriter(path, engine='openpyxl', mode='a', if_sheet_exists="replace")
    full_results.to_excel(writer, index=False, sheet_name=int_name + " Stats")
    writer.close()

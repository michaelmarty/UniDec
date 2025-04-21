import numpy as np
import pandas as pd
import os

lipids = ["PS", "PC", "PE"]

affs = [0.7, 1.0, 1.1, 2.7, 4.1, 4.3]

R = 8.31446261815324 / 1000.
units = "(kJ/mol)"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"
T = 0 + 273.15
dgs = - R * T * np.log(affs)

print(dgs)


T = 22 + 273.15



path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\Papers\\LX-MS\\230517_LX_Supplemental_Data_3.xlsx"

df = pd.read_excel(path, sheet_name="Table S5 - Brain Lipid Results", skiprows=1)

os.chdir(os.path.dirname(path))

comps = ["Elu_vs_FT", "AChR_vs_Ctrl", "Elu_vs_FT.1"]
titles = ["Control Exchange", "AChR Assembly", "AChR LX"]

n_annulus = 45  # Ellena JF, Blazing MA, McNamee MG. Lipid-protein interactions in reconstituted membranes containing acetylcholine receptor. Biochemistry. 1983;22(24):5523-35. doi: 10.1021/bi00293a012.
n_disc = 200  # CharmmGUI
n_empty = 240

bulkratio = (n_disc - n_annulus) / n_disc
print(bulkratio)
outdf = pd.DataFrame()
outdf["Name"] = df["Name"].values

exit()
for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    fc = df[fcname].values

    fc2 = 2 ** fc
    # print(fc2)
    if i != 0:
        fc2 = (1 / (1 - bulkratio)) * fc2 - (bulkratio / (1 - bulkratio))
        # print(fc2)
        bneg = fc2 < 0
        if np.any(bneg):
            # print("Negative Values:", fc2[bneg])
            pass
        fc2[bneg] = 1e-12
        fc2 = - R * T * np.log(fc2)
        fc2[bneg] = np.inf
    else:
        fc2 = - R * T * np.log(fc2)

    data = fc2

    outname = titles[i] + " Delta G " + units
    outdf[outname] = data

with pd.ExcelWriter(path, mode='a', if_sheet_exists="replace") as writer:
    outdf.to_excel(writer, sheet_name="Table S7 - Brain Thermodynamics", index=False)


df = pd.read_excel(path, sheet_name="Table S4 - Ecoli Lipid Results", skiprows=1)


comps = ["FT_vs_Elu", "AmtB_vs_Ctrl", "AqpZ_vs_Ctrl", "FT_vs_Elu.1", "FT_vs_Elu.2"]
titles = ["Control Exchange", "AmtB Assembly", "AqpZ Assembly", "AmtB LX", "AqpZ LX"]

n_annulus = 75  # Marty, 2016 Ange. Chem.
n_disc = 150  # Keener, 2019 JACS
n_empty = 240

bulkratio = (n_disc - n_annulus) / n_disc
print(bulkratio)

bulkratio = (n_disc - n_annulus) / n_disc
print(bulkratio)
outdf = pd.DataFrame()
outdf["Name"] = df["Name"].values

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    fc = df[fcname].values

    fc2 = 2 ** fc
    # print(fc2)
    if i != 0:
        fc2 = (1 / (1 - bulkratio)) * fc2 - (bulkratio / (1 - bulkratio))
        # print(fc2)
        bneg = fc2 < 0
        if np.any(bneg):
            # print("Negative Values:", fc2[bneg])
            pass
        fc2[bneg] = 1e-12
        fc2 = - R * T * np.log(fc2)
        fc2[bneg] = np.inf
    else:
        fc2 = - R * T * np.log(fc2)

    data = fc2

    outname = titles[i] + " Delta G " + units
    outdf[outname] = data

with pd.ExcelWriter(path, mode='a', if_sheet_exists="replace") as writer:
    outdf.to_excel(writer, sheet_name="Table S6 - Ecoli Thermodynamics", index=False)

print(outdf)


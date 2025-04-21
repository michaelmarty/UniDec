import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib



from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 11
rcParams['font.sans-serif'] = "Arial"

df = pd.read_excel("Z:\\Group Share\\Deseree\\FPOP Paper 2\\DMPG.xlsx")


def grouped_barplot(df, cat,subcat, val , err):
    u = df[cat].unique()
    x = np.arange(len(u))
    subx = df[subcat].unique()
    offsets = (np.arange(len(subx))-np.arange(len(subx)).mean())/(len(subx)+1.)
    width= np.diff(offsets).mean()
    for i,gr in enumerate(subx):
        dfg = df[df[subcat] == gr]
        plt.bar(x + offsets[i], dfg[val].values, width=width,
                label="{} {}".format(subcat, gr), yerr=dfg[err].values)
    plt.xlabel(cat)
    plt.ylabel(val)
    plt.xticks(x, u)
    plt.legend()
    plt.savefig("Z:\\Group Share\\Deseree\\FPOP Paper 2\\DMPGtest.pdf")
    plt.show()


cat = "Stoichiometry"
subcat = "Ratio"
val = "Mean"
err = "Stdev"

#call the function with df from the question
grouped_barplot(df, cat, subcat, val, err )

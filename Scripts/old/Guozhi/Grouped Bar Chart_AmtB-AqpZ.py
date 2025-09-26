import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams.update({'errorbar.capsize': 3})
df = pd.read_excel("Z:\\Group Share\\Guozhi Zhang\\2D_LC\\20230226LX_AmtB_AqpZ\\230309_data_analysis\\Python\\230418_group-Barchart\\230418_AmtB_AqpZ_LX.xlsx")


def grouped_barplot(df, cat,subcat, val , err):
    u = df[cat].unique()
    x = np.arange(len(u))
    subx = df[subcat].unique()
    offsets = (np.arange(len(subx))-np.arange(len(subx)).mean())/(len(subx)+1.)
    width= np.diff(offsets).mean()
    for i,gr in enumerate(subx):
        dfg = df[df[subcat] == gr]
        plt.bar(x+offsets[i], dfg[val].values, width=width,
                label="{} {}".format(subcat, gr), yerr=dfg[err].values)
    plt.xlabel(cat)
    plt.ylabel(val)
    plt.xticks(x, u)
    plt.legend(fontsize=5)
    plt.ylim(-1,1)
    plt.axhline(y = 0.07, color = 'black', linestyle = '--')
    plt.axhline(y=-0.07, color='black', linestyle='--')
    plt.axhline(y=0, color='black', linestyle='-')
    plt.show()
plt.rcParams.update({'font.size': 15})


cat = "MP"
subcat = "Class"
val = "Tail Unsat. Diff."
err = "T.Stdev"

#call the function with df from the question
grouped_barplot(df, cat, subcat, val, err )
#plt.savefig()

import os
from Lipidomics_Stats_Tests import *
from statsmodels.stats.multitest import multipletests

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
path = "NA_Brain_Classes.xlsx"

# Load the data
df = pd.read_excel(path, sheet_name="Class Comparison Pvals")

for i in range(1, 10):
    pvals = df.to_numpy()[:, i].flatten()
    pvals = np.sort(pvals)

    alpha = 0.1

    results = multipletests(pvals, alpha=alpha, method="bonferroni")
    b1 = results[0]

    x = np.arange(len(pvals))
    y = alpha / (len(x) - x)
    c = np.sum(1/(x+1))
    y = x * alpha / (len(x) * c)
    #b1 = pvals < y

    if np.sum(b1) == 0:
        print("None Significant")
    else:
        max_pos = np.amax(x[b1])
        print(i, max_pos, pvals[max_pos])


'''
plt.plot(x, pvals)
plt.plot(x, y, color="red")
#plt.xlim(0,10)
#plt.ylim(0,0.1)
plt.show()'''
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("WxAgg")

directory = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\CurrentData"
os.chdir(directory)

filename = "FINAL_HEK_widewindow_Isolist.csv"
filename = "CleanedIsoList_HEK_1.csv"

df = pd.read_csv(filename)
print(len(df))

plt.figure(figsize=(6, 4))
for i, row in df.iterrows():
    rt = row["Retention Time (min)"]
    rtwindow = row["RT Window (min)"]
    mz = row["m/z"]
    polarity = row["Polarity"]

    if polarity == "Positive":
        color = "red"
    else:
        color = "blue"

    plt.hlines(y=mz, xmin=rt-rtwindow/2, xmax=rt+rtwindow/2, color=color, alpha=0.5)
    plt.scatter(rt, mz, color=color, edgecolor="black", zorder=5)

# Add legend with positive red and negative blue
plt.scatter([], [], color="red", edgecolor="black", label="Positive", zorder=5)
plt.scatter([], [], color="blue", edgecolor="black", label="Negative", zorder=5)
plt.legend(title="Polarity")
plt.title("PRM Isolation Targets")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Retention Time (min)")
# Make this italic
plt.ylabel("m/z", fontstyle="italic")
plt.tight_layout()

plt.savefig("IsolationListPlot.png", dpi=600, transparent=True)
plt.savefig("IsolationListPlot.pdf")
plt.show()


import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
import pandas as pd

def get_scan_number(path):
    imp = ImporterFactory.create_importer(path)
    return imp.get_max_scan()

topdir = r"Z:\Group Share\Annika\IQX vs Stellar DDA\Raw Data"
os.chdir(topdir)
files = ud.match_files(topdir, "*.raw")
if False:
    rows = []
    for f in files:
        if "blank" in f.lower():
            continue
        if "pos" in f.lower():
            polarity = "Positive"
        elif "neg" in f.lower():
            polarity = "Negative"
        else:
            polarity = "Unknown"

        if "ITIT" in f:
            type = "ITIT"
        elif "OTOT" in f:
            type = "OTOT"
        elif "OTIT" in f:
            type = "OTIT"
        elif "UntDDA" in f:
            type = "Stellar"
        else:
            type = "Unknown"

        scan_number = get_scan_number(f)
        max_time = ImporterFactory.create_importer(f).get_max_time()

        rows.append({"File": os.path.basename(f), "Scan Number": scan_number, "Max Time": max_time, "Polarity": polarity, "Type": type})

    df = pd.DataFrame(rows)
    df.to_csv("Scan_Numbers.csv", index=False)
else:
    df = pd.read_csv("Scan_Numbers.csv")

rows = df.to_dict(orient="records")

unique_types = set(row["Type"] for row in rows)
# Sort by the order of ITIT, OTOT, OTIT, Stellar
type_order = ["OTOT", "OTIT", "ITIT", "Stellar"]
unique_types = sorted(unique_types, key=lambda x: type_order.index(x) if x in type_order else len(type_order))
for t in unique_types:
    type_rows = [row for row in rows if row["Type"] == t]
    pos_rows = [row for row in type_rows if row["Polarity"] == "Positive"]
    neg_rows = [row for row in type_rows if row["Polarity"] == "Negative"]
    pos_avg = np.mean([row["Scan Number"] for row in pos_rows]) if pos_rows else 0
    neg_avg = np.mean([row["Scan Number"] for row in neg_rows]) if neg_rows else 0
    print(f"{t} - Positive Avg Scan Number: {pos_avg}, Negative Avg Scan Number: {neg_avg}")
    pos_std = np.std([row["Scan Number"] for row in pos_rows]) if pos_rows else 0
    neg_std = np.std([row["Scan Number"] for row in neg_rows]) if neg_rows else 0

    total_avg = pos_avg + neg_avg
    total_std = np.sqrt(pos_std**2 + neg_std**2)
    plt.bar(t, total_avg, yerr=total_std, label=t)

plt.xlabel("Type")
plt.ylabel("Number of Scans")
plt.title("Distribution of Scan Numbers by Type")
plt.legend()
plt.savefig("Scan_Number_Distribution.png", dpi=300)
plt.savefig("Scan_Number_Distribution.pdf")
plt.show()





import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import rcParams

# Set backend to Agg
mpl.use('WxAgg')
rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['font.family'] = 'Arial'

file = "Z:\\Group Share\\JGP\\IsoDec_figs\\uvpd_etd_roc_data.csv"

roc_df = pd.read_csv(file)


engines = ["8bit-Python", "4bit-Python",  "IsoDec-THRASH"]
engine_labels = ["IsoDec 8-bin", "IsoDec 4-bin", "IsoDec-THRASH"]
# colors = ['b', 'r', 'orange', 'm', 'y', 'k', 'olive']
colors = ['darkblue', "skyblue", "orange"]

recalls = [[] for i in range(len(engines))]
precisions = [[] for i in range(len(engines))]
f1_scores = [[] for i in range(len(engines))]


for i, engine in enumerate(engines):
    recalls[i] = roc_df[engine + "-Recall"].values
    precisions[i] = roc_df[engine + "-Precision"].values
    f1_scores[i] = (2 * (precisions[i] * recalls[i]) / (precisions[i] + recalls[i]))

#Find the top F1 score and the precision and recall at that point
for i in range(len(engines)):
    best_f1 = np.max(f1_scores[i])
    best_f1_index = np.argmax(f1_scores[i])
    print(engine_labels[i], "Best F1:", best_f1, "Recall:", recalls[i][best_f1_index], "Precision:", precisions[i][best_f1_index])

topFD_precision = 0.42505
topFD_recall = 0.69099

flashdeconv_precision = 0.632
flashdeconv_recall = 0.238

eTHRASH_precision = 0.349
eTHRASH_recall = 0.721

id8_precision = 0.649
id8_recall = 0.730

plt.figure(figsize=(7,4))
# Create ROC plot
for i, engine in enumerate(engines):
    plt.plot(recalls[i], precisions[i], linestyle='-', color=colors[i], marker='', label=engine_labels[i])
    #Also label the best F1 score
    # best_f1 = np.max(f1_scores[i])
    # best_f1_index = np.argmax(f1_scores[i])
    # plt.scatter(recalls[i][best_f1_index], precisions[i][best_f1_index], color='k')


plt.scatter([flashdeconv_recall], [flashdeconv_precision], color='yellow', edgecolor='k', label="FLASHDeconv", clip_on=False)
plt.scatter([topFD_recall], [topFD_precision], color='r', edgecolor='k', label="TopFD")
plt.scatter([eTHRASH_recall], [eTHRASH_precision], color='g', edgecolor='k', label="e-THRASH")
plt.scatter([id8_recall], [id8_precision], color='darkblue', edgecolor='k', label="ID-8")

# plt.plot(0.0266,1,marker='o', color='black', markersize=5)
# plt.text(0.0266,1.05, "A", fontsize=16)
#
# plt.plot(0.5344,0.8407,marker='o', color='black', markersize=5)
# plt.text(0.5344,0.8907, "B", fontsize=16)
#
# plt.plot(0.731,0.650,marker='o', color='black', markersize=5)
# plt.text(0.731,0.700, "C", fontsize=16)

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

plt.legend(loc=[0.01,0.01], fontsize=12, frameon=False, ncol=2)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Recall", fontsize=18)
plt.ylabel("Precision", fontsize=18)
plt.xlim(0, 1.0)
plt.ylim(0, 1.0)
plt.xticks(np.arange(0,1.1,0.5), fontsize=14)
plt.yticks(np.arange(0,1.1,0.5), fontsize=14)
plt.tight_layout()
if True:
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\uvpd_etd_roc_fromfile.pdf", dpi=300, transparent=True)
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\uvpd_etd_roc_fromfile.png", dpi=300, transparent=True)
plt.show()

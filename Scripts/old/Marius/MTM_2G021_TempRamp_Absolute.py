import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file ="C:\\Users\\marty\\University of Arizona\martylab - Documents\Papers\Thermal Unfolding of AAVs\Excel sheet\\2G021 no normalization.xlsx"

df = pd.read_excel(file)
t = df["temperature"].to_numpy()

things = ["E", "P", "F"]
for i, thing in enumerate(things):
    keys = [thing+" 1", thing+" 2", thing+" 3"]

    plt.subplot(1, 3, i+1)
    es = df[keys].to_numpy()
    esnorm = []
    for e in es.transpose():
        norme = e/e[0]
        esnorm.append(norme)
        plt.plot(norme)
    esnorm = np.array(esnorm)
plt.show()





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def count_fac(df, years):
    starts = df["Start"].to_numpy()
    ends = df["End"].to_numpy()
    b1 = ends=="Current"
    ends[b1]=np.ones(len(df))[b1]*np.amax(years)
    tots = []
    for y in years:
        b1 = starts <= y
        b2 = ends >= y
        b3 = b1*b2
        tots.append(np.sum(b3))
    return np.array(tots)

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnFac.xlsx"

df = pd.read_excel(filepath)

years = [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023]
faculty = np.zeros_like(years)
totals = np.zeros_like(years)

unis = np.unique(df["Institution"])

per_year_data = []

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

for uni in unis:
    df2 = df[df["Institution"]==uni]
    ytots = []
    for y in years:
        tot = np.sum(df2[str(y)+"_docs"])
        ytots.append(tot)
    ytots = np.array(ytots)
    counts = count_fac(df2, years)
    ax1.plot(years, ytots, label=uni)
    ax2.plot(years, ytots/counts, label=uni)
ax1.set_title("Papers per year")
ax2.legend()
ax2.set_title("Paper per fac per year")
plt.show()


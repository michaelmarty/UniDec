import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def count_fac(df, years):
    starts = df["Start"].to_numpy()
    ends = df["End"].to_numpy()
    b1 = ends == "Current"
    ends[b1] = np.ones(len(df))[b1] * np.amax(years)
    tots = []
    for y in years:
        b1 = starts <= y
        b2 = ends >= y
        b3 = b1 * b2
        tots.append(np.sum(b3))
    return np.array(tots)


filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnalChemFac.xlsx"

df = pd.read_excel(filepath)

years = [2017, 2018, 2019, 2020, 2021, 2022, 2023]
faculty = np.zeros_like(years)
totals = np.zeros_like(years)

uniarray=df["Affiliation"].to_numpy().astype(str)
unis = np.unique(uniarray)

per_year_data = []

data = []
totarray = []

print(len(unis))
for uni in unis:
    df2 = df[df["Affiliation"] == uni]
    ytots = []
    for y in years:
        tot = np.sum(df2[str(y) + "_docs"])
        ytots.append(tot)
    ytots = np.array(ytots)
    counts = len(df2)
    unitot = np.sum(ytots)
    data.append([counts, unitot, unitot/counts/len(years)])
    totarray.append(ytots)



sindex = np.argsort(np.array(data)[:, -1])[::-1]
data = np.array(data)[sindex]
totarray = np.array(totarray)[sindex]
unis = unis[sindex]

# Combine all this data into a single dataframe and export to excel
df2 = pd.DataFrame(data, columns=["Faculty", "Total", "Papers per Fac per Year"])
df2["University"] = unis
df2["Papers per year"] = data[:, 1] / len(years)
# reorder with uni first
df2 = df2[["University", "Faculty", "Total", "Papers per Fac per Year", "Papers per year"]]
df2.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnalChemFacStats.xlsx", index=False)



ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
for i, d in enumerate(data[:15]):
    ax1.plot(years, totarray[i], label=unis[i])
    ax2.plot(years, np.array(totarray[i]) / float(d[0]), label=unis[i])

ax1.set_title("Papers per year")
ax2.legend()
ax2.set_title("Paper per fac per year")
plt.show()

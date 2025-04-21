import pybliometrics.scopus as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\FacultyDatabase2.xlsx"
    df = pd.read_excel(filepath)
    years = np.array([2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023])
    # years = np.array([2018, 2019, 2020])

    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    for index, row in df.iterrows():
        if row["End"] == "Current":
            ydat = []
            name = row["Last Name"]
            for i, y in enumerate(years):
                p = row[str(y) + "_docs"]
                #if y == 2022:
                #    p *= 2
                ydat.append(p)
            ydat = np.array(ydat)
            try:
                fit = np.polyfit(years - years[0], ydat / np.amax(ydat), 1)
            except:
                fit = [-0.5, 0]
            print(name, fit)

            scaling = 10
            if fit[0] >= 0:
                ax1.plot(years, ydat, label=name, linewidth = fit[0]*scaling)
            else:
                ax2.plot(years, ydat, label=name, linewidth = -fit[0]*scaling)
    ax1.legend()
    ax2.legend()
    ax1.set_xlabel("Year")
    ax2.set_xlabel("Year")
    ax1.set_ylabel("Total Pubs")
    ax1.set_title("Positive Slope")
    ax2.set_title("Negative Slope")
    plt.show()

import pybliometrics.scopus as sc
import pandas as pd
import numpy as np
from AnFacGetID import *

if __name__ == '__main__':

    filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnalChemFac.xlsx"

    df = pd.read_excel(filepath)

    years = [2017, 2018, 2019, 2020, 2021, 2022, 2023]
    faculty = np.zeros_like(years) + len(df)
    totals = np.zeros_like(years)

    for index, row in df.iterrows():
        id = row["Author ID"]
        city = row["City"]
        uni = row["Affiliation"]
        affid = row["Affiliation ID"]
        if type(id) == type("str") or type(id) == type(1):
            print(id, index)
            au = sc.AuthorRetrieval(id)

            affdat = au.affiliation_current[0]
            if affdat.parent_preferred_name is not None:
                unifix = affdat.parent_preferred_name
            else:
                unifix = affdat.preferred_name
            affidfix = affdat.id
            df.loc[index, "Affiliation"] = unifix
            df.loc[index, "Affiliation ID"] = affidfix

            if uni != unifix:
                print(uni, unifix, au.affiliation_current)

            docs = pd.DataFrame(au.get_documents(refresh=10))
            docs['year'] = docs['coverDate'].str[:4]
            counts = docs['year'].value_counts()
            citedby = np.sum(docs["citedby_count"].to_numpy(dtype=int))
            print(citedby)
            df.loc[index, "scopus_citedby"] = citedby
            df.loc[index, "scopus_documents"] = len(docs)
            countarray = []
            tot = 0
            for i, y in enumerate(years):
                c = 0
                for j, d in docs.iterrows():
                    if d['year'] == str(y):
                        if check_affil(d, city, uni, affid):
                            c += 1
                            tot += 1
                            totals[i] += 1
                countarray.append(c)
            df.loc[index, "Total"] = tot
            for i, y in enumerate(years):
                df.loc[index, str(y) + "_docs"] = countarray[i]

    print(faculty)
    print(totals)

    try:
        avg = np.array(totals) / np.array(faculty)
        print(avg)
    except:
        pass
    df.to_excel(filepath, index=False)

    import matplotlib.pyplot as plt

    plt.plot(years, avg)
    plt.xlabel("Year")
    plt.ylabel("Avg. # Papers/Fac")
    plt.show()

import pybliometrics.scopus as sc
import pandas as pd
import numpy as np


def check_affil(doc):
    cities = doc["affiliation_city"]
    unis = doc["affilname"]
    fails = 0
    try:
        if "Tucson" in cities:
            return True
    except:
        fails += 1
        pass
    try:
        if "Arizona" in unis:
            return True
    except:
        fails += 1
        pass
    try:
        if "60010065" in doc['afid'] or "60010065" in doc['author_afids']:
            return True
    except:
        fails += 1
        pass
    if fails == 3:
        return True
    return False


if __name__ == '__main__':

    filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\FacultyDatabase2.xlsx"
    # filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\CTandStaffDatabase.xlsx"
    # filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\PostdocDatabase.xlsx"
    #filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase2.xlsx"

    df = pd.read_excel(filepath)

    years = [2017, 2018, 2019, 2020, 2021, 2022]
    faculty = np.zeros_like(years)
    totals = np.zeros_like(years)

    for index, row in df.iterrows():
        id = row["scopus_id"]
        start = float(row["Start"])
        end = row["End"]
        if end == "Current":
            end= 2023
        else:
            end = float(end)
        print(end)
        for i, y in enumerate(years):
            if y >= start and y <= end:
                faculty[i] += 1

        if type(id) == type("str"):
            id = id.split('-')[2]
            print(id)
            au = sc.AuthorRetrieval(id)
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
                        if check_affil(d):
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

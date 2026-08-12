import pybliometrics.scopus as sc
import pandas as pd
import numpy as np
import os

# http://api.elsevier.com/content/abstract/doi/10.1016/S0014-5793(01)03313-0?apiKey=07070f5eff2d979e7792ebcc2096daf6

sc.init(keys=["07070f5eff2d979e7792ebcc2096daf6"])
print(sc.get_keys())

def get_affil_name(id):
    if id != "":
        aff = sc.AffiliationRetrieval(id)
        return aff.affiliation_name
    else:
        return None

def get_coauthors(id, startyear=2020, endyear=None):
    au = sc.AuthorRetrieval(id)
    docs = pd.DataFrame(au.get_documents(refresh=10))
    docs['year'] = docs['coverDate'].str[:4]

    coauthorlist = []
    affilist = []

    for i, row in docs.iterrows():
        year = int(row['year'])
        if startyear is not None:
            if year < startyear:
                continue
        if endyear is not None:
            if year > endyear:
                continue
        authornames = row['author_names'].split(";")
        try:
            authoraffs = row['author_afids'].split(";")
        except Exception as e:
            print("ERROR:", e, row)
            continue

        affnames = [get_affil_name(affid) for affid in authoraffs]
        #print(affnames)

        for a in authornames:
            coauthorlist.append(a)

        for a in affnames:
            affilist.append(a)


    uniquenames, uniqueindexes = np.unique(coauthorlist, return_index=True)
    uniqueaffs = np.array(affilist)[uniqueindexes]

    sortindex = np.argsort(uniquenames)

    uniquenames = uniquenames[sortindex]
    uniqueaffs = uniqueaffs[sortindex]

    #print(uniquenames)
    #print(uniqueaffs)

    df = pd.DataFrame()
    df["Name"] = uniquenames
    df["Affiliation"] = uniqueaffs
    #return coauthors
    return df


if __name__ == '__main__':
    id = "9-s2.0-53984653400"
    #id = "9-s2.0-8302292300"
    #id = "9-s2.0-56909034900"
    #id = "9-s2.0-55444380400"
    #id = "9-s2.0-8765172500"
    #id = "9-s2.0-51161148500"
    #id = "9-s2.0-57191446043"
    #id = "9-s2.0-13008328400"
    coauthors = get_coauthors(id, startyear=2021, endyear=2026)
    print(coauthors)
    path = r"C:\Users\mm96978\OneDrive - The University of Texas at Austin\Career"
    os.chdir(path)
    coauthors.to_excel("MMcoauthors2021-2026.xlsx", index=False)

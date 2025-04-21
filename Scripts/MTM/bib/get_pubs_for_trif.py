import pybliometrics.scopus as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import os

def get_affil_name(id):
    if id != "":
        aff = sc.AffiliationRetrieval(id)
        return aff.affiliation_name
    else:
        return None

def get_papers(id, startyear=2019, endyear=None, fy=True):
    au = sc.AuthorRetrieval(id)
    docs = pd.DataFrame(au.get_documents(refresh=10))
    docs['year'] = docs['coverDate'].str[:4]

    publists = []


    for i, row in docs.iterrows():
        year = int(row['year'])
        month = int(row['coverDate'][5:7])
        if startyear is not None:
            if year < startyear:
                continue
            elif year == startyear and month < 7 and fy:
                continue
        if endyear is not None:
            if year > endyear:
                continue
            elif year == endyear and month > 6 and fy:
                continue

        doi = row['doi']

        title = row['title']
        journal = row['publicationName']

        publists.append([title, doi, "y", journal, month, year])

    df = pd.DataFrame(publists, columns=["Publication Title","DOI #","Academic Peer Reviewed Journal (Enter Y or N)",
                                         "Journal Title", "Published Month", "Published Year"])

    return df


if __name__ == '__main__':
    id = "9-s2.0-53984653400"
    papers = get_papers(id, startyear=2022, endyear=2023, fy=True)
    print(papers)
    #path = "C:\\Users\\marty\\OneDrive - University of Arizona\\Grants\\Collaborations\\Michael Johnson"
    #os.chdir(path)
    papers.to_excel("trifpapersFY23.xlsx", index=False)

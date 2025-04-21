from scholarly import scholarly as sc
import pandas as pd

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\FacultyDatabase.xlsx"

df = pd.read_excel(filepath)

for index, row in df.iterrows():
    name = row["First Name"] + " " + row["Last Name"]
    sq = sc.search_author(name)
    match = False
    i = 0
    id = ""
    citedby=""
    while match is False and i < 5:
        try:
            result = next(sq)
        except:
            id = ""
            citedby = ""
            break
        affiliation = result['affiliation']
        email = result['email_domain']
        if "arizona" in affiliation.lower() or "arizona" in email.lower():

            id = result['scholar_id']
            citedby = result['citedby']
            match = True
        i += 1
    df.loc[index, "scholar_id"] = id
    df.loc[index, "citedby"] = citedby

df.to_excel(filepath, index=False)

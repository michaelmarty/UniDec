import pybliometrics.scopus as sc
import pandas as pd
import time

def check_affil(doc):
    cities = doc["affiliation_city"]
    unis = doc["affilname"]
    fails=0
    try:
        if "Tucson" in cities:
            return True
    except:
        fails+=1
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

    return False

def find_affil(docs, years):
    for i, y in enumerate(years):
        for j, d in docs.iterrows():
            if d['year'] == str(y):
                if check_affil(d):
                    return True
    return False


filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\FacultyDatabase.xlsx"
filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\CTandStaffDatabase.xlsx"
#filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\PostdocDatabase.xlsx"
#filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase2.xlsx"

df = pd.read_excel(filepath)

for index, row in df.iterrows():
    sid = row["scopus_id"]
    if not pd.isna(sid):
        continue

    first_name = row["First Name"]
    last_name = row["Last Name"]
    name = 'AUTHFIRST(' + first_name + ") and AUTHLAST(" + last_name + ")"
    print(name)
    sq = sc.AuthorSearch(name).authors
    match = False
    i = 0
    id = ""
    numdocuments = ""
    if sq is not None and len(sq) > 0:
        for s in sq:
            try:
                aff = s.affiliation.lower()
            except:
                aff = ""
            if "arizona" in aff:
                id = s.eid
                numdocuments = s.documents
                print(s, numdocuments)
                match = True
                break
    if not match:
        time.sleep(0.1)
        first_name = row["First Name"]
        first_name = first_name.split()[0]
        last_name = row["Last Name"]
        last_name = last_name.split()[0]
        name = 'AUTHFIRST(' + first_name + ") and AUTHLAST(" + last_name + ")"
        print(name)
        sq = sc.AuthorSearch(name).authors
        if sq is not None and len(sq) > 0:
            for s in sq:
                try:
                    aff = s.affiliation.lower()
                except:
                    aff = ""
                if "arizona" in aff:
                    id = s.eid
                    numdocuments = s.documents
                    print(s, numdocuments)
                    match = True
                    break
        if not match:
            years = [2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021]
            if sq is not None and len(sq) > 0:
                for s in sq:
                    newid = s.eid
                    newid = newid.split('-')[2]
                    au = sc.AuthorRetrieval(newid)
                    docs = pd.DataFrame(au.get_documents(refresh=10))
                    docs['year'] = docs['coverDate'].str[:4]
                    match = find_affil(docs, years)
                    if match:
                        id = s.eid
                        numdocuments = s.documents
                        print(s, numdocuments)
                        match = True
                        break
    df.loc[index, "scopus_id"] = id
    df.loc[index, "scopus_documents"] = numdocuments
    time.sleep(0.1)

df.to_excel(filepath, index=False)

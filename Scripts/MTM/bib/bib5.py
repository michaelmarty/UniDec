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

    filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\FacultyDatabase.xlsx"
    # filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\CTandStaffDatabase.xlsx"
    # filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\PostdocDatabase.xlsx"
    # filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase2.xlsx"

    df = pd.read_excel(filepath,engine='openpyxl')

    years = [2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022]
    journal_list = []
    journal_ids = []
    docids = []
    affs = []
    aff_countries = []
    titles=[]

    for index, row in df.iterrows():
        id = row["scopus_id"]
        start = float(row["Start"])
        end = float(row["End"])

        if type(id) == type("str"):
            id = id.split('-')[2]
            print(id)
            au = sc.AuthorRetrieval(id)
            docs = pd.DataFrame(au.get_documents(refresh=10))
            docs['year'] = docs['coverDate'].str[:4]
            month = docs['coverDate'].str[5:7]


            for i, y in enumerate(years):
                for j, d in docs.iterrows():
                    if d['year'] == str(y):
                        if check_affil(d):
                            journal = d['publicationName']
                            jid = d["issn"]
                            eid = d["eid"]
                            title = d["title"]
                            if eid not in docids:
                                docids.append(eid)
                                journal_list.append(journal)
                                journal_ids.append(jid)
                                titles.append(title)
                                try:
                                    affils = np.unique(d['affilname'].split(';'))
                                    countries = np.unique(d['affiliation_country'].split(';'))
                                except:
                                    affils = []
                                    countries = []
                                for a in affils:
                                    if a != 'The University of Arizona':
                                        affs.append(a)
                                for c in countries:
                                    if c != 'United States':
                                        aff_countries.append(c)

    df2 = pd.DataFrame()
    df2['Titles'] = pd.Series(titles)
    df2.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\Titles.xlsx",
                 index=False)

    df2 = pd.DataFrame()
    df2['Journals'] = pd.Series(journal_list)
    df2["ISSNs"] = pd.Series(journal_ids)
    df2.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\Journals.xlsx",
                 index=False)

    hist = df2['Journals'].value_counts()
    df3 = pd.DataFrame(hist)
    df3.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\JournalCounts.xlsx",
                 index=True)
    hist = hist[hist > 3]

    df4 = pd.DataFrame()
    df4['Affs'] = pd.Series(affs)
    df4.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\Affiliations.xlsx",
                 index=False)
    hist = df4['Affs'].value_counts()
    df5 = pd.DataFrame(hist)
    df5.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AffCounts.xlsx",
                 index=True)
    hist = hist[hist > 3]

    df6 = pd.DataFrame()
    df6['Countries'] = pd.Series(aff_countries)
    df6.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\Countries.xlsx",
                 index=False)
    hist = df6['Countries'].value_counts()
    df7 = pd.DataFrame(hist)
    df7.to_excel("C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\CountryCounts.xlsx",
                 index=True)

    hist = hist[hist > 0]

    import re
    tcloud = ""
    for t in titles:
        t = re.sub("<.*?>", "", str(t))
        tcloud = tcloud + " " + t.lower()

    import matplotlib.pyplot as plt
    from wordcloud import WordCloud
    wc = WordCloud().generate(tcloud)
    plt.imshow(wc, interpolation="bilinear")
    plt.axis("off")
    plt.show()


    exit()



    x = np.arange(len(hist))
    values = hist.to_numpy()
    print(values)
    journals = hist.keys().to_list()
    print(journals)
    rect = [0.1, 0.6, 0.8, 0.3]
    fig = plt.figure()
    ax = plt.axes(rect)
    ax.bar(x, values)
    ax.xaxis.set_major_locator(plt.FixedLocator(x))
    ax.xaxis.set_ticklabels(journals, rotation=90)
    # plt.xlabel("Year")
    plt.ylabel("Number of Unique Papers (2014-2021)")
    plt.tight_layout()
    plt.show()

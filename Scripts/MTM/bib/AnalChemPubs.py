import pandas as pd
from pybliometrics.scopus import SerialSearch, ScopusSearch
s = SerialSearch(query={"title": "Journal of the American Society for Mass Spectrometry"})
df = pd.DataFrame(pd.DataFrame(s.results))

j = df.iloc[0]
print(j)
acid = j["source-id"]
issn = j["issn"]
eissn = j["eIssn"]
print(acid, issn, eissn)


s2 = ScopusSearch("ISSN("+ issn + ") AND PUBYEAR AFT 2017")
df2 = pd.DataFrame(pd.DataFrame(s2.results))

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\JASMS.xlsx"
df2.to_excel(filepath)
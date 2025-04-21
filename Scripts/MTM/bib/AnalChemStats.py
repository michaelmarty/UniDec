import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnalChem.xlsx"

df = pd.read_excel(file)

print(df.head())

print(df.columns)

# remove amp; from affilname text
df["affilname"] = df["affilname"].str.replace("amp;", "")

df = df.dropna(subset = ["affilname"])
df["lastcountry"] = df["affiliation_country"].str.split(";").str[-1]
# Select only country is US
df = df[df["lastcountry"] == "United States"]

# Select only subtypeDescription is "Article"
df = df[df["subtype"] == "ar"]

# Strip out last affilname
df["lastaff"] = df["affilname"].str.split(";").str[-1]

# Strip out last affid
df["lastaffid"] = df["afid"].str.split(";").str[-1]

# Last City
df["lastcity"] = df["affiliation_city"].str.split(";").str[-1]

# Strip out last author_id
df["lastauthor_id"] = df["author_ids"].str.split(";").str[-1]

# Strip out last affilname
df["lastauthor"] = df["author_names"].str.split(";").str[-1]

# Get unique last author indexes
uindexes = np.unique(df["lastauthor_id"], return_index = True)[1]
df2 = df.iloc[uindexes]

# Drop all but specific columns
df2 = df2[["lastauthor", "lastauthor_id", "lastaff", "lastaffid", "lastcity"]]

#rename columns
df2 = df2.rename(columns = {"lastauthor": "Author", "lastauthor_id": "Author ID",
                            "lastaff": "Affiliation", "lastaffid": "Affiliation ID", "lastcity": "City"})

# Drop affiliations that don't have [School, College, University] in the name
df2 = df2[df2["Affiliation"].str.contains("School|College|University|Institute|Department|Center")]

# Sort by affiliation
df2 = df2.sort_values(by = ["Affiliation"])



print(df2)

# Save to excel
filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\AnalChemFac.xlsx"
df2.to_excel(filepath)


exit()
# Create affiliation histogram
affils = df["lastaff"].value_counts()
print(affils)

# filter affils above 10
affils = affils[affils > 15]

plt.bar(affils.index, affils.values)
# rotate tick labels 90
plt.xticks(rotation = 90)
plt.tight_layout()
plt.show()